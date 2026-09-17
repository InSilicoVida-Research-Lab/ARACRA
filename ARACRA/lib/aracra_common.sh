#!/usr/bin/env bash
# =============================================================================
#  ARACRA Pipeline — lib/aracra_common.sh
#
#  Single source of truth for values that setup.sh, run_app.sh and the
#  Streamlit app all need to agree on: environment name, default port,
#  database/work locations, pinned tool versions, and the contents of .env.
#
#  Source it, don't execute it:
#      source "$(dirname "${BASH_SOURCE[0]}")/lib/aracra_common.sh"
# =============================================================================

# ── Identity ─────────────────────────────────────────────────────────────────
ARACRA_ENV_NAME="${ARACRA_ENV_NAME:-test_ARACRA}"
# Legacy environment names still accepted by run_app.sh, newest first.
ARACRA_ENV_FALLBACKS=(test_ARACRA rnaseq_pipeline2 rnaseq_pipeline)

# ── Ports and locations ──────────────────────────────────────────────────────
# 8501 is Streamlit's own default and the port quoted throughout the README.
ARACRA_DEFAULT_PORT="${ARACRA_DEFAULT_PORT:-8501}"
ARACRA_WORK_ROOT="${ARACRA_WORK_ROOT:-${HOME}/aracra_star_work}"
ARACRA_DB_DIR="${ARACRA_DB_DIR:-${HOME}/databases}"

# ── Version policy ───────────────────────────────────────────────────────────
# Deliberately light. Two kinds of constraint only:
#
#   HARD CEILING  — used where a newer version is known to break the pipeline.
#                   Only Nextflow qualifies: main.nf uses 24.x DSL2 semantics.
#   SOFT FLOOR    — ">=" on everything else. Stops an ancient build sneaking in,
#                   but never forces a downgrade on someone whose stack is newer
#                   than this file. A downgrade is its own failure mode, and for
#                   an end user it is a worse one than a warning.
#
# Nothing here is an exact pin, so a solve will not start failing just because
# a channel dropped an old build. If the pinned solve fails anyway, setup.sh
# retries unpinned and says so — the install always completes.
ARACRA_NEXTFLOW_SPEC="${ARACRA_NEXTFLOW_SPEC:->=24.04,<25}"
# Salmon < 1.10 hits the GLIBC_2.34 failure in the troubleshooting section.
ARACRA_SALMON_SPEC="${ARACRA_SALMON_SPEC:->=1.10}"

# Soft floors: package -> spec. Kept in one place so there is a single file to
# edit when a floor needs raising.
ARACRA_TOOL_SPECS=(
    "openjdk=17"
    "star>=2.7.10"
    "hisat2>=2.2.1"
    "bowtie2>=2.5"
    "samtools>=1.17"
    "fastp>=0.23"
    "subread>=2.0.3"
    "rseqc>=5.0"
    "multiqc>=1.19"
    "sra-tools>=3.0"
)

# Minimum versions the app's tool check warns below. Same numbers as the floors
# above, in a form Python can read, so the GUI and the installer cannot disagree.
ARACRA_MIN_VERSIONS="star=2.7.10 hisat2=2.2.1 samtools=1.17 fastp=0.23 salmon=1.10 subread=2.0.3 multiqc=1.19 nextflow=24.04"

# ── Derived reference paths (given a DB_DIR) ─────────────────────────────────
# Usage: aracra_ref_paths "$DB_DIR"  → sets ARACRA_REF_* in the caller's shell.
aracra_ref_paths() {
    local db="${1:-$ARACRA_DB_DIR}"
    ARACRA_REF_DIR="${db}/hg38_reference"
    ARACRA_ANN_DIR="${db}/annotations"
    ARACRA_SCREEN_DIR="${db}/fastq_screen_genomes"

    ARACRA_REF_STAR_INDEX="${ARACRA_REF_DIR}/star_index_hg38"
    ARACRA_REF_HISAT2_INDEX="${ARACRA_REF_DIR}/hisat2_index_hg38/genome_tran"
    ARACRA_REF_SALMON_INDEX="${ARACRA_REF_DIR}/salmon_index_hg38"
    ARACRA_REF_GTF="${ARACRA_REF_DIR}/gencode.v44.primary_assembly.annotation.gtf"
    ARACRA_REF_TX2GENE="${ARACRA_REF_DIR}/tx2gene.csv"
    ARACRA_REF_BED="${ARACRA_ANN_DIR}/hg38_RefSeq.bed"
    ARACRA_REF_HK_BED="${ARACRA_ANN_DIR}/hg38_housekeeping.bed"
    ARACRA_REF_REFFLAT="${ARACRA_ANN_DIR}/hg38_refFlat.txt"
    ARACRA_REF_SCREEN_CONF="${ARACRA_SCREEN_DIR}/FastQ_Screen_Genomes/fastq_screen.conf"
}

# ── System discovery ─────────────────────────────────────────────────────────
aracra_cpu_cores() {
    nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4
}

aracra_ram_gb_int() {
    local kb
    kb=$(awk '/MemTotal/{print $2}' /proc/meminfo 2>/dev/null || echo 0)
    awk "BEGIN{printf \"%d\", ${kb}/1048576}"
}

# STAR needs ~32 GB to hold the hg38 index; below that recommend HISAT2.
aracra_recommended_aligner() {
    local ram="${1:-$(aracra_ram_gb_int)}"
    if [ "$ram" -lt 32 ]; then echo "hisat2"; else echo "star"; fi
}

# ── Conda discovery ──────────────────────────────────────────────────────────
# Echoes the conda base path, or nothing if none found.
aracra_find_conda() {
    local p
    for p in "${HOME}/miniforge3" "${HOME}/miniconda3" "${HOME}/anaconda3"; do
        if [ -f "${p}/etc/profile.d/conda.sh" ]; then echo "$p"; return 0; fi
    done
    p="$(conda info --base 2>/dev/null || true)"
    if [ -n "$p" ] && [ -f "${p}/etc/profile.d/conda.sh" ]; then echo "$p"; return 0; fi
    return 1
}

# ── .env writer ──────────────────────────────────────────────────────────────
# THE only place .env is written. Both setup.sh and run_app.sh call this, so
# the key set can never drift between them. aracra_star_app.py reads:
#   CONDA_BASE ENV_NAME ENV_PATH ENV_BIN WORK_DIR OUT_DIR LOG_FILE
#   STAR_INDEX HISAT2_INDEX SALMON_INDEX GTF_PATH BED_PATH HK_BED_PATH
#   REFFLAT_PATH SCREEN_CONF_PATH RECOMMENDED_ALIGNER
#
# Usage: aracra_write_env <pipeline_dir> <conda_base> <env_path> <db_dir> [origin]
aracra_write_env() {
    local pipeline_dir="$1" conda_base="$2" env_path="$3"
    local db_dir="${4:-$ARACRA_DB_DIR}" origin="${5:-setup.sh}"

    aracra_ref_paths "$db_dir"

    local cores ram aligner
    cores="$(aracra_cpu_cores)"
    ram="$(aracra_ram_gb_int)"
    aligner="$(aracra_recommended_aligner "$ram")"

    cat > "${pipeline_dir}/.env" <<EOF
# ARACRA Pipeline — auto-generated by ${origin} on $(date)
# System: ${ram} GB RAM | ${cores} CPUs | Recommended aligner: ${aligner}
CONDA_BASE=${conda_base}
ENV_NAME=$(basename "$env_path")
ENV_PATH=${env_path}
ENV_BIN=${env_path}/bin
WORK_DIR=${ARACRA_WORK_ROOT}/work
OUT_DIR=${ARACRA_WORK_ROOT}/results
LOG_FILE=${ARACRA_WORK_ROOT}/pipeline.log
DB_DIR=${db_dir}
STAR_INDEX=${ARACRA_REF_STAR_INDEX_OVERRIDE-$ARACRA_REF_STAR_INDEX}
HISAT2_INDEX=${ARACRA_REF_HISAT2_INDEX}
SALMON_INDEX=${ARACRA_REF_SALMON_INDEX}
GTF_PATH=${ARACRA_REF_GTF}
BED_PATH=${ARACRA_REF_BED}
HK_BED_PATH=${ARACRA_REF_HK_BED}
REFFLAT_PATH=${ARACRA_REF_REFFLAT}
SCREEN_CONF_PATH=${ARACRA_REF_SCREEN_CONF}
TX2GENE_PATH=${ARACRA_REF_TX2GENE}
TOTAL_RAM_GB=${ram}
CPU_CORES=${cores}
RECOMMENDED_ALIGNER=${aligner}
HALF_CORES=$(( cores / 2 > 0 ? cores / 2 : 1 ))
APP_PORT=${ARACRA_DEFAULT_PORT}
EOF
}


# ── Version record ───────────────────────────────────────────────────────────
# Writes ONE human-readable file listing what actually got installed. Not a
# lockfile and not used to rebuild anything — it exists so that "it worked last
# month" can be checked against "it fails now" without guesswork.
aracra_write_versions() {
    local pipeline_dir="$1" env_path="$2"
    local out="${pipeline_dir}/aracra_versions.txt"
    local bin="${env_path}/bin"

    {
        echo "ARACRA environment record — $(date)"
        echo "host    : $(uname -srm)"
        echo "env     : ${env_path}"
        echo ""
        echo "[tools]"
        for t in nextflow STAR hisat2 samtools fastp featureCounts salmon \
                 multiqc qualimap fastq_screen prefetch Rscript streamlit; do
            if [ -x "${bin}/${t}" ]; then
                printf "  %-14s %s\n" "$t" \
                    "$("${bin}/${t}" --version 2>&1 | head -1 | cut -c1-60)"
            else
                printf "  %-14s (not installed)\n" "$t"
            fi
        done
        echo ""
        echo "[R]"
        "${bin}/Rscript" -e '
            cat("  R            ", as.character(getRversion()), "\n")
            if (requireNamespace("BiocManager", quietly=TRUE))
                cat("  Bioconductor ", as.character(BiocManager::version()), "\n")
            for (p in c("DRomics","DESeq2","edgeR","sva","rtracklayer","org.Hs.eg.db"))
                if (requireNamespace(p, quietly=TRUE))
                    cat(sprintf("  %-13s %s\n", p, as.character(packageVersion(p))))
        ' 2>/dev/null
    } > "$out"
    echo "$out"
}
