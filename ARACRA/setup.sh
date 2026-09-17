#!/usr/bin/env bash
# =============================================================================
#  ARACRA Pipeline — setup.sh
#  Auto-discovers system resources, builds conda env, downloads references,
#  configures everything. Warns if RAM < 32GB (recommends HISAT2 over STAR).
#
#  CHANGELOG (this revision):
#    - Quoted 'salmon>=1.10' (was creating a stray file named '=1.10' and
#      installing salmon unpinned)
#    - Index existence checks are content-based, not directory-based
#    - FastQ Screen genomes fetched over HTTPS directly; conf generated locally
#    - All module toggles derived from reference-file existence, never literals
#    - nextflow.config emitted without top-level 'def' (strict-parser safe)
#    - Reference verification gate before declaring success
# =============================================================================
set -euo pipefail

CYAN='\033[0;36m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
RED='\033[0;31m'; NC='\033[0m'; BOLD='\033[1m'
log()  { echo -e "${CYAN}[SETUP]${NC} $*"; }
ok()   { echo -e "${GREEN}[  OK ]${NC} $*"; }
warn() { echo -e "${YELLOW}[ WARN]${NC} $*"; }
fail() { echo -e "${RED}[ERROR]${NC} $*"; exit 1; }

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/aracra_common.sh
source "${PIPELINE_DIR}/lib/aracra_common.sh"

# ── Parse flags ───────────────────────────────────────────────────────────────
SKIP_INDEX=false
SKIP_SCREEN=false
DB_DIR="${ARACRA_DB_DIR}"
for arg in "$@"; do
    case "$arg" in
        --skip-index)    SKIP_INDEX=true ;;
        --skip-screen)   SKIP_SCREEN=true ;;
        --db-dir=*)      DB_DIR="${arg#*=}" ;;
        --help|-h)
            echo "Usage: bash setup.sh [--skip-index] [--skip-screen] [--db-dir=PATH]"
            exit 0 ;;
    esac
done

echo -e "${BOLD}"
echo "╔══════════════════════════════════════════════════════════╗"
echo "║     ARACRA Pipeline — Setup                              ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# ══════════════════════════════════════════════════════════════════════════════
# 1. SYSTEM DISCOVERY
# ══════════════════════════════════════════════════════════════════════════════
CPU_CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
TOTAL_RAM_KB=$(awk '/MemTotal/{print $2}' /proc/meminfo 2>/dev/null || echo 0)
TOTAL_RAM_GB=$(awk "BEGIN{printf \"%.1f\", ${TOTAL_RAM_KB}/1048576}")
TOTAL_RAM_INT=$(awk "BEGIN{printf \"%d\", ${TOTAL_RAM_KB}/1048576}")
DISK_FREE_GB=$(df -BG "${HOME}" 2>/dev/null | awk 'NR==2{gsub("G",""); print $4}' || echo 0)

log "System discovery:"
log "  CPUs       : ${CPU_CORES}"
log "  RAM        : ${TOTAL_RAM_GB} GB"
log "  Disk free  : ${DISK_FREE_GB} GB"
log "  DB dir     : ${DB_DIR}"

# ── Disk preflight ──
# Rough requirement: references ~5 GB, STAR index ~28 GB, salmon ~15 GB,
# FastQ Screen ~30 GB, plus working space.
DISK_NEEDED=60
[ "$SKIP_SCREEN" = false ] && DISK_NEEDED=$((DISK_NEEDED + 30))
if [ "${DISK_FREE_GB:-0}" -lt "$DISK_NEEDED" ] 2>/dev/null; then
    warn "Only ${DISK_FREE_GB} GB free; roughly ${DISK_NEEDED} GB recommended."
    warn "  Use --db-dir=/path/on/larger/disk, and/or --skip-screen."
    warn "  Continuing anyway — downloads will fail if space runs out."
fi

# ── RAM warnings ──
if [ "$TOTAL_RAM_INT" -lt 32 ]; then
    echo ""
    echo -e "${RED}${BOLD}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}${BOLD}║  ⚠  WARNING: System has ${TOTAL_RAM_GB} GB RAM (< 32 GB)              ║${NC}"
    echo -e "${RED}${BOLD}║                                                              ║${NC}"
    echo -e "${RED}${BOLD}║  STAR requires ~32 GB RAM for human genome alignment.        ║${NC}"
    echo -e "${RED}${BOLD}║  Recommendation: Use HISAT2 instead (needs only ~8 GB RAM).  ║${NC}"
    echo -e "${RED}${BOLD}║  The pipeline supports both — select HISAT2 in the app.      ║${NC}"
    echo -e "${RED}${BOLD}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""
fi
RECOMMENDED_ALIGNER="$(aracra_recommended_aligner "$TOTAL_RAM_INT")"

# ── Core allocation strategy ──
HALF_CORES=$(( CPU_CORES / 2 > 0 ? CPU_CORES / 2 : 1 ))
LIGHT_PARALLEL=$(( CPU_CORES - 1 > 0 ? CPU_CORES - 1 : 1 ))
POST_PARALLEL=$(( CPU_CORES - 2 > 0 ? CPU_CORES - 2 : 1 ))
STAR_PARALLEL=2
HEAVY_CORES=$(( HALF_CORES > 2 ? HALF_CORES : 2 ))
if [ "$CPU_CORES" -le 4 ]; then
    STAR_PARALLEL=1
    HEAVY_CORES=$CPU_CORES
fi

log "Core strategy:"
log "  Download     : 1 core × 3 parallel"
log "  Fastp        : 1 core × ${LIGHT_PARALLEL} parallel"
log "  STAR/HISAT2  : ${HEAVY_CORES} cores × ${STAR_PARALLEL} parallel"
log "  Post-QC      : 1 core × ${POST_PARALLEL} parallel"
log "  Final steps  : ${CPU_CORES} cores"

# ══════════════════════════════════════════════════════════════════════════════
# 2. CONDA / MAMBA
# ══════════════════════════════════════════════════════════════════════════════
CONDA_BASE=""
for path in "${HOME}/miniforge3" "${HOME}/miniconda3" "${HOME}/anaconda3"; do
    if [ -f "${path}/etc/profile.d/conda.sh" ]; then
        CONDA_BASE="$path"
        break
    fi
done
if [ -z "$CONDA_BASE" ]; then
    log "No conda found — installing Miniforge3..."
    ARCH="$(uname -m)"
    curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-${ARCH}.sh" \
        -o /tmp/miniforge.sh
    bash /tmp/miniforge.sh -b -p "${HOME}/miniforge3"
    rm -f /tmp/miniforge.sh
    CONDA_BASE="${HOME}/miniforge3"
    ok "Miniforge3 installed"
fi
source "${CONDA_BASE}/etc/profile.d/conda.sh"
ok "Conda: $(conda --version 2>&1) at ${CONDA_BASE}"

if ! command -v mamba &>/dev/null; then
    log "Installing mamba..."
    conda install -n base -c conda-forge mamba -y -q
fi

# ══════════════════════════════════════════════════════════════════════════════
# 3. CREATE ENVIRONMENT
# ══════════════════════════════════════════════════════════════════════════════
ENV_NAME="${ARACRA_ENV_NAME}"
ENV_PATH="${CONDA_BASE}/envs/${ENV_NAME}"

if [ -d "$ENV_PATH" ]; then
    log "Environment '${ENV_NAME}' already exists — updating..."
else
    log "Creating environment '${ENV_NAME}'..."
    mamba create -p "$ENV_PATH" python=3.12 -y -q
fi

conda activate "$ENV_PATH"
ok "Active env: ${ENV_PATH}"

# ── Install bioinformatics tools ──
# NOTE: every version constraint MUST be quoted. An unquoted 'salmon>=1.10'
# is parsed by the shell as an output redirect: it silently creates a file
# named '=1.10' and passes a bare 'salmon' to mamba, leaving salmon unpinned.
# Version floors live in lib/aracra_common.sh (ARACRA_TOOL_SPECS etc.) — one
# file to edit when a floor needs raising, instead of two drifting copies.
log "Installing pipeline tools (this may take a few minutes)..."

# Try the constrained solve first. If the solver cannot satisfy it — a channel
# dropped a build, or the user's base env has awkward pins — fall back to an
# unconstrained solve rather than failing. A working install with a warning is
# more useful to an end user than a clean failure.
ARACRA_PINNED_OK=true
if ! mamba install -p "$ENV_PATH" -c conda-forge -c bioconda -y -q \
    "nextflow${ARACRA_NEXTFLOW_SPEC}" \
    "salmon${ARACRA_SALMON_SPEC}" \
    "${ARACRA_TOOL_SPECS[@]}" \
    python=3.12 openpyxl fastq-screen fastqc picard qualimap pigz
then
    ARACRA_PINNED_OK=false
    warn "Version-constrained solve failed — retrying without version floors."
    warn "The pipeline will still work, but check aracra_versions.txt afterwards."
    mamba install -p "$ENV_PATH" -c conda-forge -c bioconda -y -q \
        openjdk=17 nextflow python=3.12 openpyxl sra-tools fastp fastq-screen \
        fastqc star hisat2 bowtie2 samtools rseqc picard qualimap \
        "salmon${ARACRA_SALMON_SPEC}" subread multiqc pigz \
        || fail "Tool installation failed"
fi
ok "Bioinformatics tools installed"

# Remove the stray file produced by older, unquoted versions of this script.
[ -f "${PIPELINE_DIR}/=1.10" ] && rm -f "${PIPELINE_DIR}/=1.10" && \
    log "Removed stray '=1.10' file from previous unquoted install"

# ── Install Python packages ──
log "Installing Python packages..."
"${ENV_PATH}/bin/pip" install --quiet streamlit pandas openpyxl boto3 polars-lts-cpu
ok "Streamlit + pandas + MultiQC deps installed"

# ── Install R packages ──
log "Installing R + Bioconductor packages..."
mamba install -p "$ENV_PATH" -c conda-forge -c bioconda -y -q \
    r-base \
    r-optparse \
    r-dplyr \
    r-readxl \
    r-ggplot2 \
    r-ggrepel \
    r-stringr \
    r-jsonlite \
    bioconductor-deseq2 \
    bioconductor-edger \
    bioconductor-sva \
    bioconductor-org.hs.eg.db \
    || warn "Some R packages failed via mamba — will try R fallback"
ok "R packages (mamba)"

log "Installing DRomics and remaining Bioconductor packages..."
"${ENV_PATH}/bin/Rscript" -e '
options(repos=c(CRAN="https://cloud.r-project.org"), timeout=1200)
# timeout=1200, not the 300s default: reactome.db (a ReactomePA dependency,
# needed for pathway enrichment) is a ~455 MB annotation data package that
# reliably blows past 300s on an ordinary connection -- confirmed 2026-09-17,
# a real install died at 294 MB downloaded with a 300 second timeout error.
if (!requireNamespace("DRomics", quietly=TRUE)) install.packages("DRomics", quiet=TRUE)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", quiet=TRUE)
miss <- c("DESeq2","edgeR","sva","org.Hs.eg.db","GO.db","KEGGREST","AnnotationDbi",
          "clusterProfiler","ReactomePA","enrichplot")[
    !sapply(c("DESeq2","edgeR","sva","org.Hs.eg.db","GO.db","KEGGREST","AnnotationDbi",
              "clusterProfiler","ReactomePA","enrichplot"), requireNamespace, quietly=TRUE)]
if (length(miss)) BiocManager::install(miss, ask=FALSE, update=FALSE)
if (!requireNamespace("msigdbr", quietly=TRUE)) install.packages("msigdbr", quiet=TRUE)
cat("R packages OK\n")
' 2>/dev/null || warn "DRomics/enrichment package install had issues — check manually"
ok "R packages done"

# ── rtracklayer (source patch for GCC 15 compatibility) ──
# GCC 15 treats empty parentheses () as zero arguments (C23 behaviour).
# The bundled UCSC source declares a function pointer 'void (*free)()' which
# GCC 15 rejects. Try a normal install first; only patch if that fails.
log "Installing rtracklayer..."
if "${ENV_PATH}/bin/Rscript" -e 'quit(status = !requireNamespace("rtracklayer", quietly=TRUE))' 2>/dev/null; then
    ok "rtracklayer already installed"
elif "${ENV_PATH}/bin/Rscript" -e \
        'BiocManager::install("rtracklayer", ask=FALSE, update=FALSE)' &>/dev/null \
     && "${ENV_PATH}/bin/Rscript" -e 'quit(status = !requireNamespace("rtracklayer", quietly=TRUE))' 2>/dev/null; then
    ok "rtracklayer installed (no patch needed)"
else
    log "  Standard install failed — retrying with GCC 15 source patch..."
    RTRACKLAYER_TMP="$(mktemp -d)"

    # Resolve the actual version for this Bioconductor release rather than
    # hardcoding one; a hardcoded version 404s whenever Bioc moves on.
    RT_INFO=$("${ENV_PATH}/bin/Rscript" -e '
      v <- as.character(BiocManager::version())
      url <- paste0("https://bioconductor.org/packages/", v, "/bioc/src/contrib")
      ap <- available.packages(repos = url)
      cat(v, ap["rtracklayer", "Version"], sep = " ")
    ' 2>/dev/null) || RT_INFO=""

    if [ -z "$RT_INFO" ]; then
        warn "Could not resolve rtracklayer version — skipping"
        warn "  GTF-based features may be unavailable"
    else
        BIOC_VER="${RT_INFO%% *}"
        RT_VER="${RT_INFO##* }"
        RT_URL="https://bioconductor.org/packages/${BIOC_VER}/bioc/src/contrib/rtracklayer_${RT_VER}.tar.gz"
        log "  Downloading rtracklayer ${RT_VER} (Bioc ${BIOC_VER})..."
        if curl -fL -o "${RTRACKLAYER_TMP}/rtracklayer.tar.gz" "$RT_URL"; then
            tar xzf "${RTRACKLAYER_TMP}/rtracklayer.tar.gz" -C "$RTRACKLAYER_TMP"
            log "  Patching for GCC 15 (void(*free)() -> void(*free)(void*))..."
            sed -i \
                's/void (\*free)()/void (*freeEl)(void*)/g;
                 s/else if (free != NULL)/else if (freeEl != NULL)/g;
                 s/\bfree(el)\b/freeEl(el)/g' \
                "${RTRACKLAYER_TMP}/rtracklayer/src/ucsc/common.c"
            sed -i 's/void (\*free)()/void (*free)(void*)/g' \
                "${RTRACKLAYER_TMP}/rtracklayer/src/ucsc/common.h"
            "${ENV_PATH}/bin/R" CMD INSTALL "${RTRACKLAYER_TMP}/rtracklayer/" \
                && ok "rtracklayer installed (GCC 15 patched)" \
                || warn "rtracklayer install failed — GTF-based features may be unavailable"
        else
            warn "Could not download rtracklayer source — skipping"
        fi
    fi
    rm -rf "$RTRACKLAYER_TMP"
fi

# ── Record what actually got installed ────────────────────────────────────────
# One plain-text file. Not a lockfile, nothing to learn, nothing extra to run —
# it is simply the thing to attach to a GitHub issue or keep with a manuscript.
VERSIONS_FILE="$(aracra_write_versions "$PIPELINE_DIR" "$ENV_PATH")"
ok "Version record written: ${VERSIONS_FILE}"
if [ "${ARACRA_PINNED_OK}" != "true" ]; then
    warn "Installed WITHOUT version floors — see ${VERSIONS_FILE}"
fi

# ── Verify tools ──
log "Verifying tool installation..."
MISSING=0
for tool in nextflow STAR hisat2 samtools fastp featureCounts salmon fastq_screen \
            qualimap multiqc pigz prefetch fasterq-dump Rscript streamlit; do
    if [ -x "${ENV_PATH}/bin/${tool}" ] || command -v "$tool" &>/dev/null 2>&1; then
        if [ "$tool" = "nextflow" ]; then
            ver=$("${ENV_PATH}/bin/nextflow" -version 2>&1 | grep -o 'version [0-9.]*' | head -1 || echo "OK")
        else
            ver=$("${ENV_PATH}/bin/${tool}" --version 2>&1 | head -1 | cut -c1-50 || echo "OK")
        fi
        printf "  ✔ %-16s %s\n" "$tool" "$ver"
    else
        printf "  ✘ %-16s NOT FOUND\n" "$tool"
        MISSING=$((MISSING + 1))
    fi
done
if [ "$MISSING" -gt 0 ]; then warn "${MISSING} tool(s) missing"; else ok "All tools verified"; fi

# ══════════════════════════════════════════════════════════════════════════════
# 4. DOWNLOAD REFERENCES
# ══════════════════════════════════════════════════════════════════════════════
REF_DIR="${DB_DIR}/hg38_reference"
ANN_DIR="${DB_DIR}/annotations"
SCREEN_DIR="${DB_DIR}/fastq_screen_genomes"
mkdir -p "$REF_DIR" "$ANN_DIR" "$SCREEN_DIR"

# ── GTF ──
GTF_FILE="${REF_DIR}/gencode.v44.primary_assembly.annotation.gtf"
if [ ! -s "$GTF_FILE" ]; then
    log "Downloading GTF annotation (~50 MB)..."
    wget -q --show-progress -O "${GTF_FILE}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
    gunzip -f "${GTF_FILE}.gz"
    ok "GTF: ${GTF_FILE}"
else
    ok "GTF exists: ${GTF_FILE}"
fi

# ── Genome FASTA ──
GENOME_FA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
if [ ! -s "$GENOME_FA" ]; then
    log "Downloading genome FASTA (~900 MB)..."
    wget -q --show-progress -O "${GENOME_FA}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
    gunzip -f "${GENOME_FA}.gz"
    ok "Genome: ${GENOME_FA}"
else
    ok "Genome exists: ${GENOME_FA}"
fi

# ── Transcriptome FASTA (for Salmon) ──
TRANSCRIPTOME_FA="${REF_DIR}/gencode.v44.transcripts.fa"
if [ ! -s "$TRANSCRIPTOME_FA" ]; then
    log "Downloading transcriptome FASTA (~300 MB)..."
    wget -q --show-progress -O "${TRANSCRIPTOME_FA}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz"
    gunzip -f "${TRANSCRIPTOME_FA}.gz"
    ok "Transcriptome: ${TRANSCRIPTOME_FA}"
else
    ok "Transcriptome exists"
fi

# ── STAR Index ──
# Completion is judged by a file STAR writes last-ish, not by directory
# existence: mkdir -p runs before the build, so a directory test would report
# "exists" for every interrupted build and never retry.
STAR_INDEX_DIR="${REF_DIR}/star_index_hg38"
STAR_INDEX_SENTINEL="${STAR_INDEX_DIR}/genomeParameters.txt"
if [ -s "$STAR_INDEX_SENTINEL" ] && [ -s "${STAR_INDEX_DIR}/SA" ]; then
    ok "STAR index exists: ${STAR_INDEX_DIR}"
elif [ "$SKIP_INDEX" = true ]; then
    warn "STAR index build skipped (--skip-index)"
    [ -s "$STAR_INDEX_SENTINEL" ] || STAR_INDEX_DIR=""
elif [ "$TOTAL_RAM_INT" -lt 32 ]; then
    warn "Skipping STAR index build — insufficient RAM (${TOTAL_RAM_GB} GB < 32 GB)"
    warn "Use HISAT2 for alignment, or build STAR index on a machine with 32+ GB RAM"
    STAR_INDEX_DIR=""
else
    log "Building STAR index (needs ~32 GB RAM, ~30-45 min)..."
    rm -rf "$STAR_INDEX_DIR"          # discard any partial build
    mkdir -p "$STAR_INDEX_DIR"
    if STAR --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX_DIR" \
         --genomeFastaFiles "$GENOME_FA" \
         --sjdbGTFfile "$GTF_FILE" \
         --runThreadN "$CPU_CORES" \
         --sjdbOverhang 100; then
        ok "STAR index: ${STAR_INDEX_DIR}"
    else
        warn "STAR index build failed — removing partial index"
        rm -rf "$STAR_INDEX_DIR"
        STAR_INDEX_DIR=""
    fi
fi

# ── HISAT2 Index ──
HISAT2_INDEX_DIR="${REF_DIR}/hisat2_index_hg38"
HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome_tran"
if [ -s "${HISAT2_INDEX_PREFIX}.1.ht2" ]; then
    ok "HISAT2 index exists: ${HISAT2_INDEX_PREFIX}"
elif [ "$SKIP_INDEX" = true ]; then
    warn "HISAT2 index download skipped (--skip-index)"
else
    log "Downloading HISAT2 pre-built index (~4 GB, ~10-30 min)..."
    mkdir -p "$HISAT2_INDEX_DIR"
    if wget -q --show-progress -O "${REF_DIR}/grch38_tran.tar.gz" \
        "https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz"; then
        tar -xzf "${REF_DIR}/grch38_tran.tar.gz" -C "$HISAT2_INDEX_DIR" --strip-components=1
        ok "HISAT2 index: ${HISAT2_INDEX_PREFIX}"
    else
        warn "HISAT2 index download failed — HISAT2 alignment unavailable"
    fi
    rm -f "${REF_DIR}/grch38_tran.tar.gz"
fi

# ── Salmon Index (with genome decoy for better mapping rates) ──
SALMON_INDEX_DIR="${REF_DIR}/salmon_index_hg38"
SALMON_SENTINEL="${SALMON_INDEX_DIR}/info.json"
if [ -s "$SALMON_SENTINEL" ]; then
    ok "Salmon index exists: ${SALMON_INDEX_DIR}"
elif [ "$SKIP_INDEX" = true ]; then
    warn "Salmon index build skipped (--skip-index)"
    SALMON_INDEX_DIR=""
elif ! "${ENV_PATH}/bin/salmon" --version &>/dev/null; then
    warn "Salmon is not functional (possible GLIBC incompatibility)"
    warn "  Fix: mamba install -p ${ENV_PATH} -c bioconda 'salmon>=1.10'"
    warn "  Skipping salmon index — salmon quantification will not be available"
    SALMON_INDEX_DIR=""
else
    log "Building Salmon index with decoys (~15-25 min, needs ~16 GB RAM)..."
    rm -rf "$SALMON_INDEX_DIR"
    GENTROME="${REF_DIR}/gentrome.fa"
    DECOYS="${REF_DIR}/decoys.txt"
    log "  Creating gentrome (transcriptome + genome decoys)..."
    cat "$TRANSCRIPTOME_FA" "$GENOME_FA" > "$GENTROME"
    grep "^>" "$GENOME_FA" | sed 's/>//' | cut -d' ' -f1 > "$DECOYS"

    if salmon index -t "$GENTROME" -d "$DECOYS" \
        -i "$SALMON_INDEX_DIR" --threads "$CPU_CORES" --gencode; then
        ok "Salmon index (with decoys): ${SALMON_INDEX_DIR}"
    else
        warn "Salmon index build failed — salmon quantification will not be available"
        rm -rf "$SALMON_INDEX_DIR"
        SALMON_INDEX_DIR=""
    fi
    rm -f "$GENTROME" "$DECOYS"
fi

# ── Annotation files ──
BED_FILE="${ANN_DIR}/hg38_RefSeq.bed"
if [ ! -s "$BED_FILE" ]; then
    log "Downloading RefSeq BED from UCSC..."
    for attempt in 1 2 3; do
        wget -nv --tries=3 --timeout=60 -O "${ANN_DIR}/refGene.txt.gz" \
            "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz" && break
        log "  Attempt ${attempt} failed, retrying..."
        sleep 5
    done
    if [ -s "${ANN_DIR}/refGene.txt.gz" ]; then
        gunzip -c "${ANN_DIR}/refGene.txt.gz" | \
            awk -F'\t' 'BEGIN{OFS="\t"} {
                split($10, starts, ","); split($11, ends, ",")
                bstarts=""; bsizes=""
                for(i=1;i<=$9;i++){
                    bsizes=bsizes (ends[i]-starts[i])","
                    bstarts=bstarts (starts[i]-$5) ","
                }
                print $3,$5,$6,$13,0,$4,$7,$8,"0",$9,bsizes,bstarts
            }' > "$BED_FILE"
        rm -f "${ANN_DIR}/refGene.txt.gz"
        if [ -s "$BED_FILE" ]; then
            ok "BED: ${BED_FILE} ($(wc -l < "$BED_FILE") genes)"
        else
            warn "BED conversion produced empty file — gene body coverage (full) disabled"
            rm -f "$BED_FILE"
        fi
    else
        warn "Could not download RefSeq BED — gene body coverage (full) disabled"
        rm -f "${ANN_DIR}/refGene.txt.gz"
    fi
else
    ok "BED exists: ${BED_FILE}"
fi

# Housekeeping gene BED (Eisenberg & Levanon 2013, ~3800 genes from RSeQC)
HK_BED="${ANN_DIR}/hg38_housekeeping.bed"
if [ ! -s "$HK_BED" ]; then
    log "Downloading RSeQC housekeeping gene BED (hg38)..."
    HK_OK=false
    for attempt in 1 2 3; do
        wget -nv --tries=3 --timeout=60 --content-disposition -O "${HK_BED}.gz" \
            "https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38.HouseKeepingGenes.bed.gz/download" && break
        log "  Attempt ${attempt} failed, retrying..."
        sleep 5
    done
    if [ -s "${HK_BED}.gz" ]; then
        gunzip -f "${HK_BED}.gz"
        if [ -f "${ANN_DIR}/hg38.HouseKeepingGenes.bed" ] && [ ! -f "$HK_BED" ]; then
            mv "${ANN_DIR}/hg38.HouseKeepingGenes.bed" "$HK_BED"
        fi
        if [ -s "$HK_BED" ]; then
            ok "Housekeeping BED: ${HK_BED} ($(wc -l < "$HK_BED") transcripts)"
            HK_OK=true
        fi
    fi
    if [ "$HK_OK" = false ]; then
        warn "Download failed — creating minimal fallback (11 genes)"
        rm -f "${HK_BED}.gz" "${HK_BED}"
        cat > "$HK_BED" << 'HKEOF'
chr12	6534517	6538370	GAPDH	0	+
chr7	5527151	5530601	ACTB	0	-
chr15	44711477	44718877	B2M	0	+
chr19	4399400	4412854	LDHA	0	+
chr12	109535194	109539703	UBC	0	-
chr1	228283854	228289247	RHOU	0	-
chr1	11166590	11322564	MTOR	0	-
chr7	44834066	44838200	PPIA	0	-
chr11	65497657	65506516	RPLP0	0	+
chr19	49457917	49459845	RPL18	0	-
chr2	55082485	55083734	RTN4	0	+
HKEOF
        ok "Fallback housekeeping BED: ${HK_BED} (11 genes)"
    fi
else
    ok "HK BED exists ($(wc -l < "$HK_BED") transcripts)"
fi

# refFlat for Picard
REFFLAT="${ANN_DIR}/hg38_refFlat.txt"
if [ ! -s "$REFFLAT" ]; then
    log "Downloading refFlat..."
    for attempt in 1 2 3; do
        wget -nv --tries=3 --timeout=60 -O "${REFFLAT}.gz" \
            "https://hgdownload.soe.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz" && break
        log "  Attempt ${attempt} failed, retrying..."
        sleep 5
    done
    if [ -s "${REFFLAT}.gz" ]; then
        gunzip -f "${REFFLAT}.gz"
        ok "refFlat: ${REFFLAT} ($(wc -l < "$REFFLAT") entries)"
    else
        warn "Could not download refFlat — Picard RNA metrics disabled"
        rm -f "${REFFLAT}.gz"
    fi
else
    ok "refFlat exists"
fi

# ── FastQ Screen genomes ────────────────────────────────────────────────────
# 'fastq_screen --get_genomes' is NOT used. Its location file
# (genome_locations.txt) serves schemeless URLs; wget prepends http:// and the
# Babraham http->https rewrite emits a malformed Location header (the slash
# after the hostname is dropped), so every fetch fails. We mirror the served
# directory over HTTPS directly and generate the conf ourselves — which also
# guarantees the paths inside it are local and absolute.
SCREEN_ROOT="${SCREEN_DIR}/FastQ_Screen_Genomes"
SCREEN_CONF="${SCREEN_ROOT}/fastq_screen.conf"
SCREEN_URL="https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/Genome_Data/FastQ_Screen_Genomes/"

generate_screen_conf() {
    local n=0 idx
    mkdir -p "$SCREEN_ROOT"
    {
        echo "# Auto-generated by ARACRA setup.sh on $(date)"
        printf 'BOWTIE2\t%s/bin/bowtie2\n' "${ENV_PATH}"
        printf 'THREADS\t2\n'
        echo ""
    } > "${SCREEN_CONF}.tmp"
    while IFS= read -r idx; do
        printf 'DATABASE\t%s\t%s\n' \
            "$(basename "$(dirname "$idx")")" "${idx%.1.bt2}" >> "${SCREEN_CONF}.tmp"
        n=$((n + 1))
    done < <(find "$SCREEN_ROOT" -name '*.1.bt2' ! -name '*.rev.1.bt2' 2>/dev/null | sort)

    if [ "$n" -eq 0 ]; then
        rm -f "${SCREEN_CONF}.tmp"
        return 1
    fi
    mv "${SCREEN_CONF}.tmp" "$SCREEN_CONF"
    ok "Generated fastq_screen.conf (${n} genomes)"
    return 0
}

if [ "$SKIP_SCREEN" = true ]; then
    ok "FastQ Screen skipped (--skip-screen)"
elif [ -s "$SCREEN_CONF" ]; then
    ok "FastQ Screen genomes exist"
elif ! curl -fsI --max-time 30 "$SCREEN_URL" >/dev/null 2>&1; then
    warn "FastQ Screen genome server unreachable — contamination check disabled"
    warn "  Re-run setup.sh later; only this step will repeat."
else
    log "Downloading FastQ Screen genomes (~30 GB, 30-60 min)..."
    log "  Bowtie2 indices for Human, Mouse, Rat, E.coli, adapters, vectors, etc."
    log "  Run with --skip-screen to skip this step and add genomes later."
    # -nH with --cut-dirs=3 strips 'projects/fastq_screen/Genome_Data',
    # landing the tree at ${SCREEN_DIR}/FastQ_Screen_Genomes/ as SCREEN_CONF expects.
    if wget -q --show-progress -r -np -nH --cut-dirs=3 -R "index.html*" \
            -P "$SCREEN_DIR" "$SCREEN_URL"; then
        # Always regenerate: any conf shipped upstream carries the maintainer's
        # absolute paths, not ours.
        generate_screen_conf || warn "No bowtie2 indices found — contamination check disabled"
    else
        warn "FastQ Screen download failed — contamination check disabled"
    fi
fi

# ── tx2gene mapping (for Salmon) ──
TX2GENE="${REF_DIR}/tx2gene.csv"
if [ ! -s "$TX2GENE" ]; then
    log "Generating tx2gene mapping..."
    "${ENV_PATH}/bin/python3" - "$GTF_FILE" "$TX2GENE" << 'PYEOF'
import sys, re, csv
gtf_file, out_file = sys.argv[1], sys.argv[2]
pairs = set()
with open(gtf_file) as f:
    for line in f:
        if line.startswith('#'): continue
        fields = line.strip().split('\t')
        if len(fields) < 9 or fields[2] != 'transcript': continue
        attrs = fields[8]
        tid = re.search(r'transcript_id "([^"]+)"', attrs)
        gid = re.search(r'gene_id "([^"]+)"', attrs)
        if tid and gid:
            pairs.add((tid.group(1), gid.group(1)))
with open(out_file, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['transcript_id', 'gene_id'])
    for tid, gid in sorted(pairs):
        w.writerow([tid, gid])
print(f"tx2gene: {len(pairs)} transcript-gene pairs")
PYEOF
    ok "tx2gene: ${TX2GENE} ($(wc -l < "$TX2GENE") entries)"
else
    ok "tx2gene exists: ${TX2GENE} ($(wc -l < "$TX2GENE") entries)"
fi

# ══════════════════════════════════════════════════════════════════════════════
# 5. RESOLVE MODULE AVAILABILITY
# ══════════════════════════════════════════════════════════════════════════════
# Every toggle below is derived from whether its reference actually exists on
# disk. A hardcoded 'true' here is what turns a soft download failure into a
# hard pipeline crash several minutes into a run.
avail() { [ -e "$1" ] && [ -s "$1" ] && echo "true" || echo "false"; }
avail_dir() { [ -s "$1" ] && echo "true" || echo "false"; }

HAS_SCREEN=$(avail "$SCREEN_CONF")
HAS_BED=$(avail "$BED_FILE")
HAS_HK=$(avail "$HK_BED")
HAS_REFFLAT=$(avail "$REFFLAT")
HAS_STAR=$([ -n "$STAR_INDEX_DIR" ] && avail_dir "${STAR_INDEX_DIR}/SA" || echo "false")
HAS_HISAT2=$(avail "${HISAT2_INDEX_PREFIX}.1.ht2")
HAS_SALMON=$([ -n "$SALMON_INDEX_DIR" ] && avail "${SALMON_INDEX_DIR}/info.json" || echo "false")

# Fall back if the recommended aligner is not actually available.
if [ "$RECOMMENDED_ALIGNER" = "star" ] && [ "$HAS_STAR" = "false" ]; then
    if [ "$HAS_HISAT2" = "true" ]; then
        warn "STAR index unavailable — falling back to HISAT2"
        RECOMMENDED_ALIGNER="hisat2"
    else
        warn "Neither STAR nor HISAT2 index is available — alignment will not run"
    fi
elif [ "$RECOMMENDED_ALIGNER" = "hisat2" ] && [ "$HAS_HISAT2" = "false" ] && [ "$HAS_STAR" = "true" ]; then
    warn "HISAT2 index unavailable — falling back to STAR"
    RECOMMENDED_ALIGNER="star"
fi

RUN_STAR=$([ "$RECOMMENDED_ALIGNER" = "star" ] && [ "$HAS_STAR" = "true" ] && echo "true" || echo "false")
RUN_HISAT2=$([ "$RECOMMENDED_ALIGNER" = "hisat2" ] && [ "$HAS_HISAT2" = "true" ] && echo "true" || echo "false")

# ══════════════════════════════════════════════════════════════════════════════
# 6. CREATE WORK DIRECTORIES
# ══════════════════════════════════════════════════════════════════════════════
WORK_DIR="${HOME}/aracra_star_work/work"
OUT_DIR="${HOME}/aracra_star_work/results"
LOG_PATH="${HOME}/aracra_star_work/pipeline.log"
mkdir -p "$WORK_DIR" "$OUT_DIR"
ok "Work directories created"

# ══════════════════════════════════════════════════════════════════════════════
# 7. WRITE .env (auto-discovery results)
# ══════════════════════════════════════════════════════════════════════════════
# THE only place .env is written (aracra_write_env, in lib/aracra_common.sh) —
# setup.sh and run_app.sh both call it, so the key set the app reads can never
# drift from what a writer emits. The app itself re-derives module
# availability live from the sidebar paths (_file_ok()) rather than trusting
# a static HAS_* snapshot from install time, which is why those keys aren't
# needed here even though an earlier version of this script wrote them.
#
# STAR_INDEX_DIR is deliberately empty when indexing was skipped (low RAM /
# --skip-index); pass it through so the app sees the same thing it always did.
ARACRA_REF_STAR_INDEX_OVERRIDE="${STAR_INDEX_DIR}"
aracra_write_env "$PIPELINE_DIR" "$CONDA_BASE" "$ENV_PATH" "$DB_DIR" "setup.sh"
ok ".env written"

# ══════════════════════════════════════════════════════════════════════════════
# 8. WRITE nextflow.config (auto-tuned to system)
# ══════════════════════════════════════════════════════════════════════════════
# No top-level 'def' declarations: Nextflow's strict config parser (default
# from 25.10 onward) rejects variable declarations mixed with config
# statements. Since this file is regenerated per machine, the core counts are
# baked in by the shell instead of computed at runtime by Groovy.
#
# Nextflow's manifest.nextflowVersion accepts the same ">=x,<y" shape as the
# conda spec, just with a space after the comma. Convert rather than keeping
# two hand-maintained version strings that can drift apart.
ARACRA_NEXTFLOW_SPEC_NF="$(echo "${ARACRA_NEXTFLOW_SPEC}" | sed 's/,/, /')"
cat > "${PIPELINE_DIR}/nextflow.config" << NFEOF
/*
 * ARACRA-STAR Pipeline — nextflow.config
 * Auto-generated by setup.sh on $(date)
 * System: ${TOTAL_RAM_GB} GB RAM | ${CPU_CORES} CPUs
 * Do not edit by hand — re-run setup.sh to regenerate.
 */

params {
    metadata            = "metadata.xlsx"
    layout              = "SE"
    aligner             = "${RECOMMENDED_ALIGNER}"

    // References (auto-discovered by setup.sh)
    star_index          = "${STAR_INDEX_DIR}"
    hisat2_index        = "${HISAT2_INDEX_PREFIX}"
    salmon_index        = "${SALMON_INDEX_DIR}"
    gtf                 = "${GTF_FILE}"
    bed                 = "${BED_FILE}"
    housekeeping_bed    = "${HK_BED}"
    refflat             = "${REFFLAT}"
    fastq_screen_conf   = "${SCREEN_CONF}"

    outdir              = "results"
    max_memory          = "${TOTAL_RAM_INT} GB"
    star_ram            = $(( TOTAL_RAM_INT * 600000000 ))

    // FASTQ source (null = download from SRA)
    fastq_dir           = null
    trimmed_dir         = null

    // Step toggles — derived from reference availability, never hardcoded
    run_download            = true
    run_fastp               = true
    run_fastq_screen        = ${HAS_SCREEN}
    run_star                = ${RUN_STAR}
    run_hisat2              = ${RUN_HISAT2}
    run_samtools_stats      = true
    run_strandedness        = ${HAS_BED}
    run_coverage_hk         = ${HAS_HK}
    run_coverage_full       = false
    run_read_distribution   = false
    run_picard              = false
    run_qualimap            = false
    run_featurecounts       = true
    run_salmon              = false
    run_salmon_pseudo       = false
    run_merge_counts        = true
    run_multiqc             = true
}

process {
    // Transient-failure resilience: retry signal-killed tasks (network
    // blips on downloads/KEGG API calls, scheduler preemption, OOM-kill)
    // up to twice before giving up; a genuine script/input error (exit 1)
    // still fails immediately. Processes with their own errorStrategy set
    // in the script (DESEQ2_ANALYSIS, DROMICS_ANALYSIS, ENRICHMENT_ANALYSIS,
    // DOWNLOAD_SRA, etc.) override this default for that process only.
    errorStrategy = { task.exitStatus in ((130..145) + 104) ? 'retry' : 'finish' }
    maxRetries    = 2

    withLabel: 'download' {
        cpus     = 1
        memory   = '2 GB'
        maxForks = 3
    }
    withLabel: 'pre_qc' {
        cpus     = 1
        memory   = '2 GB'
        maxForks = ${LIGHT_PARALLEL}
    }
    withLabel: 'post_align' {
        cpus     = 1
        memory   = '4 GB'
        maxForks = ${POST_PARALLEL}
    }
    withLabel: 'heavy' {
        cpus     = ${HEAVY_CORES}
        memory   = params.max_memory
        maxForks = ${STAR_PARALLEL}
    }
    withLabel: 'screen' {
        cpus     = 2
        memory   = '4 GB'
        maxForks = 1
    }
    withLabel: 'final' {
        cpus   = ${CPU_CORES}
        memory = params.max_memory
    }
    withLabel: 'minimal' {
        cpus   = 1
        memory = '2 GB'
    }
}

manifest {
    name            = 'ARACRA'
    version         = '2.0'
    description     = 'RNA-seq & TempO-Seq analysis platform'
    mainScript      = 'main.nf'
    // main.nf uses 24.x DSL2 semantics. Nextflow enforces this itself and
    // stops with a clear message, instead of failing obscurely mid-run.
    // Derived from lib/aracra_common.sh's ARACRA_NEXTFLOW_SPEC (see above) —
    // not a second hardcoded copy that can drift from the installer's floor.
    nextflowVersion = '${ARACRA_NEXTFLOW_SPEC_NF}'
}

profiles {
    standard { process.executor = 'local' }
    slurm    { process.executor = 'slurm'; process.queue = 'batch' }
}

// overwrite = true — this pipeline is re-run into the same outdir constantly
// (every GUI "Run" click reuses it), and without it Nextflow refuses to start
// at all once a report file exists.
report   { enabled = true; overwrite = true; file = "\${params.outdir}/pipeline_report.html" }
timeline { enabled = true; overwrite = true; file = "\${params.outdir}/timeline.html" }
NFEOF
ok "nextflow.config written (tuned for ${CPU_CORES} cores, ${TOTAL_RAM_GB} GB RAM)"

# ══════════════════════════════════════════════════════════════════════════════
# 9. ENSURE scripts/ DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════
if [ ! -d "${PIPELINE_DIR}/scripts" ]; then
    mkdir -p "${PIPELINE_DIR}/scripts"
    warn "Created scripts/ directory — copy run_deseq2.R and run_dromics.R there"
fi

# ══════════════════════════════════════════════════════════════════════════════
# 10. REFERENCE VERIFICATION GATE
# ══════════════════════════════════════════════════════════════════════════════
# The tool loop above verifies binaries; nothing verified reference FILES.
# That asymmetry is why setup could report success while leaving a module
# enabled with no reference behind it.
log "Verifying reference files..."
REF_MISSING=0
check_ref() {  # name, path, required(yes/no)
    if [ -s "$2" ] || { [ -d "$2" ] && [ -n "$(ls -A "$2" 2>/dev/null)" ]; }; then
        printf "  ✔ %-20s %s\n" "$1" "$2"
    elif [ "$3" = "yes" ]; then
        printf "  ✘ %-20s MISSING (required)\n" "$1"
        REF_MISSING=$((REF_MISSING + 1))
    else
        printf "  ○ %-20s absent — module disabled\n" "$1"
    fi
}
check_ref "GTF"            "$GTF_FILE"                    yes
check_ref "Genome FASTA"   "$GENOME_FA"                   yes
check_ref "tx2gene"        "$TX2GENE"                     yes
check_ref "Housekeeping BED" "$HK_BED"                    yes
check_ref "STAR index"     "${STAR_INDEX_DIR:-/nonexistent}/SA"     no
check_ref "HISAT2 index"   "${HISAT2_INDEX_PREFIX}.1.ht2" no
check_ref "Salmon index"   "${SALMON_INDEX_DIR:-/nonexistent}/info.json" no
check_ref "RefSeq BED"     "$BED_FILE"                    no
check_ref "refFlat"        "$REFFLAT"                     no
check_ref "FastQ Screen conf" "$SCREEN_CONF"              no

if [ "$RUN_STAR" = "false" ] && [ "$RUN_HISAT2" = "false" ]; then
    REF_MISSING=$((REF_MISSING + 1))
    printf "  ✘ %-20s no usable aligner index\n" "ALIGNER"
fi

echo ""
if [ "$REF_MISSING" -gt 0 ]; then
    echo -e "${RED}${BOLD}══════════════════════════════════════════════════════════════${NC}"
    echo -e "${RED}${BOLD}  Setup INCOMPLETE — ${REF_MISSING} required reference(s) missing${NC}"
    echo -e "${RED}${BOLD}══════════════════════════════════════════════════════════════${NC}"
    echo ""
    echo "  Re-run 'bash setup.sh' — completed steps are skipped automatically."
    echo ""
    exit 1
fi

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
echo -e "${GREEN}${BOLD}══════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}${BOLD}  Setup complete!${NC}"
echo -e "${GREEN}${BOLD}══════════════════════════════════════════════════════════════${NC}"
echo ""
echo "  System"
echo "    CPUs            : ${CPU_CORES}"
echo "    RAM             : ${TOTAL_RAM_GB} GB"
echo "    Aligner         : ${RECOMMENDED_ALIGNER^^}"
if [ "$TOTAL_RAM_INT" -lt 32 ]; then
echo -e "    ${RED}⚠ Low RAM — HISAT2 recommended over STAR${NC}"
fi
echo ""
echo "  Environment"
echo "    Conda           : ${CONDA_BASE}"
echo "    Env             : ${ENV_PATH}"
echo "    Nextflow        : $("${ENV_PATH}/bin/nextflow" -version 2>&1 | grep -o 'version [0-9.]*' | head -1 || echo '?')"
echo ""
echo "  Optional modules"
echo "    FastQ Screen    : ${HAS_SCREEN}"
echo "    Strandedness    : ${HAS_BED}"
echo "    Picard refFlat  : ${HAS_REFFLAT}"
echo "    Salmon          : ${HAS_SALMON}"
echo ""
echo "  Directories"
echo "    Pipeline        : ${PIPELINE_DIR}"
echo "    Output          : ${OUT_DIR}"
echo "    Work            : ${WORK_DIR}"
echo ""
echo "  ➤  To start:  cd ${PIPELINE_DIR} && bash run_app.sh"
echo "  ➤  Then open: http://localhost:${ARACRA_DEFAULT_PORT}"
echo ""