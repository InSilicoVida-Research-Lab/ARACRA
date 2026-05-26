#!/usr/bin/env bash
# =============================================================================
#  ARACRA Pipeline — setup.sh 
#  Auto-discovers system resources, builds conda env, downloads references,
#  configures everything. Warns if RAM < 32GB (recommends HISAT2 over STAR).
# =============================================================================
set -euo pipefail

CYAN='\033[0;36m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
RED='\033[0;31m'; NC='\033[0m'; BOLD='\033[1m'
log()  { echo -e "${CYAN}[SETUP]${NC} $*"; }
ok()   { echo -e "${GREEN}[  OK ]${NC} $*"; }
warn() { echo -e "${YELLOW}[ WARN]${NC} $*"; }
fail() { echo -e "${RED}[ERROR]${NC} $*"; exit 1; }

# ── Error trap ────────────────────────────────────────────────────────────────
# With `set -e`, any unhandled failure aborts the script. This trap makes the
# abort informative: it reports the line number and the command that failed.
# Steps that are allowed to fail must handle it explicitly (|| warn ...).
on_error() {
    local exit_code=$?
    local line_no=$1
    echo -e "${RED}[ERROR]${NC} Setup failed at line ${line_no} (exit ${exit_code})." >&2
    echo -e "${RED}[ERROR]${NC} Failing command: ${BASH_COMMAND}" >&2
    echo -e "${RED}[ERROR]${NC} Fix the issue above and re-run — completed steps will be skipped." >&2
    exit "$exit_code"
}
trap 'on_error ${LINENO}' ERR

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Download helper: fetch a URL, then verify the result is real ──────────────
# Guards against the silent-corruption failure mode: wget can "succeed" (exit 0)
# while returning a truncated file or an HTML error page. After download we
# verify the file exists and exceeds a minimum plausible size (bytes).
# Usage: fetch_verified <url> <output_path> <min_bytes> [description]
fetch_verified() {
    local url="$1" out="$2" min_bytes="$3" desc="${4:-file}"
    local tmp="${out}.partial"
    rm -f "$tmp"
    if ! wget -q --show-progress --tries=3 --timeout=120 -O "$tmp" "$url"; then
        rm -f "$tmp"
        fail "Download failed (${desc}): ${url}"
    fi
    local size
    size=$(stat -c%s "$tmp" 2>/dev/null || echo 0)
    if [ "$size" -lt "$min_bytes" ]; then
        rm -f "$tmp"
        fail "Downloaded ${desc} is only ${size} bytes (expected >= ${min_bytes}) — likely an error page, not real data: ${url}"
    fi
    mv "$tmp" "$out"
}

# ── Gzip integrity check ──────────────────────────────────────────────────────
# Verifies a .gz file is a valid, complete gzip archive before it is used.
verify_gzip() {
    local f="$1" desc="${2:-archive}"
    gzip -t "$f" 2>/dev/null || fail "Corrupt or incomplete gzip (${desc}): ${f}"
}

# ── Parse flags ───────────────────────────────────────────────────────────────
SKIP_INDEX=false
# FastQ Screen genome builds are opt-in: ~2.5 GB of FASTA downloads plus
# multi-hour bowtie2 index builds for human/mouse/rat. A plain `bash setup.sh`
# skips it; `bash setup.sh --with-screen` runs it as part of the same command.
WITH_SCREEN=false
DB_DIR="${HOME}/databases"
for arg in "$@"; do
    case "$arg" in
        --skip-index)    SKIP_INDEX=true ;;
        --with-screen)   WITH_SCREEN=true ;;
        --skip-screen)   WITH_SCREEN=false ;;   # alias of the default; kept for compatibility
        --db-dir=*)      DB_DIR="${arg#*=}" ;;
        --help|-h)
            echo "Usage: bash setup.sh [--skip-index] [--with-screen] [--db-dir=PATH]"
            echo ""
            echo "  --skip-index    Skip building STAR/HISAT2/Salmon indexes"
            echo "  --with-screen   Also download FASTA and build FastQ Screen bowtie2"
            echo "                  indexes (Human, Mouse, Rat, E.coli, PhiX, Adapters,"
            echo "                  Vectors). Adds several hours — runs in the same command."
            echo "  --db-dir=PATH   Database directory (default: \${HOME}/databases)"
            exit 0 ;;
    esac
done

echo -e "${BOLD}"
echo "╔══════════════════════════════════════════════════════════╗"
echo "║     ARACRA Pipeline — Setup                 ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo -e "${NC}"

# ══════════════════════════════════════════════════════════════════════════════
# 1. SYSTEM DISCOVERY
# ══════════════════════════════════════════════════════════════════════════════
CPU_CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
TOTAL_RAM_KB=$(awk '/MemTotal/{print $2}' /proc/meminfo 2>/dev/null || echo 0)
TOTAL_RAM_GB=$(awk "BEGIN{printf \"%.1f\", ${TOTAL_RAM_KB}/1048576}")
TOTAL_RAM_INT=$(awk "BEGIN{printf \"%d\", ${TOTAL_RAM_KB}/1048576}")
DISK_FREE_GB=$(df -BG "${HOME}" 2>/dev/null | awk 'NR==2{gsub("G",""); print $4}' || echo "?")

log "System discovery:"
log "  CPUs       : ${CPU_CORES}"
log "  RAM        : ${TOTAL_RAM_GB} GB"
log "  Disk free  : ${DISK_FREE_GB} GB"
log "  DB dir     : ${DB_DIR}"

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
    RECOMMENDED_ALIGNER="hisat2"
else
    RECOMMENDED_ALIGNER="star"
fi

# ── Core allocation strategy ──
HALF_CORES=$(( CPU_CORES / 2 > 0 ? CPU_CORES / 2 : 1 ))
LIGHT_PARALLEL=$(( CPU_CORES - 1 > 0 ? CPU_CORES - 1 : 1 ))
POST_PARALLEL=$(( CPU_CORES - 2 > 0 ? CPU_CORES - 2 : 1 ))
STAR_CORES=$HALF_CORES
STAR_PARALLEL=2
# If very few cores, run STAR single
if [ "$CPU_CORES" -le 4 ]; then
    STAR_CORES=$CPU_CORES
    STAR_PARALLEL=1
fi

log "Core strategy:"
log "  Download     : 1 core × 3 parallel"
log "  Fastp        : 1 core × ${LIGHT_PARALLEL} parallel"
log "  STAR/HISAT2  : ${STAR_CORES} cores × ${STAR_PARALLEL} parallel"
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

# Ensure mamba
if ! command -v mamba &>/dev/null; then
    log "Installing mamba..."
    conda install -n base -c conda-forge mamba -y -q
fi

# ══════════════════════════════════════════════════════════════════════════════
# 3. CREATE ENVIRONMENT
# ══════════════════════════════════════════════════════════════════════════════
ENV_NAME="test_ARACRA"
ENV_PATH="${CONDA_BASE}/envs/${ENV_NAME}"

if [ -d "$ENV_PATH" ]; then
    log "Environment '${ENV_NAME}' already exists — updating..."
    INSTALL_CMD="mamba install -p ${ENV_PATH}"
else
    log "Creating environment '${ENV_NAME}'..."
    mamba create -p "$ENV_PATH" python=3.12 -y -q
    INSTALL_CMD="mamba install -p ${ENV_PATH}"
fi

conda activate "$ENV_PATH"
ok "Active env: ${ENV_PATH}"

# ── Install bioinformatics tools ──
log "Installing pipeline tools (this may take a few minutes)..."
mamba install -p "$ENV_PATH" -c conda-forge -c bioconda -y -q \
    openjdk=17 \
    'nextflow>=24.04,<25' \
    python=3.12 \
    openpyxl \
    sra-tools \
    fastp \
    fastq-screen \
    fastqc \
    star \
    hisat2 \
    bowtie2 \
    samtools \
    rseqc \
    picard \
    qualimap \
    salmon>=1.10 \
    subread \
    multiqc \
    pigz \
    || fail "Tool installation failed"
ok "Bioinformatics tools installed"

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

# DRomics + Bioconductor packages from R
log "Installing DRomics and remaining Bioconductor packages..."
"${ENV_PATH}/bin/Rscript" -e '
options(repos=c(CRAN="https://cloud.r-project.org"))
if (!requireNamespace("DRomics", quietly=TRUE)) install.packages("DRomics", quiet=TRUE)
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", quiet=TRUE)
miss <- c("DESeq2","edgeR","sva","org.Hs.eg.db","GO.db","KEGGREST","AnnotationDbi")[
    !sapply(c("DESeq2","edgeR","sva","org.Hs.eg.db","GO.db","KEGGREST","AnnotationDbi"), requireNamespace, quietly=TRUE)]
if (length(miss)) BiocManager::install(miss, ask=FALSE, update=FALSE)
cat("R packages OK\n")
' 2>/dev/null || warn "DRomics install had issues — check manually"
ok "R packages done"

# ── Install rtracklayer ──
# Strategy: prefer the prebuilt bioconda package (works on linux-64 and
# linux-aarch64 — the binary is already compiled). If that is unavailable,
# build from source.
#
# GCC 15 defaults to the C23 standard, where an empty parameter list '()'
# means "zero arguments" instead of "unspecified". The bundled UCSC source
# in rtracklayer declares function pointers like 'void (*free)()', which
# C23 then rejects. Rather than text-patching the source (fragile — the
# exact spelling changes between package versions), we compile this one
# package under the older gnu17 standard via a temporary ~/.R/Makevars.
# This requires no source edits at all and is version-independent.
log "Installing rtracklayer..."
if "${ENV_PATH}/bin/Rscript" -e 'requireNamespace("rtracklayer", quietly=TRUE)' 2>/dev/null | grep -q TRUE; then
    ok "rtracklayer already installed"
else
    # --- Attempt 1: prebuilt bioconda binary (no compilation needed) ---
    log "  Trying prebuilt bioconda package..."
    if mamba install -y -q -p "$ENV_PATH" -c bioconda -c conda-forge \
        bioconductor-rtracklayer >/dev/null 2>&1 \
        && "${ENV_PATH}/bin/Rscript" -e 'requireNamespace("rtracklayer", quietly=TRUE)' \
            2>/dev/null | grep -q TRUE; then
        ok "rtracklayer installed (prebuilt bioconda binary)"
    else
        # --- Attempt 2: source build under the gnu17 C standard ---
        log "  Prebuilt package unavailable — building from source (gnu17 standard)..."

        # Temporarily inject -std=gnu17 via ~/.R/Makevars so R CMD INSTALL
        # compiles all C sources for this package under pre-C23 rules.
        # Back up any existing Makevars and restore it afterwards.
        MAKEVARS_DIR="${HOME}/.R"
        MAKEVARS_FILE="${MAKEVARS_DIR}/Makevars"
        MAKEVARS_BAK=""
        mkdir -p "$MAKEVARS_DIR"
        if [ -f "$MAKEVARS_FILE" ]; then
            MAKEVARS_BAK="${MAKEVARS_FILE}.aracra_bak_$$"
            cp "$MAKEVARS_FILE" "$MAKEVARS_BAK"
        fi
        # CFLAGS covers .c sources; the *_STD vars are honoured by newer R toolchains.
        cat >> "$MAKEVARS_FILE" <<'MAKEVARS_EOF'
CFLAGS += -std=gnu17
C_STD = gnu17
CXXSTD = gnu++17
MAKEVARS_EOF

        # restore_makevars: undo our temporary Makevars change on any exit path
        restore_makevars() {
            if [ -n "$MAKEVARS_BAK" ]; then
                mv "$MAKEVARS_BAK" "$MAKEVARS_FILE"
            else
                rm -f "$MAKEVARS_FILE"
            fi
        }

        # Run the installer. BiocManager resolves the correct version for the
        # installed Bioconductor release itself. Note: we deliberately do NOT
        # use the installer's exit code to decide success — post-install steps
        # (vignette building, HTML index updates) can yield a nonzero status
        # even when the package itself installed fine ("* DONE"). The only
        # reliable test is whether the package can actually be loaded.
        log "  Compiling and installing rtracklayer (this may take a few minutes)..."
        "${ENV_PATH}/bin/Rscript" -e \
            'BiocManager::install("rtracklayer", update=FALSE, ask=FALSE)' || true

        # Restore the user's Makevars BEFORE verifying, so the check runs in a
        # clean environment unaffected by our temporary compiler flags.
        restore_makevars

        # Sole arbiter of success: can rtracklayer be loaded?
        if "${ENV_PATH}/bin/Rscript" -e \
                'if (requireNamespace("rtracklayer", quietly=TRUE)) cat("RTOK")' \
                2>/dev/null | grep -q RTOK; then
            ok "rtracklayer installed (built from source, gnu17 standard)"
        else
            fail "rtracklayer installation failed — see output above"
        fi
    fi
fi

# ── Verify tools ──
log "Verifying tool installation..."
MISSING=0
for tool in nextflow STAR hisat2 samtools fastp featureCounts salmon fastq_screen \
            qualimap multiqc pigz prefetch fasterq-dump Rscript streamlit; do
    if command -v "$tool" &>/dev/null 2>&1 || [ -x "${ENV_PATH}/bin/${tool}" ]; then
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
[ "$MISSING" -gt 0 ] && warn "${MISSING} tool(s) missing" || ok "All tools verified"

# ══════════════════════════════════════════════════════════════════════════════
# 4. DOWNLOAD REFERENCES
# ══════════════════════════════════════════════════════════════════════════════
REF_DIR="${DB_DIR}/hg38_reference"
ANN_DIR="${DB_DIR}/annotations"
SCREEN_DIR="${DB_DIR}/fastq_screen_genomes"
mkdir -p "$REF_DIR" "$ANN_DIR" "$SCREEN_DIR"

# ── GTF ──
GTF_FILE="${REF_DIR}/gencode.v44.primary_assembly.annotation.gtf"
if [ ! -f "$GTF_FILE" ] || [ ! -s "$GTF_FILE" ]; then
    log "Downloading GTF annotation (~50 MB)..."
    fetch_verified \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" \
        "${GTF_FILE}.gz" 10000000 "GTF annotation"
    verify_gzip "${GTF_FILE}.gz" "GTF annotation"
    gunzip -f "${GTF_FILE}.gz"
    [ -s "$GTF_FILE" ] || fail "GTF is empty after extraction: ${GTF_FILE}"
    ok "GTF: ${GTF_FILE}"
else
    ok "GTF exists: ${GTF_FILE}"
fi

# ── Genome FASTA ──
GENOME_FA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
if [ ! -f "$GENOME_FA" ] || [ ! -s "$GENOME_FA" ]; then
    log "Downloading genome FASTA (~900 MB)..."
    fetch_verified \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz" \
        "${GENOME_FA}.gz" 500000000 "genome FASTA"
    verify_gzip "${GENOME_FA}.gz" "genome FASTA"
    gunzip -f "${GENOME_FA}.gz"
    [ -s "$GENOME_FA" ] || fail "Genome FASTA is empty after extraction: ${GENOME_FA}"
    ok "Genome: ${GENOME_FA}"
else
    ok "Genome exists: ${GENOME_FA}"
fi

# ── Transcriptome FASTA (for Salmon) ──
TRANSCRIPTOME_FA="${REF_DIR}/gencode.v44.transcripts.fa"
if [ ! -f "$TRANSCRIPTOME_FA" ] || [ ! -s "$TRANSCRIPTOME_FA" ]; then
    log "Downloading transcriptome FASTA (~300 MB)..."
    fetch_verified \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz" \
        "${TRANSCRIPTOME_FA}.gz" 30000000 "transcriptome FASTA"
    verify_gzip "${TRANSCRIPTOME_FA}.gz" "transcriptome FASTA"
    gunzip -f "${TRANSCRIPTOME_FA}.gz"
    [ -s "$TRANSCRIPTOME_FA" ] || fail "Transcriptome FASTA is empty after extraction: ${TRANSCRIPTOME_FA}"
    ok "Transcriptome: ${TRANSCRIPTOME_FA}"
else
    ok "Transcriptome exists"
fi

# ── STAR Index ──
# Validity check: a complete STAR index always contains 'SAindex'. A bare
# directory (e.g. from a crashed build) is treated as invalid and rebuilt.
STAR_INDEX_DIR="${REF_DIR}/star_index_hg38"
star_index_valid() { [ -s "${STAR_INDEX_DIR}/SAindex" ] && [ -s "${STAR_INDEX_DIR}/Genome" ]; }
if [ "$SKIP_INDEX" = false ] && ! star_index_valid && [ "$TOTAL_RAM_INT" -ge 32 ]; then
    if [ -d "$STAR_INDEX_DIR" ]; then
        warn "STAR index directory exists but is incomplete — rebuilding from scratch"
        rm -rf "$STAR_INDEX_DIR"
    fi
    log "Building STAR index (needs ~32 GB RAM, ~30-45 min)..."
    mkdir -p "$STAR_INDEX_DIR"
    STAR --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX_DIR" \
         --genomeFastaFiles "$GENOME_FA" \
         --sjdbGTFfile "$GTF_FILE" \
         --runThreadN "$CPU_CORES" \
         --sjdbOverhang 100
    star_index_valid || fail "STAR index build finished but index is incomplete: ${STAR_INDEX_DIR}"
    ok "STAR index: ${STAR_INDEX_DIR}"
elif [ "$TOTAL_RAM_INT" -lt 32 ] && ! star_index_valid; then
    warn "Skipping STAR index build — insufficient RAM (${TOTAL_RAM_GB} GB < 32 GB)"
    warn "Use HISAT2 for alignment, or build STAR index on a machine with 32+ GB RAM"
    STAR_INDEX_DIR=""
else
    ok "STAR index exists and is valid: ${STAR_INDEX_DIR}"
fi

# ── HISAT2 Index ──
# Validity check: the pre-built grch38_tran index ships 8 .ht2 files; the
# last one ('.8.ht2') is the sentinel for a complete extraction.
HISAT2_INDEX_DIR="${REF_DIR}/hisat2_index_hg38"
HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome_tran"
hisat2_index_valid() { [ -s "${HISAT2_INDEX_PREFIX}.1.ht2" ] && [ -s "${HISAT2_INDEX_PREFIX}.8.ht2" ]; }
if [ "$SKIP_INDEX" = false ] && ! hisat2_index_valid; then
    if [ -d "$HISAT2_INDEX_DIR" ]; then
        warn "HISAT2 index directory exists but is incomplete — re-downloading"
        rm -rf "$HISAT2_INDEX_DIR"
    fi
    log "Downloading HISAT2 pre-built index (~4 GB, ~10-30 min)..."
    mkdir -p "$HISAT2_INDEX_DIR"
    fetch_verified \
        "https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz" \
        "${REF_DIR}/grch38_tran.tar.gz" 3000000000 "HISAT2 index archive"
    verify_gzip "${REF_DIR}/grch38_tran.tar.gz" "HISAT2 index archive"
    tar -xzf "${REF_DIR}/grch38_tran.tar.gz" -C "$HISAT2_INDEX_DIR" --strip-components=1
    rm -f "${REF_DIR}/grch38_tran.tar.gz"
    hisat2_index_valid || fail "HISAT2 index extracted but is incomplete: ${HISAT2_INDEX_PREFIX}"
    ok "HISAT2 index: ${HISAT2_INDEX_PREFIX}"
else
    ok "HISAT2 index exists and is valid: ${HISAT2_INDEX_PREFIX}"
fi

# ── Salmon Index (with genome decoy for better mapping rates) ──
# Validity check: 'info.json' is written at the end of every successful build
# (all salmon versions). The second sentinel is format-dependent: salmon >=1.11
# uses the SSHash index format ('sshash.bin'); older salmon used the dense
# pufferfish format ('pos.bin'). Accept either so the check is version-robust.
SALMON_INDEX_DIR="${REF_DIR}/salmon_index_hg38"
salmon_index_valid() {
    [ -s "${SALMON_INDEX_DIR}/info.json" ] || return 1
    [ -s "${SALMON_INDEX_DIR}/sshash.bin" ] || [ -s "${SALMON_INDEX_DIR}/pos.bin" ]
}
if [ "$SKIP_INDEX" = false ] && ! salmon_index_valid; then
    if [ -d "$SALMON_INDEX_DIR" ]; then
        warn "Salmon index directory exists but is incomplete — rebuilding from scratch"
        rm -rf "$SALMON_INDEX_DIR"
    fi
    # Verify salmon works before attempting index build
    if ! salmon --version &>/dev/null; then
        warn "Salmon is not functional (possible GLIBC incompatibility)"
        warn "  Fix: conda install -c bioconda 'salmon>=1.10'"
        warn "  Skipping salmon index — salmon quantification will not be available"
    else
        log "Building Salmon index with decoys (~15-25 min, needs ~16 GB RAM)..."

        # Create gentrome (transcriptome + genome as decoy)
        GENTROME="${REF_DIR}/gentrome.fa"
        DECOYS="${REF_DIR}/decoys.txt"
        log "  Creating gentrome (transcriptome + genome decoys)..."
        cat "$TRANSCRIPTOME_FA" "$GENOME_FA" > "$GENTROME"
        grep "^>" "$GENOME_FA" | sed 's/>//' | cut -d' ' -f1 > "$DECOYS"

        # salmon index is allowed to fail without aborting the whole script;
        # temporarily lift the ERR trap for this one explicitly-handled command.
        salmon_rc=0
        salmon index -t "$GENTROME" -d "$DECOYS" \
            -i "$SALMON_INDEX_DIR" --threads "$CPU_CORES" --gencode || salmon_rc=$?
        if [ "$salmon_rc" -eq 0 ] && salmon_index_valid; then
            ok "Salmon index (with decoys): ${SALMON_INDEX_DIR}"
        else
            warn "Salmon index build failed or incomplete — salmon quantification will not be available"
            rm -rf "$SALMON_INDEX_DIR"
        fi
        rm -f "$GENTROME" "$DECOYS"
    fi
else
    ok "Salmon index exists and is valid: ${SALMON_INDEX_DIR}"
fi

# ── Annotation files ──
# RefSeq BED for RSeQC (from UCSC — reliable source)
BED_FILE="${ANN_DIR}/hg38_RefSeq.bed"
if [ ! -f "$BED_FILE" ] || [ ! -s "$BED_FILE" ]; then
    log "Downloading RefSeq BED from UCSC..."
    for attempt in 1 2 3; do
        wget -q --tries=3 --timeout=60 -O "${ANN_DIR}/refGene.txt.gz" \
            "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz" && break
        log "  Attempt ${attempt} failed, retrying..."
        sleep 5
    done
    if [ -f "${ANN_DIR}/refGene.txt.gz" ] && [ -s "${ANN_DIR}/refGene.txt.gz" ]; then
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
            warn "BED conversion produced empty file — gene body coverage (full) may not work"
            rm -f "$BED_FILE"
        fi
    else
        warn "Could not download RefSeq BED — gene body coverage (full) may not work"
        rm -f "${ANN_DIR}/refGene.txt.gz"
    fi
else
    ok "BED exists: ${BED_FILE}"
fi

# Housekeeping gene BED (Eisenberg & Levanon 2013, ~3800 genes from RSeQC)
HK_BED="${ANN_DIR}/hg38_housekeeping.bed"
if [ ! -f "$HK_BED" ]; then
    log "Downloading RSeQC housekeeping gene BED (hg38)..."
    HK_OK=false
    for attempt in 1 2 3; do
        wget -q --tries=3 --timeout=60 --content-disposition -O "${HK_BED}.gz" \
            "https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38.HouseKeepingGenes.bed.gz/download" && break
        log "  Attempt ${attempt} failed, retrying..."
        sleep 5
    done
    if [ -f "${HK_BED}.gz" ] && [ -s "${HK_BED}.gz" ]; then
        gunzip -f "${HK_BED}.gz"
        # gunzip may produce the original filename instead of our target name
        if [ -f "${ANN_DIR}/hg38.HouseKeepingGenes.bed" ] && [ ! -f "$HK_BED" ]; then
            mv "${ANN_DIR}/hg38.HouseKeepingGenes.bed" "$HK_BED"
        fi
        if [ -f "$HK_BED" ] && [ -s "$HK_BED" ]; then
            N_HK=$(wc -l < "$HK_BED")
            ok "Housekeeping BED: ${HK_BED} (${N_HK} transcripts)"
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
if [ ! -f "$REFFLAT" ] || [ ! -s "$REFFLAT" ]; then
    log "Downloading refFlat..."
    for attempt in 1 2 3; do
        wget -q --tries=3 --timeout=60 -O "${REFFLAT}.gz" \
            "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refFlat.txt.gz" && break
        log "  Attempt ${attempt} failed, retrying..."
        sleep 5
    done
    if [ -f "${REFFLAT}.gz" ] && [ -s "${REFFLAT}.gz" ]; then
        gunzip -f "${REFFLAT}.gz"
        ok "refFlat: ${REFFLAT} ($(wc -l < "$REFFLAT") entries)"
    else
        warn "Could not download refFlat — Picard RNA metrics may not work"
        rm -f "${REFFLAT}.gz"
    fi
else
    ok "refFlat exists"
fi

# ── FastQ Screen genomes (opt-in: --with-screen) ──
# `fastq_screen --get_genomes` is broken in v0.16.0 (the current bioconda
# release): an HTTP->HTTPS redirect drops a path slash, producing an
# unresolvable host. Instead of relying on it, we build the screen panel
# ourselves: download genome FASTA from stable infrastructure hosts (Ensembl,
# NCBI), build bowtie2 indexes locally, and generate fastq_screen.conf.
#
# Each genome is built independently and validity-checked, so an interrupted
# run resumes where it left off rather than rebuilding everything.
SCREEN_GENOMES_DIR="${SCREEN_DIR}/FastQ_Screen_Genomes"
SCREEN_CONF="${SCREEN_GENOMES_DIR}/fastq_screen.conf"
# The panel is complete only when the conf lists all 7 expected genomes.
# Checking merely "exists with >=1 DATABASE line" would let a stale partial
# conf from an interrupted run masquerade as a finished panel — so we count.
SCREEN_EXPECTED_GENOMES=7
screen_conf_valid() {
    [ -s "$SCREEN_CONF" ] || return 1
    local n
    n=$(grep -c '^DATABASE' "$SCREEN_CONF" 2>/dev/null || echo 0)
    [ "$n" -ge "$SCREEN_EXPECTED_GENOMES" ]
}

if [ "$WITH_SCREEN" = false ]; then
    ok "FastQ Screen skipped (run with --with-screen to build the genome panel)"
elif screen_conf_valid; then
    ok "FastQ Screen genomes already built: ${SCREEN_GENOMES_DIR}"
else
    log "Building FastQ Screen genome panel (--with-screen)..."
    log "  Downloads ~2.5 GB of FASTA and builds bowtie2 indexes."
    log "  Mammalian index builds are slow — expect several hours total."

    # A bowtie2 index is complete when both the forward '.1.bt2' and the
    # reverse '.rev.1.bt2' files exist and are non-empty.
    bt2_index_valid() {
        local prefix="$1"
        { [ -s "${prefix}.1.bt2" ] && [ -s "${prefix}.rev.1.bt2" ]; } || \
        { [ -s "${prefix}.1.bt2l" ] && [ -s "${prefix}.rev.1.bt2l" ]; }
    }

    # build_screen_genome <name> <fasta_url> <min_bytes> <gzipped:true|false>
    # Downloads FASTA into a per-genome dir, builds the bowtie2 index, and
    # leaves a clean prefix at <SCREEN_GENOMES_DIR>/<name>/<name>. Each genome
    # failure is non-fatal: a warn is emitted and the genome is left out of
    # the conf, so one bad download doesn't abort the whole panel.
    build_screen_genome() {
        local name="$1" url="$2" min_bytes="$3" gzipped="$4"
        local gdir="${SCREEN_GENOMES_DIR}/${name}"
        local prefix="${gdir}/${name}"
        local fasta="${gdir}/${name}.fa"

        if bt2_index_valid "$prefix"; then
            ok "  ${name}: bowtie2 index already built"
            return 0
        fi
        mkdir -p "$gdir"
        log "  ${name}: downloading FASTA..."

        local rc=0
        if [ "$gzipped" = "true" ]; then
            wget -q --show-progress --tries=3 --timeout=180 \
                -O "${fasta}.gz" "$url" || rc=$?
            if [ "$rc" -ne 0 ]; then
                warn "  ${name}: download failed — skipping this genome"; rm -f "${fasta}.gz"; return 1
            fi
            local size; size=$(stat -c%s "${fasta}.gz" 2>/dev/null || echo 0)
            if [ "$size" -lt "$min_bytes" ]; then
                warn "  ${name}: download too small (${size} B) — skipping"; rm -f "${fasta}.gz"; return 1
            fi
            gzip -t "${fasta}.gz" 2>/dev/null || { warn "  ${name}: corrupt gzip — skipping"; rm -f "${fasta}.gz"; return 1; }
            gunzip -f "${fasta}.gz"
        else
            wget -q --show-progress --tries=3 --timeout=180 \
                -O "$fasta" "$url" || rc=$?
            if [ "$rc" -ne 0 ]; then
                warn "  ${name}: download failed — skipping this genome"; rm -f "$fasta"; return 1
            fi
            local size; size=$(stat -c%s "$fasta" 2>/dev/null || echo 0)
            if [ "$size" -lt "$min_bytes" ]; then
                warn "  ${name}: download too small (${size} B) — skipping"; rm -f "$fasta"; return 1
            fi
        fi
        [ -s "$fasta" ] || { warn "  ${name}: FASTA empty after extraction — skipping"; return 1; }

        log "  ${name}: building bowtie2 index (${CPU_CORES} threads)..."
        if bowtie2-build --threads "$CPU_CORES" "$fasta" "$prefix" >/dev/null 2>&1 \
           && bt2_index_valid "$prefix"; then
            rm -f "$fasta"   # FASTA no longer needed once the index exists
            ok "  ${name}: bowtie2 index built"
            return 0
        else
            warn "  ${name}: bowtie2-build failed — skipping this genome"
            return 1
        fi
    }

    mkdir -p "$SCREEN_GENOMES_DIR"

    # Adapters: short, well-known Illumina sequences — written directly,
    # no download needed. bowtie2 needs them as a FASTA to index.
    ADAPTERS_DIR="${SCREEN_GENOMES_DIR}/Adapters"
    ADAPTERS_FA="${ADAPTERS_DIR}/Adapters.fa"
    if ! bt2_index_valid "${ADAPTERS_DIR}/Adapters"; then
        mkdir -p "$ADAPTERS_DIR"
        cat > "$ADAPTERS_FA" << 'ADAPTEREOF'
>Illumina_Universal_Adapter
AGATCGGAAGAG
>Illumina_TruSeq_Adapter_Read1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Illumina_TruSeq_Adapter_Read2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>Nextera_Transposase_Adapter
CTGTCTCTTATACACATCT
>Illumina_Small_RNA_3p_Adapter
TGGAATTCTCGGGTGCCAAGG
>Illumina_Small_RNA_5p_Adapter
GTTCAGAGTTCTACAGTCCGACGATC
>PolyA
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
>PolyG
GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
ADAPTEREOF
        log "  Adapters: building bowtie2 index..."
        if bowtie2-build --threads "$CPU_CORES" "$ADAPTERS_FA" \
                "${ADAPTERS_DIR}/Adapters" >/dev/null 2>&1; then
            ok "  Adapters: bowtie2 index built"
        else
            warn "  Adapters: bowtie2-build failed — skipping"
        fi
    else
        ok "  Adapters: bowtie2 index already built"
    fi

    # Genomes: each call is independently resumable and non-fatal on failure.
    # min_bytes is a conservative floor to catch error pages / truncation.
    build_screen_genome "Human" \
        "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
        500000000 true || true
    build_screen_genome "Mouse" \
        "https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz" \
        500000000 true || true
    build_screen_genome "Rat" \
        "https://ftp.ensembl.org/pub/release-111/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz" \
        500000000 true || true
    build_screen_genome "Ecoli" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz" \
        500000 true || true
    build_screen_genome "PhiX" \
        "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" \
        1000 true || true
    build_screen_genome "Vectors" \
        "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec" \
        100000 false || true

    # Generate fastq_screen.conf from whatever genomes built successfully.
    log "  Writing fastq_screen.conf..."
    {
        echo "# FastQ Screen configuration — auto-generated by ARACRA setup.sh"
        echo "# $(date)"
        echo ""
        echo "THREADS	${CPU_CORES}"
        echo ""
        for g in Human Mouse Rat Ecoli PhiX Vectors Adapters; do
            prefix="${SCREEN_GENOMES_DIR}/${g}/${g}"
            if bt2_index_valid "$prefix"; then
                printf 'DATABASE\t%s\t%s\n' "$g" "$prefix"
            fi
        done
    } > "$SCREEN_CONF"

    n_db=$(grep -c "^DATABASE" "$SCREEN_CONF" 2>/dev/null || echo 0)
    if [ "$n_db" -gt 0 ]; then
        ok "FastQ Screen panel ready: ${n_db} genome(s) — ${SCREEN_CONF}"
    else
        warn "FastQ Screen: no genomes built successfully — check warnings above"
        warn "  The pipeline will run, but set run_fastq_screen=false in nextflow.config"
    fi
fi

# ── tx2gene mapping (for Salmon — using Python for cross-system compatibility) ──
TX2GENE="${REF_DIR}/tx2gene.csv"
if [ ! -f "$TX2GENE" ] || [ ! -s "$TX2GENE" ]; then
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
# 5. CREATE WORK DIRECTORIES
# ══════════════════════════════════════════════════════════════════════════════
WORK_DIR="${HOME}/aracra_star_work/work"
OUT_DIR="${HOME}/aracra_star_work/results"
LOG_PATH="${HOME}/aracra_star_work/pipeline.log"
mkdir -p "$WORK_DIR" "$OUT_DIR"
ok "Work directories created"

# ══════════════════════════════════════════════════════════════════════════════
# 6. WRITE .env (auto-discovery results)
# ══════════════════════════════════════════════════════════════════════════════
cat > "${PIPELINE_DIR}/.env" << EOF
# ARACRA-STAR Pipeline — auto-generated $(date)
# System: ${TOTAL_RAM_GB} GB RAM | ${CPU_CORES} CPUs | Recommended aligner: ${RECOMMENDED_ALIGNER}
CONDA_BASE=${CONDA_BASE}
ENV_NAME=${ENV_NAME}
ENV_PATH=${ENV_PATH}
ENV_BIN=${ENV_PATH}/bin
WORK_DIR=${WORK_DIR}
OUT_DIR=${OUT_DIR}
LOG_FILE=${LOG_PATH}
DB_DIR=${DB_DIR}
STAR_INDEX=${STAR_INDEX_DIR}
HISAT2_INDEX=${HISAT2_INDEX_PREFIX}
SALMON_INDEX=${SALMON_INDEX_DIR}
GTF_PATH=${GTF_FILE}
BED_PATH=${BED_FILE}
HK_BED_PATH=${HK_BED}
REFFLAT_PATH=${REFFLAT}
SCREEN_CONF_PATH=${SCREEN_CONF}
TX2GENE_PATH=${TX2GENE}
TOTAL_RAM_GB=${TOTAL_RAM_GB}
CPU_CORES=${CPU_CORES}
RECOMMENDED_ALIGNER=${RECOMMENDED_ALIGNER}
HALF_CORES=${HALF_CORES}
EOF
ok ".env written"

# ══════════════════════════════════════════════════════════════════════════════
# 7. WRITE nextflow.config (auto-tuned to system)
# ══════════════════════════════════════════════════════════════════════════════
cat > "${PIPELINE_DIR}/nextflow.config" << NFEOF
/*
 * ARACRA-STAR Pipeline — nextflow.config
 * Auto-generated by setup.sh on $(date)
 * System: ${TOTAL_RAM_GB} GB RAM | ${CPU_CORES} CPUs
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

    // ── Resource detection ────────────────────────────────────────────────
    // CPU count is detected at pipeline-launch time by Groovy, so this adapts
    // automatically to whatever machine runs the pipeline.
    max_cpus            = Runtime.runtime.availableProcessors()

    // RAM: detected by setup.sh (reliably, via /proc/meminfo) and set here as
    // the default. It is a normal param, so moving to a different machine just
    // needs an override:  nextflow run ... --max_memory_gb 64
    // (We deliberately avoid JVM reflection for RAM — a bad cast there would
    //  break the entire config parse, and JVM internals vary across builds.)
    max_memory_gb       = ${TOTAL_RAM_INT}
    max_memory          = "\${params.max_memory_gb} GB"
    star_ram            = params.max_memory_gb * 600000000L

    // Derived parallelism — runtime-dynamic, follows max_cpus
    light_parallel      = Math.max(1, params.max_cpus - 1)
    post_parallel       = Math.max(1, params.max_cpus - 2)
    heavy_cpus          = Math.max(2, params.max_cpus.intdiv(2))
    heavy_parallel      = params.max_cpus <= 4 ? 1 : 2

    // FASTQ source (null = download from SRA)
    fastq_dir           = null
    trimmed_dir         = null

    // Step toggles (auto-tuned: aligner=${RECOMMENDED_ALIGNER})
    run_download            = true
    run_fastp               = true
    run_fastq_screen        = true
    run_star                = $([ "$RECOMMENDED_ALIGNER" = "star" ] && echo "true" || echo "false")
    run_hisat2              = $([ "$RECOMMENDED_ALIGNER" = "hisat2" ] && echo "true" || echo "false")
    run_samtools_stats      = true
    run_strandedness        = true
    run_coverage_hk         = true
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
    withLabel: 'download' {
        cpus     = 1
        memory   = '2 GB'
        maxForks = 3
    }
    withLabel: 'pre_qc' {
        cpus     = 1
        memory   = '2 GB'
        maxForks = params.light_parallel
    }
    withLabel: 'post_align' {
        cpus     = 1
        memory   = '4 GB'
        maxForks = params.post_parallel
    }
    withLabel: 'heavy' {
        cpus     = params.heavy_cpus
        memory   = params.max_memory
        maxForks = params.heavy_parallel
    }
    withLabel: 'screen' {
        cpus     = 2
        memory   = '4 GB'
        maxForks = 1
    }
    withLabel: 'final' {
        cpus   = params.max_cpus
        memory = params.max_memory
    }
    withLabel: 'minimal' {
        cpus   = 1
        memory = '2 GB'
    }
}

profiles {
    standard { process.executor = 'local' }
    slurm    { process.executor = 'slurm'; process.queue = 'batch' }
}

report   { enabled = true; overwrite = true; file = "\${params.outdir}/pipeline_report.html" }
timeline { enabled = true; overwrite = true; file = "\${params.outdir}/timeline.html" }
NFEOF
ok "nextflow.config written (tuned for ${CPU_CORES} cores, ${TOTAL_RAM_GB} GB RAM)"

# ══════════════════════════════════════════════════════════════════════════════
# 8. ENSURE scripts/ DIRECTORY
# ══════════════════════════════════════════════════════════════════════════════
if [ ! -d "${PIPELINE_DIR}/scripts" ]; then
    mkdir -p "${PIPELINE_DIR}/scripts"
    warn "Created scripts/ directory — copy run_deseq2.R and run_dromics.R there"
fi

# ══════════════════════════════════════════════════════════════════════════════
# SUMMARY
# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo -e "${GREEN}${BOLD}══════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}${BOLD}  Setup complete!${NC}"
echo -e "${GREEN}${BOLD}══════════════════════════════════════════════════════════════${NC}"
echo ""
echo "  System"
echo "    CPUs            : ${CPU_CORES}"
echo "    RAM             : ${TOTAL_RAM_GB} GB"
echo "    Recommended     : ${RECOMMENDED_ALIGNER^^}"
if [ "$TOTAL_RAM_INT" -lt 32 ]; then
echo -e "    ${RED}⚠ Low RAM — HISAT2 recommended over STAR${NC}"
fi
echo ""
echo "  Environment"
echo "    Conda           : ${CONDA_BASE}"
echo "    Env             : ${ENV_PATH}"
echo ""
echo "  References"
echo "    GTF             : ${GTF_FILE}"
echo "    STAR index      : ${STAR_INDEX_DIR:-SKIPPED (low RAM)}"
echo "    HISAT2 index    : ${HISAT2_INDEX_PREFIX}"
echo "    Salmon index    : ${SALMON_INDEX_DIR}"
echo ""
echo "  Directories"
echo "    Pipeline        : ${PIPELINE_DIR}"
echo "    Output          : ${OUT_DIR}"
echo "    Work            : ${WORK_DIR}"
echo ""
echo "  ➤  To start:  cd ${PIPELINE_DIR} && bash run_app.sh"
echo ""
