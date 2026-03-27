#!/usr/bin/env bash
# =============================================================================
#  ARACRA-STAR Pipeline — setup.sh (v3)
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

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Parse flags ───────────────────────────────────────────────────────────────
SKIP_INDEX=false
SKIP_SCREEN=false
DB_DIR="${HOME}/databases"
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
echo "║     ARACRA-STAR Pipeline — Setup  (v3)                  ║"
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
    curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
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
    nextflow \
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

# ── Install rtracklayer (requires source patch for GCC 15 compatibility) ──
# GCC 15 treats empty parentheses () as zero arguments (C23 behaviour).
# The UCSC source in rtracklayer uses a function pointer named 'free' with
# signature void(*free)() which GCC 15 rejects. We patch the source before
# compiling to fix the signature to void(*free)(void*).
log "Installing rtracklayer (applying GCC 15 source patch)..."
if ! "${ENV_PATH}/bin/Rscript" -e 'requireNamespace("rtracklayer", quietly=TRUE)' 2>/dev/null | grep -q TRUE; then
    RTRACKLAYER_TMP="/tmp/rtracklayer_build_$$"
    mkdir -p "$RTRACKLAYER_TMP"

    # Get the correct version for the installed Bioconductor
    BIOC_VER=$("${ENV_PATH}/bin/Rscript" -e 'cat(as.character(BiocManager::version()))' 2>/dev/null)
    RT_URL="https://bioconductor.org/packages/${BIOC_VER}/bioc/src/contrib/rtracklayer_1.66.0.tar.gz"

    log "  Downloading rtracklayer source..."
    curl -L -o "${RTRACKLAYER_TMP}/rtracklayer.tar.gz" "$RT_URL" \
        || fail "Failed to download rtracklayer source"

    tar xzf "${RTRACKLAYER_TMP}/rtracklayer.tar.gz" -C "$RTRACKLAYER_TMP"

    log "  Patching for GCC 15 (void(*free)() → void(*free)(void*))..."
    # Rename parameter 'free' to 'freeEl' to avoid collision with stdlib macro
    sed -i \
        's/void (\*free)()/void (*freeEl)(void*)/g;
         s/else if (free != NULL)/else if (freeEl != NULL)/g;
         s/\bfree(el)\b/freeEl(el)/g' \
        "${RTRACKLAYER_TMP}/rtracklayer/src/ucsc/common.c"

    # Also fix the header declaration
    sed -i 's/void (\*free)()/void (*free)(void*)/g' \
        "${RTRACKLAYER_TMP}/rtracklayer/src/ucsc/common.h"

    log "  Compiling and installing rtracklayer..."
    "${ENV_PATH}/bin/R" CMD INSTALL "${RTRACKLAYER_TMP}/rtracklayer/" \
        || fail "rtracklayer installation failed — see output above"

    rm -rf "$RTRACKLAYER_TMP"
    ok "rtracklayer installed (GCC 15 patched)"
else
    ok "rtracklayer already installed"
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
if [ ! -f "$GTF_FILE" ]; then
    log "Downloading GTF annotation (~50 MB)..."
    wget -q --show-progress -O "${GTF_FILE}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz"
    gunzip "${GTF_FILE}.gz"
    ok "GTF: ${GTF_FILE}"
else
    ok "GTF exists: ${GTF_FILE}"
fi

# ── Genome FASTA ──
GENOME_FA="${REF_DIR}/GRCh38.primary_assembly.genome.fa"
if [ ! -f "$GENOME_FA" ]; then
    log "Downloading genome FASTA (~900 MB)..."
    wget -q --show-progress -O "${GENOME_FA}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz"
    gunzip "${GENOME_FA}.gz"
    ok "Genome: ${GENOME_FA}"
else
    ok "Genome exists: ${GENOME_FA}"
fi

# ── Transcriptome FASTA (for Salmon) ──
TRANSCRIPTOME_FA="${REF_DIR}/gencode.v44.transcripts.fa"
if [ ! -f "$TRANSCRIPTOME_FA" ]; then
    log "Downloading transcriptome FASTA (~300 MB)..."
    wget -q --show-progress -O "${TRANSCRIPTOME_FA}.gz" \
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz"
    gunzip "${TRANSCRIPTOME_FA}.gz"
    ok "Transcriptome: ${TRANSCRIPTOME_FA}"
else
    ok "Transcriptome exists"
fi

# ── STAR Index ──
STAR_INDEX_DIR="${REF_DIR}/star_index_hg38"
if [ "$SKIP_INDEX" = false ] && [ ! -d "$STAR_INDEX_DIR" ] && [ "$TOTAL_RAM_INT" -ge 32 ]; then
    log "Building STAR index (needs ~32 GB RAM, ~30-45 min)..."
    mkdir -p "$STAR_INDEX_DIR"
    STAR --runMode genomeGenerate \
         --genomeDir "$STAR_INDEX_DIR" \
         --genomeFastaFiles "$GENOME_FA" \
         --sjdbGTFfile "$GTF_FILE" \
         --runThreadN "$CPU_CORES" \
         --sjdbOverhang 100
    ok "STAR index: ${STAR_INDEX_DIR}"
elif [ "$TOTAL_RAM_INT" -lt 32 ] && [ ! -d "$STAR_INDEX_DIR" ]; then
    warn "Skipping STAR index build — insufficient RAM (${TOTAL_RAM_GB} GB < 32 GB)"
    warn "Use HISAT2 for alignment, or build STAR index on a machine with 32+ GB RAM"
    STAR_INDEX_DIR=""
else
    ok "STAR index exists: ${STAR_INDEX_DIR}"
fi

# ── HISAT2 Index ──
HISAT2_INDEX_DIR="${REF_DIR}/hisat2_index_hg38"
HISAT2_INDEX_PREFIX="${HISAT2_INDEX_DIR}/genome_tran"
if [ "$SKIP_INDEX" = false ] && [ ! -f "${HISAT2_INDEX_PREFIX}.1.ht2" ]; then
    log "Downloading HISAT2 pre-built index (~4 GB, ~10-30 min)..."
    mkdir -p "$HISAT2_INDEX_DIR"
    wget -q --show-progress -O "${REF_DIR}/grch38_tran.tar.gz" \
        "https://genome-idx.s3.amazonaws.com/hisat/grch38_tran.tar.gz"
    tar -xzf "${REF_DIR}/grch38_tran.tar.gz" -C "$HISAT2_INDEX_DIR" --strip-components=1
    rm -f "${REF_DIR}/grch38_tran.tar.gz"
    ok "HISAT2 index: ${HISAT2_INDEX_PREFIX}"
else
    ok "HISAT2 index exists: ${HISAT2_INDEX_PREFIX}"
fi

# ── Salmon Index (with genome decoy for better mapping rates) ──
SALMON_INDEX_DIR="${REF_DIR}/salmon_index_hg38"
if [ "$SKIP_INDEX" = false ] && [ ! -d "$SALMON_INDEX_DIR" ]; then
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

        if salmon index -t "$GENTROME" -d "$DECOYS" \
            -i "$SALMON_INDEX_DIR" --threads "$CPU_CORES" --gencode; then
            ok "Salmon index (with decoys): ${SALMON_INDEX_DIR}"
        else
            warn "Salmon index build failed — salmon quantification will not be available"
            rm -rf "$SALMON_INDEX_DIR"
        fi
        rm -f "$GENTROME"
    fi
else
    ok "Salmon index exists: ${SALMON_INDEX_DIR}"
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

# ── FastQ Screen genomes ──
SCREEN_CONF="${SCREEN_DIR}/FastQ_Screen_Genomes/fastq_screen.conf"
if [ "$SKIP_SCREEN" = false ] && [ ! -f "$SCREEN_CONF" ]; then
    if command -v fastq_screen &>/dev/null; then
        log "Downloading FastQ Screen genomes (~30-60 min, ~30 GB)..."
        log "  This downloads bowtie2 indices for Human, Mouse, Rat, E.coli, etc."
        log "  Run with --skip-screen to skip this step and add genomes later."
        fastq_screen --get_genomes --outdir "$SCREEN_DIR" \
            || warn "FastQ Screen genome download failed — contamination check will be skipped"
    else
        warn "fastq_screen not found in PATH — skipping genome download"
        warn "  Install: conda install -c bioconda fastq-screen"
        warn "  Then re-run setup.sh, or run: fastq_screen --get_genomes --outdir ${SCREEN_DIR}"
    fi
elif [ "$SKIP_SCREEN" = true ]; then
    ok "FastQ Screen skipped (--skip-screen)"
else
    ok "FastQ Screen genomes exist"
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
    max_memory          = "${TOTAL_RAM_INT} GB"
    star_ram            = $(( TOTAL_RAM_INT * 600000000 ))

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

def _totalCores = Runtime.runtime.availableProcessors()

process {
    withLabel: 'download' {
        cpus     = 1
        memory   = '2 GB'
        maxForks = 3
    }
    withLabel: 'pre_qc' {
        cpus     = 1
        memory   = '2 GB'
        maxForks = Math.max(1, _totalCores - 1)
    }
    withLabel: 'post_align' {
        cpus     = 1
        memory   = '4 GB'
        maxForks = Math.max(1, _totalCores - 2)
    }
    withLabel: 'heavy' {
        cpus     = Math.max(2, (int)(_totalCores / 2))
        memory   = params.max_memory
        maxForks = ${STAR_PARALLEL}
    }
    withLabel: 'screen' {
        cpus     = 2
        memory   = '4 GB'
        maxForks = 1
    }
    withLabel: 'final' {
        cpus   = _totalCores
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

report   { enabled = true; file = "\${params.outdir}/pipeline_report.html" }
timeline { enabled = true; file = "\${params.outdir}/timeline.html" }
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
