#!/usr/bin/env bash
# =============================================================================
#  ARACRA-STAR Pipeline — run_app.sh
#  Launches Streamlit app with test_ARCRA conda env
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PORT="${1:-8502}"

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║       ARACRA— Launching App                       ║"
echo "╚══════════════════════════════════════════════════════════╝"

# Try to detect conda
CONDA_BASE="${HOME}/miniforge3"
[ ! -d "$CONDA_BASE" ] && CONDA_BASE="${HOME}/miniconda3"
[ ! -d "$CONDA_BASE" ] && CONDA_BASE="$(conda info --base 2>/dev/null || echo "")"

if [ -z "$CONDA_BASE" ] || [ ! -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
    echo "ERROR: Cannot find conda. Install miniforge3 or set CONDA_BASE."
    exit 1
fi

source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Activate rnaseq_pipeline env — try common names
ENV_NAME=""
for try_name in test_ARACRA rnaseq_pipeline2 rnaseq_pipeline; do
    if conda activate "$try_name" 2>/dev/null; then
        ENV_NAME="$try_name"
        echo "  ✔ Activated: ${ENV_NAME}"
        break
    fi
done
if [ -z "$ENV_NAME" ]; then
    echo "ERROR: Environment 'test_ARACRA' (or rnaseq_pipeline2) not found."
    echo "  Run setup.sh first."
    exit 1
fi

# Check streamlit
if ! command -v streamlit &>/dev/null; then
    echo "  Installing streamlit..."
    pip install streamlit pandas openpyxl --quiet
fi

# Write .env if not exists
ENV_PATH="${CONDA_BASE}/envs/${ENV_NAME}"
if [ ! -f "${PIPELINE_DIR}/.env" ]; then
    cat > "${PIPELINE_DIR}/.env" <<EOF
CONDA_BASE=${CONDA_BASE}
ENV_NAME=${ENV_NAME}
ENV_PATH=${ENV_PATH}
ENV_BIN=${ENV_PATH}/bin
WORK_DIR=${HOME}/aracra_star_work/work
OUT_DIR=${HOME}/aracra_star_work/results
LOG_FILE=${HOME}/aracra_star_work/pipeline.log
STAR_INDEX=${HOME}/databases/hg38_reference/star_index_hg38
HISAT2_INDEX=${HOME}/databases/hg38_reference/hisat2_index_hg38/genome_tran
SALMON_INDEX=${HOME}/databases/hg38_reference/salmon_index_hg38
GTF_PATH=${HOME}/databases/hg38_reference/gencode.v44.primary_assembly.annotation.gtf
BED_PATH=${HOME}/databases/annotations/hg38_RefSeq.bed
HK_BED_PATH=${HOME}/databases/annotations/hg38_housekeeping.bed
REFFLAT_PATH=${HOME}/databases/annotations/hg38_refFlat.txt
SCREEN_CONF_PATH=${HOME}/databases/fastq_screen_genomes/FastQ_Screen_Genomes/fastq_screen.conf
EOF
    echo "  ✔ .env created"
fi

mkdir -p "${HOME}/aracra_star_work/work" "${HOME}/aracra_star_work/results"

# Quick tool check
echo ""
echo "  Tool check:"
for tool in nextflow STAR hisat2 samtools fastp featureCounts salmon qualimap multiqc streamlit; do
    if command -v "$tool" &>/dev/null; then
        echo "    ✔ ${tool}"
    else
        echo "    ✘ ${tool}"
    fi
done

# Stop old instance
OLD_PID=$(lsof -ti tcp:"${PORT}" 2>/dev/null || true)
if [ -n "$OLD_PID" ]; then
    echo "  Stopping old instance (PID ${OLD_PID})..."
    kill "$OLD_PID" 2>/dev/null || true
    sleep 1
fi

# Launch
echo ""
APP_LOG="${HOME}/aracra_star_work/streamlit.log"
nohup streamlit run "${PIPELINE_DIR}/aracra_star_app.py" \
    --server.port "$PORT" \
    --server.headless true \
    --server.address 0.0.0.0 \
    --server.fileWatcherType none \
    > "$APP_LOG" 2>&1 &

APP_PID=$!
echo "$APP_PID" > "${PIPELINE_DIR}/.app.pid"
sleep 2

echo "  ✔ App running (PID ${APP_PID})"
echo ""
echo "  Open browser:"
echo "  ➤ http://localhost:${PORT}"
echo "  ➤ http://$(hostname -I 2>/dev/null | awk '{print $1}' || echo 'your-ip'):${PORT}"
echo ""
echo "  Stop: kill \$(cat ${PIPELINE_DIR}/.app.pid)"
echo "  Log:  tail -f ${APP_LOG}"
echo ""
