#!/usr/bin/env bash
# =============================================================================
#  ARACRA Pipeline — run_app.sh
#  Launches the Streamlit app inside the ARACRA conda env.
# =============================================================================
set -euo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib/aracra_common.sh
source "${PIPELINE_DIR}/lib/aracra_common.sh"

PORT="${1:-$ARACRA_DEFAULT_PORT}"

echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║       ARACRA— Launching App                       ║"
echo "╚══════════════════════════════════════════════════════════╝"

# Try to detect conda
CONDA_BASE="$(aracra_find_conda || true)"

if [ -z "$CONDA_BASE" ]; then
    echo "ERROR: Cannot find conda. Run 'bash setup.sh' first, or install Miniforge3."
    exit 1
fi

source "${CONDA_BASE}/etc/profile.d/conda.sh"

# Activate the ARACRA env — try the canonical name, then legacy names
ENV_NAME=""
for try_name in "${ARACRA_ENV_FALLBACKS[@]}"; do
    if conda activate "$try_name" 2>/dev/null; then
        ENV_NAME="$try_name"
        echo "  ✔ Activated: ${ENV_NAME}"
        break
    fi
done
if [ -z "$ENV_NAME" ]; then
    echo "ERROR: Environment '${ARACRA_ENV_NAME}' not found."
    echo "  Run setup.sh first."
    exit 1
fi

# Check streamlit
if ! command -v streamlit &>/dev/null; then
    echo "  Installing streamlit..."
    pip install streamlit pandas openpyxl --quiet
fi

# Write .env if missing. setup.sh is the normal author; this is a fallback that
# uses the SAME writer, so the key set can never drift between the two scripts.
ENV_PATH="${CONDA_BASE}/envs/${ENV_NAME}"
if [ ! -f "${PIPELINE_DIR}/.env" ]; then
    echo "  ! .env not found — writing defaults (run setup.sh for a tuned config)"
    aracra_write_env "$PIPELINE_DIR" "$CONDA_BASE" "$ENV_PATH" "$ARACRA_DB_DIR" "run_app.sh"
    echo "  ✔ .env created"
fi

mkdir -p "${ARACRA_WORK_ROOT}/work" "${ARACRA_WORK_ROOT}/results"

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
APP_LOG="${ARACRA_WORK_ROOT}/streamlit.log"
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
