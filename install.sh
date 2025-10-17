#!/bin/bash

# ==============================================================================
# ARACRA - Automated RNA-seq and Comprehensive 'omics Analysis
# Installation Script for Linux
#
# This script creates a Conda environment with all necessary dependencies.
#
# Prerequisites:
#   - A Linux-based operating system.
#   - Conda (Miniconda or Anaconda) must be installed. Mamba is recommended for
#     faster installation.
#
# Usage:
#   1. Make the script executable: chmod +x install.sh
#   2. Run the script: ./install.sh
# ==============================================================================

set -e # Exit immediately if a command exits with a non-zero status.

# --- Configuration ---
ENV_NAME="aracra_env"
PYTHON_VERSION="3.9"

# --- Helper Functions ---
# For colored output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
    exit 1
}

# --- Main Script ---

info "Starting ARACRA environment setup..."

# 1. Check for Conda or Mamba
# We prefer mamba for its speed, but fall back to conda.
if command -v mamba &> /dev/null; then
    INSTALLER="mamba"
    info "Mamba detected. Using Mamba for a faster installation."
elif command -v conda &> /dev/null; then
    INSTALLER="conda"
    info "Conda detected. Using Conda for installation."
    warn "Installation may be slow. For faster setup, consider installing Mamba: 'conda install -n base -c conda-forge mamba'"
else
    error "Conda/Mamba not found. Please install Miniconda or Anaconda first."
    echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# 2. Check if the environment already exists
if $INSTALLER env list | grep -q "^${ENV_NAME}\s"; then
    warn "Conda environment '${ENV_NAME}' already exists."
    read -p "Do you want to remove the existing environment and reinstall? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        info "Removing existing environment '${ENV_NAME}'..."
        $INSTALLER env remove --name "${ENV_NAME}"
    else
        info "Skipping installation. To use the application, activate the environment:"
        info "conda activate ${ENV_NAME}"
        exit 0
    fi
fi


# 3. Create the Conda environment and install all dependencies
info "Creating Conda environment '${ENV_NAME}'... This may take several minutes."

# List of all dependencies for a single install command
# Channels are prioritized: conda-forge > bioconda > defaults
# This ensures package compatibility.
$INSTALLER create -n "${ENV_NAME}" --yes \
    -c conda-forge \
    -c bioconda \
    -c defaults \
    "python=${PYTHON_VERSION}" \
    \
    # Bioinformatics Command-Line Tools
    "sra-tools" \
    "hisat2" \
    "samtools" \
    "subread" \
    "trimmomatic" \
    "openjdk" \
    \
    # Python Packages
    "pip" \
    "streamlit" \
    "pandas" \
    "numpy" \
    "psutil" \
    "openpyxl" \
    \
    # R and R Packages
    "r-base" \
    "r-essentials" \
    "r-readxl" \
    "r-stringr" \
    "r-ggrepel" \
    "r-dromics" \
    "r-jsonlite" \
    "bioconductor-deseq2" \
    "bioconductor-edger" \
    "bioconductor-sva" \
    "bioconductor-org.hs.eg.db"

if [ $? -ne 0 ]; then
    error "Failed to create the Conda environment. Please check the error messages above."
fi

# 4. Final instructions
info "Installation complete!"
echo
info "--- NEXT STEPS ---"
echo
echo -e "1. Activate the new environment:"
echo -e "   ${YELLOW}conda activate ${ENV_NAME}${NC}"
echo
echo -e "2. Run the ARACRA application:"
echo -e "   ${YELLOW}streamlit run rnaseq_pipeline_gui.py${NC}"
echo
echo -e "3. ${RED}IMPORTANT:${NC} First-Time Configuration:"
echo -e "   Before running the pipeline, you must download necessary data files and configure their paths in the GUI sidebar:"
echo -e "   - ${YELLOW}Trimmomatic JAR Path:${NC} The JAR is inside your new conda env. You can find it with:"
echo -e "     (conda activate ${ENV_NAME} && find \$CONDA_PREFIX -name \"trimmomatic*jar\")"
echo -e "   - ${YELLOW}Adapter Sequences:${NC} Download common Illumina adapter sequences if you don't have them."
echo -e "   - ${YELLOW}HISAT2 Index Path:${NC} Download a pre-built index for your reference genome (e.g., GRCh38)."
echo -e "   - ${YELLOW}Reference GTF:${NC} Download the corresponding GTF annotation file for your genome."
echo
info "Setup finished successfully."
