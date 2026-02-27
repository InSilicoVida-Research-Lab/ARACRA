# 🧬 ARACRA-STAR Pipeline

**Automated RNA-seq & Chemical Risk Assessment Pipeline**

A complete Nextflow pipeline with Streamlit GUI for RNA-seq and TempO-Seq analysis:
SRA → fastp → STAR/HISAT2 → QC → featureCounts/Salmon → DESeq2 → DRomics/BMD

---

## Features

- **Dual aligner support**: STAR (fast, needs 32 GB RAM) or HISAT2 (memory-efficient, 8 GB RAM)
- **Quantification choice**: featureCounts, Salmon, or both
- **Comprehensive QC**: fastp, FastQ Screen, RSeQC strandedness, gene body coverage, Picard RNA metrics, MultiQC
- **Dose-response analysis**: DESeq2 differential expression + DRomics benchmark dose (BMD)
- **Direct analysis mode**: Upload count matrix → skip alignment → DESeq2/DRomics instantly
- **Web GUI**: Streamlit app with live log, progress tracking, results browser
- **Detached execution**: Pipeline runs in background — safe to close browser/terminal
- **Auto-tuned**: `setup.sh` detects CPU cores, RAM, and configures optimal parallelism

---

## System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU       | 4 cores | 8+ cores |
| RAM       | 8 GB (HISAT2 only) | 32+ GB (STAR + HISAT2) |
| Disk      | 100 GB free | 250+ GB free |
| OS        | Ubuntu 20.04+ / WSL2 | Ubuntu 22.04+ |

> **RAM < 32 GB?** Use HISAT2 — the pipeline will warn you automatically and recommend it.

---

## Installation

### Step 0: Windows Users — Enable WSL2 + Ubuntu

> **Linux / macOS users**: Skip to [Step 1](#step-1-install-miniforge).

#### 0.1 Enable WSL2

Open **PowerShell as Administrator** (right-click → "Run as administrator") and run:

```powershell
wsl --install
```

This installs WSL2 and Ubuntu automatically. **Restart your computer** when prompted.

#### 0.2 Set Up Ubuntu

After restart, Ubuntu will open automatically and ask you to create a username and password:

```
Enter new UNIX username: saurav
New password: ********
```

> Remember this password — you'll need it for `sudo` commands.

#### 0.3 Update Ubuntu

```bash
sudo apt update && sudo apt upgrade -y
```

#### 0.4 Install Essential Build Tools

```bash
sudo apt install -y build-essential wget curl git unzip gfortran \
    libxml2-dev libcurl4-openssl-dev libssl-dev libfontconfig1-dev \
    libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
    libharfbuzz-dev libfribidi-dev liblapack-dev libblas-dev bc
```

#### 0.5 Increase WSL Memory (Important!)

By default, WSL2 only uses half your RAM. To allow STAR to work, create/edit the WSL config:

Open **PowerShell** (not Ubuntu) and run:

```powershell
notepad "$env:USERPROFILE\.wslconfig"
```

Add these lines (adjust `memory` to your total RAM minus 4 GB):

```ini
[wsl2]
memory=28GB
processors=8
swap=8GB
```

Then restart WSL:

```powershell
wsl --shutdown
```

Re-open Ubuntu from the Start menu.

#### 0.6 Accessing Windows Files from WSL

Your Windows drives are available under `/mnt/`:

```bash
# Your Windows Desktop
ls /mnt/c/Users/YourName/Desktop/

# Copy files from Windows to Linux home
cp /mnt/c/Users/YourName/Downloads/metadata.xlsx ~/
```

> **Important**: Always run the pipeline from the Linux filesystem (`~/`), not from `/mnt/c/`. Running from Windows drives is extremely slow.

---

### Step 1: Install Miniforge

> Skip if you already have `conda` or `mamba` installed.

```bash
curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o miniforge.sh
bash miniforge.sh -b -p ~/miniforge3
rm miniforge.sh

# Activate conda
source ~/miniforge3/etc/profile.d/conda.sh
conda init bash
```

Close and reopen terminal, then verify:

```bash
conda --version
# Should print: conda 25.x.x
```

---

### Step 2: Download the Pipeline

```bash
cd ~
mkdir -p ARACRA_STAR && cd ARACRA_STAR
```

Place all pipeline files in this directory:

```
~/ARACRA_STAR/
├── setup.sh
├── main.nf
├── aracra_star_app.py
├── run_app.sh
├── scripts/
│   ├── run_deseq2.R
│   └── run_dromics.R
```

---

### Step 3: Run Setup

```bash
cd ~/ARACRA_STAR
bash setup.sh
```

This will automatically:

1. **Detect your system** — CPUs, RAM, disk space
2. **Create conda environment** — all bioinformatics tools (STAR, HISAT2, Salmon, fastp, samtools, featureCounts, R, etc.)
3. **Download references** (~15 GB total):
   - Human genome FASTA (GRCh38)
   - GTF annotation (GENCODE v44)
   - STAR index (built if RAM ≥ 32 GB)
   - HISAT2 index (pre-built, always downloaded)
   - Salmon index
   - FastQ Screen contamination genomes
   - RSeQC/Picard annotation files
4. **Write configuration** — `.env` and `nextflow.config` tuned to your machine

> **First run takes 1–2 hours** (mostly downloading references). Subsequent runs skip existing files.

#### Setup Options

```bash
# Skip index building (if you already have them)
bash setup.sh --skip-index

# Skip FastQ Screen genomes download
bash setup.sh --skip-screen

# Use a custom database directory
bash setup.sh --db-dir=/data/references
```

#### Verify Setup

After setup completes, you should see:

```
══════════════════════════════════════════════════════════════
  Setup complete!
══════════════════════════════════════════════════════════════

  System
    CPUs            : 8
    RAM             : 62.6 GB
    Recommended     : STAR

  References
    GTF             : ~/databases/hg38_reference/gencode.v44...
    STAR index      : ~/databases/hg38_reference/star_index_hg38
    HISAT2 index    : ~/databases/hg38_reference/hisat2_index_hg38/genome_tran
    Salmon index    : ~/databases/hg38_reference/salmon_index_hg38

  ➤  To start:  cd ~/ARACRA_STAR && bash run_app.sh
```

---

### Step 4: Launch the App

```bash
cd ~/ARACRA_STAR
bash run_app.sh
```

Open your browser:

```
http://localhost:8502
```

> **Windows WSL users**: This URL works directly in your Windows browser (Edge, Chrome, etc.).

---

## Usage

### Full Pipeline (SRA → DESeq2/DRomics)

1. Open the **🚀 Full Pipeline** tab
2. Upload your **metadata Excel file** (must contain SRR accession column)
3. Select **aligner** (STAR or HISAT2), **quantification** (featureCounts/Salmon), and **QC options**
4. Set **Treatment** and **Control** groups
5. Click **🚀 Start Full Pipeline**

#### Required Metadata Format

Your Excel/CSV must have these columns:

| Sample_Name | Treatment | Dose | Batch | Type |
|-------------|-----------|------|-------|------|
| SRR29684827 | Bisphenol A | 0.1 | 1 | treatment |
| SRR29684828 | Bisphenol A | 1.0 | 1 | treatment |
| SRR29684829 | DMSO | 0 | 1 | control |

### Direct Analysis (Count Matrix → DESeq2/DRomics)

1. Open the **⚡ Direct Analysis** tab
2. Upload a **count matrix** (featureCounts `.out`, `.csv`, or `.tsv`)
3. Upload **metadata** (same format as above)
4. Click **⚡ Start Direct Analysis**

### Command Line Usage

```bash
# Activate environment
conda activate rnaseq_pipeline

# Full pipeline with STAR
nextflow run main.nf --metadata metadata.xlsx --aligner star

# Full pipeline with HISAT2
nextflow run main.nf --metadata metadata.xlsx --aligner hisat2

# Direct analysis (skip alignment)
nextflow run main.nf --count_matrix counts.csv --metadata metadata.csv \
    --treatment "Bisphenol A" --control "DMSO" --run_name my_analysis

# Resume after interruption
nextflow run main.nf --metadata metadata.xlsx -resume
```

---

## Pipeline Architecture

```
┌─────────────────────────────────────────────────────┐
│ Phase 1: Download & QC (parallel)                    │
│   SRA Download (3 parallel) → fastp (7 parallel)    │
│   FastQ Screen (independent, 2 cores)               │
├─────────────────────────────────────────────────────┤
│ Phase 2: Alignment                                   │
│   STAR (shared memory, 4 cores × 2 parallel)        │
│   — OR —                                            │
│   HISAT2 (4 cores × 2 parallel, ~8 GB RAM)          │
├─────────────────────────────────────────────────────┤
│ Phase 3: Post-alignment QC (6 parallel)              │
│   samtools stats │ strandedness │ coverage │ Picard  │
├─────────────────────────────────────────────────────┤
│ Phase 4: Quantification (all cores)                  │
│   featureCounts (all BAMs) and/or Salmon             │
├─────────────────────────────────────────────────────┤
│ Phase 5: Analysis                                    │
│   DESeq2 → Volcano/PCA │ DRomics → BMD curves      │
│   MultiQC report                                     │
└─────────────────────────────────────────────────────┘
```

### STAR vs HISAT2

| Feature | STAR | HISAT2 |
|---------|------|--------|
| RAM needed | ~32 GB | ~8 GB |
| Speed | Faster | Slightly slower |
| Shared memory | Yes (loads once) | No (per-sample) |
| Splice-aware | Yes | Yes |
| When to use | 32+ GB RAM systems | Low-memory systems |

---

## Output Structure

```
results/
├── data/                    # Raw FASTQ files
├── trimmed/                 # Trimmed FASTQ files
├── qc/
│   ├── fastp/               # fastp HTML/JSON reports
│   ├── fastq_screen/        # Contamination reports
│   ├── samtools/            # flagstat, idxstats
│   ├── strandedness/        # Library strandedness
│   ├── genebody_coverage/   # RSeQC coverage plots
│   └── picard/              # Picard RNA metrics
├── star/ (or hisat2/)       # BAM files + alignment logs
├── featurecounts/           # Raw gene counts
├── salmon/                  # Salmon quantification
├── counts/
│   ├── gene_count_matrix.csv    # ← Main count matrix
│   └── sample_metadata.csv
├── multiqc/
│   └── multiqc_report.html      # ← QC summary
├── runs/
│   └── BPA_vs_DMSO_20260226/
│       ├── deg_results/
│       │   ├── All_Results.csv
│       │   ├── Custom_Filtered_DEGs.csv
│       │   ├── Volcano_Plot.png
│       │   ├── Final_PCA_Plot.png
│       │   └── analysis_summary.json
│       └── dromics_results/
│           ├── bmd_results.csv
│           ├── dose_response_curves.png
│           ├── bmd_distribution.png
│           └── dromics_summary.json
├── pipeline_report.html
└── timeline.html
```

---

## Troubleshooting

### STAR fails with "out of memory"

Your system doesn't have enough RAM for STAR. Switch to HISAT2:

```bash
nextflow run main.nf --metadata metadata.xlsx --aligner hisat2
```

Or select HISAT2 in the Streamlit app.

### Pipeline stuck at STAR_LOAD_GENOME

A previous crashed run left orphaned shared memory. Clean it up:

```bash
STAR --genomeDir ~/databases/hg38_reference/star_index_hg38 --genomeLoad Remove
```

> The pipeline now does this automatically before loading, but manual cleanup may be needed after a hard kill.

### "No common samples" error in DESeq2

The `Sample_Name` column in your metadata must match the SRR IDs in the count matrix columns. Check:

```bash
# See what's in the count matrix
head -1 results/counts/gene_count_matrix.csv

# See what's in metadata
head results/metadata/metadata.csv
```

### WSL: "Cannot allocate memory"

Increase WSL memory limit — edit `~/.wslconfig` on Windows (see Step 0.5 above).

### Pipeline reports fail to render

This warning is cosmetic and doesn't affect results:

```
WARN: Failed to render execution report
```

The pipeline completed successfully — check your results folder.

### Streamlit app not loading

```bash
# Check if it's running
cat ~/ARACRA_STAR/.app.pid
ps aux | grep streamlit

# Restart
kill $(cat ~/ARACRA_STAR/.app.pid) 2>/dev/null
bash run_app.sh
```

### Conda environment issues

```bash
# Remove and recreate
conda deactivate
conda env remove -n rnaseq_pipeline
bash setup.sh
```

---

## Tool Versions

| Tool | Version | Purpose |
|------|---------|---------|
| Nextflow | 25.10+ | Workflow manager |
| STAR | 2.7.11b | Splice-aware aligner |
| HISAT2 | 2.2.1 | Memory-efficient aligner |
| fastp | 1.0+ | Read trimming & QC |
| samtools | 1.22+ | BAM processing |
| featureCounts | 2.0+ | Gene-level counting |
| Salmon | 1.10+ | Transcript quantification |
| FastQ Screen | 0.16+ | Contamination screening |
| RSeQC | 5.0+ | RNA-seq QC |
| Picard | 3.0+ | RNA metrics |
| MultiQC | 1.33+ | QC report aggregation |
| DESeq2 | 1.42+ | Differential expression |
| DRomics | 2.5+ | Dose-response modeling |
| R | 4.4+ | Statistical computing |
| Streamlit | 1.50+ | Web GUI |

---

## Quick Reference

```bash
# Setup (first time only)
bash setup.sh

# Launch app
bash run_app.sh

# Stop app
kill $(cat .app.pid)

# Run from command line
conda activate rnaseq_pipeline
nextflow run main.nf --metadata metadata.xlsx

# Resume interrupted run
nextflow run main.nf --metadata metadata.xlsx -resume

# Clean work directory
rm -rf ~/aracra_star_work/work/*

# Clean STAR shared memory
STAR --genomeDir ~/databases/hg38_reference/star_index_hg38 --genomeLoad Remove
```
