# ARACRA: Automated RNA-seq and Comprehensive 'omics Analysis

ARACRA is a Streamlit-based graphical user interface that integrates a complete RNA-seq processing pipeline with downstream differential expression (DEG) and dose-response (DRomics) analysis.

## Features

-   **End-to-End RNA-seq Pipeline**: Fetches data from SRA, performs QC, trimming, alignment, and generates a count matrix.
-   **Differential Expression**: Implements robust preprocessing and uses DESeq2 for robust DEG analysis.
-   **Dose-Response Modeling**: Utilizes the DRomics R package to model dose-dependent effects and calculate Benchmark Doses (BMD).
-   **Integrated GUI**: A user-friendly interface to configure and run all analyses, view results, and download plots/tables.

## Setup and Installation (Linux)

This project uses Conda to manage its complex dependencies across Python, R, and command-line bioinformatics tools.

### Prerequisites

1.  **Linux Operating System**: Tested on Ubuntu 20.04/22.04.
2.  **Git**: To clone the repository (`sudo apt-get install git`).
3.  **Conda**: You must have Anaconda or Miniconda installed. We highly recommend installing **Mambaforge**, which includes the faster `mamba` package manager.
    -   [Install Mambaforge (Recommended)](https://github.com/conda-forge/mambaforge#installation)
    -   [Install Miniconda](https://docs.conda.io/en/latest/miniconda.html)

### Installation Steps

1.  **Clone the Repository**
    ```bash
    git clone https://github.com/your-username/your-repo-name.git
    cd your-repo-name
    ```

2.  **Run the Installation Script**
    The `install.sh` script will create a self-contained Conda environment named `aracra_env` with all required software.

    ```bash
    chmod +x install.sh
    ./install.sh
    ```
    This process can take several minutes.

### Running the Application

1.  **Activate the Conda Environment**
    You must activate the environment in every new terminal session before running the app.
    ```bash
    conda activate aracra_env
    ```

2.  **Launch the Streamlit App**
    ```bash
    streamlit run rnaseq_pipeline_gui.py
    ```
    Your web browser should open with the ARACRA interface.

### First-Time Configuration (Crucial!)

Before running the RNA-seq pipeline for the first time, you need to download some reference files and configure their paths in the application's sidebar.

1.  **Trimmomatic JAR Path**:
    The installation script installs Trimmomatic, but the GUI needs the explicit path to its `.jar` file. Find it by running:
    ```bash
    # First, activate the environment
    conda activate aracra_env

    # Then, find the path
    find $CONDA_PREFIX -name "trimmomatic*jar"
    ```
    Copy the output path (e.g., `/path/to/mambaforge/envs/aracra_env/share/trimmomatic-0.39-2/trimmomatic.jar`) and paste it into the "Trimmomatic JAR Path" field in the GUI sidebar.

2.  **Adapter Sequences**:
    You need a FASTA file containing adapter sequences for Trimmomatic. If you don't have one, you can download common Illumina adapters.

3.  **HISAT2 Genome Index**:
    Download a pre-built HISAT2 index for your reference genome (e.g., human GRCh38).
    -   Go to the [HISAT2 Index Download Page](https://daehwankimlab.github.io/hisat2/download/#grch38-and-grcm38).
    -   Download the index archive (e.g., `genome_tran.tar.gz`).
    -   Extract it: `tar -zxvf genome_tran.tar.gz`.
    -   In the GUI's "HISAT2 Index Path" field, provide the path to the index basename (e.g., `/path/to/your/indices/grch38/genome_tran`).

4.  **Reference GTF File**:
    Download the corresponding gene annotation file (GTF format) for your reference genome from a source like GENCODE, Ensembl, or UCSC.

Once these paths are correctly set in the sidebar, you are ready to run the pipeline!
