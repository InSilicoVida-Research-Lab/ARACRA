# ARACRA: RNA-seq Pipeline & Analysis Suite - Windows Installation Guide

This repository contains the ARACRA Streamlit application and an automated setup script to install all necessary dependencies and tools on a Windows 10/11 machine using the Windows Subsystem for Linux (WSL2).

The setup process is designed to be as simple as possible. It will automatically:
1.  Install the Windows Subsystem for Linux (WSL2) and Ubuntu.
2.  Install a Linux version of the Miniconda package manager.
3.  Create a dedicated `aracra` Conda environment with all required bioinformatics tools, Python libraries, and R packages.

---

## Prerequisites

-   Windows 10 (version 2004 or higher) or Windows 11.
-   Administrator privileges on the machine.
-   A stable internet connection for downloading ~5-10 GB of software.

---

## Installation

The installation is a **two-stage process**. You will run the setup script once, reboot your computer, and then run the exact same script a second time.

### Step 1: Download the Project Files

-   Click the green `<> Code` button on this GitHub page.
-   Select **Download ZIP**.
-   Unzip the contents to a simple path on your computer, for example: `C:\ARACRA_Project`.

### Step 2: Run the Automated Setup

1.  Navigate to the project folder (e.g., `C:\ARACRA_Project`).
2.  Right-click the **`setup.bat`** file and choose **"Run as administrator"**.

#### First Run: Installing WSL2
-   The script will detect that WSL2 is not installed and will begin the installation.
-   It will prompt you to **restart your computer**. Press Enter to allow it to restart.
-   After your computer reboots, an **Ubuntu terminal window will open**. It will ask you to create a **username** and a **password**. These are for your new Linux environment. Please complete this step.

#### Second Run: Installing Software
1.  After setting up your Ubuntu user, go back to the project folder.
2.  Once again, right-click the **`setup.bat`** file and choose **"Run as administrator"**.
3.  This time, the script will detect that WSL2 is ready. It will proceed to install Miniconda and all the software for the pipeline. This is a long process (15-40 minutes) that requires no user interaction.
4.  When the script is finished, it will display a green **"SETUP COMPLETE!"** message.

### Step 3: Prepare Your Reference Data Files

The setup installs the *tools*, but you must provide the *data files* for the pipeline to work.

1.  Inside your project folder (`C:\ARACRA_Project`), create a new sub-folder named **`data`**.
2.  **Download and prepare the following files**, placing them in your new `data` folder:
    -   **HISAT2 Genome Index:** Download a pre-built index for your reference genome (e.g., GRCh38 for human) from the [HISAT2 index page](https://daehwankimlab.github.io/hisat2/download/#grch38-and-grcm39). Unzip it into a new folder like `C:\ARACRA_Project\data\grch38_index\`.
    -   **Reference Genome Annotation (GTF):** Download the GTF file for your reference genome from a source like [Ensembl](https://ftp.ensembl.org/pub/current_gtf/homo_sapiens/) or [GENCODE](https://www.gencodegenes.org/human/). Place the file in `C:\ARACRA_Project\data\`.
    -   **Adapter Sequences:** Place your Trimmomatic adapter file (e.g., `TruSeq3-PE.fa`) in `C:\ARACRA_Project\data\`.

---

## How to Run the Application

After the setup is complete, starting the application is simple:

1.  Navigate to the project folder (`C:\ARACRA_Project`).
2.  **Double-click the `launch-app.bat` file.**

A command window will open, and your default web browser will automatically open with the ARACRA application running.

### **IMPORTANT: First-Time Configuration**

When you run the app for the first time, you **must** configure the file paths in the GUI. The application is running inside Linux, so you must provide **Linux-style paths**.

Use the table below as a guide, assuming your project is in `C:\ARACRA_Project`.

| Field in UI         | What to Enter in the UI (Linux Path)                                                |
| ------------------- | ----------------------------------------------------------------------------------- |
| **Base Output Directory** | `/mnt/c/ARACRA_Project/output`                                                      |
| **Trimmomatic JAR Path**  | `/home/YOUR_LINUX_USER/miniconda3/envs/aracra/share/trimmomatic-0.39-2/trimmomatic.jar` |
| **Adapter Sequences**   | `/mnt/c/ARACRA_Project/data/TruSeq3-PE.fa`                                          |
| **HISAT2 Index Path**   | `/mnt/c/ARACRA_Project/data/grch38_index/genome_tran` (or your index prefix)      |
| **Reference GTF**     | `/mnt/c/ARACRA_Project/data/Your_GTF_File.gtf`                                    |

*Replace `YOUR_LINUX_USER` with the Linux username you created during setup.*

---

## File Descriptions

-   `app.py` / `rnaseq_pipeline_core.py`: The core Python source code for the Streamlit application.
-   `setup.bat`: The main script to run for the entire automated installation.
-   `setup.ps1`: The powerful PowerShell script called by `setup.bat` that contains the setup logic.
-   `launch-app.bat`: The simple script to double-click to run the application after setup is complete.