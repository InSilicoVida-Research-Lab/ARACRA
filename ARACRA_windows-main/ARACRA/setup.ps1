# ARACRA Automated Setup Script for Windows + WSL2
# This script must be run as an Administrator.

# Helper function for colored output to make messages clearer
function Write-Host-Color {
    param(
        [string]$Message,
        [string]$Color
    )
    Write-Host $Message -ForegroundColor $Color
}

Write-Host-Color "=====================================================" "Cyan"
Write-Host-Color "     ARACRA Automated Environment Setup Script       " "Cyan"
Write-Host-Color "=====================================================" "Cyan"
Write-Host ""

# --- Part 0: Check for Administrator Privileges ---
$currentUser = New-Object Security.Principal.WindowsPrincipal $([Security.Principal.WindowsIdentity]::GetCurrent())
if (-Not $currentUser.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)) {
    Write-Host-Color "ERROR: This script must be run as an Administrator." "Red"
    Write-Host-Color "Please right-click the 'setup.bat' file and choose 'Run as administrator'." "Red"
    Read-Host "Press Enter to exit"
    exit
}

# --- Part 1: Check for and Install WSL2 and Ubuntu ---
Write-Host-Color "[STEP 1] Checking for Windows Subsystem for Linux (WSL2) installation..." "Yellow"
$wsl_status = wsl --status 2>&1
if ($LASTEXITCODE -ne 0) {
    Write-Host-Color "WSL is not installed. Beginning installation..." "Green"
    wsl --install -d Ubuntu
    Write-Host-Color "------------------- IMPORTANT ACTION REQUIRED -------------------" "Red"
    Write-Host-Color "1. Your computer will now RESTART to complete the WSL installation." "Red"
    Write-Host-Color "2. After restarting, an Ubuntu window will open. Please set your desired username and password." "Red"
    Write-Host-Color "3. Once that setup is done, please run this 'setup.bat' script AGAIN to continue." "Red"
    Read-Host "Press Enter to restart your computer."
    Restart-Computer -Force
    exit
} else {
    Write-Host-Color "WSL2 is already installed. Proceeding to the next step." "Green"
}


# --- Part 2: This script will now create the CORRECT environment.yml file for you ---
$ProjectRoot = $PSScriptRoot
$EnvFile = Join-Path $ProjectRoot "environment.yml"
$WslProjectRoot = "/mnt/c" + ($ProjectRoot -replace 'C:', '' -replace '\\', '/')

Write-Host-Color "[STEP 2] Creating the Conda environment definition file (environment.yml)..." "Yellow"

@'
name: aracra
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.9
  - pip
  - r-base >=4.2,<4.3
  - sra-tools
  - hisat2
  - samtools
  - subread
  - trimmomatic
  - r-dplyr
  - r-readxl
  - r-ggplot2
  - r-ggrepel
  - r-stringr
  - r-jsonlite
  - bioconductor-deseq2
  - bioconductor-edger
  - bioconductor-sva
  - bioconductor-org.hs.eg.db
  - r-ggfortify
  - r-remotes
  - pip:
    - streamlit
    - pandas
    - numpy
    - psutil
    - openpyxl
'@ | Set-Content -Path $EnvFile -Encoding UTF8

Write-Host-Color "environment.yml file created successfully." "Green"
Write-Host-Color "[STEP 3] Setting up the ARACRA environment inside Ubuntu..." "Yellow"
Write-Host-Color "This will install the main environment and then install DRomics from GitHub." "Yellow"

$wslScriptBlock = @"
    set -e
    echo '--- Updating the 'aracra' conda environment with all dependencies ---'
    ~/miniconda3/bin/conda env update --file ${WslProjectRoot}/environment.yml --prune

    echo '--- Installing DRomics from GitHub using the remotes package ---'
~/miniconda3/bin/conda run -n aracra R -e \"remotes::install_github('lbbe-software/DRomics')\"
    echo '--- The ARACRA environment is now fully configured! ---'
"@

$wslScriptBlock = $wslScriptBlock -replace "`r", ""
wsl -- bash -c $wslScriptBlock

# --- Part 3: Final Confirmation ---
if ($LASTEXITCODE -eq 0) {
    Write-Host-Color "=====================================================" "Cyan"
    Write-Host-Color "    SETUP COMPLETE! ARACRA is ready to launch.     " "Cyan"
    Write-Host-Color "=====================================================" "Cyan"
    Write-Host-Color "You can now close this window." "Green"
    Write-Host-Color "Use the 'launch-app.bat' file to start the application anytime." "Green"
} else {
    Write-Host-Color "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" "Red"
    Write-Host-Color "        An error occurred during the WSL setup.      " "Red"
    Write-Host-Color "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" "Red"
    Write-Host-Color "Please review the messages above to diagnose the issue." "Red"
}

Read-Host "Press Enter to exit"