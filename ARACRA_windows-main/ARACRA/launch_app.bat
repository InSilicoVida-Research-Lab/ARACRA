@echo off
ECHO Starting ARACRA Application via WSL2...
ECHO.

REM 1. Get the current directory and convert to lowercase drive
set "WIN_PATH=%~dp0"
set "DRIVE=%WIN_PATH:~0,1%"
if "%DRIVE%"=="C" set "DRIVE=c"
if "%DRIVE%"=="D" set "DRIVE=d"

REM 2. Create the clean Linux path
set "PROC_PATH=%WIN_PATH:~3,-1%"
set "WSL_PATH=/mnt/%DRIVE%/%PROC_PATH:\=/%"

REM 3. The Adaptable Launch:
REM We use 'cd' first, then run everything else inside the safety of the target folder.
REM This avoids the "Unexpected Token" error from your Windows System PATH.
wsl cd "%WSL_PATH%" ^&^& export PATH="$HOME/miniconda3/envs/aracra/bin:$PATH" ^&^& streamlit run rnaseq_pipeline_gui.py

pause