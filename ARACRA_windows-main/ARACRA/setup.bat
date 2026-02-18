@echo off
ECHO #######################################################
ECHO # ARACRA Automated Setup Launcher                     #
ECHO # This will request Administrator permissions.        #
ECHO #######################################################
ECHO.
PowerShell -NoProfile -ExecutionPolicy Bypass -Command "& {Start-Process PowerShell -ArgumentList '-NoProfile -ExecutionPolicy Bypass -File ""%~dp0setup.ps1""' -Verb RunAs}"
pause