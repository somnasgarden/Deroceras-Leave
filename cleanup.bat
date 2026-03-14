@echo off
cd /d "%~dp0"
del /F /Q HEAD config hooks objects refs .Rhistory .bash_profile .bashrc .gitconfig .gitmodules .mcp.json .profile .ripgreprc .zprofile .zshrc 2>nul
del /F /Q .idea .vscode 2>nul
rmdir /S /Q results\01_genome results\02_methylome results\03_transcriptomics 2>nul
rmdir /S /Q figures\fig01 figures\fig02 figures\fig03 figures\fig04 figures\fig05 figures\fig06 figures\fig07 figures\fig08 figures\fig09 figures\fig10 2>nul
echo Done.
