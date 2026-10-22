@echo off
setlocal
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" > nul 2>&1
cd /d "%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -Command "& {.\scripts\build_debug.ps1}" 2>&1
