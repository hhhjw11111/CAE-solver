@echo off
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
cd /d "C:\Users\Administrator\Desktop\gemini??\Console3\solver_v2"
powershell -ExecutionPolicy Bypass -File .\scripts\build_debug.ps1
