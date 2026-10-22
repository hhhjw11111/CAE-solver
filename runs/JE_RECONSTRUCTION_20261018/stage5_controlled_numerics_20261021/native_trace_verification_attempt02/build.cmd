@echo off
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" >nul
if errorlevel 1 exit /b %errorlevel%
set "PATH=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin;%PATH%"
set "LIB=C:\Program Files (x86)\Intel\oneAPI\compiler\latest\lib;%LIB%"
"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin\ifx.exe" @ifx.rsp
exit /b %errorlevel%
