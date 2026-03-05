@echo off
echo Minimalist GNU for Windows with Fortran compiler
set PATH=C:\bin\mingw64\bin;%PATH%
rem echo %PATH%
echo working directory:
cd

echo.
mingw32-make Test_HaltaFall

echo.
pause
