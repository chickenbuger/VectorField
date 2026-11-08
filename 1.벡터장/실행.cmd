@echo off
chcp 65001 >nul
pushd "%~dp0"
"%~dp0실행파일\VectorField.exe" --capture-dir "%~dp0캡처" %*
set "appExitCode=%errorlevel%"
popd
if not "%appExitCode%"=="0" pause
exit /b %appExitCode%