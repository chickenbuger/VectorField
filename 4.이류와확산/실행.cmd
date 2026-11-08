@echo off
chcp 65001 >nul
if not exist "%~dp0실행파일\VectorField.exe" (
    call "%~dp0..\run.cmd" --stage 4 --capture-dir "%~dp0캡처" %*
    exit /b
)
pushd "%~dp0"
"%~dp0실행파일\VectorField.exe" --stage 4 --capture-dir "%~dp0캡처" %*
set "appExitCode=%errorlevel%"
popd
if not "%appExitCode%"=="0" pause
exit /b %appExitCode%
