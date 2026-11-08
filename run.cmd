@echo off
if not exist "%~dp0x64\Release\opgnel_project_smoke.exe" (
    call "%~dp0build.cmd"
    if errorlevel 1 (
        pause
        exit /b 1
    )
)
pushd "%~dp0x64\Release"
"%~dp0x64\Release\opgnel_project_smoke.exe" %*
set "appExitCode=%errorlevel%"
popd
if not "%appExitCode%"=="0" pause
exit /b %appExitCode%
