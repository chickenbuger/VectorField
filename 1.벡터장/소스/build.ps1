param(
    [ValidateSet('Debug', 'Release')]
    [string]$Configuration = 'Release',
    [ValidateSet('x64', 'Win32')]
    [string]$Platform = 'x64'
)

$ErrorActionPreference = 'Stop'

function Invoke-Checked {
    param([string]$Program, [string[]]$ArgumentList)
    & $Program @ArgumentList
    if ($LASTEXITCODE -ne 0) {
        throw "$Program failed with exit code $LASTEXITCODE."
    }
}

$repoRoot = $PSScriptRoot
$freeglutCommit = '566cd459201ae29003ef5686c06506d8646e55c6'
$vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'

Push-Location -LiteralPath $repoRoot
try {
    if (-not (Test-Path -LiteralPath $vswhere)) {
        throw 'Install Visual Studio 2022 with Desktop development with C++ first.'
    }
    $vsRoot = & $vswhere -latest -products '*' -version '[17.0,18.0)' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if ($LASTEXITCODE -ne 0 -or -not $vsRoot) {
        throw 'Visual Studio 2022 C++ tools were not found. Add Desktop development with C++ and a Windows SDK.'
    }
    $msbuild = Join-Path $vsRoot 'MSBuild\Current\Bin\MSBuild.exe'
    $cmake = Join-Path $vsRoot 'Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe'
    if (-not (Test-Path -LiteralPath $cmake)) {
        throw 'Add the C++ CMake tools for Windows component in Visual Studio Installer.'
    }
    if (-not (Test-Path -LiteralPath $msbuild)) {
        throw 'MSBuild was not found in the Visual Studio installation.'
    }
    $git = (Get-Command git -ErrorAction Stop).Source
    $depsRoot = Join-Path $repoRoot '.deps'
    $freeglutSource = Join-Path $depsRoot 'freeglut-src'
    $freeglutBuild = Join-Path $depsRoot "freeglut-build-$Platform"
    $freeglutInstall = Join-Path $depsRoot "freeglut-$Platform"
    New-Item -ItemType Directory -Path $depsRoot -Force | Out-Null

    if (-not (Test-Path -LiteralPath $freeglutSource)) {
        Write-Host 'Downloading freeglut 3.6.0...'
        Invoke-Checked -Program $git -ArgumentList @('clone', '--depth', '1', '--branch', 'v3.6.0', 'https://github.com/freeglut/freeglut.git', $freeglutSource)
    }
    $actualCommit = & $git -C $freeglutSource rev-parse HEAD
    if ($LASTEXITCODE -ne 0 -or $actualCommit -ne $freeglutCommit) {
        throw "Unexpected freeglut source revision. Expected $freeglutCommit in $freeglutSource."
    }

    Write-Host "Building freeglut (Release, $Platform)..."
    Invoke-Checked -Program $cmake -ArgumentList @(
        '-S', $freeglutSource, '-B', $freeglutBuild,
        '-G', 'Visual Studio 17 2022', '-A', $Platform,
        "-DCMAKE_GENERATOR_INSTANCE=$vsRoot",
        "-DCMAKE_INSTALL_PREFIX=$freeglutInstall",
        '-DFREEGLUT_BUILD_DEMOS=OFF', '-DFREEGLUT_BUILD_SHARED_LIBS=ON',
        '-DFREEGLUT_BUILD_STATIC_LIBS=OFF', '-DFREEGLUT_REPLACE_GLUT=OFF'
    )
    Invoke-Checked -Program $cmake -ArgumentList @('--build', $freeglutBuild, '--config', 'Release', '--parallel')
    Invoke-Checked -Program $cmake -ArgumentList @('--install', $freeglutBuild, '--config', 'Release')

    $solutionPlatform = if ($Platform -eq 'Win32') { 'x86' } else { 'x64' }
    Write-Host "Building VectorField ($Configuration, $Platform)..."
    Invoke-Checked -Program $msbuild -ArgumentList @(
        (Join-Path $repoRoot 'opgnel_project_smoke.sln'), '/m', '/nologo', '/verbosity:minimal',
        "/p:Configuration=$Configuration", "/p:Platform=$solutionPlatform"
    )
    $outputDirectory = if ($Platform -eq 'Win32') { $Configuration } else { "x64\$Configuration" }
    $executable = Join-Path $repoRoot "$outputDirectory\opgnel_project_smoke.exe"
    if (-not (Test-Path -LiteralPath $executable)) {
        throw "Build completed without the expected executable: $executable"
    }
    Write-Host "Ready: $executable"
}
catch {
    [Console]::Error.WriteLine($_.Exception.Message)
    exit 1
}
finally {
    Pop-Location
}
