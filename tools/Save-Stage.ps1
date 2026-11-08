param(
    [Parameter(Mandatory=$true)][int]$Number,
    [Parameter(Mandatory=$true)][string]$Folder
)
$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot
$stageRoot = Join-Path $repoRoot $Folder
$checkpointRoot = Join-Path $repoRoot "build\checkpoint-$Number"
New-Item -ItemType Directory -Path "$stageRoot\실행파일","$stageRoot\소스","$stageRoot\캡처", "$checkpointRoot\opgnel_project_smoke", "$checkpointRoot\tests" -Force | Out-Null
Copy-Item -LiteralPath "$repoRoot\x64\Release\opgnel_project_smoke.exe" -Destination "$stageRoot\실행파일\VectorField.exe" -Force
Copy-Item -LiteralPath "$repoRoot\x64\Release\freeglut.dll" -Destination "$stageRoot\실행파일\freeglut.dll" -Force
Copy-Item -LiteralPath "$repoRoot\.deps\freeglut-src\COPYING" -Destination "$stageRoot\실행파일\freeglut-LICENSE.txt" -Force
Get-ChildItem -LiteralPath "$repoRoot\opgnel_project_smoke" -File | Where-Object { $_.Extension -in '.h','.cpp','.vcxproj','.filters' } | ForEach-Object {
    Copy-Item -LiteralPath $_.FullName -Destination "$checkpointRoot\opgnel_project_smoke\$($_.Name)" -Force
}
New-Item -ItemType Directory -Path "$checkpointRoot\opgnel_project_smoke\include" -Force | Out-Null
Copy-Item -Path "$repoRoot\opgnel_project_smoke\include\*" -Destination "$checkpointRoot\opgnel_project_smoke\include" -Recurse -Force
Get-ChildItem -LiteralPath "$repoRoot\tests" -Filter '*.h' -File | ForEach-Object { Copy-Item -LiteralPath $_.FullName -Destination "$checkpointRoot\tests\$($_.Name)" -Force }
foreach ($fileName in @('build.cmd','build.ps1','run.cmd','opgnel_project_smoke.sln')) {
    Copy-Item -LiteralPath "$repoRoot\$fileName" -Destination "$checkpointRoot\$fileName" -Force
}
Copy-Item -LiteralPath "$repoRoot\opgnel_project_smoke\FieldLab.cpp" -Destination "$stageRoot\소스\FieldLab.cpp" -Force
foreach ($header in @('FieldMath.h','DynamicsMath.h','FluidGrid.h','WaterParticles.h','TwoMaterialSimulation.h','SharpInterfaceSimulation.h')) {
    if (Test-Path -LiteralPath "$repoRoot\opgnel_project_smoke\$header") { Copy-Item -LiteralPath "$repoRoot\opgnel_project_smoke\$header" -Destination "$stageRoot\소스\$header" -Force }
}
Compress-Archive -Path "$checkpointRoot\*" -DestinationPath "$stageRoot\소스\전체소스.zip" -Force
$launcher = (Get-Content -LiteralPath (Join-Path $PSScriptRoot 'Stage-Launcher.cmd.template') -Raw -Encoding UTF8).Replace('@STAGE@', [string]$Number)
$launcher = $launcher -replace '\r?\n', "`r`n"
[IO.File]::WriteAllText("$stageRoot\실행.cmd",$launcher,(New-Object Text.UTF8Encoding($false)))
$readme = "이 단계 완성 시점의 소스입니다. 전체소스.zip을 별도 폴더에 풀고 build.cmd를 실행하면 다시 빌드할 수 있습니다. FieldLab.cpp와 계산 헤더는 압축을 풀지 않고도 읽을 수 있게 함께 보관했습니다.`r`n"
[IO.File]::WriteAllText("$stageRoot\소스\읽어주세요.txt",$readme,(New-Object Text.UTF8Encoding($false)))
New-Item -ItemType Directory -Path "$repoRoot\배포" -Force | Out-Null
$archiveName = ('{0:D2}_{1}.zip' -f $Number, ($Folder -replace '^\d+\.', ''))
Compress-Archive -LiteralPath $stageRoot -DestinationPath "$repoRoot\배포\$archiveName" -Force
Write-Output "Saved stage $Number : $stageRoot"
