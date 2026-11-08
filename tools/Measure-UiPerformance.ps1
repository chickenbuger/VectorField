param([int]$Frames=120,[int]$Repeats=3)
$ErrorActionPreference='Stop'
$taskRepo=Split-Path -Parent $PSScriptRoot
$taskExe=Join-Path $taskRepo 'x64\Release\opgnel_project_smoke.exe'
$taskOutput=Join-Path $taskRepo '9.정밀계면\최적화\UI측정'
New-Item -ItemType Directory -Path $taskOutput -Force | Out-Null
function Number([string]$Value){return [double]::Parse($Value,[Globalization.CultureInfo]::InvariantCulture)}
function Quantile([double[]]$Values,[double]$Probability){
    $taskSorted=@($Values|Sort-Object)
    return $taskSorted[[Math]::Max(0,[Math]::Ceiling($Probability*$taskSorted.Count)-1)]
}
$taskCases=@(
    @{Name='dam_reference_ui';Args=@('--legacy-frame-loop','--legacy-grid-render')},
    @{Name='dam_optimized_ui';Args=@()},
    @{Name='flat_legacy_draw';Args=@('--hydrostatic','--pressure-view','--legacy-frame-loop','--legacy-grid-render')},
    @{Name='flat_batched_draw';Args=@('--hydrostatic','--pressure-view','--legacy-frame-loop')}
)
$taskSummary=@()
for($taskRound=1;$taskRound -le $Repeats;$taskRound++){
    $taskOrder=@($taskCases)
    if($taskRound%2 -eq 0){[Array]::Reverse($taskOrder)}
    foreach($taskCase in $taskOrder){
        $taskRaw=Join-Path $taskOutput ($taskCase.Name+'_'+$taskRound+'.csv')
        $taskArgs=@('--stage','9','--profile-frames',[string]$Frames,'--profile-log',$taskRaw)+$taskCase.Args
        $taskNativeOutput=& $taskExe @taskArgs
        $taskExit=$LASTEXITCODE
        $taskNativeOutput|Out-File -LiteralPath ($taskRaw+'.stdout.txt') -Encoding utf8
        if($taskExit -ne 0){throw ('UI benchmark failed: '+$taskCase.Name)}
        $taskRows=@(Import-Csv -LiteralPath $taskRaw)
        if($taskRows.Count -ne $Frames){throw 'Unexpected profile frame count'}
        # First 20 frames are warm-up; full-interval simulation throughput is separate.
        $taskMeasured=@($taskRows|Select-Object -Skip 20)
        $taskTimer=@($taskMeasured|ForEach-Object{Number $_.timer_cpu_ms})
        $taskDraw=@($taskMeasured|ForEach-Object{Number $_.draw_cpu_ms})
        $taskIntervals=@($taskMeasured|ForEach-Object{Number $_.frame_interval_ms})
        $taskWall=(@($taskRows|Select-Object -Skip 1|ForEach-Object{Number $_.frame_interval_ms})|Measure-Object -Sum).Sum/1000
        $taskSim=(Number $taskRows[-1].simulation_seconds)-(Number $taskRows[0].simulation_seconds)
        $taskSummary+=[PSCustomObject]@{
            case=$taskCase.Name;trial=$taskRound;frames=$Frames;warmup_frames=20
            timer_p50_ms=(Quantile $taskTimer .5);timer_p95_ms=(Quantile $taskTimer .95);timer_max_ms=($taskTimer|Measure-Object -Maximum).Maximum
            draw_p50_ms=(Quantile $taskDraw .5);draw_p95_ms=(Quantile $taskDraw .95)
            frame_interval_p50_ms=(Quantile $taskIntervals .5);frame_interval_p95_ms=(Quantile $taskIntervals .95)
            measured_wall_seconds=$taskWall;simulation_seconds=$taskSim;simulation_to_wall_ratio=$taskSim/$taskWall
            discarded_wall_seconds=(Number $taskRows[-1].discarded_wall_seconds)-(Number $taskRows[0].discarded_wall_seconds)
        }
        Write-Output ('Measured '+$taskCase.Name+' trial '+$taskRound)
    }
}
$taskSummary|Export-Csv -LiteralPath (Join-Path $taskOutput 'summary.csv') -NoTypeInformation -Encoding utf8
$taskSummary|Format-Table case,trial,timer_p95_ms,draw_p50_ms,frame_interval_p95_ms,simulation_to_wall_ratio -AutoSize
