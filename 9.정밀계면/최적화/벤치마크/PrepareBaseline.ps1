param([string]$RepositoryRoot, [string]$OutputDirectory)
$ErrorActionPreference = 'Stop'
$source = (& git -C $RepositoryRoot show '75cf82b:opgnel_project_smoke/SharpInterfaceSimulation.h') -join "`n"
if ($LASTEXITCODE -ne 0) { throw 'Cannot read immutable baseline commit 75cf82b.' }
$source += "`n"
$encoding = New-Object Text.UTF8Encoding($false)
New-Item -ItemType Directory -Path $OutputDirectory -Force | Out-Null
[IO.File]::WriteAllText((Join-Path $OutputDirectory 'SharpInterfaceBaseline.original.h'), $source, $encoding)
$baseline = $source.Replace('FIELDLAB_SHARP_INTERFACE_SIMULATION_H','FIELDLAB_BASELINE_INTERFACE_SIMULATION_H').Replace('SharpInterfaceSimulation','BaselineSharpInterfaceSimulation')
[IO.File]::WriteAllText((Join-Path $OutputDirectory 'SharpInterfaceBaseline.h'), $baseline, $encoding)
$profile = $source.Replace('FIELDLAB_SHARP_INTERFACE_SIMULATION_H','FIELDLAB_PROFILE_INTERFACE_SIMULATION_H').Replace('SharpInterfaceSimulation','ProfiledSharpInterfaceSimulation')
$profile = $profile.Replace('#include <functional>', "#include <functional>`n#include <chrono>")
$profile = $profile.Replace('    int Size() const', @'
    struct Profile { double snapshots=0, phi=0, reinitialize=0, area=0, coefficients=0, velocity=0, pressure=0, dye=0; std::size_t substeps=0; } profile;
    int Size() const
'@)
$profile = $profile.Replace('        const auto oldU=u_,oldV=v_,oldPhi=phi_,oldDye=dye_;', @'
        using Clock=std::chrono::steady_clock;
        auto phaseStart=Clock::now();
        const auto elapsed=[&](){return std::chrono::duration<double,std::milli>(Clock::now()-phaseStart).count();};
        const auto oldU=u_,oldV=v_,oldPhi=phi_,oldDye=dye_;
        profile.snapshots+=elapsed();phaseStart=Clock::now();++profile.substeps;
'@)
$profile = $profile.Replace('        ScalarBoundary(phi_);Reinitialize();CorrectArea();BuildCoefficients();', @'
        ScalarBoundary(phi_);profile.phi+=elapsed();phaseStart=Clock::now();
        Reinitialize();profile.reinitialize+=elapsed();phaseStart=Clock::now();
        CorrectArea();profile.area+=elapsed();phaseStart=Clock::now();
        BuildCoefficients();profile.coefficients+=elapsed();phaseStart=Clock::now();
'@)
$profile = $profile.Replace('        if(projection_)ProjectOnly(dt);', @'
        profile.velocity+=elapsed();phaseStart=Clock::now();
        if(projection_)ProjectOnly(dt);
'@)
$profile = $profile.Replace('        // Passive gas tracer, with material clipping separate from pressure.', @'
        profile.pressure+=elapsed();phaseStart=Clock::now();
        // Passive gas tracer, with material clipping separate from pressure.
'@)
$profile = $profile.Replace("        ScalarBoundary(dye_);`n    }`n};", "        ScalarBoundary(dye_);profile.dye+=elapsed();`n    }`n};")
if (-not $profile.Contains('profile.dye+=elapsed()')) { throw 'Profiling instrumentation insertion failed.' }
[IO.File]::WriteAllText((Join-Path $OutputDirectory 'SharpInterfaceProfile.h'), $profile, $encoding)
