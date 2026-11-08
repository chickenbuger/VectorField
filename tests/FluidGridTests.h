#ifndef FIELDLAB_FLUID_GRID_TESTS_H
#define FIELDLAB_FLUID_GRID_TESTS_H

#include "../opgnel_project_smoke/FluidGrid.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace fieldlab {
namespace fluid_grid_test_detail {

inline bool Near(double a, double b, double tolerance = 1e-9) {
    return std::isfinite(a) && std::abs(a - b) <= tolerance;
}

struct DyeStats {
    double sum = 0.0;
    double maximum = 0.0;
    Vec2 center{};
};

inline DyeStats Measure(const FluidGrid& grid) {
    DyeStats stats;
    for (int j = 1; j <= grid.Size(); ++j) {
        for (int i = 1; i <= grid.Size(); ++i) {
            const double density = grid.Dye(i, j);
            stats.sum += density;
            stats.maximum = (std::max)(stats.maximum, density);
            stats.center.x += density * (i - 0.5) / grid.Size();
            stats.center.y += density * (j - 0.5) / grid.Size();
        }
    }
    if (stats.sum > 0.0) stats.center = stats.center / stats.sum;
    return stats;
}

inline bool FiniteAndBounded(const FluidGrid& grid) {
    for (int j = 1; j <= grid.Size(); ++j) {
        for (int i = 1; i <= grid.Size(); ++i) {
            if (!std::isfinite(grid.Dye(i, j)) || grid.Dye(i, j) < 0.0
                || grid.Dye(i, j) > 10.000001 || !std::isfinite(grid.U(i, j))
                || !std::isfinite(grid.V(i, j)) || std::abs(grid.U(i, j)) > 100.0
                || std::abs(grid.V(i, j)) > 100.0) return false;
        }
    }
    return std::isfinite(grid.DivergenceRms());
}

inline bool ClosedWalls(const FluidGrid& grid) {
    for (int i = 0; i <= grid.Size(); ++i) {
        const double t = double(i) / grid.Size();
        if (!Near(grid.SampleVelocity({ 0.0, t }).x, 0.0)
            || !Near(grid.SampleVelocity({ 1.0, t }).x, 0.0)
            || !Near(grid.SampleVelocity({ t, 0.0 }).y, 0.0)
            || !Near(grid.SampleVelocity({ t, 1.0 }).y, 0.0)) return false;
    }
    return true;
}

} // namespace fluid_grid_test_detail

inline int RunFluidGridTests() {
    using namespace fluid_grid_test_detail;
    int passed = 0, failed = 0;
    const auto check = [&](const std::string& name, bool result) {
        std::cout << (result ? "PASS " : "FAIL ") << name << '\n';
        result ? ++passed : ++failed;
    };

    FluidGrid empty(32);
    empty.Step(1.0 / 60.0);
    check("fluid zero state stays zero", Near(empty.DivergenceRms(), 0.0)
        && Near(Measure(empty).sum, 0.0) && Near(empty.SampleVelocity({ 0.5, 0.5 }).Length(), 0.0));
    check("fluid minimum resolution is safe", FluidGrid(1).Size() == 8);

    FluidGrid uniform(32);
    uniform.FillDye(0.75);
    uniform.SetPresetVelocity(FieldType::Uniform, 0.6);
    uniform.StepDye(0.05, 0.04);
    bool preserved = true;
    for (int j = 1; j <= uniform.Size(); ++j) {
        for (int i = 1; i <= uniform.Size(); ++i) preserved &= Near(uniform.Dye(i, j), 0.75);
    }
    check("constant dye survives diffusion and backtracing", preserved);
    check("uniform velocity interpolates away from walls",
        Near(uniform.SampleVelocity({ 0.34, 0.62 }).x, 0.6)
        && Near(uniform.SampleVelocity({ 0.34, 0.62 }).y, 0.15));
    check("wall normal velocity is zero", ClosedWalls(uniform));

    FluidGrid interpolation(32);
    interpolation.SetPresetVelocity(FieldType::Rotation, 0.7);
    const Vec2 position{ 0.424, 0.583 };
    const Vec2 sampled = interpolation.SampleVelocity(position);
    const Vec2 exact = Evaluate(FieldType::Rotation, position * 4.0 - Vec2{ 2.0, 2.0 }, 0.7);
    check("staggered bilinear velocity reproduces a linear field",
        Near(sampled.x, exact.x) && Near(sampled.y, exact.y));
    check("cell velocity access matches center interpolation",
        Near(interpolation.U(11, 17), interpolation.SampleVelocity({ 10.5 / 32, 16.5 / 32 }).x)
        && Near(interpolation.V(11, 17), interpolation.SampleVelocity({ 10.5 / 32, 16.5 / 32 }).y));
    interpolation.AddDye(0.42, 0.57, 2.0, 0.12);
    const int i = 13, j = 18;
    const double sx = 0.3, sy = 0.7;
    const Vec2 dyePosition{ (i + sx - 0.5) / 32, (j + sy - 0.5) / 32 };
    const double expectedDye = (1.0 - sy) * ((1.0 - sx) * interpolation.Dye(i, j) + sx * interpolation.Dye(i + 1, j))
        + sy * ((1.0 - sx) * interpolation.Dye(i, j + 1) + sx * interpolation.Dye(i + 1, j + 1));
    check("dye sample uses four surrounding cell weights", Near(interpolation.SampleDye(dyePosition), expectedDye));
    check("out of bounds velocity samples clamp safely", std::isfinite(interpolation.SampleVelocity({ -5.0, 9.0 }).Length()));

    FluidGrid transport(64);
    transport.SetPresetVelocity(FieldType::Uniform, 1.0);
    transport.AddDye(0.32, 0.43, 2.0, 0.09);
    const DyeStats start = Measure(transport);
    for (int step = 0; step < 20; ++step) transport.StepDye(0.05, 0.0);
    const DyeStats transported = Measure(transport);
    const Vec2 actualDisplacement = transported.center - start.center;
    check("dye centroid follows uniform flow in physical units",
        Near(actualDisplacement.x, 0.25, 0.012) && Near(actualDisplacement.y, 0.0625, 0.012));
    check("semi Lagrangian advection does not create a higher dye peak",
        transported.maximum <= start.maximum + 1e-10 && FiniteAndBounded(transport));
    std::cout << std::setprecision(10) << "Fluid transport displacement: "
        << actualDisplacement.x << ", " << actualDisplacement.y << " (expected 0.25, 0.0625).\n";

    FluidGrid diffusion(48);
    diffusion.AddDye(0.51, 0.49, 3.0, 0.09);
    const DyeStats beforeDiffusion = Measure(diffusion);
    for (int step = 0; step < 12; ++step) diffusion.StepDye(0.05, 0.04);
    const DyeStats afterDiffusion = Measure(diffusion);
    check("implicit diffusion lowers a compact dye peak", afterDiffusion.maximum < 0.75 * beforeDiffusion.maximum);
    check("closed-wall diffusion approximately preserves dye sum",
        std::abs(afterDiffusion.sum - beforeDiffusion.sum) < beforeDiffusion.sum * 1e-5);
    std::cout << "Fluid diffusion peak: " << beforeDiffusion.maximum << " -> " << afterDiffusion.maximum
        << "; relative dye sum error: " << std::abs(afterDiffusion.sum - beforeDiffusion.sum) / beforeDiffusion.sum << ".\n";

    FluidGrid projection(64);
    projection.AddVelocity(0.46, 0.54, 2.0, -1.3, 0.16);
    projection.Project();
    check("MAC projection reduces divergence by at least 95 percent",
        projection.divergenceBeforeProjection > 0.01
        && projection.divergenceAfterProjection < 0.05 * projection.divergenceBeforeProjection);
    check("projection diagnostic matches the measured final field",
        Near(projection.divergenceAfterProjection, projection.DivergenceRms()));
    check("projection keeps walls impermeable", ClosedWalls(projection));
    std::cout << "Fluid pressure projection RMS divergence: " << projection.divergenceBeforeProjection
        << " -> " << projection.divergenceAfterProjection << "; max: " << projection.maxDivergenceBefore
        << " -> " << projection.maxDivergenceAfter << ".\n";

    FluidGrid compressing(32);
    compressing.SetPresetVelocity(FieldType::Source, 0.8);
    const double sourceBefore = compressing.DivergenceRms();
    compressing.Project();
    check("pressure removes a global source field in a sealed box",
        compressing.DivergenceRms() < sourceBefore * 0.1);
    std::cout << "Fluid source field RMS divergence: " << sourceBefore << " -> " << compressing.DivergenceRms() << ".\n";

    FluidGrid invalid(16);
    const double nan = std::numeric_limits<double>::quiet_NaN();
    invalid.AddDye(nan, 0.5, 1.0);
    invalid.AddVelocity(0.5, 0.5, nan, 1.0);
    invalid.Step(nan);
    invalid.Step(-1.0);
    check("invalid sources and time steps leave finite zero state", FiniteAndBounded(invalid)
        && Near(Measure(invalid).sum, 0.0) && Near(invalid.DivergenceRms(), 0.0));

    FluidGrid longRun(32);
    const auto timerStart = std::chrono::steady_clock::now();
    bool allFinite = true, allClosed = true;
    for (int step = 0; step < 300; ++step) {
        longRun.AddDye(0.5, 0.12, 0.08, 0.06);
        longRun.AddVelocity(0.5, 0.12, 0.02 * std::sin(step * 0.03), 0.04, 0.07);
        longRun.Step(1.0 / 60.0, 0.00015, 0.00002, true, true, true);
        allFinite &= FiniteAndBounded(longRun);
        allClosed &= ClosedWalls(longRun);
    }
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - timerStart).count();
    check("300 buoyant nozzle steps stay finite and density bounded", allFinite && Measure(longRun).sum > 0.0);
    check("300 buoyant nozzle steps keep closed walls", allClosed);
    std::cout << "Fluid 300-step scenario (32x32): " << seconds << " s; final divergence RMS "
        << longRun.DivergenceRms() << "; dye peak " << Measure(longRun).maximum << ".\n";

    longRun.Clear();
    check("clear resets fields and diagnostics", Near(longRun.DivergenceRms(), 0.0)
        && Near(Measure(longRun).sum, 0.0) && Near(longRun.maxDivergenceBefore, 0.0));
    std::cout << "FluidGrid tests: " << passed << " passed, " << failed << " failed, "
        << (passed + failed) << " total.\n";
    return failed == 0 ? 0 : 1;
}

} // namespace fieldlab
#endif
