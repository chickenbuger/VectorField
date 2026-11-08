#ifndef FIELDLAB_DYNAMICS_TESTS_H
#define FIELDLAB_DYNAMICS_TESTS_H

#include "../opgnel_project_smoke/DynamicsMath.h"
#include <iostream>
#include <limits>
#include <string>

namespace fieldlab {
namespace dynamics_test_detail {

inline bool Near(double actual, double expected, double tolerance = 1e-10) {
    return std::isfinite(actual) && std::abs(actual - expected) <= tolerance;
}

inline bool NearVector(Vec2 actual, Vec2 expected, double tolerance = 1e-10) {
    return Near(actual.x, expected.x, tolerance) && Near(actual.y, expected.y, tolerance);
}

} // namespace dynamics_test_detail

inline int RunDynamicsTests() {
    using namespace dynamics_test_detail;
    int passed = 0;
    int failed = 0;
    const auto check = [&](const std::string& name, bool result) {
        std::cout << (result ? "PASS " : "FAIL ") << name << '\n';
        if (result) { ++passed; }
        else { ++failed; }
    };

    const std::vector<GravitySource> one{ { { 0.0, 0.0 }, 1.0 } };
    const std::vector<GravitySource> zeroMass{ { { 0.0, 0.0 }, 0.0 } };
    check("gravity with no sources is zero", NearVector(GravityField({ 0.7, -0.2 }, {}), {}));
    check("softened gravity at source center is finite zero", NearVector(GravityField({}, one), {}));
    check("center remains finite with very small epsilon", NearVector(GravityField({}, one, 0.4, 1e-200), {}));
    check("zero source mass produces zero acceleration", NearVector(GravityField({ 1.0, 2.0 }, zeroMass), {}));
    check("zero gravitational constant produces zero acceleration", NearVector(GravityField({ 1.0, 2.0 }, one, 0.0), {}));
    const Vec2 right = GravityField({ 1.0, 0.0 }, one);
    const Vec2 left = GravityField({ -1.0, 0.0 }, one);
    check("gravity points toward the source", right.x < 0.0 && Near(right.y, 0.0));
    check("single source has reflection symmetry", NearVector(left, right * -1.0));
    check("single source is rotationally symmetric",
        NearVector(GravityField({ 0.0, 1.0 }, one), { right.y, right.x }));
    check("gravity magnitude decreases with distance",
        GravityField({ 2.0, 0.0 }, one).Length() < right.Length());
    check("inverse-square falloff approaches ratio four",
        Near(GravityField({ 1.0, 0.0 }, one, 1.0, 1e-6).Length()
            / GravityField({ 2.0, 0.0 }, one, 1.0, 1e-6).Length(), 4.0, 1e-9));
    check("source mass scales acceleration linearly",
        NearVector(GravityField({ 0.7, -0.4 }, { { {}, 3.0 } }), GravityField({ 0.7, -0.4 }, one) * 3.0));
    check("G scales acceleration linearly",
        NearVector(GravityField({ 0.7, -0.4 }, one, 1.2), GravityField({ 0.7, -0.4 }, one) * 3.0));
    check("negative G reverses acceleration",
        NearVector(GravityField({ 0.7, -0.4 }, one, -0.4), GravityField({ 0.7, -0.4 }, one) * -1.0));
    check("coincident duplicate sources equal doubled mass",
        NearVector(GravityField({ 0.8, 0.2 }, { { {}, 1.0 }, { {}, 1.0 } }),
            GravityField({ 0.8, 0.2 }, { { {}, 2.0 } })));
    const std::vector<GravitySource> pair{ { { -1.0, 0.0 }, 2.0 }, { { 1.0, 0.0 }, 2.0 } };
    check("equal sources cancel at the midpoint", NearVector(GravityField({}, pair), {}));
    check("equal sources cancel horizontal acceleration above midpoint", Near(GravityField({ 0.0, 1.0 }, pair).x, 0.0));
    const Vec2 probe{ 0.3, -0.6 };
    const Vec2 separate = GravityField(probe, { pair[0] }) + GravityField(probe, { pair[1] });
    check("gravity obeys vector superposition", NearVector(GravityField(probe, pair), separate));
    check("source order does not change the result",
        NearVector(GravityField(probe, pair), GravityField(probe, { pair[1], pair[0] })));
    const std::vector<GravitySource> several{
        { { -1.0, 0.4 }, 2.0 }, { { 0.8, -0.3 }, 0.7 }, { { 0.2, 1.0 }, 3.1 }, { { -0.5, -0.8 }, 1.6 }
    };
    std::vector<GravitySource> reversed = several;
    std::reverse(reversed.begin(), reversed.end());
    check("multiple source ordering changes only negligible roundoff",
        NearVector(GravityField(probe, several), GravityField(probe, reversed), 1e-12));
    check("repeated gravity evaluation has no accumulated state",
        NearVector(GravityField(probe, several), GravityField(probe, several)));
    const Vec2 shift{ 4.0, -7.0 };
    std::vector<GravitySource> shifted = several;
    for (GravitySource& source : shifted) source.position = source.position + shift;
    check("gravity is invariant under translating all positions",
        NearVector(GravityField(probe, several), GravityField(probe + shift, shifted), 1e-12));
    check("larger softening reduces near-source acceleration",
        GravityField({ 0.05, 0.0 }, one, 0.4, 0.2).Length()
        < GravityField({ 0.05, 0.0 }, one, 0.4, 0.1).Length());
    check("coincident source does not hide other sources",
        NearVector(GravityField({}, { { {}, 5.0 }, { { 1.0, 0.0 }, 1.0 } }),
            GravityField({}, { { { 1.0, 0.0 }, 1.0 } })));

    const auto uniform = [](Vec2) { return Vec2{ 1.0, 0.25 }; };
    check("RK2 integrates uniform velocity exactly",
        NearVector(AdvectRK2({ 0.5, -1.0 }, 2.0, uniform), { 2.5, -0.5 }));
    check("RK2 field overload matches callable overload",
        NearVector(AdvectRK2(probe, 0.04, FieldType::Rotation, 1.5),
            AdvectRK2(probe, 0.04, [](Vec2 p) { return Evaluate(FieldType::Rotation, p, 1.5); })));
    check("RK2 traces backward in a uniform field",
        NearVector(AdvectRK2({ 0.5, -1.0 }, -2.0, uniform), { -1.5, -1.5 }));
    check("RK2 zero velocity preserves position",
        NearVector(AdvectRK2(probe, 10.0, [](Vec2) { return Vec2{}; }), probe));
    int samples = 0;
    const Vec2 zeroStep = AdvectRK2(probe, 0.0, [&](Vec2) { ++samples; return Vec2{ 1.0, 1.0 }; });
    check("RK2 zero dt preserves position", NearVector(zeroStep, probe));
    check("RK2 zero dt performs no velocity samples", samples == 0);
    samples = 0;
    (void)AdvectRK2(probe, 0.01, [&](Vec2 p) { ++samples; return Evaluate(FieldType::Rotation, p); });
    check("RK2 takes two velocity samples per step", samples == 2);
    check("RK2 fixed point at rotational origin stays zero",
        NearVector(AdvectRK2({}, 1.0, FieldType::Rotation), {}));

    Vec2 rk{ 1.0, 0.0 };
    Vec2 euler = rk;
    for (int step = 0; step < 1000; ++step) {
        rk = AdvectRK2(rk, 0.01, FieldType::Rotation);
        euler = euler + Evaluate(FieldType::Rotation, euler) * 0.01;
    }
    const Vec2 analytic{ std::cos(10.0), std::sin(10.0) };
    std::cout << "Rotation diagnostics (1000 x dt=0.01): RK2 radius error="
        << std::abs(rk.Length() - 1.0) << ", RK2 position error=" << (rk - analytic).Length()
        << ", Euler radius error=" << std::abs(euler.Length() - 1.0)
        << ", Euler position error=" << (euler - analytic).Length() << '\n';
    check("RK2 rotation radius drift after 1000 steps is below 2e-6", std::abs(rk.Length() - 1.0) < 2e-6);
    check("RK2 rotation position matches analytic orbit within 2e-4", (rk - analytic).Length() < 2e-4);
    check("RK2 rotation radius improves on Euler by over 1000x",
        std::abs(rk.Length() - 1.0) * 1000.0 < std::abs(euler.Length() - 1.0));
    check("RK2 rotation position improves on Euler by over 100x",
        (rk - analytic).Length() * 100.0 < (euler - analytic).Length());
    check("RK2 zero strength preserves position", NearVector(AdvectRK2(probe, 1.0, FieldType::Source, 0.0), probe));
    check("RK2 short source-field step matches exponential growth",
        NearVector(AdvectRK2({ 1.0, -0.5 }, 0.01, FieldType::Source),
            Vec2{ 1.0, -0.5 } * std::exp(0.01), 2e-7));
    const auto sourceError = [](double dt, int steps) {
        Vec2 p{ 1.0, 0.0 };
        for (int step = 0; step < steps; ++step) p = AdvectRK2(p, dt, FieldType::Source);
        return std::abs(p.x - std::exp(1.0));
    };
    check("halving RK2 step reduces source-field error by at least 3.5x",
        sourceError(0.05, 20) * 3.5 < sourceError(0.1, 10));
    const Vec2 byGravity = AdvectRK2(probe, 0.001, [&](Vec2 p) { return GravityField(p, one); });
    check("RK2 accepts a softened gravity-driven velocity callback", dynamics_detail::Finite(byGravity));

    bool invalidSoftening = false;
    try { (void)GravityField({}, one, 0.4, 0.0); }
    catch (const std::invalid_argument&) { invalidSoftening = true; }
    check("gravity rejects zero softening", invalidSoftening);
    bool invalidStep = false;
    try { (void)AdvectRK2(probe, std::numeric_limits<double>::quiet_NaN(), uniform); }
    catch (const std::invalid_argument&) { invalidStep = true; }
    check("RK2 rejects nonfinite dt", invalidStep);
    bool invalidVelocity = false;
    try { (void)AdvectRK2(probe, 0.1, [](Vec2) { return Vec2{ std::numeric_limits<double>::infinity(), 0.0 }; }); }
    catch (const std::overflow_error&) { invalidVelocity = true; }
    check("RK2 rejects nonfinite sampled velocity", invalidVelocity);

    std::cout << "Dynamics tests: " << passed << " passed, " << failed << " failed, "
        << (passed + failed) << " total.\n";
    return failed == 0 ? 0 : 1;
}

} // namespace fieldlab

#endif
