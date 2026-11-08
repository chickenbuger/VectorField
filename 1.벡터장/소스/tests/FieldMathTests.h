#ifndef FIELDLAB_FIELD_MATH_TESTS_H
#define FIELDLAB_FIELD_MATH_TESTS_H

#include "../opgnel_project_smoke/FieldMath.h"
#include <iostream>
#include <limits>
#include <string>

namespace fieldlab {
namespace field_math_test_detail {

inline bool Near(double actual, double expected, double tolerance = 1e-9) {
    return std::isfinite(actual) && std::abs(actual - expected) <= tolerance;
}

inline bool NearVector(Vec2 actual, Vec2 expected, double tolerance = 1e-9) {
    return Near(actual.x, expected.x, tolerance) && Near(actual.y, expected.y, tolerance);
}

} // namespace field_math_test_detail

inline int RunFieldMathTests() {
    using namespace field_math_test_detail;
    int passed = 0;
    int failed = 0;
    const auto check = [&](const std::string& name, bool result) {
        std::cout << (result ? "PASS " : "FAIL ") << name << '\n';
        if (result) { ++passed; }
        else { ++failed; }
    };

    check("default vector is zero", NearVector(Vec2{}, { 0.0, 0.0 }));
    check("3-4 vector magnitude is 5", Near(Vec2{ 3.0, 4.0 }.Length(), 5.0));
    check("vector addition", NearVector(Vec2{ 2.0, -3.0 } + Vec2{ -1.0, 5.0 }, { 1.0, 2.0 }));
    check("vector subtraction", NearVector(Vec2{ 2.0, -3.0 } - Vec2{ -1.0, 5.0 }, { 3.0, -8.0 }));
    check("negative scalar multiplication", NearVector(Vec2{ 2.0, -3.0 } * -2.0, { -4.0, 6.0 }));
    check("left scalar multiplication", NearVector(0.5 * Vec2{ 2.0, -3.0 }, { 1.0, -1.5 }));
    check("scalar division", NearVector(Vec2{ 6.0, -8.0 } / 2.0, { 3.0, -4.0 }));
    check("zero normalization stays finite zero", NearVector(Normalize({}), {}));
    check("3-4 vector normalization", NearVector(Normalize({ 3.0, 4.0 }), { 0.6, 0.8 }));
    check("negative vector normalization", NearVector(Normalize({ -3.0, -4.0 }), { -0.6, -0.8 }));
    check("large vector normalization avoids overflow", NearVector(Normalize({ 3e200, 4e200 }), { 0.6, 0.8 }));
    check("tiny vector normalization avoids underflow", NearVector(Normalize({ 3e-200, 4e-200 }), { 0.6, 0.8 }));
    check("nonfinite normalization has finite fallback", NearVector(Normalize({ std::numeric_limits<double>::infinity(), 1.0 }), {}));

    check("uniform field is position independent",
        NearVector(Evaluate(FieldType::Uniform, { -9.0, 7.0 }), { 1.0, 0.25 }));
    check("rotation at (2,-3)", NearVector(Evaluate(FieldType::Rotation, { 2.0, -3.0 }), { 3.0, 2.0 }));
    check("source at (2,-3)", NearVector(Evaluate(FieldType::Source, { 2.0, -3.0 }), { 2.0, -3.0 }));
    check("saddle at (2,-3)", NearVector(Evaluate(FieldType::Saddle, { 2.0, -3.0 }), { 2.0, 3.0 }));
    check("rotation origin is zero", NearVector(Evaluate(FieldType::Rotation, {}), {}));
    check("source origin is zero", NearVector(Evaluate(FieldType::Source, {}), {}));
    check("saddle origin is zero", NearVector(Evaluate(FieldType::Saddle, {}), {}));
    check("rotation velocity is perpendicular to radius", []() {
        const Vec2 p{ 1.25, -0.75 };
        const Vec2 velocity = Evaluate(FieldType::Rotation, p);
        return Near(p.x * velocity.x + p.y * velocity.y, 0.0);
    }());
    check("source velocity points outward", []() {
        const Vec2 p{ -1.25, 0.75 };
        const Vec2 velocity = Evaluate(FieldType::Source, p);
        return p.x * velocity.x + p.y * velocity.y > 0.0;
    }());

    const FieldType types[] = { FieldType::Uniform, FieldType::Rotation, FieldType::Source, FieldType::Saddle };
    const char* names[] = { "uniform", "rotation", "source", "saddle" };
    const double expectedDivergence[] = { 0.0, 0.0, 2.0, 0.0 };
    const double expectedCurl[] = { 0.0, 2.0, 0.0, 0.0 };
    const double steps[] = { 1e-2, 1e-4, 1e-6 };
    for (int index = 0; index < 4; ++index) {
        const std::string label = names[index];
        const Vec2 point{ 0.7, -1.2 };
        check(label + " strength scales linearly",
            NearVector(Evaluate(types[index], point, 2.5), Evaluate(types[index], point) * 2.5));
        check(label + " negative strength reverses field",
            NearVector(Evaluate(types[index], point, -1.0), Evaluate(types[index], point) * -1.0));
        check(label + " zero strength gives zero vector", NearVector(Evaluate(types[index], point, 0.0), {}));
        for (int stepIndex = 0; stepIndex < 3; ++stepIndex) {
            const Differential derivative = Derivatives(types[index], point, 1.0, steps[stepIndex]);
            check(label + " divergence at h=" + std::to_string(steps[stepIndex]),
                Near(derivative.divergence, expectedDivergence[index], 1e-7));
            check(label + " curl at h=" + std::to_string(steps[stepIndex]),
                Near(derivative.curl, expectedCurl[index], 1e-7));
        }
        const Differential atOrigin = Derivatives(types[index], {});
        check(label + " derivatives at origin",
            Near(atOrigin.divergence, expectedDivergence[index]) && Near(atOrigin.curl, expectedCurl[index]));
        const Differential scaled = Derivatives(types[index], point, -2.0);
        check(label + " derivatives scale with strength",
            Near(scaled.divergence, -2.0 * expectedDivergence[index]) && Near(scaled.curl, -2.0 * expectedCurl[index]));
    }

    bool rejectedZeroStep = false;
    try { (void)Derivatives(FieldType::Source, {}, 1.0, 0.0); }
    catch (const std::invalid_argument&) { rejectedZeroStep = true; }
    check("central differences reject zero step", rejectedZeroStep);
    bool rejectedNegativeStep = false;
    try { (void)Derivatives(FieldType::Source, {}, 1.0, -1e-4); }
    catch (const std::invalid_argument&) { rejectedNegativeStep = true; }
    check("central differences reject negative step", rejectedNegativeStep);
    bool rejectedNonfiniteStep = false;
    try { (void)Derivatives(FieldType::Source, {}, 1.0, std::numeric_limits<double>::quiet_NaN()); }
    catch (const std::invalid_argument&) { rejectedNonfiniteStep = true; }
    check("central differences reject nonfinite step", rejectedNonfiniteStep);

    std::cout << "FieldMath tests: " << passed << " passed, " << failed << " failed, "
        << (passed + failed) << " total.\n";
    return failed == 0 ? 0 : 1;
}

} // namespace fieldlab

#endif
