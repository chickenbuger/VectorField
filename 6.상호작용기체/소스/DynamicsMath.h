#ifndef FIELDLAB_DYNAMICS_MATH_H
#define FIELDLAB_DYNAMICS_MATH_H

#include "FieldMath.h"
#include <vector>

namespace fieldlab {

struct GravitySource {
    Vec2 position;
    double mass = 1.0;
};

namespace dynamics_detail {

inline bool Finite(Vec2 value) {
    return std::isfinite(value.x) && std::isfinite(value.y);
}

// Neumaier compensated summation reduces roundoff when sources cancel.
struct Sum {
    double value = 0.0;
    double correction = 0.0;
    void Add(double next) {
        const double combined = value + next;
        correction += std::abs(value) >= std::abs(next)
            ? (value - combined) + next : (next - combined) + value;
        value = combined;
    }
    double Result() const { return value + correction; }
};

} // namespace dynamics_detail

// Acceleration, not force: the probe particle's mass is not part of the sum.
// r points FROM the probe TO each source. Softening makes the center finite.
inline Vec2 GravityField(Vec2 position, const std::vector<GravitySource>& sources,
    double gravitationalConstant = 0.4, double epsilon = 0.12) {
    if (!dynamics_detail::Finite(position) || !std::isfinite(gravitationalConstant)
        || !(epsilon > 0.0) || !std::isfinite(epsilon)) {
        throw std::invalid_argument("Gravity requires finite coordinates/G and positive finite epsilon.");
    }
    dynamics_detail::Sum sumX;
    dynamics_detail::Sum sumY;
    for (const GravitySource& source : sources) {
        if (!dynamics_detail::Finite(source.position) || !std::isfinite(source.mass)) {
            throw std::invalid_argument("Gravity source coordinates and mass must be finite.");
        }
        const Vec2 r = source.position - position;
        // A coincident softened source contributes exactly zero acceleration.
        if ((r.x == 0.0 && r.y == 0.0) || source.mass == 0.0 || gravitationalConstant == 0.0) {
            continue;
        }
        const double radius = std::hypot(r.Length(), epsilon);
        const double scale = (gravitationalConstant * source.mass / radius) / radius;
        const Vec2 contribution = (r / radius) * scale;
        if (!dynamics_detail::Finite(contribution)) {
            throw std::overflow_error("Gravity acceleration exceeds the finite numeric range.");
        }
        sumX.Add(contribution.x);
        sumY.Add(contribution.y);
    }
    const Vec2 result{ sumX.Result(), sumY.Result() };
    if (!dynamics_detail::Finite(result)) {
        throw std::overflow_error("The combined gravity acceleration exceeds the finite numeric range.");
    }
    return result;
}

// Midpoint Runge-Kutta integration of dx/dt = velocity(x).
// Negative dt is allowed for backward tracing; zero dt performs no sampling.
template <typename Velocity>
inline Vec2 AdvectRK2(Vec2 position, double dt, Velocity velocity) {
    if (!dynamics_detail::Finite(position) || !std::isfinite(dt)) {
        throw std::invalid_argument("RK2 requires a finite position and time step.");
    }
    if (dt == 0.0) return position;
    const Vec2 first = velocity(position);
    if (!dynamics_detail::Finite(first)) {
        throw std::overflow_error("RK2 sampled a nonfinite velocity.");
    }
    const Vec2 midpoint = position + first * (0.5 * dt);
    if (!dynamics_detail::Finite(midpoint)) {
        throw std::overflow_error("RK2 midpoint exceeds the finite numeric range.");
    }
    const Vec2 second = velocity(midpoint);
    const Vec2 result = position + second * dt;
    if (!dynamics_detail::Finite(second) || !dynamics_detail::Finite(result)) {
        throw std::overflow_error("RK2 result exceeds the finite numeric range.");
    }
    return result;
}

inline Vec2 AdvectRK2(Vec2 position, double dt, FieldType type, double strength = 1.0) {
    return AdvectRK2(position, dt, [type, strength](Vec2 p) {
        return Evaluate(type, p, strength);
    });
}

} // namespace fieldlab

#endif
