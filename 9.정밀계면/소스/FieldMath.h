#ifndef FIELDLAB_FIELD_MATH_H
#define FIELDLAB_FIELD_MATH_H

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fieldlab {

struct Vec2 {
    double x = 0.0;
    double y = 0.0;

    double Length() const { return std::hypot(x, y); }
    Vec2 operator+(const Vec2& other) const { return { x + other.x, y + other.y }; }
    Vec2 operator-(const Vec2& other) const { return { x - other.x, y - other.y }; }
    Vec2 operator*(double scalar) const { return { x * scalar, y * scalar }; }
    Vec2 operator/(double scalar) const { return { x / scalar, y / scalar }; }
};

inline Vec2 operator*(double scalar, const Vec2& vector) {
    return vector * scalar;
}

// Scaling before normalization also avoids overflow/underflow in x*x + y*y.
inline Vec2 Normalize(Vec2 vector) {
    const double scale = (std::max)(std::abs(vector.x), std::abs(vector.y));
    if (scale == 0.0 || !std::isfinite(vector.x) || !std::isfinite(vector.y)) {
        return {};
    }
    const Vec2 scaled = vector / scale;
    return scaled / scaled.Length();
}

enum class FieldType { Uniform, Rotation, Source, Saddle };

inline Vec2 Evaluate(FieldType type, Vec2 position, double strength = 1.0) {
    switch (type) {
    case FieldType::Uniform:
        return Vec2{ 1.0, 0.25 } * strength;
    case FieldType::Rotation:
        return Vec2{ -position.y, position.x } * strength;
    case FieldType::Source:
        return position * strength;
    case FieldType::Saddle:
        return Vec2{ position.x, -position.y } * strength;
    }
    return {};
}

struct Differential {
    double divergence = 0.0;
    double curl = 0.0;
};

// For F=(u,v), divergence=du/dx+dv/dy and 2D curl=dv/dx-du/dy.
inline Differential Derivatives(FieldType type, Vec2 position,
    double strength = 1.0, double h = 1e-4) {
    if (!(h > 0.0) || !std::isfinite(h)) {
        throw std::invalid_argument("Central-difference step h must be positive and finite.");
    }
    const Vec2 dx = (Evaluate(type, position + Vec2{ h, 0.0 }, strength)
        - Evaluate(type, position - Vec2{ h, 0.0 }, strength)) / (2.0 * h);
    const Vec2 dy = (Evaluate(type, position + Vec2{ 0.0, h }, strength)
        - Evaluate(type, position - Vec2{ 0.0, h }, strength)) / (2.0 * h);
    return { dx.x + dy.y, dx.y - dy.x };
}

} // namespace fieldlab

#endif
