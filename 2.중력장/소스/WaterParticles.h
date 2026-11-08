#ifndef FIELDLAB_WATER_PARTICLES_H
#define FIELDLAB_WATER_PARTICLES_H

#include "FieldMath.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace fieldlab {

struct WaterParticle {
    Vec2 position;
    Vec2 velocity;
};

// A two-dimensional, equal-mass Position Based Fluids solver. Density is
// estimated with a normalized 2D poly6 kernel, and its exact gradient is used
// for the Jacobi constraint solve. Free-surface density deficits are allowed;
// positive compression is projected, rather than pulling sparse particles
// together with negative pressure.
class WaterParticles {
public:
    WaterParticles() {
        gridWidth_ = static_cast<int>(std::ceil(4.0 / smoothingRadius_));
        gridHead_.resize(static_cast<std::size_t>(gridWidth_ * gridWidth_), -1);
        const int reach = static_cast<int>(std::ceil(smoothingRadius_ / spacing_));
        double referenceKernelSum = 0.0;
        for (int y = -reach; y <= reach; ++y) {
            for (int x = -reach; x <= reach; ++x) {
                referenceKernelSum += Kernel(Vec2{ x * spacing_, y * spacing_ });
            }
        }
        // Calibrating against a filled lattice avoids a resolution-dependent
        // error between nominal rest density and the discrete density estimate.
        particleMass_ = nominalRestDensity_ / referenceKernelSum;
        Reset();
    }

    const std::vector<WaterParticle>& Particles() const { return particles_; }
    std::size_t Count() const { return particles_.size(); }
    std::size_t MaximumParticles() const { return maximumParticles_; }
    double ParticleRadius() const { return particleRadius_; }
    double ParticleSpacing() const { return spacing_; }
    double SmoothingRadius() const { return smoothingRadius_; }
    double ParticleMass() const { return particleMass_; }
    double RestDensity() const { return restDensity_; }
    double Viscosity() const { return viscosity_; }
    Vec2 Gravity() const { return gravity_; }
    double FixedTimeStep() const { return fixedTimeStep_; }
    double SimulatedTime() const { return simulatedTime_; }
    int ConstraintIterations() const { return constraintIterations_; }

    // These are mean positive relative density excesses, not absolute errors.
    // A low-density free-surface particle must not be treated as compressed.
    double AverageDensityError() const { return averageDensityError_; }
    double LastPredictedDensityError() const { return lastPredictedDensityError_; }
    double LastSolvedDensityError() const { return lastSolvedDensityError_; }

    void SetViscosity(double viscosity) {
        RequireFinite(viscosity, "Viscosity must be finite.");
        viscosity_ = Clamp(viscosity, 0.0, 2.0);
    }

    void SetRestDensity(double restDensity) {
        RequireFinite(restDensity, "Rest density must be finite.");
        if (!(restDensity > 0.0)) {
            throw std::invalid_argument("Rest density must be positive.");
        }
        // Particle mass stays fixed: changing rest density changes the target
        // occupied area, instead of silently cancelling the material change.
        restDensity_ = Clamp(restDensity, 100.0, 5000.0);
        RefreshDensityMetrics();
    }

    void SetGravity(Vec2 gravity) {
        if (!Finite(gravity)) {
            throw std::invalid_argument("Gravity must be finite.");
        }
        gravity_ = gravity;
    }

    void SetConstraintIterations(int iterations) {
        constraintIterations_ = (std::max)(1, (std::min)(12, iterations));
    }

    void Reset() {
        particles_.clear();
        particles_.reserve(maximumParticles_);
        // A 22 by 22 dam, initially on the left side of the container.
        for (int y = 0; y < 22; ++y) {
            for (int x = 0; x < 22; ++x) {
                particles_.push_back({ { -1.72 + x * spacing_, -1.75 + y * spacing_ }, {} });
            }
        }
        accumulator_ = 0.0;
        simulatedTime_ = 0.0;
        ResizeStorage();
        RefreshDensityMetrics();
    }

    // Explicit state injection supports deterministic validation and saved
    // initial conditions. All input is checked before modifying the solver.
    void SetParticleState(const std::vector<WaterParticle>& particles) {
        if (particles.size() > maximumParticles_) {
            throw std::invalid_argument("Too many water particles.");
        }
        for (const WaterParticle& particle : particles) {
            if (!Finite(particle.position) || !Finite(particle.velocity)) {
                throw std::invalid_argument("Water particle state must be finite.");
            }
        }
        particles_ = particles;
        for (WaterParticle& particle : particles_) {
            particle.position = ConstrainToContainer(particle.position);
            particle.velocity = LimitSpeed(particle.velocity);
        }
        accumulator_ = 0.0;
        simulatedTime_ = 0.0;
        ResizeStorage();
        RefreshDensityMetrics();
    }

    // Elapsed time is accumulated into exact 1/120-second substeps. A stalled
    // interactive frame contributes at most 0.05 s, preventing a catch-up spiral.
    void Step(double dt) {
        if (!(dt > 0.0) || !std::isfinite(dt)) {
            return;
        }
        accumulator_ += (std::min)(dt, 0.05);
        int substeps = 0;
        while (accumulator_ + 1e-12 >= fixedTimeStep_ && substeps < 6) {
            Substep(fixedTimeStep_);
            accumulator_ = (std::max)(0.0, accumulator_ - fixedTimeStep_);
            ++substeps;
        }
    }

    // delta is a world-space pointer displacement, not a force with SI units.
    // A compact radial weight distributes a bounded velocity impulse.
    void Stir(Vec2 position, Vec2 delta) {
        if (!Finite(position) || !Finite(delta)) {
            return;
        }
        const double radius = 0.45;
        const double deltaLength = delta.Length();
        if (deltaLength > 0.25) {
            delta = delta * (0.25 / deltaLength);
        }
        for (WaterParticle& particle : particles_) {
            const double distance = (particle.position - position).Length();
            if (distance < radius) {
                const double weight = 1.0 - distance / radius;
                particle.velocity = LimitSpeed(particle.velocity + delta * (35.0 * weight));
            }
        }
    }

    void AddBlock(Vec2 center) {
        if (!Finite(center)) {
            return;
        }
        const double minimumSeparation = 0.72 * spacing_;
        for (int y = 0; y < 8 && particles_.size() < maximumParticles_; ++y) {
            for (int x = 0; x < 8 && particles_.size() < maximumParticles_; ++x) {
                const Vec2 position = ConstrainToContainer(center
                    + Vec2{ (x - 3.5) * spacing_, (y - 3.5) * spacing_ });
                bool available = true;
                for (const WaterParticle& particle : particles_) {
                    if ((particle.position - position).Length() < minimumSeparation) {
                        available = false;
                        break;
                    }
                }
                if (available) {
                    particles_.push_back({ position, {} });
                }
            }
        }
        ResizeStorage();
        RefreshDensityMetrics();
    }

    bool IsFinite() const {
        if (!std::isfinite(accumulator_) || !std::isfinite(simulatedTime_)
            || !std::isfinite(averageDensityError_)
            || !std::isfinite(lastPredictedDensityError_)
            || !std::isfinite(lastSolvedDensityError_)) {
            return false;
        }
        for (const WaterParticle& particle : particles_) {
            if (!Finite(particle.position) || !Finite(particle.velocity)) {
                return false;
            }
        }
        return true;
    }

private:
    static double Clamp(double value, double low, double high) {
        return (std::max)(low, (std::min)(high, value));
    }

    static bool Finite(Vec2 value) {
        return std::isfinite(value.x) && std::isfinite(value.y);
    }

    static double SquaredLength(Vec2 value) {
        return value.x * value.x + value.y * value.y;
    }

    static void RequireFinite(double value, const char* message) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(message);
        }
    }

    Vec2 ConstrainToContainer(Vec2 position) const {
        const double low = -2.0 + particleRadius_;
        const double high = 2.0 - particleRadius_;
        return { Clamp(position.x, low, high), Clamp(position.y, low, high) };
    }

    Vec2 LimitSpeed(Vec2 velocity) const {
        const double speed = velocity.Length();
        return speed > maximumSpeed_ ? velocity * (maximumSpeed_ / speed) : velocity;
    }

    double Kernel(Vec2 displacement) const {
        const double h2 = smoothingRadius_ * smoothingRadius_;
        const double q = h2 - SquaredLength(displacement);
        if (!(q > 0.0)) {
            return 0.0;
        }
        const double h4 = h2 * h2;
        const double coefficient = 4.0 / (pi_ * h4 * h4);
        return coefficient * q * q * q;
    }

    Vec2 KernelGradient(Vec2 displacement) const {
        const double h2 = smoothingRadius_ * smoothingRadius_;
        const double q = h2 - SquaredLength(displacement);
        if (!(q > 0.0)) {
            return {};
        }
        const double h4 = h2 * h2;
        return displacement * (-24.0 * q * q / (pi_ * h4 * h4));
    }

    int CellCoordinate(double value) const {
        const int coordinate = static_cast<int>(std::floor((value + 2.0) / smoothingRadius_));
        return (std::max)(0, (std::min)(gridWidth_ - 1, coordinate));
    }

    void ResizeStorage() {
        const std::size_t size = particles_.size();
        oldPositions_.resize(size);
        predicted_.resize(size);
        candidate_.resize(size);
        corrections_.resize(size);
        velocityBuffer_.resize(size);
        densities_.resize(size);
        multipliers_.resize(size);
        gridNext_.resize(size);
        neighbors_.resize(size);
        for (std::vector<int>& neighbors : neighbors_) {
            if (neighbors.capacity() < 48) {
                neighbors.reserve(48);
            }
        }
    }

    void BuildNeighbors(const std::vector<Vec2>& positions) {
        std::fill(gridHead_.begin(), gridHead_.end(), -1);
        for (std::size_t i = 0; i < positions.size(); ++i) {
            const int cell = CellCoordinate(positions[i].x)
                + gridWidth_ * CellCoordinate(positions[i].y);
            gridNext_[i] = gridHead_[static_cast<std::size_t>(cell)];
            gridHead_[static_cast<std::size_t>(cell)] = static_cast<int>(i);
        }
        const double radiusSquared = smoothingRadius_ * smoothingRadius_;
        for (std::size_t i = 0; i < positions.size(); ++i) {
            std::vector<int>& neighbors = neighbors_[i];
            neighbors.clear();
            const int cx = CellCoordinate(positions[i].x);
            const int cy = CellCoordinate(positions[i].y);
            for (int y = (std::max)(0, cy - 1); y <= (std::min)(gridWidth_ - 1, cy + 1); ++y) {
                for (int x = (std::max)(0, cx - 1); x <= (std::min)(gridWidth_ - 1, cx + 1); ++x) {
                    int j = gridHead_[static_cast<std::size_t>(x + gridWidth_ * y)];
                    while (j >= 0) {
                        if (SquaredLength(positions[i] - positions[static_cast<std::size_t>(j)]) < radiusSquared) {
                            neighbors.push_back(j);
                        }
                        j = gridNext_[static_cast<std::size_t>(j)];
                    }
                }
            }
        }
    }

    double ComputeDensities(const std::vector<Vec2>& positions) {
        BuildNeighbors(positions);
        double error = 0.0;
        for (std::size_t i = 0; i < positions.size(); ++i) {
            double density = 0.0;
            for (int index : neighbors_[i]) {
                density += particleMass_ * Kernel(positions[i] - positions[static_cast<std::size_t>(index)]);
            }
            densities_[i] = density;
            error += (std::max)(0.0, density / restDensity_ - 1.0);
        }
        return positions.empty() ? 0.0 : error / static_cast<double>(positions.size());
    }

    void RefreshDensityMetrics() {
        for (std::size_t i = 0; i < particles_.size(); ++i) {
            predicted_[i] = particles_[i].position;
        }
        averageDensityError_ = ComputeDensities(predicted_);
        lastPredictedDensityError_ = averageDensityError_;
        lastSolvedDensityError_ = averageDensityError_;
    }

    Vec2 SeparationDirection(std::size_t i, std::size_t j) const {
        const std::size_t first = (std::min)(i, j);
        const std::size_t second = (std::max)(i, j);
        const double angle = static_cast<double>((first * 73856093u) ^ (second * 19349663u)) * 0.0001;
        const Vec2 direction{ std::cos(angle), std::sin(angle) };
        return i < j ? direction : direction * -1.0;
    }

    void ComputeMultipliers() {
        const double massOverDensity = particleMass_ / restDensity_;
        for (std::size_t i = 0; i < particles_.size(); ++i) {
            const double constraint = (std::max)(0.0, densities_[i] / restDensity_ - 1.0);
            Vec2 ownGradient{};
            double gradientNormSquared = 0.0;
            for (int index : neighbors_[i]) {
                const std::size_t j = static_cast<std::size_t>(index);
                if (j == i) {
                    continue;
                }
                const Vec2 gradient = KernelGradient(predicted_[i] - predicted_[j]) * massOverDensity;
                ownGradient = ownGradient + gradient;
                gradientNormSquared += SquaredLength(gradient);
            }
            gradientNormSquared += SquaredLength(ownGradient);
            multipliers_[i] = -constraint / (gradientNormSquared + 1e-6);
        }
    }

    double ComputeCorrections() {
        const double massOverDensity = particleMass_ / restDensity_;
        const double contactDistance = 0.72 * spacing_;
        const double maximumCorrection = 0.18 * smoothingRadius_;
        double largestCorrection = 0.0;
        for (std::size_t i = 0; i < particles_.size(); ++i) {
            Vec2 correction{};
            for (int index : neighbors_[i]) {
                const std::size_t j = static_cast<std::size_t>(index);
                if (j == i) {
                    continue;
                }
                const Vec2 displacement = predicted_[i] - predicted_[j];
                correction = correction + KernelGradient(displacement)
                    * (massOverDensity * (multipliers_[i] + multipliers_[j]));
                const double distance = displacement.Length();
                if (distance < contactDistance) {
                    // This short-range contact safeguard resolves coincident
                    // particles, whose radial kernel gradient is exactly zero.
                    const Vec2 direction = distance > 1e-10
                        ? displacement / distance : SeparationDirection(i, j);
                    correction = correction + direction * (0.5 * (contactDistance - distance));
                }
            }
            const double length = correction.Length();
            if (length > maximumCorrection) {
                correction = correction * (maximumCorrection / length);
            }
            corrections_[i] = correction;
            largestCorrection = (std::max)(largestCorrection, correction.Length());
        }
        return largestCorrection;
    }

    void SolveDensity() {
        double currentError = ComputeDensities(predicted_);
        lastPredictedDensityError_ = currentError;
        for (int iteration = 0; iteration < constraintIterations_; ++iteration) {
            ComputeMultipliers();
            if (ComputeCorrections() < 1e-10) {
                break;
            }
            bool accepted = false;
            double scale = 1.0;
            // A Jacobi step can over-correct a crowded boundary neighborhood.
            // Backtracking keeps the measured mean compression non-increasing.
            for (int attempt = 0; attempt < 5; ++attempt) {
                for (std::size_t i = 0; i < predicted_.size(); ++i) {
                    candidate_[i] = ConstrainToContainer(predicted_[i] + corrections_[i] * scale);
                }
                const double candidateError = ComputeDensities(candidate_);
                if (candidateError <= currentError + 1e-10) {
                    predicted_.swap(candidate_);
                    currentError = candidateError;
                    accepted = true;
                    break;
                }
                scale *= 0.5;
            }
            if (!accepted) {
                // Restore neighbor/density data to the positions retained.
                currentError = ComputeDensities(predicted_);
                break;
            }
        }
        lastSolvedDensityError_ = ComputeDensities(predicted_);
        averageDensityError_ = lastSolvedDensityError_;
    }

    Vec2 BoundaryVelocity(Vec2 position, Vec2 velocity) const {
        const double low = -2.0 + particleRadius_;
        const double high = 2.0 - particleRadius_;
        const double tolerance = 1e-8;
        if ((position.x <= low + tolerance && velocity.x < 0.0)
            || (position.x >= high - tolerance && velocity.x > 0.0)) {
            velocity.x *= -0.02;
        }
        if ((position.y <= low + tolerance && velocity.y < 0.0)
            || (position.y >= high - tolerance && velocity.y > 0.0)) {
            velocity.y *= -0.02;
        }
        if (position.y <= low + tolerance) {
            velocity.x *= 0.99;
        }
        return LimitSpeed(velocity);
    }

    void ApplyViscosity(double dt) {
        const double blend = 1.0 - std::exp(-viscosity_ * dt
            / (smoothingRadius_ * smoothingRadius_));
        for (std::size_t i = 0; i < particles_.size(); ++i) {
            const Vec2 velocity = particles_[i].velocity;
            Vec2 smoothing{};
            for (int index : neighbors_[i]) {
                const std::size_t j = static_cast<std::size_t>(index);
                if (j == i) {
                    continue;
                }
                // A symmetric weight bounded by either particle's density
                // makes this a momentum-preserving convex XSPH smoothing step.
                const double denominator = (std::max)(densities_[i], densities_[j]);
                if (denominator > 1e-12) {
                    const double weight = particleMass_ * Kernel(predicted_[i] - predicted_[j]) / denominator;
                    smoothing = smoothing + (particles_[j].velocity - velocity) * weight;
                }
            }
            velocityBuffer_[i] = BoundaryVelocity(predicted_[i], velocity + smoothing * blend);
        }
        for (std::size_t i = 0; i < particles_.size(); ++i) {
            particles_[i].velocity = velocityBuffer_[i];
        }
    }

    void Substep(double dt) {
        for (std::size_t i = 0; i < particles_.size(); ++i) {
            oldPositions_[i] = particles_[i].position;
            particles_[i].velocity = LimitSpeed(particles_[i].velocity + gravity_ * dt);
            predicted_[i] = ConstrainToContainer(oldPositions_[i] + particles_[i].velocity * dt);
        }
        SolveDensity();
        for (std::size_t i = 0; i < particles_.size(); ++i) {
            particles_[i].velocity = BoundaryVelocity(predicted_[i],
                (predicted_[i] - oldPositions_[i]) / dt);
        }
        ApplyViscosity(dt);
        for (std::size_t i = 0; i < particles_.size(); ++i) {
            particles_[i].position = predicted_[i];
        }
        simulatedTime_ += dt;
    }

    const double pi_ = 3.14159265358979323846;
    const double nominalRestDensity_ = 1000.0;
    const double smoothingRadius_ = 0.18;
    const double spacing_ = 0.078;
    const double particleRadius_ = 0.033;
    const double fixedTimeStep_ = 1.0 / 120.0;
    const double maximumSpeed_ = 12.0;
    const std::size_t maximumParticles_ = 1200;
    double particleMass_ = 1.0;
    double restDensity_ = 1000.0;
    double viscosity_ = 0.025;
    Vec2 gravity_{ 0.0, -9.81 };
    int constraintIterations_ = 6;
    int gridWidth_ = 1;
    double accumulator_ = 0.0;
    double simulatedTime_ = 0.0;
    double averageDensityError_ = 0.0;
    double lastPredictedDensityError_ = 0.0;
    double lastSolvedDensityError_ = 0.0;
    std::vector<WaterParticle> particles_;
    std::vector<Vec2> oldPositions_;
    std::vector<Vec2> predicted_;
    std::vector<Vec2> candidate_;
    std::vector<Vec2> corrections_;
    std::vector<Vec2> velocityBuffer_;
    std::vector<double> densities_;
    std::vector<double> multipliers_;
    std::vector<int> gridHead_;
    std::vector<int> gridNext_;
    std::vector<std::vector<int>> neighbors_;
};

} // namespace fieldlab

#endif
