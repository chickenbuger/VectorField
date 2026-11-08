#ifndef FIELDLAB_TWO_MATERIAL_SIMULATION_H
#define FIELDLAB_TWO_MATERIAL_SIMULATION_H

#include "FluidGrid.h"
#include "WaterParticles.h"
#include <limits>

namespace fieldlab {

enum class CouplingMode { Off = 0, OneWay = 1, TwoWay = 2 };

// A regularized immersed-drag coupling of the existing MAC carrier fluid and
// PBF liquid. The gas mesh extends through the liquid. Occupancy excludes its
// visual tracer, but is not a cut-cell or sharp-interface pressure boundary.
// The exchange suboperator conserves discrete momentum and dissipates kinetic
// energy. Projection, gravity, emitters and PBF walls remain separate operators.
class TwoMaterialSimulation {
public:
    explicit TwoMaterialSimulation(int n = 48)
        : gas_(n), occupancy_(static_cast<std::size_t>(gas_.Size() * gas_.Size()), 0.0) { Reset(); }

    FluidGrid& Gas() { return gas_; }
    const FluidGrid& Gas() const { return gas_; }
    WaterParticles& Liquid() { return liquid_; }
    const WaterParticles& Liquid() const { return liquid_; }
    CouplingMode Mode() const { return mode_; }
    bool JetEnabled() const { return jetEnabled_; }
    double GasDensity() const { return gasDensity_; }
    double LiquidDensity() const { return liquid_.RestDensity(); }
    double DragRate() const { return dragRate_; }
    double CouplingStrength() const { return dragRate_; }
    double JetStrength() const { return jetStrength_; }
    double Time() const { return time_; }
    double FixedTimeStep() const { return 1.0 / 120.0; }
    double LiquidMass() const { return liquid_.Count() * liquid_.ParticleMass(); }
    double ReferenceGasMass() const { return 16.0 * gasDensity_; }
    double Occupancy(int i, int j) const {
        i = (std::max)(1, (std::min)(gas_.Size(), i));
        j = (std::max)(1, (std::min)(gas_.Size(), j));
        return occupancy_[static_cast<std::size_t>((j - 1) * gas_.Size() + i - 1)];
    }
    Vec2 WaterMomentum() const { return liquid_.Momentum(); }
    Vec2 GasMomentum() const { return gas_.DiscreteMomentum(gasDensity_); }
    Vec2 TotalMomentum() const { return WaterMomentum() + GasMomentum(); }
    double TotalKineticEnergy() const {
        return liquid_.KineticEnergy() + gas_.DiscreteKineticEnergy(gasDensity_);
    }

    // All three vectors below are independently measured from actual states
    // immediately before/after ExchangeOnly, rather than assigned +J/-J ledgers.
    Vec2 couplingWaterImpulse{};
    Vec2 couplingGasImpulse{};
    Vec2 couplingResidual{};
    double maxCouplingResidual = 0.0;
    double cumulativeExchangeMagnitude = 0.0;
    double lastExchangeMagnitude = 0.0;
    double lastExchangeEnergyChange = 0.0;
    double removedDyeTracer = 0.0;

    void SetCouplingMode(CouplingMode mode) {
        if (mode != CouplingMode::Off && mode != CouplingMode::OneWay && mode != CouplingMode::TwoWay)
            throw std::invalid_argument("Unknown material coupling mode.");
        mode_ = mode;
    }
    void SetJetEnabled(bool enabled) { jetEnabled_ = enabled; }
    void SetGasDensity(double density) {
        if (!std::isfinite(density) || density <= 0.0)
            throw std::invalid_argument("Gas reference density must be finite and positive.");
        gasDensity_ = (std::max)(0.1, (std::min)(1000.0, density));
    }
    void SetDragRate(double rate) {
        if (!std::isfinite(rate) || rate < 0.0)
            throw std::invalid_argument("Drag rate must be finite and nonnegative.");
        dragRate_ = (std::min)(300.0, rate);
    }
    void SetCouplingStrength(double rate) { SetDragRate(rate); }
    void SetJetStrength(double strength) {
        if (!std::isfinite(strength) || strength < 0.0)
            throw std::invalid_argument("Jet strength must be finite and nonnegative.");
        jetStrength_ = (std::min)(4.0, strength);
    }

    void Reset() {
        gas_.Clear();
        liquid_.SetRestDensity(1000.0);
        liquid_.SetViscosity(0.025);
        liquid_.SetGravity({0.0, -9.81});
        std::vector<WaterParticle> block;
        for (int j = 0; j < 20; ++j) for (int i = 0; i < 20; ++i)
            block.push_back({{-0.75 + i * liquid_.ParticleSpacing(), -1.82 + j * liquid_.ParticleSpacing()}, {}});
        liquid_.SetParticleState(block);
        mode_ = CouplingMode::TwoWay;
        gasDensity_ = 80.0; // Pedagogical density ratio 12.5, not physical room-air density.
        dragRate_ = 45.0; jetStrength_ = 1.0; jetEnabled_ = true;
        time_ = accumulator_ = 0.0;
        couplingWaterImpulse = couplingGasImpulse = couplingResidual = {};
        maxCouplingResidual = cumulativeExchangeMagnitude = lastExchangeMagnitude = 0.0;
        lastExchangeEnergyChange = removedDyeTracer = 0.0;
        RebuildOccupancy();
    }

    // Smooth 2D volume fraction from the same normalized poly6 kernel used in
    // the liquid density estimate. Clipping makes this a bounded phase indicator;
    // exact water mass remains particle count * fixed particle mass.
    void RebuildOccupancy() {
        std::fill(occupancy_.begin(), occupancy_.end(), 0.0);
        const int n = gas_.Size();
        const double cell = gas_.CellWidth();
        const double radius = liquid_.SmoothingRadius();
        const double radius2 = radius * radius;
        const double radius4 = radius2 * radius2;
        const double coefficient = 4.0 / (3.14159265358979323846 * radius4 * radius4);
        const double volume = liquid_.ParticleMass() / liquid_.RestDensity();
        for (const auto& particle : liquid_.Particles()) {
            const int minI = (std::max)(1, static_cast<int>(std::floor((particle.position.x - radius + 2.0) / cell)) + 1);
            const int maxI = (std::min)(n, static_cast<int>(std::floor((particle.position.x + radius + 2.0) / cell)) + 1);
            const int minJ = (std::max)(1, static_cast<int>(std::floor((particle.position.y - radius + 2.0) / cell)) + 1);
            const int maxJ = (std::min)(n, static_cast<int>(std::floor((particle.position.y + radius + 2.0) / cell)) + 1);
            for (int j = minJ; j <= maxJ; ++j) for (int i = minI; i <= maxI; ++i) {
                const Vec2 cellPosition{-2.0 + (i - 0.5) * cell, -2.0 + (j - 0.5) * cell};
                const Vec2 d = cellPosition - particle.position;
                const double q = radius2 - d.x * d.x - d.y * d.y;
                if (q > 0.0) occupancy_[static_cast<std::size_t>((j - 1) * n + i - 1)] += volume * coefficient * q * q * q;
            }
        }
        for (double& value : occupancy_) value = (std::max)(0.0, (std::min)(1.0, value));
    }

    // dt is the exchange duration, not a frame count. Density and the gas face
    // masses are fixed during this operation. Large dt still has alpha in [0,1].
    void ExchangeOnly(double dt) {
        couplingWaterImpulse = couplingGasImpulse = couplingResidual = {};
        lastExchangeMagnitude = lastExchangeEnergyChange = 0.0;
        if (!(dt > 0.0) || !std::isfinite(dt) || mode_ == CouplingMode::Off || dragRate_ == 0.0) return;
        RebuildOccupancy();
        const Vec2 waterBefore = WaterMomentum(), gasBefore = GasMomentum();
        const double energyBefore = TotalKineticEnergy();
        const double particleMass = liquid_.ParticleMass();
        const double gasFaceMass = gasDensity_ * gas_.CellWidth() * gas_.CellWidth();
        for (std::size_t index = 0; index < liquid_.Count(); ++index) {
            const Vec2 normalized = (liquid_.Particles()[index].position + Vec2{2.0,2.0}) / 4.0;
            // A regularized exposed-surface weight, with 5% interior coupling.
            // It is frozen for this particle's component exchange.
            const double exposure = 0.05 + 0.95 * (1.0 - SampleOccupancy(normalized));
            const double alpha = -std::expm1(-dragRate_ * dt * exposure);
            for (int component = 0; component <= 1; ++component) {
                const auto stencil = gas_.InteriorVelocityStencil(normalized, component);
                if (stencil.empty()) continue;
                double gasVelocity = 0.0, inverseGasMass = 0.0;
                for (const auto& face : stencil) {
                    gasVelocity += face.weight * gas_.MacFaceVelocity(component, face.i, face.j);
                    inverseGasMass += face.weight * face.weight / gasFaceMass;
                }
                const Vec2 velocity = liquid_.Particles()[index].velocity;
                const double liquidVelocity = component == 0 ? velocity.x : velocity.y;
                const double inverseMass = 1.0 / particleMass + (mode_ == CouplingMode::TwoWay ? inverseGasMass : 0.0);
                const double impulse = (gasVelocity - liquidVelocity) * alpha / inverseMass;
                const Vec2 actual = liquid_.AddVelocityIncrement(index, component == 0
                    ? Vec2{impulse / particleMass, 0.0} : Vec2{0.0, impulse / particleMass});
                const double appliedImpulse = particleMass * (component == 0 ? actual.x : actual.y);
                lastExchangeMagnitude += std::abs(appliedImpulse);
                if (mode_ == CouplingMode::TwoWay) for (const auto& face : stencil)
                    gas_.AddMacFaceVelocity(component, face.i, face.j, -appliedImpulse * face.weight / gasFaceMass);
            }
        }
        gas_.RefreshBoundaries();
        couplingWaterImpulse = WaterMomentum() - waterBefore;
        couplingGasImpulse = GasMomentum() - gasBefore;
        couplingResidual = couplingWaterImpulse + couplingGasImpulse;
        lastExchangeEnergyChange = TotalKineticEnergy() - energyBefore;
        maxCouplingResidual = (std::max)(maxCouplingResidual, couplingResidual.Length());
        cumulativeExchangeMagnitude += lastExchangeMagnitude;
    }

    void Step(double dt) {
        if (!(dt > 0.0) || !std::isfinite(dt)) return;
        accumulator_ += (std::min)(0.05, dt);
        const double fixed = FixedTimeStep();
        int count = 0;
        while (accumulator_ + 1e-12 >= fixed && count < 6) {
            if (jetEnabled_) {
                // A horizontal/upward gas jet meets the left side of the water.
                const double nozzleY = 0.19 + 0.012 * std::sin(time_ * 1.3);
                gas_.AddDye(0.10, nozzleY, 5.0 * fixed * jetStrength_, 0.085);
                gas_.AddVelocity(0.10, nozzleY, 75.0 * fixed * jetStrength_, 18.0 * fixed * jetStrength_, 0.10);
            }
            // Dye is passive: no concentration buoyancy is used, so the gas
            // velocity in Off and OneWay benchmarks has exactly the same inputs.
            gas_.Step(fixed, 0.00015, 0.00002, true, true, false);
            ExchangeOnly(fixed);
            liquid_.Step(fixed);
            // Project after immersed impulses. Its momentum change is outside
            // the exchange-conservation measurement and can include box reaction.
            gas_.Project();
            RebuildOccupancy();
            removedDyeTracer += gas_.ExcludeDyeFromLiquid(occupancy_);
            time_ += fixed; accumulator_ = (std::max)(0.0, accumulator_ - fixed); ++count;
        }
    }

    bool IsFinite() const {
        if (!liquid_.IsFinite() || !std::isfinite(time_) || !std::isfinite(accumulator_)
            || !std::isfinite(maxCouplingResidual) || !std::isfinite(TotalKineticEnergy())) return false;
        for (double value : occupancy_) if (!std::isfinite(value) || value < 0.0 || value > 1.0) return false;
        for (int j = 1; j <= gas_.Size(); ++j) for (int i = 1; i <= gas_.Size(); ++i) {
            if (!std::isfinite(gas_.U(i,j)) || !std::isfinite(gas_.V(i,j))
                || !std::isfinite(gas_.Dye(i,j)) || gas_.Dye(i,j) < 0.0 || gas_.Dye(i,j) > 10.000001) return false;
        }
        return true;
    }

private:
    double SampleOccupancy(Vec2 p) const {
        const int n = gas_.Size();
        const double x = (std::max)(0.0,(std::min)(double(n-1),p.x*n-0.5));
        const double y = (std::max)(0.0,(std::min)(double(n-1),p.y*n-0.5));
        const int i = static_cast<int>(std::floor(x)), j = static_cast<int>(std::floor(y));
        const int right = (std::min)(n-1,i+1), top = (std::min)(n-1,j+1);
        const double sx=x-i, sy=y-j;
        return (1.0-sy)*((1.0-sx)*Occupancy(i+1,j+1)+sx*Occupancy(right+1,j+1))
            + sy*((1.0-sx)*Occupancy(i+1,top+1)+sx*Occupancy(right+1,top+1));
    }
    FluidGrid gas_;
    WaterParticles liquid_;
    std::vector<double> occupancy_;
    CouplingMode mode_ = CouplingMode::TwoWay;
    double gasDensity_ = 80.0;
    double dragRate_ = 45.0;
    double jetStrength_ = 1.0;
    bool jetEnabled_ = true;
    double time_ = 0.0;
    double accumulator_ = 0.0;
};

} // namespace fieldlab
#endif
