#ifndef FIELDLAB_WATER_PARTICLE_TESTS_H
#define FIELDLAB_WATER_PARTICLE_TESTS_H

#include "../opgnel_project_smoke/WaterParticles.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace fieldlab {
namespace water_particle_test_detail {

inline bool Near(double value, double expected, double tolerance = 1e-9) {
    return std::isfinite(value) && std::abs(value - expected) <= tolerance;
}

inline double Energy(const WaterParticles& water) {
    double energy = 0.0;
    for (const WaterParticle& particle : water.Particles()) {
        energy += 0.5 * water.ParticleMass()
            * (particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);
    }
    return energy;
}

inline bool InsideContainer(const WaterParticles& water) {
    const double limit = 2.0 - water.ParticleRadius() + 1e-9;
    for (const WaterParticle& particle : water.Particles()) {
        if (std::abs(particle.position.x) > limit || std::abs(particle.position.y) > limit) {
            return false;
        }
    }
    return true;
}

} // namespace water_particle_test_detail

inline int RunWaterParticleTests() {
    using namespace water_particle_test_detail;
    int passed = 0;
    int failed = 0;
    const auto check = [&](const std::string& name, bool success) {
        std::cout << (success ? "PASS " : "FAIL ") << name << '\n';
        if (success) { ++passed; }
        else { ++failed; }
    };

    WaterParticles water;
    const double dt = water.FixedTimeStep();
    check("water starts with 484 particles", water.Count() == 484);
    check("water initial state is finite and inside the container", water.IsFinite() && InsideContainer(water));
    check("water particle mass is positive and finite", std::isfinite(water.ParticleMass()) && water.ParticleMass() > 0.0);
    check("discrete rest lattice has no positive density excess", water.AverageDensityError() < 1e-10);
    check("fixed water timestep is 1/120 second", Near(dt, 1.0 / 120.0));

    water.SetParticleState({ { { 0.0, 0.5 }, {} } });
    water.SetViscosity(0.0);
    water.Step(dt);
    check("gravity changes isolated particle velocity", Near(water.Particles()[0].velocity.y, -9.81 * dt));
    check("semi-implicit gravity prediction changes position", Near(water.Particles()[0].position.y, 0.5 - 9.81 * dt * dt));
    check("isolated free-surface particle is not pulled to a target density", Near(water.AverageDensityError(), 0.0));

    WaterParticles splitStep;
    WaterParticles wholeStep;
    splitStep.SetParticleState({ { { 0.0, 0.5 }, { 0.2, 0.0 } } });
    wholeStep.SetParticleState(splitStep.Particles());
    splitStep.Step(dt * 0.5);
    check("half timestep is accumulated before integration", Near(splitStep.SimulatedTime(), 0.0));
    splitStep.Step(dt * 0.5);
    wholeStep.Step(dt);
    check("two half frame intervals equal one fixed substep",
        Near(splitStep.Particles()[0].position.x, wholeStep.Particles()[0].position.x)
        && Near(splitStep.Particles()[0].position.y, wholeStep.Particles()[0].position.y));
    const double timeBeforeInvalidStep = water.SimulatedTime();
    water.Step(0.0);
    water.Step(-1.0);
    water.Step(std::numeric_limits<double>::quiet_NaN());
    check("invalid frame intervals do not advance time", Near(water.SimulatedTime(), timeBeforeInvalidStep));

    water.SetParticleState({ { { 1.96, 0.0 }, { 12.0, 0.0 } }, { { 0.0, -1.96 }, { 0.0, -12.0 } } });
    water.Step(dt);
    check("wall and floor block particle penetration", InsideContainer(water));
    check("wall collision removes outward normal velocity", water.Particles()[0].velocity.x <= 0.0);
    check("floor collision removes downward normal velocity", water.Particles()[1].velocity.y >= 0.0);
    water.SetParticleState({ { { -100.0, 100.0 }, {} } });
    check("injected initial positions are constrained to the container", InsideContainer(water));

    water.SetGravity({});
    water.SetParticleState({ { { 0.0, 0.0 }, {} } });
    water.Step(dt);
    check("zero gravity and zero velocity preserve an isolated particle position",
        Near(water.Particles()[0].position.x, 0.0) && Near(water.Particles()[0].position.y, 0.0));
    water.SetParticleState({ { { 0.0, 0.0 }, {} }, { { 0.0, 0.0 }, {} } });
    water.Step(dt);
    check("coincident particles receive finite separating contact corrections",
        water.IsFinite() && (water.Particles()[0].position - water.Particles()[1].position).Length() > 1e-4);

    std::vector<WaterParticle> compressed;
    for (int y = 0; y < 9; ++y) {
        for (int x = 0; x < 9; ++x) {
            compressed.push_back({ { (x - 4) * 0.043, (y - 4) * 0.043 }, {} });
        }
    }
    water.SetParticleState(compressed);
    water.Step(dt);
    std::cout << "METRIC water compressed-block mean excess: predicted="
        << water.LastPredictedDensityError() << " solved=" << water.LastSolvedDensityError() << '\n';
    check("compressed block has nonzero predicted density excess", water.LastPredictedDensityError() > 0.1);
    check("PBF correction reduces predicted mean density excess",
        water.LastSolvedDensityError() < 0.8 * water.LastPredictedDensityError());
    check("compression metric equals the committed solved metric",
        Near(water.AverageDensityError(), water.LastSolvedDensityError()));
    check("compressed-block correction remains finite", water.IsFinite());

    WaterParticles inviscid;
    WaterParticles viscous;
    const std::vector<WaterParticle> pair{
        { { -0.06, 0.0 }, { 1.0, 0.25 } },
        { { 0.06, 0.0 }, { -1.0, -0.25 } }
    };
    inviscid.SetGravity({});
    viscous.SetGravity({});
    inviscid.SetViscosity(0.0);
    viscous.SetViscosity(0.8);
    inviscid.SetParticleState(pair);
    viscous.SetParticleState(pair);
    inviscid.Step(dt);
    viscous.Step(dt);
    check("XSPH viscosity reduces relative kinetic energy without gravity",
        Energy(viscous) < Energy(inviscid));
    check("symmetric viscosity preserves isolated pair linear momentum",
        Near(viscous.Particles()[0].velocity.x + viscous.Particles()[1].velocity.x, 0.0)
        && Near(viscous.Particles()[0].velocity.y + viscous.Particles()[1].velocity.y, 0.0));
    viscous.SetViscosity(100.0);
    check("viscosity is clamped to the stable material range", Near(viscous.Viscosity(), 2.0));
    viscous.SetViscosity(-1.0);
    check("negative viscosity is not applied", Near(viscous.Viscosity(), 0.0));

    water.SetParticleState({ { { 0.0, 0.0 }, {} }, { { 1.0, 1.0 }, {} } });
    water.Stir({ 0.0, 0.0 }, { 0.05, 0.01 });
    check("stir input adds velocity to nearby fluid", water.Particles()[0].velocity.Length() > 0.0);
    check("stir input is spatially local", Near(water.Particles()[1].velocity.Length(), 0.0));
    water.Stir({ 0.0, 0.0 }, { 1e100, 1e100 });
    check("large stir input remains bounded and finite", water.IsFinite() && water.Particles()[0].velocity.Length() <= 12.0 + 1e-9);

    water.Reset();
    const std::size_t countBeforeAdd = water.Count();
    water.AddBlock({ 1.0, 1.0 });
    check("add block creates additional liquid particles", water.Count() > countBeforeAdd && water.Count() <= countBeforeAdd + 64);
    check("added block respects container and capacity", InsideContainer(water) && water.Count() <= water.MaximumParticles());
    const double massBeforeDensityChange = water.ParticleMass();
    water.SetRestDensity(1200.0);
    check("rest density is independently configurable", Near(water.RestDensity(), 1200.0));
    check("changing rest density does not silently rescale particle mass", Near(water.ParticleMass(), massBeforeDensityChange));
    water.SetRestDensity(1000.0);

    bool rejectedNonfiniteState = false;
    const std::size_t countBeforeInvalidState = water.Count();
    try {
        water.SetParticleState({ { { std::numeric_limits<double>::infinity(), 0.0 }, {} } });
    }
    catch (const std::invalid_argument&) { rejectedNonfiniteState = true; }
    check("nonfinite injected state is rejected without modifying particles",
        rejectedNonfiniteState && water.Count() == countBeforeInvalidState);
    bool rejectedNonfiniteGravity = false;
    try { water.SetGravity({ 0.0, std::numeric_limits<double>::quiet_NaN() }); }
    catch (const std::invalid_argument&) { rejectedNonfiniteGravity = true; }
    check("nonfinite gravity is rejected", rejectedNonfiniteGravity);

    water.SetParticleState({});
    water.Step(dt);
    check("empty fluid has finite zero compression", water.IsFinite() && Near(water.AverageDensityError(), 0.0));

    water.SetGravity({ 0.0, -9.81 });
    water.Reset();
    bool finiteThroughout = true;
    bool insideThroughout = true;
    bool correctionNeverIncreasedCompression = true;
    double peakPredicted = 0.0;
    double peakSolved = 0.0;
    for (int step = 0; step < 600; ++step) {
        if (step == 150) {
            water.Stir({ -0.9, -1.3 }, { 0.2, 0.12 });
        }
        water.Step(dt);
        finiteThroughout = finiteThroughout && water.IsFinite();
        insideThroughout = insideThroughout && InsideContainer(water);
        correctionNeverIncreasedCompression = correctionNeverIncreasedCompression
            && water.LastSolvedDensityError() <= water.LastPredictedDensityError() + 2e-10;
        peakPredicted = (std::max)(peakPredicted, water.LastPredictedDensityError());
        peakSolved = (std::max)(peakSolved, water.LastSolvedDensityError());
    }
    std::cout << "METRIC water 600 fixed steps: peak predicted mean excess=" << peakPredicted
        << " peak solved mean excess=" << peakSolved << '\n';
    check("484-particle dam and stir stay finite for 600 substeps", finiteThroughout);
    check("484-particle dam and stir stay inside the tank for 600 substeps", insideThroughout);
    check("density projection never increases mean positive compression", correctionNeverIncreasedCompression);
    check("fixed-step simulated time matches five seconds", Near(water.SimulatedTime(), 5.0, 1e-9));
    check("closed liquid simulation conserves particle count", water.Count() == 484);
    check("settled liquid mean compression remains below five percent", water.AverageDensityError() < 0.05);

    std::cout << "Water particle tests: " << passed << " passed, " << failed << " failed.\n";
    return failed == 0 ? 0 : 1;
}

} // namespace fieldlab

#endif
