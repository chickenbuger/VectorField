#ifndef FIELDLAB_FLUID_GRID_H
#define FIELDLAB_FLUID_GRID_H

#include "FieldMath.h"
#include <vector>

namespace fieldlab {

// A length-4 square, cell-centered dye and staggered (MAC) face velocities.
// Public positions are normalized to [0,1]^2; velocities are length/second.
class FluidGrid {
public:
    explicit FluidGrid(int n = 64)
        : n_((std::max)(8, (std::min)(n, 256))), stride_(n_ + 2),
          h_(4.0 / n_), u_(Count(), 0.0), v_(Count(), 0.0),
          dye_(Count(), 0.0), u0_(Count(), 0.0), v0_(Count(), 0.0),
          dye0_(Count(), 0.0), pressure_(Count(), 0.0), divergence_(Count(), 0.0) {}

    int Size() const { return n_; }
    double CellWidth() const { return h_; }
    double Dye(int i, int j) const { return dye_[Index(Cell(i), Cell(j))]; }
    double U(int i, int j) const {
        i = Cell(i); j = Cell(j);
        return 0.5 * (u_[Index(i - 1, j)] + u_[Index(i, j)]);
    }
    double V(int i, int j) const {
        i = Cell(i); j = Cell(j);
        return 0.5 * (v_[Index(i, j - 1)] + v_[Index(i, j)]);
    }

    double divergenceBeforeProjection = 0.0;
    double divergenceAfterProjection = 0.0;
    double maxDivergenceBefore = 0.0;
    double maxDivergenceAfter = 0.0;

    void Clear() {
        for (auto* values : { &u_, &v_, &dye_, &u0_, &v0_, &dye0_, &pressure_, &divergence_ }) {
            std::fill(values->begin(), values->end(), 0.0);
        }
        divergenceBeforeProjection = divergenceAfterProjection = 0.0;
        maxDivergenceBefore = maxDivergenceAfter = 0.0;
    }

    void FillDye(double value) {
        value = std::isfinite(value) ? Clamp(value, 0.0, 10.0) : 0.0;
        std::fill(dye_.begin(), dye_.end(), value);
    }

    void SeedDemo() {
        Clear();
        SetPresetVelocity(FieldType::Rotation, 0.65);
        AddDye(0.50, 0.27, 3.0, 0.10);
        AddDye(0.30, 0.62, 2.0, 0.08);
        AddDye(0.68, 0.65, 2.5, 0.07);
    }

    void AddDye(double x, double y, double amount, double radius = 0.04) {
        if (!ValidSource(x, y, amount, radius) || amount <= 0.0) return;
        x = Clamp(x, 0.0, 1.0); y = Clamp(y, 0.0, 1.0);
        radius = (std::max)(radius, 0.75 / n_);
        for (int j = 1; j <= n_; ++j) {
            for (int i = 1; i <= n_; ++i) {
                const double weight = SplatWeight((i - 0.5) / n_, (j - 0.5) / n_, x, y, radius);
                dye_[Index(i, j)] = Clamp(dye_[Index(i, j)] + amount * weight, 0.0, 10.0);
            }
        }
        Boundary(dye_, 0);
    }

    // Adds an instantaneous velocity increment. Callers multiply forces by dt.
    void AddVelocity(double x, double y, double forceX, double forceY, double radius = 0.06) {
        if (!ValidSource(x, y, forceX, radius) || !std::isfinite(forceY)) return;
        x = Clamp(x, 0.0, 1.0); y = Clamp(y, 0.0, 1.0);
        radius = (std::max)(radius, 0.75 / n_);
        for (int j = 1; j <= n_; ++j) {
            for (int i = 1; i < n_; ++i) {
                const double weight = SplatWeight(double(i) / n_, (j - 0.5) / n_, x, y, radius);
                u_[Index(i, j)] = Clamp(u_[Index(i, j)] + forceX * weight, -20.0, 20.0);
            }
        }
        for (int j = 1; j < n_; ++j) {
            for (int i = 1; i <= n_; ++i) {
                const double weight = SplatWeight((i - 0.5) / n_, double(j) / n_, x, y, radius);
                v_[Index(i, j)] = Clamp(v_[Index(i, j)] + forceY * weight, -20.0, 20.0);
            }
        }
        Boundary(u_, 1); Boundary(v_, 2);
    }

    void SetPresetVelocity(FieldType type, double strength = 1.0) {
        if (!std::isfinite(strength)) return;
        for (int j = 1; j <= n_; ++j) {
            for (int i = 0; i <= n_; ++i) {
                const Vec2 p{ 4.0 * i / n_ - 2.0, 4.0 * (j - 0.5) / n_ - 2.0 };
                u_[Index(i, j)] = Clamp(Evaluate(type, p, strength).x, -20.0, 20.0);
            }
        }
        for (int j = 0; j <= n_; ++j) {
            for (int i = 1; i <= n_; ++i) {
                const Vec2 p{ 4.0 * (i - 0.5) / n_ - 2.0, 4.0 * j / n_ - 2.0 };
                v_[Index(i, j)] = Clamp(Evaluate(type, p, strength).y, -20.0, 20.0);
            }
        }
        Boundary(u_, 1); Boundary(v_, 2);
    }

    Vec2 SampleVelocity(Vec2 normalized) const { return SampleVelocityFrom(u_, v_, normalized); }
    double SampleDye(Vec2 normalized) const {
        if (!Finite(normalized)) return 0.0;
        return Bilinear(dye_, Clamp(normalized.x, 0.0, 1.0) * n_ + 0.5,
            Clamp(normalized.y, 0.0, 1.0) * n_ + 0.5);
    }

    void StepDye(double dt, double diffusion = 0.00002) {
        if (!(dt > 0.0) || !std::isfinite(dt)) return;
        dt = (std::min)(dt, 0.1);
        dye0_ = dye_;
        Diffuse(dye_, dye0_, dt, diffusion, 0);
        dye0_ = dye_;
        for (int j = 1; j <= n_; ++j) {
            for (int i = 1; i <= n_; ++i) {
                const Vec2 p{ (i - 0.5) / n_, (j - 0.5) / n_ };
                const Vec2 back = p - SampleVelocity(p) * (dt / 4.0);
                dye_[Index(i, j)] = (std::max)(0.0,
                    Bilinear(dye0_, Clamp(back.x, 0.0, 1.0) * n_ + 0.5,
                        Clamp(back.y, 0.0, 1.0) * n_ + 0.5));
            }
        }
        Boundary(dye_, 0);
    }

    void Step(double dt, double viscosity = 0.00015, double diffusion = 0.00002,
        bool project = true, bool advectVelocity = true, bool buoyancy = false) {
        if (!(dt > 0.0) || !std::isfinite(dt)) return;
        dt = (std::min)(dt, 0.1);
        if (buoyancy) {
            for (int j = 1; j < n_; ++j) {
                for (int i = 1; i <= n_; ++i) {
                    const double density = 0.5 * (dye_[Index(i, j)] + dye_[Index(i, j + 1)]);
                    v_[Index(i, j)] = Clamp(v_[Index(i, j)] + dt * 1.2 * density, -20.0, 20.0);
                }
            }
        }
        u0_ = u_; v0_ = v_;
        Diffuse(u_, u0_, dt, viscosity, 1);
        Diffuse(v_, v0_, dt, viscosity, 2);
        if (project) Project();
        if (advectVelocity) {
            u0_ = u_; v0_ = v_;
            for (int j = 1; j <= n_; ++j) {
                for (int i = 1; i < n_; ++i) {
                    const Vec2 p{ double(i) / n_, (j - 0.5) / n_ };
                    const Vec2 back = p - SampleVelocityFrom(u0_, v0_, p) * (dt / 4.0);
                    u_[Index(i, j)] = SampleVelocityFrom(u0_, v0_, back).x;
                }
            }
            for (int j = 1; j < n_; ++j) {
                for (int i = 1; i <= n_; ++i) {
                    const Vec2 p{ (i - 0.5) / n_, double(j) / n_ };
                    const Vec2 back = p - SampleVelocityFrom(u0_, v0_, p) * (dt / 4.0);
                    v_[Index(i, j)] = SampleVelocityFrom(u0_, v0_, back).y;
                }
            }
            Boundary(u_, 1); Boundary(v_, 2);
            if (project) Project();
        }
        if (!project) {
            divergenceBeforeProjection = divergenceAfterProjection = DivergenceRms();
            maxDivergenceBefore = maxDivergenceAfter = DivergenceMaximum();
        }
        StepDye(dt, diffusion);
    }

    // Solves for a pressure impulse q=dt*p/rho, so the velocity update is u-=grad(q).
    void Project() {
        Boundary(u_, 1); Boundary(v_, 2);
        divergenceBeforeProjection = DivergenceRms();
        maxDivergenceBefore = DivergenceMaximum();
        std::fill(pressure_.begin(), pressure_.end(), 0.0);
        for (int j = 1; j <= n_; ++j) {
            for (int i = 1; i <= n_; ++i) divergence_[Index(i, j)] = CellDivergence(i, j);
        }
        // Neumann pressure walls. Their missing neighbor also leaves the diagonal.
        // SOR accelerates the five-point Poisson solve without a checkerboard mode.
        const double relaxation = 1.7;
        for (int iteration = 0; iteration < 100; ++iteration) {
            for (int j = 1; j <= n_; ++j) {
                for (int i = 1; i <= n_; ++i) {
                    double sum = 0.0;
                    int neighbors = 0;
                    if (i > 1) { sum += pressure_[Index(i - 1, j)]; ++neighbors; }
                    if (i < n_) { sum += pressure_[Index(i + 1, j)]; ++neighbors; }
                    if (j > 1) { sum += pressure_[Index(i, j - 1)]; ++neighbors; }
                    if (j < n_) { sum += pressure_[Index(i, j + 1)]; ++neighbors; }
                    const double target = (sum - h_ * h_ * divergence_[Index(i, j)]) / neighbors;
                    double& pressure = pressure_[Index(i, j)];
                    pressure += relaxation * (target - pressure);
                }
            }
        }
        for (int j = 1; j <= n_; ++j) {
            for (int i = 1; i < n_; ++i) {
                u_[Index(i, j)] -= (pressure_[Index(i + 1, j)] - pressure_[Index(i, j)]) / h_;
            }
        }
        for (int j = 1; j < n_; ++j) {
            for (int i = 1; i <= n_; ++i) {
                v_[Index(i, j)] -= (pressure_[Index(i, j + 1)] - pressure_[Index(i, j)]) / h_;
            }
        }
        Boundary(u_, 1); Boundary(v_, 2);
        divergenceAfterProjection = DivergenceRms();
        maxDivergenceAfter = DivergenceMaximum();
    }

    double DivergenceRms() const {
        double sum = 0.0;
        for (int j = 1; j <= n_; ++j) {
            for (int i = 1; i <= n_; ++i) {
                const double value = CellDivergence(i, j);
                sum += value * value;
            }
        }
        return std::sqrt(sum / (n_ * n_));
    }

private:
    int n_, stride_;
    double h_;
    std::vector<double> u_, v_, dye_, u0_, v0_, dye0_, pressure_, divergence_;

    std::size_t Count() const { return std::size_t(stride_) * stride_; }
    std::size_t Index(int i, int j) const { return std::size_t(i) + std::size_t(stride_) * j; }
    int Cell(int value) const { return (std::max)(1, (std::min)(value, n_)); }
    static double Clamp(double x, double a, double b) { return (std::max)(a, (std::min)(b, x)); }
    static bool Finite(Vec2 p) { return std::isfinite(p.x) && std::isfinite(p.y); }
    static bool ValidSource(double x, double y, double amount, double radius) {
        return std::isfinite(x) && std::isfinite(y) && std::isfinite(amount)
            && std::isfinite(radius) && radius > 0.0;
    }
    static double SplatWeight(double px, double py, double x, double y, double radius) {
        const double distance2 = ((px - x) * (px - x) + (py - y) * (py - y)) / (radius * radius);
        const double falloff = (std::max)(0.0, 1.0 - distance2);
        return falloff * falloff;
    }

    double Bilinear(const std::vector<double>& values, double x, double y) const {
        x = Clamp(x, 0.0, double(n_ + 1)); y = Clamp(y, 0.0, double(n_ + 1));
        const int i = int(std::floor(x)), j = int(std::floor(y));
        const int right = (std::min)(i + 1, n_ + 1), top = (std::min)(j + 1, n_ + 1);
        const double sx = x - i, sy = y - j;
        return (1.0 - sy) * ((1.0 - sx) * values[Index(i, j)] + sx * values[Index(right, j)])
            + sy * ((1.0 - sx) * values[Index(i, top)] + sx * values[Index(right, top)]);
    }

    Vec2 SampleVelocityFrom(const std::vector<double>& u, const std::vector<double>& v, Vec2 p) const {
        if (!Finite(p)) return {};
        p.x = Clamp(p.x, 0.0, 1.0); p.y = Clamp(p.y, 0.0, 1.0);
        return { Bilinear(u, p.x * n_, p.y * n_ + 0.5),
            Bilinear(v, p.x * n_ + 0.5, p.y * n_) };
    }

    void Boundary(std::vector<double>& values, int component) const {
        if (component == 0) {
            for (int i = 1; i <= n_; ++i) {
                values[Index(0, i)] = values[Index(1, i)];
                values[Index(n_ + 1, i)] = values[Index(n_, i)];
                values[Index(i, 0)] = values[Index(i, 1)];
                values[Index(i, n_ + 1)] = values[Index(i, n_)];
            }
            values[Index(0, 0)] = values[Index(1, 1)];
            values[Index(0, n_ + 1)] = values[Index(1, n_)];
            values[Index(n_ + 1, 0)] = values[Index(n_, 1)];
            values[Index(n_ + 1, n_ + 1)] = values[Index(n_, n_)];
        } else if (component == 1) {
            for (int j = 0; j <= n_ + 1; ++j) {
                values[Index(0, j)] = values[Index(n_, j)] = values[Index(n_ + 1, j)] = 0.0;
            }
            for (int i = 0; i <= n_; ++i) {
                values[Index(i, 0)] = values[Index(i, 1)];
                values[Index(i, n_ + 1)] = values[Index(i, n_)];
            }
        } else {
            for (int i = 0; i <= n_ + 1; ++i) {
                values[Index(i, 0)] = values[Index(i, n_)] = values[Index(i, n_ + 1)] = 0.0;
            }
            for (int j = 0; j <= n_; ++j) {
                values[Index(0, j)] = values[Index(1, j)];
                values[Index(n_ + 1, j)] = values[Index(n_, j)];
            }
        }
    }

    void Diffuse(std::vector<double>& output, const std::vector<double>& input,
        double dt, double coefficient, int component) const {
        output = input;
        if (!(coefficient > 0.0) || !std::isfinite(coefficient)) { Boundary(output, component); return; }
        const double a = (std::min)(coefficient, 1.0) * dt / (h_ * h_);
        const int maxI = component == 1 ? n_ - 1 : n_;
        const int maxJ = component == 2 ? n_ - 1 : n_;
        Boundary(output, component);
        for (int iteration = 0; iteration < 24; ++iteration) {
            for (int j = 1; j <= maxJ; ++j) {
                for (int i = 1; i <= maxI; ++i) {
                    output[Index(i, j)] = (input[Index(i, j)] + a *
                        (output[Index(i - 1, j)] + output[Index(i + 1, j)]
                            + output[Index(i, j - 1)] + output[Index(i, j + 1)])) / (1.0 + 4.0 * a);
                }
            }
            Boundary(output, component);
        }
    }

    double CellDivergence(int i, int j) const {
        return (u_[Index(i, j)] - u_[Index(i - 1, j)]
            + v_[Index(i, j)] - v_[Index(i, j - 1)]) / h_;
    }
    double DivergenceMaximum() const {
        double result = 0.0;
        for (int j = 1; j <= n_; ++j) {
            for (int i = 1; i <= n_; ++i) result = (std::max)(result, std::abs(CellDivergence(i, j)));
        }
        return result;
    }
};

} // namespace fieldlab
#endif
