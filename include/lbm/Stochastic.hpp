#pragma once

// Random forcing that can be switched on or off independently:
//   InflowModel  the inlet velocity profile wanders (Ornstein-Uhlenbeck mode amplitudes)
//                while the total inflow stays exactly fixed
//   Stirring     a few smooth, divergence-free body-force patterns with random,
//                slowly varying amplitudes, shaped to follow the channel
// Thermal fluctuations live in the collision kernel (CollisionConfig::kT) and the
// dots' random walk lives in Tracers.

#include "lbm/Domain.hpp"
#include "lbm/Random.hpp"

#include <cstdint>
#include <functional>
#include <vector>

namespace lbm {

struct InflowNoiseParams {
    bool enabled = true;
    int modes = 3;                     // cos(m pi s), m = 1..modes
    double amplitude = 0.15;           // std of mode m is amplitude / m
    double minCorrelationTime = 6.0;   // seconds; mode m gets a time between min and max
    double maxCorrelationTime = 15.0;
    double floor = 0.3;                // smallest allowed local multiplier
    bool meander = true;               // slow transverse tilt of the inflow
    double meanderAmplitude = 0.05;    // std of transverse / streamwise speed
    double meanderCorrelationTime = 10.0;
};

// Inlet profile U(s, t) = Ubar W g(s) B(s, t) / A(t), s in [0, 1] across the inlet:
//   g(s) = 1 - (2s - 1)^6 vanishes at both banks (no-slip corners);
//   B = max(1 + sum_m a_m(t) cos(m pi s), floor);
//   A(t) makes the discrete flux over the inlet nodes exactly Ubar W.
class InflowModel {
public:
    // nodeS: cross-inlet coordinate of each inlet node; nodeWeight: the width each
    // node represents (page units); their sum is the inlet width W.
    InflowModel(const InflowNoiseParams& params, std::uint64_t seed,
                std::vector<double> nodeS, std::vector<double> nodeWeight);

    void advance(double dtSeconds);

    // Streamwise speed (along the inflow) and transverse speed at s, for mean speed ubar.
    double streamwise(double s, double ubar) const;
    double transverse(double s, double ubar) const;

    // Exact discrete flux sum_k U(s_k) w_k for mean speed ubar (diagnostic/test).
    double discreteFlux(double ubar) const;
    const InflowNoiseParams& params() const { return params_; }

private:
    double shape(double s) const;   // g(s) B(s, t)
    void renormalize();

    InflowNoiseParams params_;
    RandomStream rng_;
    std::vector<OrnsteinUhlenbeck> modes_;
    OrnsteinUhlenbeck meander_;
    std::vector<double> nodeS_, nodeWeight_;
    double width_ = 0.0;
    double area_ = 1.0;   // A(t) = sum_k g B w_k / W
};

struct StirringParams {
    bool enabled = false;
    bool periodic = false;          // smooth Fourier basis on a periodic box
    int modes = 4;
    double intensity = 0.15;          // target velocity perturbation / mean inlet speed
    double correlationTime = 4.0;     // seconds
};

// Body force F = sum_k b_k(t) curl(m psi_k), evaluated on the lattice with central
// differences of the analytic stream functions, so its discrete divergence is zero.
class Stirring {
public:
    // channel(X, Y) maps lattice coordinates to channel coordinates (xi, eta) in
    // [0, 1]^2: xi across the river from the bank, eta along it.
    using ChannelMap = std::function<void(double X, double Y, double& xi, double& eta)>;

    Stirring(const StirringParams& params, std::uint64_t seed, const Domain& domain,
             const ChannelMap& channel);

    // Advance amplitudes and rebuild the force field. forceScale is the lattice force
    // that corresponds to intensity 1 (set by the river's units).
    void advance(double dtSeconds, double forceScale);

    const std::vector<double>& fx() const { return fx_; }
    const std::vector<double>& fy() const { return fy_; }
    // Unit-rms force field of one mode (for tests).
    const std::vector<double>& modeFx(int k) const { return modeFx_[k]; }
    const std::vector<double>& modeFy(int k) const { return modeFy_[k]; }

private:
    StirringParams params_;
    RandomStream rng_;
    std::vector<OrnsteinUhlenbeck> amplitudes_;
    std::vector<std::vector<double>> modeFx_, modeFy_;
    std::vector<double> fx_, fy_;
};

}  // namespace lbm
