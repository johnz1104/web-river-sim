#pragma once

// Lattice Boltzmann time stepping on a Domain. One step() is a fused
// "pull, apply boundary rules, collide" sweep over the fluid nodes:
//   f_i(x, t+1) = f_i*(x - c_i, t)            interior links (streaming)
//   f_i(x, t+1) = boundary rule(link)          links that cross a boundary
//   f*(x, t+1)  = collide(f(x, t+1))           moment-space collision
// Populations are stored post-collision as deviations f_i - w_i rho0, which keeps
// single-precision storage accurate (only small numbers are stored).
//
// Boundary rules (link-wise; o = outgoing direction, i = opposite(o), q = wall fraction):
//   Wall            Bouzidi-Firdaouss-Lallemand interpolated bounce-back; q = 1/2 is
//                   ordinary halfway bounce-back (no-slip)
//   MovingWall      f_i = f_o* - 6 w_o rho (c_o . u_w)          (Ladd; inlets and lids)
//   PressureOutlet  f_i = -f_o* + 2 w_o rho_out [1 + 4.5 (c_o . u_w)^2 - 1.5 u_w^2]
//                   with u_w extrapolated from the two nearest nodes (anti-bounce-back)

#include "lbm/Collision.hpp"
#include "lbm/Domain.hpp"

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

namespace lbm {

class Solver {
public:
    Solver(const Domain& domain, const CollisionConfig& config, double rho0 = 1.0,
           std::uint64_t seed = 0);
    ~Solver();
    // Native fixed worker pool. Emscripten always executes on the calling thread.
    void setThreadCount(int count);
    int threadCount() const { return threadCount_; }
    void captureViscousFlux(bool enabled) {
        viscousFlux_.resize(enabled ? domain_.fluidCount() : 0);
    }
    const std::vector<std::array<double,3>>& viscousFlux() const { return viscousFlux_; }
    std::array<double,2> force(int n) const {
        return {uniformFx_+(fieldFx_ ? (*fieldFx_)[n] : 0.),
                uniformFy_+(fieldFy_ ? (*fieldFy_)[n] : 0.)};
    }
    double wallDistance(int n) const { return wallDistance_.empty() ? -1. : wallDistance_[n]; }
    // Research/test initialization: full post-collision populations.
    void initializePopulations(int n, const double* populations);

    // Equilibrium populations at a uniform or per-node state (lattice units).
    void initializeUniform(double rho, double ux, double uy);
    void initializeNode(int n, double rho, double ux, double uy);

    void step();

    // ---- configuration ----
    CollisionConfig& config() { return config_; }
    void refreshRates();   // call after changing config().tau or the model
    void setUniformForce(double Fx, double Fy) { uniformFx_ = Fx; uniformFy_ = Fy; }
    // Per-node force field (lattice units); pass nullptrs to remove it.
    void setForceField(const std::vector<double>* Fx, const std::vector<double>* Fy) {
        fieldFx_ = Fx; fieldFy_ = Fy;
    }
    // Per-node linear drag F = -k rho u, k in lattice units per step (a hidden bed depth
    // in a 2D river). Applied implicitly, so any k in [0, 1] is stable; pass nullptr to
    // remove it.
    void setDragField(const std::vector<double>* k) { dragField_ = k; }
    void setWallVelocity(int link, double ux, double uy) { wallVelocity_[link] = {ux, uy}; }
    void setOutletDensity(double rho) { rhoOut_ = rho; }
    // Backflow stabilization: where water would flow in through a pressure outlet, the
    // outlet rule uses only the tangential velocity. Off by default.
    void setOutletBackflowClamp(bool enabled) { outletBackflowClamp_ = enabled; }
    void setWallDistance(std::vector<double> d) { wallDistance_ = std::move(d); }

    // ---- state (describes the latest completed step) ----
    const Domain& domain() const { return domain_; }
    std::int64_t time() const { return time_; }
    double rho0() const { return rho0_; }
    const std::vector<double>& rho() const { return rho_; }
    const std::vector<double>& ux() const { return ux_; }
    const std::vector<double>& uy() const { return uy_; }
    double maxTauEffective() const { return maxTauEff_; }
    // Momentum-exchange force on a body id during the last step.
    std::array<double, 2> bodyForce(int body) const;
    // Full population f_i at node n (post-collision), for tests and diagnostics.
    double population(int n, int i) const { return f_[n * Q + i] + w[i] * rho0_; }

private:
    struct WorkerPool;
    std::unique_ptr<WorkerPool> pool_;
    int threadCount_ = 1;
    void stepRange(int begin, int end, int worker);
    std::vector<double> partialMaxTau_;
    std::vector<std::array<double,6>> learnedNoise_;
    double noiseDecay_=0.,noiseInnovation_=1.;
    std::vector<std::array<double,3>> viscousFlux_;
    // Link contributions are reduced in the original node/direction order.
    std::vector<double> exchange_;
    std::vector<int> bodyLinks_;
    double incoming(int n, int i, const BoundaryLink& link, int linkIndex) const;

    const Domain& domain_;
    CollisionConfig config_;
    double rho0_;
    std::uint64_t seed_;
    double rhoOut_;
    bool outletBackflowClamp_ = false;
    double baseRates_[Q];

    std::vector<real> f_, fNext_;
    std::vector<double> rho_, ux_, uy_;            // current step
    std::vector<double> rhoPrev_, uxPrev_, uyPrev_; // previous step (read by boundary rules)
    std::vector<std::array<double, 2>> wallVelocity_;
    std::vector<double> wallDistance_;
    const std::vector<double>* fieldFx_ = nullptr;
    const std::vector<double>* fieldFy_ = nullptr;
    const std::vector<double>* dragField_ = nullptr;
    double uniformFx_ = 0.0, uniformFy_ = 0.0;
    std::vector<std::array<double, 2>> bodyForce_;
    std::int64_t time_ = 0;
    double maxTauEff_ = 0.0;
};

}  // namespace lbm
