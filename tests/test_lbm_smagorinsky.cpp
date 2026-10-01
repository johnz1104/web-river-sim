// Smagorinsky LES-LBM checks.
// 1. Consistency: in a steady channel flow, the strain rate the closure infers from the
//    non-equilibrium stress, S = -Pi^neq / (2 rho c_s^2 tau_eff), matches the strain
//    from finite differences of the velocity profile.
// 2. Stability: a thin double shear layer at tau0 = 0.5 + 1e-4 stays bounded with the
//    closure; the result without it is reported for comparison.

#include "TestUtil.hpp"
#include "lbm/Solver.hpp"

#include <cmath>
#include <string>

using namespace lbm;

namespace {

void strainConsistency() {
    // Low molecular viscosity and a large constant so the eddy viscosity is a sizeable
    // fraction (~20%) of the total; the run covers a few diffusion times H^2 / nu.
    const int H = 32;
    DomainSpec spec;
    spec.nx = 4;
    spec.ny = H;
    spec.xMin = spec.xMax = SideKind::Periodic;
    Domain domain(spec);
    CollisionConfig cfg;
    cfg.tau = 0.505;
    cfg.closure = ClosureKind::Smagorinsky;
    cfg.smagorinsky = 0.3;
    cfg.oddGhostRate = 1.2;   // as in the river: damped ghosts near tau = 1/2
    Solver solver(domain, cfg, 1.0);
    const double F = 8.0 * viscosityFromTau(cfg.tau) * 0.03 / (H * H);
    solver.setUniformForce(F, 0.0);
    for (int t = 0; t < 600000; ++t) solver.step();

    double worst = 0.0, scale = 0.0;
    const int k = 1;
    for (int j = 4; j < H - 4; ++j) {
        const int n = domain.nodeAt(k, j);
        // Pre-collision populations at n: pull the post-collision values (interior node).
        double f[Q], rho = 0.0, jx = 0.0, jy = 0.0;
        for (int i = 0; i < Q; ++i) {
            f[i] = solver.population(domain.pull(n, i), i);
            rho += f[i]; jx += f[i] * cx[i]; jy += f[i] * cy[i];
        }
        const double ux = (jx + 0.5 * F) / rho, uy = jy / rho;
        double Pxx = 0.0, Pyy = 0.0, Pxy = 0.0;
        for (int i = 0; i < Q; ++i) {
            const double neq = f[i] - equilibrium(i, rho, ux, uy);
            Pxx += cx[i] * cx[i] * neq; Pyy += cy[i] * cy[i] * neq; Pxy += cx[i] * cy[i] * neq;
        }
        Pxx += ux * F;
        Pxy += 0.5 * F * uy;
        const double tauEff = smagorinskyTau(cfg.tau, cfg.smagorinsky, rho, Pxx, Pyy, Pxy);
        const double sFromStress = -Pxy / (2.0 * rho * cs2 * tauEff);
        const double sFromProfile = 0.25 * (solver.ux()[domain.nodeAt(k, j + 1)]
                                            - solver.ux()[domain.nodeAt(k, j - 1)]);
        if (j % 6 == 0)
            std::printf("  j=%2d u=%.5f  S(stress)=%.4e  S(profile)=%.4e  tau_eff=%.5f\n", j,
                        solver.ux()[n], sFromStress, sFromProfile, tauEff);
        worst = std::max(worst, std::abs(sFromStress - sFromProfile));
        scale = std::max(scale, std::abs(sFromProfile));
    }
    lbmtest::checkRange(worst / scale, 0.0, 0.02,
                        "strain from Pi^neq matches finite-difference strain (relative)");
}

// Returns the largest |u| after the run, or infinity if it blew up.
double shearLayer(ClosureKind closure) {
    const int N = 128;
    DomainSpec spec;
    spec.nx = spec.ny = N;
    spec.xMin = spec.xMax = spec.yMin = spec.yMax = SideKind::Periodic;
    Domain domain(spec);
    CollisionConfig cfg;
    cfg.model = CollisionModel::BGK;
    cfg.tau = 0.5001;
    cfg.closure = closure;
    cfg.smagorinsky = 0.17;
    Solver solver(domain, cfg, 1.0);
    const double U = 0.08, kappa = 80.0, delta = 0.05;
    for (int n = 0; n < domain.fluidCount(); ++n) {
        const double x = (domain.k(n) + 0.5) / N, y = (domain.j(n) + 0.5) / N;
        const double ux = U * (y <= 0.5 ? std::tanh(kappa * (y - 0.25)) : std::tanh(kappa * (0.75 - y)));
        const double uy = U * delta * std::sin(2.0 * pi * (x + 0.25));
        solver.initializeNode(n, 1.0, ux, uy);
    }
    double peak = 0.0;
    for (int t = 0; t < 20000; ++t) {
        solver.step();
        if (t % 500 == 0) {
            peak = 0.0;
            for (int n = 0; n < domain.fluidCount(); ++n) {
                const double s = std::hypot(solver.ux()[n], solver.uy()[n]);
                if (!std::isfinite(s)) return INFINITY;
                peak = std::max(peak, s);
            }
            if (peak > 1.0) return INFINITY;
        }
    }
    return peak;
}

}  // namespace

int main() {
    strainConsistency();
    const double laminar = shearLayer(ClosureKind::Laminar);
    const double les = shearLayer(ClosureKind::Smagorinsky);
    std::printf("double shear layer, BGK tau0 = 0.5001: laminar peak |u| = %g, "
                "Smagorinsky peak |u| = %g\n", laminar, les);
    lbmtest::checkRange(les, 0.0, 0.2, "Smagorinsky keeps the thin shear layer bounded");
    return lbmtest::finish();
}
