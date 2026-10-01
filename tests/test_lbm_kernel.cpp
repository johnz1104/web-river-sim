// Algebraic checks of the lattice, moment basis and collision kernel.

#include "TestUtil.hpp"
#include "lbm/Collision.hpp"
#include "lbm/Random.hpp"

#include <algorithm>

using namespace lbm;
using lbmtest::check;
using lbmtest::checkNear;

int main() {
    // Weighted orthogonality: sum_i w_i e_k e_l = norm_k delta_kl.
    double worst = 0.0;
    for (int k = 0; k < Q; ++k)
        for (int l = 0; l < Q; ++l) {
            double s = 0.0;
            for (int i = 0; i < Q; ++i) s += w[i] * basis[k][i] * basis[l][i];
            worst = std::max(worst, std::abs(s - (k == l ? norm[k] : 0.0)));
        }
    checkNear(worst, 0.0, 1e-14, "moment basis is weighted-orthogonal");

    // Lattice isotropy: sum w = 1, sum w c_a c_b = c_s^2 delta_ab.
    double sw = 0.0, sxx = 0.0, sxy = 0.0;
    for (int i = 0; i < Q; ++i) { sw += w[i]; sxx += w[i] * cx[i] * cx[i]; sxy += w[i] * cx[i] * cy[i]; }
    check(std::abs(sw - 1.0) < 1e-15 && std::abs(sxx - cs2) < 1e-15 && std::abs(sxy) < 1e-15,
          "weights are normalised and isotropic to second order");

    RandomStream rng(7, 0);
    double eqError = 0.0, roundTrip = 0.0, tableError = 0.0, bgkError = 0.0, massError = 0.0,
           momError = 0.0;
    for (int trial = 0; trial < 200; ++trial) {
        const double rho = 1.0 + 0.05 * (rng.uniform() - 0.5);
        const double ux = 0.1 * (rng.uniform() - 0.5), uy = 0.1 * (rng.uniform() - 0.5);

        // Analytic equilibrium moments match the transform of f^eq.
        double feq[Q], m[Q];
        for (int i = 0; i < Q; ++i) feq[i] = equilibrium(i, rho, ux, uy);
        toMoments(feq, m);
        const EquilibriumMoments eq = equilibriumMoments(rho, ux, uy);
        eqError = std::max({eqError, std::abs(m[0] - rho), std::abs(m[1] - rho * ux),
                            std::abs(m[2] - rho * uy), std::abs(m[3] - eq.bulk),
                            std::abs(m[4] - eq.normalDiff), std::abs(m[5] - eq.shear),
                            std::abs(m[6]), std::abs(m[7]), std::abs(m[8])});

        // The explicit transforms agree with the basis table, and round-trip.
        double f[Q], back[Q];
        for (int i = 0; i < Q; ++i) f[i] = feq[i] * (1.0 + 0.2 * (rng.uniform() - 0.5));
        toMoments(f, m);
        for (int k = 0; k < Q; ++k) {
            double s = 0.0;
            for (int i = 0; i < Q; ++i) s += basis[k][i] * f[i];
            tableError = std::max(tableError, std::abs(s - m[k]));
        }
        fromMoments(m, back);
        for (int i = 0; i < Q; ++i) {
            double s = 0.0;
            for (int k = 0; k < Q; ++k) s += basis[k][i] * m[k] / norm[k];
            tableError = std::max(tableError, std::abs(w[i] * s - back[i]));
            roundTrip = std::max(roundTrip, std::abs(back[i] - f[i]));
        }

        // MRT kernel with BGK rates equals BGK written on populations (with Guo forcing).
        const double tau = 0.51 + 1.5 * rng.uniform();
        const double Fx = 1e-4 * (rng.uniform() - 0.5), Fy = 1e-4 * (rng.uniform() - 0.5);
        CollisionConfig cfg;
        cfg.model = CollisionModel::BGK;
        cfg.tau = tau;
        double rates[Q];
        relaxationRates(cfg, tau, rates);
        const double rho0 = 1.0;
        double dev[Q], ref[Q];
        for (int i = 0; i < Q; ++i) { dev[i] = f[i] - w[i] * rho0; ref[i] = f[i]; }
        collide(dev, rho0, Fx, Fy, cfg, rates, -1.0, nullptr);
        collideBGKPopulations(ref, tau, Fx, Fy);
        for (int i = 0; i < Q; ++i)
            bgkError = std::max(bgkError, std::abs(dev[i] + w[i] * rho0 - ref[i]));

        // MRT with free rates conserves mass exactly and gains exactly F in momentum.
        cfg.model = CollisionModel::MRT;
        relaxationRates(cfg, tau, rates);
        double before[3] = {0, 0, 0}, after[3] = {0, 0, 0};
        for (int i = 0; i < Q; ++i) {
            dev[i] = f[i] - w[i] * rho0;
            before[0] += dev[i]; before[1] += dev[i] * cx[i]; before[2] += dev[i] * cy[i];
        }
        collide(dev, rho0, Fx, Fy, cfg, rates, -1.0, nullptr);
        for (int i = 0; i < Q; ++i) {
            after[0] += dev[i]; after[1] += dev[i] * cx[i]; after[2] += dev[i] * cy[i];
        }
        massError = std::max(massError, std::abs(after[0] - before[0]));
        momError = std::max({momError, std::abs(after[1] - before[1] - Fx),
                             std::abs(after[2] - before[2] - Fy)});
    }
    checkNear(eqError, 0.0, 1e-14, "analytic equilibrium moments match E f^eq");
    checkNear(tableError, 0.0, 1e-14, "explicit moment transforms match the basis table");
    checkNear(roundTrip, 0.0, 1e-14, "moment transform round trip");
    checkNear(bgkError, 0.0, 1e-14, "MRT kernel with BGK rates equals population-space BGK");
    checkNear(massError, 0.0, 1e-15, "MRT collision conserves mass");
    checkNear(momError, 0.0, 1e-15, "MRT collision adds exactly F to momentum");

    // TRT magic parameter: (tau - 1/2)(1/s- - 1/2) = Lambda.
    CollisionConfig trt;
    trt.model = CollisionModel::TRT;
    double rates[Q];
    relaxationRates(trt, 0.73, rates);
    checkNear((0.73 - 0.5) * (1.0 / rates[6] - 0.5), 3.0 / 16.0, 1e-14,
              "TRT odd rate satisfies the magic parameter");

    // Smagorinsky closed form satisfies its defining quadratic:
    // nu_eff = nu_0 + C^2 |S|, with |S| = sqrt(2) |Pi| / (2 rho c_s^2 tau_eff).
    const double tau0 = 0.52, C = 0.15, rho = 1.01, Pxx = 1e-4, Pyy = -3e-5, Pxy = 2e-4;
    const double te = smagorinskyTau(tau0, C, rho, Pxx, Pyy, Pxy);
    const double piNorm = std::sqrt(Pxx * Pxx + Pyy * Pyy + 2 * Pxy * Pxy);
    const double S = std::sqrt(2.0) * piNorm / (2.0 * rho * cs2 * te);
    checkNear(viscosityFromTau(te), viscosityFromTau(tau0) + C * C * S, 1e-15,
              "Smagorinsky tau_eff solves nu_eff = nu_0 + C^2 |S|");

    return lbmtest::finish();
}
