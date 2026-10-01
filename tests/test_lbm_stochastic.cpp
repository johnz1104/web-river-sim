// Checks of the four random ingredients:
//   Ornstein-Uhlenbeck statistics; exact inflow flux under noise; divergence-free
//   stirring; thermal equipartition (fluctuating LBM); tracer random-walk diffusivity.

#include "TestUtil.hpp"
#include "HomeBank.hpp"
#include "lbm/Solver.hpp"
#include "lbm/Stochastic.hpp"
#include "lbm/Tracers.hpp"

#include <cmath>
#include <string>
#include <vector>

using namespace lbm;

namespace {

void ornsteinUhlenbeck() {
    RandomStream rng(11, 0);
    OrnsteinUhlenbeck ou;
    ou.sigma = 0.7;
    ou.correlationTime = 3.0;
    const double dt = 0.05;
    const int n = 400000, lag = static_cast<int>(ou.correlationTime / dt);
    std::vector<double> a(n);
    for (int t = 0; t < n; ++t) { ou.step(dt, rng.normal()); a[t] = ou.value; }
    double mean = 0.0, var = 0.0, cov = 0.0;
    for (double v : a) mean += v;
    mean /= n;
    for (int t = 0; t < n; ++t) var += (a[t] - mean) * (a[t] - mean);
    for (int t = 0; t + lag < n; ++t) cov += (a[t] - mean) * (a[t + lag] - mean);
    var /= n;
    cov /= (n - lag);
    lbmtest::checkNear(std::sqrt(var), ou.sigma, 0.03 * ou.sigma, "OU stationary std");
    lbmtest::checkNear(cov / var, std::exp(-1.0), 0.03, "OU autocorrelation at lag T");
}

void inflowFlux() {
    std::vector<double> s, width;
    for (int k = 0; k < 126; ++k) { s.push_back((k + 0.3) / 126.0); width.push_back(8.0); }
    InflowNoiseParams p;
    p.amplitude = 0.4;   // stronger than the default, to exercise the floor
    InflowModel model(p, 5, s, width);
    const double ubar = 30.0, target = ubar * 126 * 8.0;
    double worst = 0.0, minSpeed = 1e9;
    for (int t = 0; t < 20000; ++t) {
        model.advance(0.004);
        worst = std::max(worst, std::abs(model.discreteFlux(ubar) - target) / target);
        for (double sk : s) minSpeed = std::min(minSpeed, model.streamwise(sk, ubar));
    }
    lbmtest::checkRange(worst, 0.0, 1e-12, "inflow flux stays exact under noise (relative)");
    lbmtest::checkRange(minSpeed, 0.0, 1e9, "inflow never reverses");
}

void stirringDivergence() {
    DomainSpec spec;
    spec.nx = 60;
    spec.ny = 80;
    Domain domain(spec);
    StirringParams p;
    p.enabled = true;
    p.modes = 6;
    Stirring stirring(p, 3, domain, [&](double X, double Y, double& xi, double& eta) {
        xi = X / spec.nx;
        eta = Y / spec.ny;
    });
    stirring.advance(0.01, 1.0);
    double worst = 0.0, rms = 0.0;
    for (int n = 0; n < domain.fluidCount(); ++n) {
        const auto& fx = stirring.fx();
        const auto& fy = stirring.fy();
        rms += fx[n] * fx[n] + fy[n] * fy[n];
        const int k = domain.k(n), j = domain.j(n);
        if (k < 1 || j < 1 || k > spec.nx - 2 || j > spec.ny - 2) continue;
        const double div = 0.5 * (fx[domain.nodeAt(k + 1, j)] - fx[domain.nodeAt(k - 1, j)])
                         + 0.5 * (fy[domain.nodeAt(k, j + 1)] - fy[domain.nodeAt(k, j - 1)]);
        worst = std::max(worst, std::abs(div));
    }
    rms = std::sqrt(rms / domain.fluidCount());
    lbmtest::checkRange(worst / rms, 0.0, 1e-12, "stirring force has zero discrete divergence");
}

void equipartition(CollisionModel model, const std::string& label) {
    DomainSpec spec;
    spec.nx = spec.ny = 32;
    spec.xMin = spec.xMax = spec.yMin = spec.yMax = SideKind::Periodic;
    Domain domain(spec);
    CollisionConfig cfg;
    cfg.model = model;
    cfg.tau = 0.8;
    cfg.kT = 1e-4;
    Solver solver(domain, cfg, 1.0, 99);
    for (int t = 0; t < 2000; ++t) solver.step();
    double energy = 0.0, density = 0.0;
    long samples = 0;
    for (int t = 0; t < 20000; ++t) {
        solver.step();
        if (t % 10) continue;
        for (int n = 0; n < domain.fluidCount(); ++n) {
            const double r = solver.rho()[n];
            energy += 0.5 * r * (solver.ux()[n] * solver.ux()[n] + solver.uy()[n] * solver.uy()[n]);
            density += (r - 1.0) * (r - 1.0);
            ++samples;
        }
    }
    // Per node: <rho u_x^2> = <rho u_y^2> = kT, and <drho^2> = rho kT / c_s^2.
    lbmtest::checkNear(energy / samples, cfg.kT, 0.05 * cfg.kT,
                       label + " thermal equipartition <rho u^2 / 2> = kT");
    lbmtest::checkNear(density / samples, cfg.kT / cs2, 0.05 * cfg.kT / cs2,
                       label + " density fluctuations <drho^2> = rho kT / c_s^2");
}

void randomWalk() {
    BankSpline bank(homeBankSegments());
    TracerParams p;
    p.diffusivity = 2.0;
    p.maxCount = 100000;
    Tracers tracers(p, 17, bank, kHomeWidth, kHomeHeight);
    const int n = 20000;
    const double x0 = 1250.0, y0 = 600.0;   // far (> 150 units) from every wall
    for (int i = 0; i < n; ++i) tracers.add(x0, y0);
    const double dt = 1.0 / 60.0, T = 10.0;
    for (int f = 0; f < static_cast<int>(T / dt); ++f)
        tracers.advance(dt, [](double, double, double& u, double& v) { u = v = 0.0; });
    double msdX = 0.0, msdY = 0.0;
    for (int i = 0; i < tracers.count(); ++i) {
        msdX += (tracers.x(i) - x0) * (tracers.x(i) - x0);
        msdY += (tracers.y(i) - y0) * (tracers.y(i) - y0);
    }
    msdX /= tracers.count();
    msdY /= tracers.count();
    const double expected = 2.0 * p.diffusivity * T;
    lbmtest::checkNear(msdX, expected, 0.03 * expected, "dot random walk: <dx^2> = 2 D t");
    lbmtest::checkNear(msdY, expected, 0.03 * expected, "dot random walk: <dy^2> = 2 D t");
}

}  // namespace

int main() {
    ornsteinUhlenbeck();
    inflowFlux();
    stirringDivergence();
    equipartition(CollisionModel::BGK, "BGK");
    equipartition(CollisionModel::MRT, "MRT");
    randomWalk();
    return lbmtest::finish();
}
