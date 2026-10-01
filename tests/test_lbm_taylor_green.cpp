// Taylor-Green vortex decay in a periodic box. The kinetic energy decays as
// exp(-4 nu k^2 t) (velocity as exp(-2 nu k^2 t)), which measures the viscosity the
// lattice actually produces and compares it with nu = c_s^2 (tau - 1/2).

#include "TestUtil.hpp"
#include "lbm/Solver.hpp"

#include <cmath>
#include <string>

using namespace lbm;

namespace {

double kineticEnergy(const Solver& s) {
    double e = 0.0;
    for (size_t n = 0; n < s.rho().size(); ++n)
        e += 0.5 * s.rho()[n] * (s.ux()[n] * s.ux()[n] + s.uy()[n] * s.uy()[n]);
    return e;
}

// Relative error of the measured viscosity.
double viscosityError(CollisionModel model, int N, double tau) {
    DomainSpec spec;
    spec.nx = spec.ny = N;
    spec.xMin = spec.xMax = spec.yMin = spec.yMax = SideKind::Periodic;
    Domain domain(spec);

    CollisionConfig cfg;
    cfg.model = model;
    cfg.tau = tau;
    const double rho0 = 1.0;
    Solver solver(domain, cfg, rho0);

    const double k = 2.0 * pi / N;
    const double U = 0.02 * 32.0 / N;   // diffusive scaling keeps the Mach error shrinking
    for (int n = 0; n < domain.fluidCount(); ++n) {
        const double X = domain.k(n) + 0.5, Y = domain.j(n) + 0.5;
        const double ux = -U * std::cos(k * X) * std::sin(k * Y);
        const double uy = U * std::sin(k * X) * std::cos(k * Y);
        const double p = -0.25 * rho0 * U * U * (std::cos(2 * k * X) + std::cos(2 * k * Y));
        solver.initializeNode(n, rho0 + p / cs2, ux, uy);
    }

    const double nu = viscosityFromTau(tau);
    const int total = static_cast<int>(std::log(2.0) / (4.0 * nu * k * k));
    const int t1 = total / 4;
    double e1 = 0.0;
    for (int t = 1; t <= total; ++t) {
        solver.step();
        if (t == t1) e1 = kineticEnergy(solver);
    }
    const double e2 = kineticEnergy(solver);
    const double nuMeasured = std::log(e1 / e2) / (4.0 * k * k * (total - t1));
    return std::abs(nuMeasured - nu) / nu;
}

const char* name(CollisionModel m) {
    return m == CollisionModel::BGK ? "BGK" : m == CollisionModel::TRT ? "TRT" : "MRT";
}

}  // namespace

int main() {
    for (CollisionModel model : {CollisionModel::BGK, CollisionModel::TRT, CollisionModel::MRT}) {
        const double e32 = viscosityError(model, 32, 0.8);
        const double e64 = viscosityError(model, 64, 0.8);
        lbmtest::checkRange(e64, 0.0, 0.01,
                            std::string(name(model)) + " viscosity within 1% (N=64, tau=0.8)");
        lbmtest::checkRange(e32 / std::max(e64, 1e-16), 3.0, 1e9,
                            std::string(name(model)) + " error ratio N=32/N=64 (2nd order ~ 4)");
        const double eLow = viscosityError(model, 64, 0.55);
        lbmtest::checkRange(eLow, 0.0, 0.01,
                            std::string(name(model)) + " viscosity within 1% (N=64, tau=0.55)");
    }
    return lbmtest::finish();
}
