// Body-force-driven Poiseuille flow between two halfway bounce-back walls.
// Exact solution: u(Y) = F Y (H - Y) / (2 rho nu), walls at Y = 0 and Y = H.
// TRT/MRT with Lambda = 3/16 reproduce it to round-off; BGK converges at 2nd order.

#include "TestUtil.hpp"
#include "lbm/Solver.hpp"

#include <cmath>
#include <string>

using namespace lbm;

namespace {

double relativeError(CollisionModel model, int H, double tau) {
    DomainSpec spec;
    spec.nx = 4;
    spec.ny = H;
    spec.xMin = spec.xMax = SideKind::Periodic;
    spec.yMin = spec.yMax = SideKind::Wall;
    Domain domain(spec);

    CollisionConfig cfg;
    cfg.model = model;
    cfg.tau = tau;
    Solver solver(domain, cfg, 1.0);
    const double nu = viscosityFromTau(tau);
    const double uMax = 0.05;
    const double F = 8.0 * nu * uMax / (H * H);
    solver.setUniformForce(F, 0.0);

    const int steps = static_cast<int>(40.0 * H * H / nu);
    for (int t = 0; t < steps; ++t) solver.step();

    double num = 0.0, den = 0.0;
    for (int n = 0; n < domain.fluidCount(); ++n) {
        const double Y = domain.j(n) + 0.5;
        const double exact = F * Y * (H - Y) / (2.0 * nu);
        num += (solver.ux()[n] - exact) * (solver.ux()[n] - exact);
        den += exact * exact;
    }
    return std::sqrt(num / den);
}

}  // namespace

int main() {
    lbmtest::checkRange(relativeError(CollisionModel::TRT, 16, 0.8), 0.0, 1e-9,
                        "TRT (Lambda=3/16) Poiseuille exact, tau=0.8");
    lbmtest::checkRange(relativeError(CollisionModel::TRT, 16, 1.7), 0.0, 1e-9,
                        "TRT (Lambda=3/16) Poiseuille exact, tau=1.7");
    // MRT is exact only when its free bulk/even-ghost rates equal 1/tau (it is then TRT);
    // with the default free rates the error is small and must shrink with resolution.
    const double m8 = relativeError(CollisionModel::MRT, 8, 0.8);
    const double m16 = relativeError(CollisionModel::MRT, 16, 0.8);
    const double m32 = relativeError(CollisionModel::MRT, 32, 0.8);
    std::printf("MRT tau=0.8 errors: H=8 %.3e, H=16 %.3e, H=32 %.3e\n", m8, m16, m32);
    lbmtest::checkRange(m16, 0.0, 1e-5, "MRT (free ghost rates) Poiseuille error, H=16");
    lbmtest::checkRange(m16 / m32, 3.5, 1e9, "MRT error ratio H=16/H=32 (at least 2nd order)");

    const double b8 = relativeError(CollisionModel::BGK, 8, 1.2);
    const double b16 = relativeError(CollisionModel::BGK, 16, 1.2);
    const double b32 = relativeError(CollisionModel::BGK, 32, 1.2);
    std::printf("BGK tau=1.2 errors: H=8 %.3e, H=16 %.3e, H=32 %.3e\n", b8, b16, b32);
    lbmtest::checkRange(b8 / b16, 3.5, 4.5, "BGK error ratio H=8/H=16 (2nd order)");
    lbmtest::checkRange(b16 / b32, 3.5, 4.5, "BGK error ratio H=16/H=32 (2nd order)");
    return lbmtest::finish();
}
