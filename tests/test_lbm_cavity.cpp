// Lid-driven cavity: u along the vertical centreline against Ghia, Ghia & Shin (1982),
// Table I, at Re = 100 and 400. Exercises the moving-wall rule and corner handling.

#include "TestUtil.hpp"
#include "lbm/Solver.hpp"

#include <cmath>
#include <string>

using namespace lbm;

namespace {

// (y, u/U) from Ghia et al. (1982), Table I; y measured from the stationary bottom.
const double ghiaY[] = {0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531,
                        0.5000, 0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766};
const double ghia100[] = {-0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662,
                          -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717,
                          0.73722, 0.78871, 0.84123};
const double ghia400[] = {-0.08186, -0.09266, -0.10338, -0.14612, -0.24299, -0.32726,
                          -0.17119, -0.11477, 0.02135, 0.16256, 0.29093, 0.55892,
                          0.61756, 0.68439, 0.75837};

double centrelineError(double Re, const double* reference, int N) {
    // Lattice Y runs from the bottom (YMin) to the lid (YMax).
    DomainSpec spec;
    spec.nx = spec.ny = N;
    spec.yMax = SideKind::MovingWall;
    Domain domain(spec);

    const double U = 0.1;
    CollisionConfig cfg;
    cfg.tau = tauFromViscosity(U * N / Re);
    Solver solver(domain, cfg, 1.0);
    const auto& links = domain.links();
    for (size_t li = 0; li < links.size(); ++li)
        if (links[li].kind == LinkKind::MovingWall)
            solver.setWallVelocity(static_cast<int>(li), U, 0.0);

    const long steps = static_cast<long>(60.0 * N / U);   // ~60 lid passes
    for (long t = 0; t < steps; ++t) solver.step();

    // Centreline x = N/2 lies between columns N/2 - 1 and N/2.
    auto uAt = [&](double y) {
        const double Y = y * N;
        int j0 = static_cast<int>(std::floor(Y - 0.5));
        j0 = std::max(0, std::min(N - 2, j0));
        const double t = (Y - 0.5) - j0;
        auto column = [&](int j) {
            return 0.5 * (solver.ux()[domain.nodeAt(N / 2 - 1, j)]
                          + solver.ux()[domain.nodeAt(N / 2, j)]);
        };
        return ((1.0 - t) * column(j0) + t * column(j0 + 1)) / U;
    };
    double worst = 0.0;
    for (int p = 0; p < 15; ++p) worst = std::max(worst, std::abs(uAt(ghiaY[p]) - reference[p]));
    return worst;
}

}  // namespace

int main() {
    lbmtest::checkRange(centrelineError(100.0, ghia100, 128), 0.0, 0.02,
                        "Re 100 centreline u vs Ghia (max abs error / U), N=128");
    lbmtest::checkRange(centrelineError(400.0, ghia400, 128), 0.0, 0.02,
                        "Re 400 centreline u vs Ghia (max abs error / U), N=128");
    return lbmtest::finish();
}
