// Schaefer & Turek (1996) benchmarks 2D-1 and 2D-2: a cylinder in a channel with a
// parabolic inflow and a pressure outlet. This exercises every boundary rule the river
// uses: halfway walls, Bouzidi curved walls, a moving-wall velocity inlet and the
// anti-bounce-back pressure outlet.
//   Channel 2.2 x 0.41, cylinder D = 0.1 centred at (0.2, 0.2), inflow
//   U(y) = 4 Um y (H - y) / H^2, mean U = 2 Um / 3, Re = U D / nu.
//   2D-1 (Re 20, steady):   C_D in [5.57, 5.59], C_L in [0.0104, 0.0110],
//                           dp = p(0.15, 0.2) - p(0.25, 0.2) in [0.1172, 0.1176]
//   2D-2 (Re 100, periodic): C_D,max in [3.22, 3.24], C_L,max in [0.99, 1.01],
//                           St in [0.295, 0.305]
//
// Usage: test_lbm_cylinder [D=20] [both|steady|periodic] [regression|strict]
//                          [lattice-speed-scale=1] [odd-ghost-rate=0]
// The default "regression" mode (D = 20, used by CTest) checks tolerances sized to that
// resolution: measured there, drag carries a ~1% compressibility (Ma^2) bias and the
// peak lift a ~10% resolution error. Both shrink as expected: halving the lattice speed
// brings C_D,max into its range; D = 40 brings C_L,max into its range. "strict" checks
// the published ranges and needs D >= 40 with a reduced lattice speed (about 0.5).

#include "TestUtil.hpp"
#include "lbm/Solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace lbm;

namespace {

bool strict = false;

// Strict: the published range. Regression: within `tolerance` of the range centre.
void checkBenchmark(double value, double lo, double hi, double tolerance, const std::string& what) {
    if (strict) {
        lbmtest::checkRange(value, lo, hi, what + " (benchmark range)");
    } else {
        const double centre = 0.5 * (lo + hi);
        lbmtest::checkRange(value, centre * (1.0 - tolerance), centre * (1.0 + tolerance),
                            what + " (regression tolerance)");
    }
}

struct Channel {
    int D;           // cells per diameter
    int nx, ny;
    double cx0, cy0, radius;
};

Channel makeChannel(int D) {
    return {D, 22 * D, static_cast<int>(std::lround(4.1 * D)), 2.0 * D, 2.0 * D, 0.5 * D};
}

DomainSpec channelSpec(const Channel& c) {
    DomainSpec spec;
    spec.nx = c.nx;
    spec.ny = c.ny;
    spec.phi = [c](double X, double Y) {
        return std::hypot(X - c.cx0, Y - c.cy0) - c.radius;
    };
    spec.bodyAt = [](double, double) { return 1; };
    spec.xMin = SideKind::MovingWall;      // inlet
    spec.xMax = SideKind::PressureOutlet;  // outlet
    spec.yMin = spec.yMax = SideKind::Wall;
    return spec;
}

double inflow(double Y, double H, double Um) { return 4.0 * Um * Y * (H - Y) / (H * H); }

void setInletAndInitialize(Solver& solver, const Domain& domain, const Channel& c, double Um) {
    const auto& links = domain.links();
    for (size_t li = 0; li < links.size(); ++li) {
        const BoundaryLink& L = links[li];
        if (L.kind != LinkKind::MovingWall) continue;
        const double Y = domain.j(L.node) + 0.5 + 0.5 * cy[L.dir];
        solver.setWallVelocity(static_cast<int>(li), inflow(Y, c.ny, Um), 0.0);
    }
    for (int n = 0; n < domain.fluidCount(); ++n)
        solver.initializeNode(n, 1.0, inflow(domain.j(n) + 0.5, c.ny, Um), 0.0);
}

// Density at (X, Y), interpolated in Y between rows and extrapolated linearly in X from
// the two nearest fluid nodes on the side given by `dir` (-1 upstream, +1 downstream).
double densityAtSurface(const Solver& s, const Domain& d, double X, double Y, int dir) {
    const int j0 = static_cast<int>(std::floor(Y - 0.5));
    double value = 0.0;
    for (int dj = 0; dj < 2; ++dj) {
        const int j = j0 + dj;
        const double wy = dj == 0 ? 1.0 - (Y - 0.5 - j0) : (Y - 0.5 - j0);
        int k = static_cast<int>(std::floor(X - 0.5));
        if (dir > 0) ++k;
        while (d.nodeAt(k, j) < 0) k += dir;   // first fluid node outward
        const int n1 = d.nodeAt(k, j), n2 = d.nodeAt(k + dir, j);
        const double X1 = k + 0.5, X2 = k + dir + 0.5;
        const double r1 = s.rho()[n1], r2 = s.rho()[n2];
        value += wy * (r1 + (r2 - r1) * (X - X1) / (X2 - X1));
    }
    return value;
}

void steadyCase(int D, double oddGhostRate, double speedScale) {
    const Channel c = makeChannel(D);
    Domain domain(channelSpec(c));
    const double Um = 0.05 * speedScale, U = 2.0 * Um / 3.0;
    CollisionConfig cfg;
    cfg.oddGhostRate = oddGhostRate;
    cfg.tau = tauFromViscosity(U * D / 20.0);
    Solver solver(domain, cfg, 1.0);
    setInletAndInitialize(solver, domain, c, Um);

    const auto t0 = std::chrono::steady_clock::now();
    double cd = 0.0, cl = 0.0, previous = 0.0;
    long steps = 0;
    const long maxSteps = static_cast<long>(12.0 * c.nx / U);
    while (steps < maxSteps) {
        solver.step();
        ++steps;
        if (steps % 2000 == 0) {
            const auto F = solver.bodyForce(1);
            cd = 2.0 * F[0] / (U * U * D);
            cl = 2.0 * F[1] / (U * U * D);
            if (steps > 4 * c.nx / U && std::abs(cd - previous) < 1e-6) break;
            previous = cd;
        }
    }
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const double pFront = densityAtSurface(solver, domain, c.cx0 - c.radius, c.cy0, -1);
    const double pBack = densityAtSurface(solver, domain, c.cx0 + c.radius, c.cy0, +1);
    // dp in benchmark units: lattice dp / (U_lat^2) * U_phys^2 with U_phys = 0.2.
    const double dp = cs2 * (pFront - pBack) / (U * U) * 0.04;
    std::printf("2D-1, D=%d: %ld steps, %.1f s, %.1f MLUPS, tau=%.4f\n", D, steps, seconds,
                steps * double(domain.fluidCount()) / seconds * 1e-6, cfg.tau);
    checkBenchmark(cd, 5.57, 5.59, 0.015, "2D-1 drag coefficient C_D");
    checkBenchmark(cl, 0.0104, 0.0110, 0.03, "2D-1 lift coefficient C_L");
    checkBenchmark(dp, 0.1172, 0.1176, 0.02, "2D-1 pressure difference dp");
}

void periodicCase(int D, double oddGhostRate, double speedScale) {
    const Channel c = makeChannel(D);
    Domain domain(channelSpec(c));
    const double Um = 0.075 * speedScale, U = 2.0 * Um / 3.0;
    CollisionConfig cfg;
    cfg.oddGhostRate = oddGhostRate;
    cfg.tau = tauFromViscosity(U * D / 100.0);
    Solver solver(domain, cfg, 1.0);
    setInletAndInitialize(solver, domain, c, Um);

    // A brief transverse kick breaks the symmetry so shedding starts promptly.
    solver.setUniformForce(0.0, 1e-5);
    for (int t = 0; t < static_cast<int>(D / U); ++t) solver.step();
    solver.setUniformForce(0.0, 0.0);

    // Develop, then record several shedding periods.
    const long develop = static_cast<long>(6.0 * c.nx / U);
    for (long t = 0; t < develop; ++t) solver.step();
    const long record = static_cast<long>(12.0 * D / (0.3 * U));   // about 12 periods
    std::vector<double> cdSeries, clSeries;
    for (long t = 0; t < record; ++t) {
        solver.step();
        const auto F = solver.bodyForce(1);
        cdSeries.push_back(2.0 * F[0] / (U * U * D));
        clSeries.push_back(2.0 * F[1] / (U * U * D));
    }
    // Upward zero crossings of C_L give the period.
    std::vector<double> crossings;
    for (size_t t = 1; t < clSeries.size(); ++t)
        if (clSeries[t - 1] < 0.0 && clSeries[t] >= 0.0)
            crossings.push_back(t - 1 + clSeries[t - 1] / (clSeries[t - 1] - clSeries[t]));
    double period = 0.0;
    if (crossings.size() >= 3)
        period = (crossings.back() - crossings.front()) / (crossings.size() - 1);
    const double st = period > 0 ? D / (period * U) : 0.0;
    const double cdMax = *std::max_element(cdSeries.begin(), cdSeries.end());
    const double clMax = *std::max_element(clSeries.begin(), clSeries.end());
    std::printf("2D-2, D=%d: tau=%.4f, %zu periods recorded\n", D, cfg.tau,
                crossings.size() > 0 ? crossings.size() - 1 : 0);
    checkBenchmark(st, 0.295, 0.305, 0.017, "2D-2 Strouhal number");
    checkBenchmark(cdMax, 3.22, 3.24, 0.02, "2D-2 maximum drag coefficient");
    checkBenchmark(clMax, 0.99, 1.01, 0.12, "2D-2 maximum lift coefficient");
}

}  // namespace

int main(int argc, char** argv) {
    const int D = argc > 1 ? std::atoi(argv[1]) : 20;
    const std::string which = argc > 2 ? argv[2] : "both";
    strict = argc > 3 && std::string(argv[3]) == "strict";
    const double speedScale = argc > 4 ? std::atof(argv[4]) : 1.0;   // lattice inflow speed
    const double oddGhostRate = argc > 5 ? std::atof(argv[5]) : 0.0;
    if (which != "periodic") steadyCase(D, oddGhostRate, speedScale);
    if (which != "steady") periodicCase(D, oddGhostRate, speedScale);
    return lbmtest::finish();
}
