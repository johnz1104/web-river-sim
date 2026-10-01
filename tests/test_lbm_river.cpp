// End-to-end checks of the homepage river with its default settings.
// Built twice: in double precision (60 s run, conservation and dots) and in single
// precision, the website's storage type (10 simulated minutes, stability).
// With the argument `lake`, the same checks run on the two-bank lake river
// (HomeRiver.hpp) at cell 12 and Re 3000, seen through a MacBook Air-sized view, and
// the page's designed river (left bank only, side inflow, hidden bed) is checked too.

#include "TestUtil.hpp"
#include "HomeBank.hpp"
#include "HomeRiver.hpp"
#include "lbm/River.hpp"
#include "lbm/RiverOptions.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>

using namespace lbm;

namespace {

void conservationAndDots(double seconds) {
    RiverParams params;
    params.seed = 2024;
    RiverSimulation river(homeBankSegments(), kHomeWidth, kHomeHeight, params);
    river.warmUp();
    double inSum = 0.0, outSum = 0.0, maxMach = 0.0, maxDrho = 0.0;
    bool finite = true;
    const double frame = 1.0 / 60.0;
    int samples = 0;
    for (int f = 0; f < static_cast<int>(seconds / frame); ++f) {
        river.advance(frame);
        const RiverDiagnostics d = river.diagnostics();
        finite = finite && d.finite;
        maxMach = std::max(maxMach, d.maxMach);
        maxDrho = std::max(maxDrho, d.maxDensityDeviation);
        if (d.time > 20.0) { inSum += d.inflow; outSum += d.outflow; ++samples; }
    }
    const RiverDiagnostics d = river.diagnostics();
    std::printf("river after %.0f s: inflow %.1f, outflow %.1f (time averages), "
                "max Mach %.3f, max |drho| %.4f, dots %d of %.0f expected\n", d.time,
                inSum / samples, outSum / samples, maxMach, maxDrho, d.tracers, d.expectedTracers);
    lbmtest::check(finite, "river stays finite");
    lbmtest::checkRange(std::abs(outSum - inSum) / inSum, 0.0, 0.005,
                        "time-averaged outflow matches inflow within 0.5%");
    lbmtest::checkRange(inSum / samples / d.targetFlux, 0.995, 1.005,
                        "inflow matches the requested flux within 0.5%");
    lbmtest::checkRange(maxMach, 0.0, 0.2, "peak Mach number below 0.2");
    // The steady part of |drho| is the physical pressure drop into the narrow outlet
    // (Bernoulli: drho ~ 1.5 u_out^2 ~ 1% at this lattice speed).
    lbmtest::checkRange(maxDrho, 0.0, 0.02, "density deviation below 2%");
    lbmtest::checkRange(d.tracers / d.expectedTracers, 0.9, 1.1,
                        "dot count within 10% of density x area");
}

void determinism() {
    auto run = [] {
        RiverParams params;
        params.seed = 77;
        params.warmupTime = 1.0;
        params.rampTime = 1.0;
        params.stirring.enabled = true;
        auto river = std::make_unique<RiverSimulation>(homeBankSegments(), kHomeWidth,
                                                       kHomeHeight, params);
        for (int f = 0; f < 120; ++f) river->advance(1.0 / 60.0);
        return river;
    };
    auto a = run();
    auto b = run();
    bool same = a->tracers().count() == b->tracers().count() && a->tracers().count() > 0;
    for (int i = 0; same && i < a->tracers().count(); ++i)
        same = a->tracers().x(i) == b->tracers().x(i) && a->tracers().y(i) == b->tracers().y(i);
    for (int n = 0; same && n < a->domain().fluidCount(); ++n)
        same = a->solver().ux()[n] == b->solver().ux()[n];
    lbmtest::check(same, "identical seeds give bit-identical flow and dots");
}

void longRun(double seconds) {
    RiverParams params;
    params.seed = 9;
    params.stirring.enabled = true;   // exercise every website switch at once
    RiverSimulation river(homeBankSegments(), kHomeWidth, kHomeHeight, params);
    double maxMach = 0.0;
    bool finite = true;
    for (int f = 0; f < static_cast<int>(seconds * 30); ++f) {
        river.advance(1.0 / 30.0);
        if (f % 30 == 0) {
            const RiverDiagnostics d = river.diagnostics();
            finite = finite && d.finite;
            maxMach = std::max(maxMach, d.maxMach);
            if (!finite) break;
        }
    }
    std::printf("long run: %.0f s simulated, max Mach %.3f\n", river.time(), maxMach);
    lbmtest::check(finite, "single-precision river finite over " +
                               std::to_string(static_cast<int>(seconds)) + " s");
    lbmtest::checkRange(maxMach, 0.0, 0.2, "single-precision peak Mach below 0.2");
}

// ---- the lake river ----

// A MacBook Air window (1440 x 790 CSS px) at the page scale 5/6, in page units.
constexpr double kViewRight = 1440.0 * 6.0 / 5.0, kViewBottom = 790.0 * 6.0 / 5.0;

RiverParams lakeParams(std::uint64_t seed) {
    RiverParams params;
    params.cell = 12.0;
    params.reynolds = 3000.0;   // about today's physical viscosity: the lake carries ~3x the flux
    // The lake's outlet is narrow for its flux, so a lower lattice speed keeps Mach < 0.2.
    params.latticeInletSpeed = 0.010;
    params.seed = seed;
    return params;
}

std::unique_ptr<RiverSimulation> lake(const RiverParams& params) {
    return std::make_unique<RiverSimulation>(homeRiverLeftSegments(), homeRiverRightSegments(),
                                             kRiverWidth, kRiverHeight, params);
}

void lakeConservationAndDots(double seconds) {
    auto river = lake(lakeParams(2024));
    river->setView(kViewRight, kViewBottom);
    river->warmUp();
    double inSum = 0.0, outSum = 0.0, maxMach = 0.0, maxDrho = 0.0;
    bool finite = true;
    const double frame = 1.0 / 60.0;
    int samples = 0;
    for (int f = 0; f < static_cast<int>(seconds / frame); ++f) {
        river->advance(frame);
        const RiverDiagnostics d = river->diagnostics();
        finite = finite && d.finite;
        maxMach = std::max(maxMach, d.maxMach);
        maxDrho = std::max(maxDrho, d.maxDensityDeviation);
        if (d.time > 20.0) { inSum += d.inflow; outSum += d.outflow; ++samples; }
        if (!finite) break;
    }
    const RiverDiagnostics d = river->diagnostics();
    std::printf("lake after %.0f s: inflow %.1f, outflow %.1f (time averages), max Mach %.3f, "
                "max |drho| %.4f, dots in view %d of %.0f expected\n", d.time, inSum / samples,
                outSum / samples, maxMach, maxDrho, d.tracers, d.expectedTracers);
    lbmtest::check(finite, "lake stays finite");
    lbmtest::checkRange(std::abs(outSum - inSum) / inSum, 0.0, 0.005,
                        "lake: time-averaged outflow matches inflow within 0.5%");
    // The diagnostic is a mass flux (rho u); the lake sits ~0.6% above the outlet
    // density (the Bernoulli drop that drives it through the narrow outlet).
    lbmtest::checkRange(inSum / samples / d.targetFlux, 0.99, 1.01,
                        "lake: inflow matches the requested flux within 1%");
    lbmtest::checkRange(maxMach, 0.0, 0.2, "lake: peak Mach number below 0.2");
    // The steady Bernoulli drop into the narrow outlet is ~1.5 u^2 ~ 1.3% at Mach 0.16;
    // eddies passing the outlet add brief peaks near 3% when sampled every frame.
    lbmtest::checkRange(maxDrho, 0.0, 0.04, "lake: density deviation below 4%");
    // Dots enter the view at its edges; the recirculating bay behind the upper bend
    // is rarely entered, so the count settles below density x area.
    lbmtest::checkRange(d.tracers / d.expectedTracers, 0.75, 1.1,
                        "lake: dots in view settle near density x visible area");
}

void lakeDeterminismAndView() {
    auto run = [](bool moveView) {
        RiverParams params = lakeParams(77);
        params.warmupTime = 1.0;
        params.rampTime = 1.0;
        params.stirring.enabled = true;
        auto river = lake(params);
        river->setView(kViewRight, kViewBottom);
        for (int f = 0; f < 120; ++f) {
            // A resize mid-run changes only which dots are kept, never the flow.
            if (moveView && f == 60) river->setView(0.8 * kViewRight, 0.8 * kViewBottom);
            river->advance(1.0 / 60.0);
        }
        return river;
    };
    auto a = run(false);
    auto b = run(false);
    auto c = run(true);
    bool same = a->tracers().count() == b->tracers().count() && a->tracers().count() > 0;
    for (int i = 0; same && i < a->tracers().count(); ++i)
        same = a->tracers().x(i) == b->tracers().x(i) && a->tracers().y(i) == b->tracers().y(i);
    bool sameFlow = true;
    for (int n = 0; n < a->domain().fluidCount(); ++n) {
        same = same && a->solver().ux()[n] == b->solver().ux()[n];
        sameFlow = sameFlow && a->solver().ux()[n] == c->solver().ux()[n] &&
                   a->solver().uy()[n] == c->solver().uy()[n];
    }
    lbmtest::check(same, "lake: identical seeds give bit-identical flow and dots");
    lbmtest::check(sameFlow, "lake: changing the view leaves the flow bit-identical");
}

// A stress test of the side push (PushParams): pulsing, with stirring. The push speeds
// the water up locally, so peak Mach can pass the 0.2 accuracy guard (the tuning panel
// shows it); the check here is stability.
void lakePush(double seconds) {
    RiverParams params = lakeParams(11);
    params.push.strength = 3.0;   // the default placement, at the upper bend's pocket
    params.push.pulse = 0.3;
    params.push.period = 4.0;
    params.stirring.enabled = true;
    auto river = lake(params);
    river->setView(kViewRight, kViewBottom);
    double maxMach = 0.0;
    bool finite = true;
    for (int f = 0; f < static_cast<int>(seconds * 60); ++f) {
        river->advance(1.0 / 60.0);
        const RiverDiagnostics d = river->diagnostics();
        finite = finite && d.finite;
        maxMach = std::max(maxMach, d.maxMach);
        if (!finite) break;
    }
    std::printf("lake with a pulsing push 3 and stirring: %.0f s, max Mach %.3f\n", river->time(),
                maxMach);
    lbmtest::check(finite, "lake with the side push stays finite");
    lbmtest::checkRange(maxMach, 0.0, 0.3, "lake with the side push: peak Mach below 0.3");
}

void lakeLongRun(double seconds) {
    RiverParams params = lakeParams(9);
    params.stirring.enabled = true;
    auto river = lake(params);
    double maxMach = 0.0;
    bool finite = true;
    for (int f = 0; f < static_cast<int>(seconds * 30); ++f) {
        river->advance(1.0 / 30.0);
        if (f % 30 == 0) {
            const RiverDiagnostics d = river->diagnostics();
            finite = finite && d.finite;
            maxMach = std::max(maxMach, d.maxMach);
            if (!finite) break;
        }
    }
    std::printf("lake long run: %.0f s simulated, max Mach %.3f\n", river->time(), maxMach);
    lbmtest::check(finite, "single-precision lake finite over " +
                               std::to_string(static_cast<int>(seconds)) + " s");
    lbmtest::checkRange(maxMach, 0.0, 0.2, "single-precision lake peak Mach below 0.2");
}

// ---- the page's designed river: left bank, side inflow and the hidden bed ----

// The page river: the lake's left bank cut to the 1728 x 1080 page, as config.js does.
std::vector<BankSegment> pageBank() {
    std::vector<BankSegment> bank;
    for (BankSegment s : homeRiverLeftSegments()) {
        if (s.y0 >= 1080.0) break;
        s.y1 = std::min(s.y1, 1080.0);
        bank.push_back(s);
    }
    return bank;
}

// Runs the designed river with a hidden bed; returns the time-averaged upward flux
// (page units^2 / s) across the sill line, right of the sill's tip.
double pageDepthRun(const std::string& extra, double seconds, bool& finite, double& maxMach,
                    double& balance) {
    RiverParams params;
    const std::string error = applyRiverOptions(params,
        "cell=10 re=3000 speed=38 lattice=0.007 side=1 sideTop=940 sideSpeed=60 sideAngle=40 "
        "sideOutlet=0 inflowAmplitude=0.5 inflowFloor=0.6 seed=5 depth=1 " + extra);
    lbmtest::check(error.empty(), "page river options parse: " + error);
    RiverSimulation river(pageBank(), 1728.0, 1080.0, params);
    river.setView(kViewRight, kViewBottom);
    river.warmUp();
    const DepthParams& d = params.depth;
    double flux = 0.0, inSum = 0.0, outSum = 0.0;
    int samples = 0;
    finite = true;
    maxMach = 0.0;
    for (int f = 0; f < static_cast<int>(seconds * 60); ++f) {
        river.advance(1.0 / 60.0);
        const RiverDiagnostics diag = river.diagnostics();
        finite = finite && diag.finite;
        maxMach = std::max(maxMach, diag.maxMach);
        if (!finite) break;
        if (diag.time < 12.0 || f % 6 != 0) continue;
        inSum += diag.inflow;
        outSum += diag.outflow;
        const double step = 10.0;
        for (double x = d.sillX + 50.0; x < 1728.0 + params.side.margin; x += step) {
            double u, v;
            river.velocityAt(x, d.sillY, u, v);
            flux += -v * step;   // page y grows downward
        }
        ++samples;
    }
    balance = samples > 0 ? outSum / inSum : 0.0;
    return samples > 0 ? flux / samples : 0.0;
}

void pageDepth(double seconds) {
    bool finite = false, openFinite = false;
    double maxMach = 0.0, openMach = 0.0, balance = 0.0, openBalance = 0.0;
    const double sill = pageDepthRun("", seconds, finite, maxMach, balance);
    const double open = pageDepthRun("sillWidth=0", seconds, openFinite, openMach, openBalance);
    std::printf("page river with the hidden bed: upward flux across the sill line %.0f "
                "(no sill: %.0f), max Mach %.3f, outflow / inflow %.3f\n",
                sill, open, maxMach, balance);
    lbmtest::check(finite && openFinite, "page river with the hidden bed stays finite");
    lbmtest::checkRange(maxMach, 0.0, 0.2, "page river with the hidden bed: peak Mach below 0.2");
    lbmtest::checkRange(balance, 0.98, 1.02, "page river: outflow matches inflow within 2%");
    // Measured: about 4% at the default sill (100 units wide, depth 0.02).
    lbmtest::check(open > 0.0 && sill < 0.1 * open,
                   "the sill cuts the flow across it to under 10% of the flow without it");
}

}  // namespace

int main(int argc, char** argv) {
    const bool lakeMode = argc > 1 && std::string(argv[1]) == "lake";
#ifdef LBM_SINGLE_PRECISION
    if (lakeMode) lakeLongRun(argc > 2 ? std::atof(argv[2]) : 600.0);
    else longRun(argc > 1 ? std::atof(argv[1]) : 600.0);
#else
    if (lakeMode) {
        lakeConservationAndDots(argc > 2 ? std::atof(argv[2]) : 60.0);
        lakeDeterminismAndView();
        lakePush(60.0);
        pageDepth(30.0);
    } else {
        conservationAndDots(argc > 1 ? std::atof(argv[1]) : 60.0);
        determinism();
    }
#endif
    return lbmtest::finish();
}
