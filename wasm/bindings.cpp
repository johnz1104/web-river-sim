// C API of the WebAssembly build. One river per module instance: the website runs
// one module per Web Worker (src/river/sim/worker.js). Plain extern "C" functions
// (no Embind), so JavaScript passes only numbers and pointers into module memory.
//
// river_create(left, nLeft, right, nRight, width, height, seedLo, seedHi, options)
//   left     nLeft bank segments, 6 doubles each: y0, y1, a, b, c, d (banks.js order)
//   right    nRight segments of the right bank; nRight = 0 keeps the straight right
//            wall at x = width (the one-bank river)
//   seed     64-bit seed as two 32-bit halves; it overrides any seed= option
//   options  "key=value" list, the keys of lbm/RiverOptions.hpp ("" for defaults)
//   Returns 0, or 1 with the message in river_last_error(). A failed create leaves
//   no river; nothing aborts.
//
// river_export_flow(out) writes the flow, 3 * nx * ny floats (nx, ny: diagnostics 12
// and 13): density, then x and y velocity, in lattice units, row-major; NaN density
// marks a cell that is not water. river_start_from_flow(in) starts the river from such
// a flow instead of the hidden warm-up (RiverSimulation::startFromFlow); call it
// before the first advance. It returns 0, or 1 with the message in river_last_error().
//
// river_grid(out) writes 4 doubles that place the flow's cells on the page: x0 and y0
// (page x and y of the lattice's top-left corner; cell (k, j) is centred at
// x0 + (k + 0.5) h, y0 + (j + 0.5) h), h (page units per cell) and the page units per
// second of one lattice velocity unit.
//
// river_set_view(right, bottom) shows only x <= right, y <= bottom (page units): dots
// enter where water flows in across those edges. It never changes the flow.
//
// river_dots() points at river_dot_count() dots, river_dot_stride() floats each: x, y,
// tone, age, then the dot's tail: its past positions (x, y), newest first, taken about
// every 1/30 s. river_tail_ages() points at their ages in seconds (the same for every
// dot). Page reference units (x right, y down). Both stay valid until the next call
// into the module.
//
// river_diagnostics(out) writes RIVER_DIAGNOSTIC_COUNT doubles, in this order:
//    0 time             simulated seconds, warm-up included
//    1 steps            lattice steps so far
//    2 stepsPerSecond   lattice steps per simulated second
//    3 inflow           page units^2/s through the row next to the inlet
//    4 outflow          the same through the row next to the outlet
//    5 targetFlux       requested flux: inlet speed x inlet width, once ramped
//    6 maxMach          max |u| / c_s
//    7 maxDensityDev    max |rho - rho0| / rho0
//    8 maxTauEffective
//    9 finite           1 or 0
//   10 dots             current dot count
//   11 expectedDots     dot density x river area
//   12 nx               lattice columns
//   13 ny               lattice rows
//   14 fluidNodes
//   15 tau              molecular relaxation time

#include "lbm/River.hpp"
#include "lbm/RiverOptions.hpp"

#include <cmath>
#include <cstdint>
#include <exception>
#include <memory>
#include <string>
#include <vector>

namespace {

constexpr int RIVER_DIAGNOSTIC_COUNT = 16;

std::unique_ptr<lbm::RiverSimulation> river;
std::string lastError;

int fail(const std::string& message) {
    lastError = message;
    return 1;
}

std::vector<lbm::BankSegment> segments(const double* segs, int n) {
    std::vector<lbm::BankSegment> bank(n);
    for (int i = 0; i < n; ++i) {
        const double* s = segs + 6 * i;
        bank[i] = {s[0], s[1], s[2], s[3], s[4], s[5]};
    }
    return bank;
}

}  // namespace

extern "C" {

int river_create(const double* left, int nLeft, const double* right, int nRight, double width,
                 double height, std::uint32_t seedLo, std::uint32_t seedHi,
                 const char* options) {
    river.reset();
    lastError.clear();
    if (!left || nLeft < 1) return fail("the bank needs at least one segment");
    if (nRight < 0 || (nRight > 0 && !right)) return fail("bad right bank");
    if (!(std::isfinite(width) && width > 0 && std::isfinite(height) && height > 0))
        return fail("width and height must be positive");
    try {
        lbm::RiverParams params;
        const std::string error = lbm::applyRiverOptions(params, options ? options : "");
        if (!error.empty()) return fail(error);
        params.seed = (static_cast<std::uint64_t>(seedHi) << 32) | seedLo;
        if (nRight > 0)
            river = std::make_unique<lbm::RiverSimulation>(
                segments(left, nLeft), segments(right, nRight), width, height, params);
        else
            river = std::make_unique<lbm::RiverSimulation>(segments(left, nLeft), width, height,
                                                           params);
    } catch (const std::exception& e) {
        river.reset();
        return fail(e.what());
    }
    return 0;
}

const char* river_last_error() { return lastError.c_str(); }

void river_set_view(double right, double bottom) {
    if (river && right > 0 && bottom > 0) river->setView(right, bottom);
}

void river_warm_up() {
    if (river) river->warmUp();
}

void river_export_flow(float* out) {
    if (river && out) river->exportFlow(out);
}

void river_grid(double* out) {
    if (!river || !out) return;
    const lbm::RiverUnits& u = river->units();
    out[0] = u.x0;
    out[1] = u.y0;
    out[2] = u.h;
    out[3] = u.velocityScale;
}

int river_start_from_flow(const float* in) {
    if (!river || !in) return fail("no river or no flow");
    try {
        river->startFromFlow(in);
    } catch (const std::exception& e) {
        return fail(e.what());
    }
    return 0;
}

void river_advance(double seconds) {
    if (river && seconds > 0) river->advance(seconds);
}

// Dots to draw: those in the river, then those that have just left (drawn until their
// tails have followed them out). The diagnostics count only the first.
int river_dot_count() { return river ? river->tracers().packedCount() : 0; }

// The dots' random walk on (1) or off (0); the flow is never touched.
void river_set_dot_jitter(int on) {
    if (river) river->tracers().setJitter(on != 0);
}

const float* river_dots() { return river ? river->tracers().packed().data() : nullptr; }

int river_dot_stride() { return river ? river->tracers().packedStride() : 4; }

const float* river_tail_ages() { return river ? river->tracers().tailAges().data() : nullptr; }

int river_diagnostic_count() { return RIVER_DIAGNOSTIC_COUNT; }

int river_diagnostics(double* out) {
    if (!river || !out) return 0;
    const lbm::RiverDiagnostics d = river->diagnostics();
    const lbm::RiverUnits& u = river->units();
    const double values[RIVER_DIAGNOSTIC_COUNT] = {
        d.time, static_cast<double>(d.steps), u.stepsPerSecond, d.inflow, d.outflow,
        d.targetFlux, d.maxMach, d.maxDensityDeviation, d.maxTauEffective,
        d.finite ? 1.0 : 0.0, static_cast<double>(d.tracers), d.expectedTracers,
        static_cast<double>(u.nx), static_cast<double>(u.ny),
        static_cast<double>(river->domain().fluidCount()), u.tau,
    };
    for (int i = 0; i < RIVER_DIAGNOSTIC_COUNT; ++i) out[i] = values[i];
    return RIVER_DIAGNOSTIC_COUNT;
}

void river_destroy() { river.reset(); }

}  // extern "C"
