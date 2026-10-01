#include "lbm/RiverOptions.hpp"

#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>

namespace lbm {

namespace {

// The whole string must be the number: "", "12abc" and non-finite values are errors.
bool parseNumber(const std::string& s, double& out) {
    if (s.empty() || std::isspace(static_cast<unsigned char>(s[0]))) return false;
    char* end = nullptr;
    errno = 0;
    out = std::strtod(s.c_str(), &end);
    return errno == 0 && end == s.c_str() + s.size() && std::isfinite(out);
}

bool parseUnsigned(const std::string& s, std::uint64_t& out) {
    if (s.empty() || !std::isdigit(static_cast<unsigned char>(s[0]))) return false;
    char* end = nullptr;
    errno = 0;
    const unsigned long long v = std::strtoull(s.c_str(), &end, 10);
    if (errno != 0 || end != s.c_str() + s.size()) return false;
    out = static_cast<std::uint64_t>(v);
    return true;
}

std::string bad(const std::string& key, const std::string& value, const char* expected) {
    return "bad value for " + key + ": \"" + value + "\" (expected " + expected + ")";
}

template <class Check>
std::string setNumber(const std::string& key, const std::string& value, double& target,
                      Check ok, const char* expected) {
    double v = 0.0;
    if (!parseNumber(value, v) || !ok(v)) return bad(key, value, expected);
    target = v;
    return "";
}

std::string setSwitch(const std::string& key, const std::string& value, bool& target) {
    if (value != "0" && value != "1") return bad(key, value, "0 or 1");
    target = value == "1";
    return "";
}

bool positive(double v) { return v > 0.0; }
bool nonNegative(double v) { return v >= 0.0; }
bool rate(double v) { return v > 0.0 && v < 2.0; }
bool finite(double) { return true; }   // parseNumber already rejects non-finite values

}  // namespace

std::string applyRiverOption(RiverParams& p, const std::string& key, const std::string& value) {
    if (key == "rows") {
        std::uint64_t rows = 0;
        if (!parseUnsigned(value, rows) || rows < 8 || rows > 1024)
            return bad(key, value, "an integer from 8 to 1024");
        p.rows = static_cast<int>(rows);
        return "";
    }
    if (key == "outletCells") {
        std::uint64_t cells = 0;
        if (!parseUnsigned(value, cells) || cells > 400)
            return bad(key, value, "an integer from 0 to 400");
        p.outletCells = static_cast<int>(cells);
        return "";
    }
    if (key == "seed") {
        std::uint64_t seed = 0;
        if (!parseUnsigned(value, seed)) return bad(key, value, "an unsigned integer");
        p.seed = seed;
        return "";
    }
    if (key == "model") {
        if (value == "bgk") p.collision.model = CollisionModel::BGK;
        else if (value == "trt") p.collision.model = CollisionModel::TRT;
        else if (value == "mrt") p.collision.model = CollisionModel::MRT;
        else return bad(key, value, "bgk, trt or mrt");
        return "";
    }
    if (key == "smag") {
        const std::string error =
            setNumber(key, value, p.collision.smagorinsky, nonNegative, "a number >= 0");
        if (error.empty())
            p.collision.closure =
                p.collision.smagorinsky > 0 ? ClosureKind::Smagorinsky : ClosureKind::Laminar;
        return error;
    }
    if (key == "cell")
        return setNumber(key, value, p.cell, [](double v) { return v >= 2.0 && v <= 64.0; },
                         "page units in [2, 64]");
    if (key == "re") return setNumber(key, value, p.reynolds, positive, "a number > 0");
    if (key == "speed") return setNumber(key, value, p.inletSpeed, positive, "a number > 0");
    if (key == "lattice")
        return setNumber(key, value, p.latticeInletSpeed,
                         [](double v) { return v > 0.0 && v <= 0.3; }, "a number in (0, 0.3]");
    if (key == "sb") return setNumber(key, value, p.collision.bulkRate, rate, "a rate in (0, 2)");
    if (key == "se")
        return setNumber(key, value, p.collision.evenGhostRate, rate, "a rate in (0, 2)");
    if (key == "sq")
        return setNumber(key, value, p.collision.oddGhostRate,
                         [](double v) { return v >= 0.0 && v < 2.0; }, "a rate in [0, 2)");
    if (key == "magic") return setNumber(key, value, p.collision.magic, positive, "a number > 0");
    if (key == "ramp") return setNumber(key, value, p.rampTime, nonNegative, "seconds >= 0");
    if (key == "warmup") return setNumber(key, value, p.warmupTime, nonNegative, "seconds >= 0");
    if (key == "inflowNoise") return setSwitch(key, value, p.inflow.enabled);
    if (key == "inflowAmplitude")
        return setNumber(key, value, p.inflow.amplitude, [](double v) { return v >= 0.0 && v <= 1.0; },
                         "a number in [0, 1]");
    if (key == "inflowFloor")
        return setNumber(key, value, p.inflow.floor, [](double v) { return v >= 0.0 && v <= 1.0; },
                         "a number in [0, 1]");
    if (key == "inflowModes") {
        std::uint64_t modes = 0;
        if (!parseUnsigned(value, modes) || modes < 1 || modes > 8)
            return bad(key, value, "an integer from 1 to 8");
        p.inflow.modes = static_cast<int>(modes);
        return "";
    }
    if (key == "meander") return setSwitch(key, value, p.inflow.meander);
    if (key == "stir") return setSwitch(key, value, p.stirring.enabled);
    if (key == "stirIntensity")
        return setNumber(key, value, p.stirring.intensity, nonNegative, "a number >= 0");
    if (key == "kT") return setNumber(key, value, p.thermalKT, nonNegative, "a number >= 0");
    if (key == "push")
        return setNumber(key, value, p.push.strength, nonNegative, "a number >= 0 (0 = off)");
    if (key == "pushX") return setNumber(key, value, p.push.x, finite, "page units");
    if (key == "pushY") return setNumber(key, value, p.push.y, finite, "page units");
    if (key == "pushRadius")
        return setNumber(key, value, p.push.radius, positive, "page units > 0");
    if (key == "pushAimX") return setNumber(key, value, p.push.aimX, finite, "page units");
    if (key == "pushAimY") return setNumber(key, value, p.push.aimY, finite, "page units");
    if (key == "pushPulse")
        return setNumber(key, value, p.push.pulse, [](double v) { return v >= 0.0 && v <= 1.0; },
                         "a number in [0, 1]");
    if (key == "pushPeriod")
        return setNumber(key, value, p.push.period, positive, "seconds > 0");
    if (key == "fill")
        return setNumber(key, value, p.fill.radius, nonNegative, "page units >= 0 (0 = off)");
    if (key == "fillX") return setNumber(key, value, p.fill.x, finite, "page units");
    if (key == "fillY") return setNumber(key, value, p.fill.y, finite, "page units");
    if (key == "fillHold")
        return setNumber(key, value, p.fill.hold, nonNegative, "seconds >= 0");
    if (key == "side") return setSwitch(key, value, p.side.enabled);
    if (key == "sideTop") return setNumber(key, value, p.side.top, nonNegative, "page y >= 0");
    if (key == "sideSpeed")
        return setNumber(key, value, p.side.speed, nonNegative, "page units/s >= 0");
    if (key == "sideAngle")
        return setNumber(key, value, p.side.angle, [](double v) { return v > -80.0 && v < 80.0; },
                         "degrees in (-80, 80)");
    if (key == "sideBottom") return setNumber(key, value, p.side.bottom, positive, "page y > 0");
    if (key == "outletRight")
        return setNumber(key, value, p.side.outletRight, positive, "page x > 0");
    if (key == "sideOutlet") return setSwitch(key, value, p.side.outlet);
    if (key == "topOutlet") return setSwitch(key, value, p.side.topOutlet);
    if (key == "sideOutletBottom")
        return setNumber(key, value, p.side.outletBottom, nonNegative, "page y >= 0");
    if (key == "sideMargin")
        return setNumber(key, value, p.side.margin, nonNegative, "page units >= 0");
    if (key == "boxX") return setNumber(key, value, p.box.x, finite, "page units");
    if (key == "boxY") return setNumber(key, value, p.box.y, finite, "page units");
    if (key == "boxW") return setNumber(key, value, p.box.width, nonNegative, "page units >= 0 (0 = none)");
    if (key == "boxH") return setNumber(key, value, p.box.height, nonNegative, "page units >= 0 (0 = none)");
    if (key == "cornerJet")
        return setNumber(key, value, p.cornerJet.width, nonNegative, "page units >= 0 (0 = off)");
    if (key == "cornerSpeed")
        return setNumber(key, value, p.cornerJet.speed, positive, "page units/s > 0");
    if (key == "cornerAngle")
        return setNumber(key, value, p.cornerJet.angle, [](double v) { return v > -80.0 && v < 80.0; },
                         "degrees in (-80, 80)");
    if (key == "cornerEdge")
        return setNumber(key, value, p.cornerJet.edge, nonNegative, "page units >= 0");
    if (key == "converge") return setSwitch(key, value, p.converge.enabled);
    if (key == "convergeX") return setNumber(key, value, p.converge.x, finite, "page units");
    if (key == "convergeY") return setNumber(key, value, p.converge.y, finite, "page units");
    if (key == "convergeMax")
        return setNumber(key, value, p.converge.maxAngle,
                         [](double v) { return v > 0.0 && v < 89.0; }, "degrees in (0, 89)");
    if (key == "sideConverge") return setSwitch(key, value, p.converge.side);
    if (key == "depth") return setSwitch(key, value, p.depth.enabled);
    if (key == "bedDrag")
        return setNumber(key, value, p.depth.bedDrag, nonNegative, "a rate per second >= 0");
    if (key == "depth1") return setNumber(key, value, p.depth.depth1, positive, "a depth > 0");
    if (key == "depth2") return setNumber(key, value, p.depth.depth2, positive, "a depth > 0");
    if (key == "depth3") return setNumber(key, value, p.depth.depth3, positive, "a depth > 0");
    if (key == "sillDepth")
        return setNumber(key, value, p.depth.sillDepth, positive, "a depth > 0");
    if (key == "sillX") return setNumber(key, value, p.depth.sillX, finite, "page units");
    if (key == "sillY") return setNumber(key, value, p.depth.sillY, finite, "page units");
    if (key == "sillWidth")
        return setNumber(key, value, p.depth.sillWidth, nonNegative, "page units >= 0 (0 = no sill)");
    if (key == "curveBottom")
        return setNumber(key, value, p.depth.curveBottom, finite, "page units");
    if (key == "stepWidth")
        return setNumber(key, value, p.depth.stepWidth, nonNegative, "page units >= 0");
    if (key == "spacing")
        return setNumber(key, value, p.tracers.spacing, positive, "page units > 0");
    if (key == "diffusivity")
        return setNumber(key, value, p.tracers.diffusivity, nonNegative, "a number >= 0");
    return "unknown option: " + key;
}

std::string applyRiverOptions(RiverParams& params, const std::string& list) {
    RiverParams p = params;
    size_t i = 0;
    while (i < list.size()) {
        auto separator = [](char c) { return c == ',' || std::isspace(static_cast<unsigned char>(c)); };
        while (i < list.size() && separator(list[i])) ++i;
        size_t end = i;
        while (end < list.size() && !separator(list[end])) ++end;
        if (end == i) break;
        const std::string item = list.substr(i, end - i);
        i = end;
        const size_t eq = item.find('=');
        if (eq == std::string::npos || eq == 0) return "expected key=value, got \"" + item + "\"";
        const std::string error = applyRiverOption(p, item.substr(0, eq), item.substr(eq + 1));
        if (!error.empty()) return error;
    }
    params = p;
    return "";
}

}  // namespace lbm
