#pragma once

// Minimal test helpers: each check prints its result; main returns the failure count.

#include <cmath>
#include <cstdio>
#include <string>

namespace lbmtest {

inline int& failures() {
    static int count = 0;
    return count;
}

inline void check(bool ok, const std::string& what) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
    if (!ok) ++failures();
}

inline void checkNear(double value, double expected, double tolerance, const std::string& what) {
    const bool ok = std::isfinite(value) && std::abs(value - expected) <= tolerance;
    char buffer[256];
    std::snprintf(buffer, sizeof buffer, "%s: %.6g (expected %.6g +/- %.3g)", what.c_str(),
                  value, expected, tolerance);
    check(ok, buffer);
}

inline void checkRange(double value, double lo, double hi, const std::string& what) {
    const bool ok = std::isfinite(value) && value >= lo && value <= hi;
    char buffer[256];
    std::snprintf(buffer, sizeof buffer, "%s: %.6g (expected in [%.6g, %.6g])", what.c_str(),
                  value, lo, hi);
    check(ok, buffer);
}

inline int finish() {
    std::printf("%d failure(s)\n", failures());
    return failures() == 0 ? 0 : 1;
}

}  // namespace lbmtest
