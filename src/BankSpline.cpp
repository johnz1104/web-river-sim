#include "lbm/BankSpline.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lbm {

namespace {

double evaluate(const BankSegment& s, double y) {
    const double u = y - s.y0;
    return ((s.d * u + s.c) * u + s.b) * u + s.a;
}

double evaluateSlope(const BankSegment& s, double y) {
    const double u = y - s.y0;
    return (3.0 * s.d * u + 2.0 * s.c) * u + s.b;
}

}  // namespace

BankSpline::BankSpline(std::vector<BankSegment> segments) : segments_(std::move(segments)) {
    if (segments_.empty()) throw std::invalid_argument("BankSpline: no segments");
    for (size_t i = 0; i < segments_.size(); ++i) {
        if (!(segments_[i].y1 > segments_[i].y0))
            throw std::invalid_argument("BankSpline: segments must have increasing y");
        if (i > 0 && std::abs(segments_[i].y0 - segments_[i - 1].y1) > 1e-9)
            throw std::invalid_argument("BankSpline: segments must be contiguous");
    }
    minX_ = maxX_ = x(yMin());
    const int samples = 4000;
    for (int s = 0; s <= samples; ++s) {
        const double xs = x(yMin() + (yMax() - yMin()) * s / samples);
        minX_ = std::min(minX_, xs);
        maxX_ = std::max(maxX_, xs);
    }
}

const BankSegment& BankSpline::segmentFor(double y) const {
    for (const BankSegment& s : segments_)
        if (y <= s.y1) return s;
    return segments_.back();
}

double BankSpline::x(double y) const {
    if (y < yMin()) return evaluate(segments_.front(), yMin()) + slope(yMin()) * (y - yMin());
    if (y > yMax()) return evaluate(segments_.back(), yMax()) + slope(yMax()) * (y - yMax());
    return evaluate(segmentFor(y), y);
}

double BankSpline::slope(double y) const {
    const double yc = std::min(std::max(y, yMin()), yMax());
    return evaluateSlope(segmentFor(yc), yc);
}

double BankSpline::distance(double px, double py) const {
    // The horizontal gap bounds the distance, so only |y' - py| <= gap can be closer.
    const double gap = std::abs(px - x(py));
    auto d2 = [&](double y) { const double dx = px - x(y), dy = py - y; return dx * dx + dy * dy; };
    double best = py, bestD2 = gap * gap;
    const int samples = 64;
    for (int s = 0; s <= samples; ++s) {
        const double y = py - gap + 2.0 * gap * s / samples;
        const double v = d2(y);
        if (v < bestD2) { bestD2 = v; best = y; }
    }
    // Golden-section refinement around the best sample.
    double lo = best - 2.0 * gap / samples, hi = best + 2.0 * gap / samples;
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    for (int it = 0; it < 40; ++it) {
        const double a = hi - g * (hi - lo), b = lo + g * (hi - lo);
        if (d2(a) < d2(b)) hi = b; else lo = a;
    }
    return std::sqrt(std::min(bestD2, d2(0.5 * (lo + hi))));
}

}  // namespace lbm
