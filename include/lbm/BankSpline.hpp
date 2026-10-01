#pragma once

// The page's riverbank x = f(y), rebuilt from the cubic segments that
// src/river/banks.js computes (createCubicBank). Coefficients arrive as data, so
// the drawn line and the simulated wall are the same curve.

#include <vector>

namespace lbm {

// x = a + b u + c u^2 + d u^3, u = y - y0, for y0 <= y <= y1 (page reference units).
struct BankSegment {
    double y0, y1, a, b, c, d;
};

class BankSpline {
public:
    explicit BankSpline(std::vector<BankSegment> segments);

    // Position and slope. Beyond the ends the bank continues along its end tangent,
    // so geometry queries slightly outside the page stay well defined.
    double x(double y) const;
    double slope(double y) const;

    double yMin() const { return segments_.front().y0; }
    double yMax() const { return segments_.back().y1; }
    double minX() const { return minX_; }
    double maxX() const { return maxX_; }

    // Euclidean distance from (px, py) to the curve (px to the right of the bank).
    double distance(double px, double py) const;

    const std::vector<BankSegment>& segments() const { return segments_; }

private:
    const BankSegment& segmentFor(double y) const;
    std::vector<BankSegment> segments_;
    double minX_, maxX_;
};

}  // namespace lbm
