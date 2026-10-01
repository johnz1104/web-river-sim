#pragma once

// The visible dots: passive tracers in page reference units. They follow the flow
// (RK2 midpoint), optionally take a small random walk, enter at the inlet in
// proportion to the local inflow, reflect off the banks and leave at the outlet.
// They never interact with each other or change the flow.

#include "lbm/BankSpline.hpp"
#include "lbm/Random.hpp"

#include <cstdint>
#include <functional>
#include <vector>

namespace lbm {

struct TracerParams {
    bool enabled = true;
    double spacing = 22.0;        // mean distance between dots (page units)
    double diffusivity = 1.0;     // random-walk diffusivity (page units^2 / s); 0 = off
    double maxLifetime = 300.0;   // seconds; safety net for dots stuck near a wall
    int maxCount = 6000;
    // Tails: each dot remembers where it was about every tailInterval seconds, for the
    // last tailSamples samples (about 1 s), so the page can draw a line along its recent
    // path (its length grows with the dot's speed). 0 samples = no history.
    int tailSamples = 30;
    double tailInterval = 1.0 / 30.0;
};

// One inlet column: dots enter across [x - width/2, x + width/2] at speed `speed`
// (page units per second, toward decreasing y).
struct InletColumn {
    double x, width, speed;
};

// A piece of a view edge where water flows into view: centred at (x, y), `width` long,
// with inward speed `speed`. Horizontal pieces lie on the bottom edge (water moving up);
// vertical ones on the right edge (water moving left).
struct EntrySegment {
    double x, y, width, speed;
    bool horizontal;
};

class Tracers {
public:
    using VelocityFn = std::function<void(double x, double y, double& u, double& v)>;

    Tracers(const TracerParams& params, std::uint64_t seed, const BankSpline& bank,
            double rightWall, double inletY);

    void emit(const std::vector<InletColumn>& columns, double dt);
    // Place dots at the usual density inside a disc with a soft edge (full density out
    // to 70% of the radius, fading to none at the rim), wherever `inside` holds.
    // Adds at most `limit` dots (all it can when limit < 0); returns how many it added.
    int fillDisc(double cx, double cy, double radius,
                 const std::function<bool(double, double)>& inside, int limit = -1);
    // Dots per cell of an nx x ny grid of square cells starting at (x0, y0), row-major.
    std::vector<int> countCells(double x0, double y0, double cell, int nx, int ny) const;
    // Adds up to `count` dots at random places in the box where `inside` holds.
    int fillBox(double x0, double y0, double width, double height,
                const std::function<bool(double, double)>& inside, int count);
    void emitAlong(const std::vector<EntrySegment>& segments, double dt);
    void advance(double dt, const VelocityFn& velocity);

    // A second bank replaces the straight right wall (dots reflect off its tangent).
    void setRightBank(const BankSpline* right) { right_ = right; }
    // A solid box: a dot that ends up inside is moved out through the nearest side.
    void setObstacle(double x0, double y0, double x1, double y1) {
        hasBox_ = true; boxX0_ = x0; boxY0_ = y0; boxX1_ = x1; boxY1_ = y1;
    }
    // The right wall is open above y = top (an outlet): dots crossing it there leave.
    void setRightOutlet(double top) { rightOutletTop_ = top; }
    // Only the box [.., right] x [.., bottom] is shown: dots leaving it are removed.
    void setView(double right, double bottom) {
        hasView_ = true;
        viewRight_ = right;
        viewBottom_ = bottom;
    }

    // Dots in the river (the ones that have left are not counted).
    int count() const { return static_cast<int>(x_.size()); }
    // Dots to draw: the ones in the river, then the ones that have just left (their
    // tails are still following them out).
    int packedCount() const { return static_cast<int>(x_.size() + ghosts_.size()); }
    // The dots' random walk (diffusivity) on or off, at any time; the page switches it
    // off while it draws tracer lines, so tails follow the water without wiggling.
    void setJitter(bool on) { jitter_ = on; }
    double x(int i) const { return x_[i]; }
    double y(int i) const { return y_[i]; }
    // Per dot, for rendering (packedCount() dots): x, y, tone, age, then its tailSamples
    // past positions (x, y), newest first; packedStride() floats in all.
    const std::vector<float>& packed();
    int packedStride() const { return 4 + 2 * params_.tailSamples; }
    // How old each tail sample is, in seconds, newest first (tailSamples values). Samples
    // are taken at frame times, so they are about, not exactly, tailInterval apart; the
    // page uses these ages to cut every tail at exactly the length it wants.
    const std::vector<float>& tailAges();

    void clear();
    // Test hook: place a dot directly.
    void add(double x, double y);
    const TracerParams& params() const { return params_; }

private:
    void keepInside(double& x, double& y) const;
    void removeAt(size_t i);

    TracerParams params_;
    RandomStream rng_;
    const BankSpline& bank_;
    const BankSpline* right_ = nullptr;
    double rightWall_, inletY_;
    double rightOutletTop_ = -1.0;   // < 0: the right wall is closed
    bool hasView_ = false;
    double viewRight_ = 0.0, viewBottom_ = 0.0;
    std::vector<double> x_, y_, age_;
    std::vector<float> tone_;
    std::vector<float> packed_;
    // Tail samples: tailSamples (x, y) pairs per dot in a ring shared by all dots
    // (they are all sampled at the same times); tailHead_ is the newest slot.
    std::vector<float> tail_;
    int tailHead_ = 0;
    double tailClock_ = 0.0;   // seconds since the last sample was due
    double clock_ = 0.0;       // seconds advanced in all
    std::vector<double> tailTimes_;   // clock_ when each ring slot was sampled
    std::vector<float> tailAges_;
    void recordTails();
    bool jitter_ = true;
    bool hasBox_ = false;
    double boxX0_ = 0.0, boxY0_ = 0.0, boxX1_ = 0.0, boxY1_ = 0.0;
    // Dots that have left (through the outlet or out of view) keep following the water
    // for one tail length, so their tails flow out instead of vanishing. They are kept
    // apart from the dots in the river: never counted, jittered or re-entered, and they
    // use no random numbers, so the river's dots are exactly as without them.
    struct Ghost { double x, y, age, life; float tone; };
    std::vector<Ghost> ghosts_;
    std::vector<float> ghostTail_;   // tailSamples (x, y) per ghost, in the same ring
    void addGhost(size_t i, double x, double y);
    void advanceGhosts(double dt, const VelocityFn& velocity);
};

}  // namespace lbm
