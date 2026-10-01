#include "lbm/Tracers.hpp"

#include <algorithm>
#include <cmath>

namespace lbm {

Tracers::Tracers(const TracerParams& params, std::uint64_t seed, const BankSpline& bank,
                 double rightWall, double inletY)
    : params_(params), rng_(seed, 3), bank_(bank), rightWall_(rightWall), inletY_(inletY) {}

void Tracers::clear() {
    x_.clear(); y_.clear(); age_.clear(); tone_.clear(); tail_.clear();
    ghosts_.clear(); ghostTail_.clear();
}

void Tracers::add(double x, double y) {
    x_.push_back(x);
    y_.push_back(y);
    age_.push_back(0.0);
    tone_.push_back(static_cast<float>(rng_.uniform()));
    // A new dot has been here all along: its tail starts with zero length.
    for (int k = 0; k < params_.tailSamples; ++k) {
        tail_.push_back(static_cast<float>(x));
        tail_.push_back(static_cast<float>(y));
    }
}

void Tracers::emit(const std::vector<InletColumn>& columns, double dt) {
    if (!params_.enabled) return;
    const double density = 1.0 / (params_.spacing * params_.spacing);
    for (const InletColumn& col : columns) {
        if (col.speed <= 0.0) continue;
        // Arrivals are a Poisson process with mean density x flux x dt. A deterministic
        // accumulator would make every column emit in step and draw rows of dots.
        const double mean = density * col.speed * col.width * dt;
        int arrivals = 0;
        double p = std::exp(-mean), cumulative = p;
        const double u = rng_.uniform();
        while (u > cumulative && arrivals < 64) {
            ++arrivals;
            p *= mean / arrivals;
            cumulative += p;
        }
        for (int a = 0; a < arrivals && count() < params_.maxCount; ++a) {
            // Place each arrival within the distance the water moved this frame.
            double x = col.x + (rng_.uniform() - 0.5) * col.width;
            double y = inletY_ - rng_.uniform() * col.speed * dt;
            keepInside(x, y);
            add(x, y);
        }
    }
}

void Tracers::emitAlong(const std::vector<EntrySegment>& segments, double dt) {
    if (!params_.enabled) return;
    const double density = 1.0 / (params_.spacing * params_.spacing);
    for (const EntrySegment& seg : segments) {
        if (seg.speed <= 0.0) continue;
        // Poisson arrivals, as in emit().
        const double mean = density * seg.speed * seg.width * dt;
        int arrivals = 0;
        double p = std::exp(-mean), cumulative = p;
        const double u = rng_.uniform();
        while (u > cumulative && arrivals < 64) {
            ++arrivals;
            p *= mean / arrivals;
            cumulative += p;
        }
        for (int a = 0; a < arrivals && count() < params_.maxCount; ++a) {
            // Place each arrival within the distance the water moved into view this frame.
            const double along = (rng_.uniform() - 0.5) * seg.width;
            const double inward = rng_.uniform() * seg.speed * dt;
            double x = seg.horizontal ? seg.x + along : seg.x - inward;
            double y = seg.horizontal ? seg.y - inward : seg.y + along;
            keepInside(x, y);
            add(x, y);
        }
    }
}

int Tracers::fillDisc(double cx, double cy, double radius,
                      const std::function<bool(double, double)>& inside, int limit) {
    if (!params_.enabled || !(radius > 0.0) || limit == 0) return 0;
    const double density = 1.0 / (params_.spacing * params_.spacing);
    // Candidates uniform over the disc's bounding square, kept with the soft profile.
    const int candidates = static_cast<int>(std::lround(density * 4.0 * radius * radius));
    int added = 0;
    for (int c = 0; c < candidates && count() < params_.maxCount; ++c) {
        const double x = cx + (2.0 * rng_.uniform() - 1.0) * radius;
        const double y = cy + (2.0 * rng_.uniform() - 1.0) * radius;
        const double r = std::hypot(x - cx, y - cy);
        const double keep = std::clamp((radius - r) / (0.3 * radius), 0.0, 1.0);
        if (rng_.uniform() < keep && inside(x, y)) {
            add(x, y);
            if (++added == limit) break;
        }
    }
    return added;
}

std::vector<int> Tracers::countCells(double x0, double y0, double cell, int nx, int ny) const {
    std::vector<int> counts(static_cast<size_t>(nx) * ny, 0);
    for (size_t i = 0; i < x_.size(); ++i) {
        const int cx = static_cast<int>(std::floor((x_[i] - x0) / cell));
        const int cy = static_cast<int>(std::floor((y_[i] - y0) / cell));
        if (cx >= 0 && cx < nx && cy >= 0 && cy < ny) ++counts[static_cast<size_t>(cy) * nx + cx];
    }
    return counts;
}

int Tracers::fillBox(double x0, double y0, double width, double height,
                     const std::function<bool(double, double)>& inside, int count) {
    int added = 0;
    // A bounded number of tries: a box that is mostly land may not take them all.
    for (int t = 0; t < 8 * count && added < count && this->count() < params_.maxCount; ++t) {
        const double x = x0 + rng_.uniform() * width;
        const double y = y0 + rng_.uniform() * height;
        if (inside(x, y)) {
            add(x, y);
            ++added;
        }
    }
    return added;
}

// Reflect a point off a bank's tangent line at height y (the point is (x, y)).
static void reflectOffBank(const BankSpline& bank, double xb, double& x, double& y) {
    const double slope = bank.slope(y);             // tangent (slope, 1)
    const double len2 = slope * slope + 1.0;
    const double vx = x - xb;                       // point relative to (xb, y)
    const double along = vx * slope / len2;         // projection on the tangent
    x = xb + 2.0 * along * slope - vx;
    y = y + 2.0 * along;
}

// Reflect a point that left the water back inside: off the bank's tangent line,
// off the right wall (or the right bank), and off the inlet line.
void Tracers::keepInside(double& x, double& y) const {
    if (y > inletY_) y = 2.0 * inletY_ - y;
    if (hasBox_ && x > boxX0_ && x < boxX1_ && y > boxY0_ && y < boxY1_) {
        // Out through the nearest side of the box.
        const double left = x - boxX0_, right = boxX1_ - x, top = y - boxY0_, bottom = boxY1_ - y;
        const double least = std::min({left, right, top, bottom});
        if (least == left) x = boxX0_ - 1e-6;
        else if (least == right) x = boxX1_ + 1e-6;
        else if (least == top) y = boxY0_ - 1e-6;
        else y = boxY1_ + 1e-6;
    }
    if (right_) {
        const double xr = right_->x(y);
        if (x > xr) {
            reflectOffBank(*right_, xr, x, y);
            const double xr2 = right_->x(y);
            if (x > xr2) x = xr2 - 1e-6;
        }
    } else if (x > rightWall_) {
        x = 2.0 * rightWall_ - x;
    }
    const double xb = bank_.x(y);
    if (x < xb) {
        const double slope = bank_.slope(y);           // tangent (slope, 1)
        const double len2 = slope * slope + 1.0;
        const double vx = x - xb, vy = 0.0;             // point relative to (xb, y)
        const double along = (vx * slope + vy) / len2;  // projection on the tangent
        const double tx = along * slope, ty = along;
        x = xb + 2.0 * tx - vx;
        y = y + 2.0 * ty - vy;
        // A strongly curved bank could leave the mirror image outside; clamp then.
        const double xb2 = bank_.x(y);
        if (x < xb2) x = xb2 + 1e-6;
    }
}

void Tracers::removeAt(size_t i) {
    const size_t last = x_.size() - 1;
    x_[i] = x_[last]; y_[i] = y_[last]; age_[i] = age_[last]; tone_[i] = tone_[last];
    x_.pop_back(); y_.pop_back(); age_.pop_back(); tone_.pop_back();
    const size_t n = 2 * static_cast<size_t>(params_.tailSamples);
    if (n > 0) {
        std::copy(tail_.begin() + last * n, tail_.begin() + (last + 1) * n, tail_.begin() + i * n);
        tail_.resize(last * n);
    }
}

// Every tailInterval seconds (at the first frame after it is due), each dot's position
// goes into the newest tail slot, and the slot remembers when.
void Tracers::recordTails() {
    const int k = params_.tailSamples;
    if (k <= 0 || !(params_.tailInterval > 0.0)) return;
    if (tailTimes_.empty()) tailTimes_.assign(k, 0.0);
    if (tailClock_ < params_.tailInterval) return;
    // A slow frame skips samples rather than repeating one position.
    tailClock_ = std::fmod(tailClock_, params_.tailInterval);
    {
        tailHead_ = (tailHead_ + 1) % k;
        tailTimes_[tailHead_] = clock_;
        for (size_t i = 0; i < x_.size(); ++i) {
            tail_[(i * k + tailHead_) * 2] = static_cast<float>(x_[i]);
            tail_[(i * k + tailHead_) * 2 + 1] = static_cast<float>(y_[i]);
        }
        for (size_t g = 0; g < ghosts_.size(); ++g) {
            ghostTail_[(g * k + tailHead_) * 2] = static_cast<float>(ghosts_[g].x);
            ghostTail_[(g * k + tailHead_) * 2 + 1] = static_cast<float>(ghosts_[g].y);
        }
    }
}

// Dot i has just left, now at (x, y): it carries on as a ghost for one tail length.
void Tracers::addGhost(size_t i, double x, double y) {
    const int k = params_.tailSamples;
    if (k <= 0 || !(params_.tailInterval > 0.0)) return;
    ghosts_.push_back({x, y, age_[i], k * params_.tailInterval, tone_[i]});
    const size_t n = 2 * static_cast<size_t>(k);
    ghostTail_.insert(ghostTail_.end(), tail_.begin() + i * n, tail_.begin() + (i + 1) * n);
}

void Tracers::advanceGhosts(double dt, const VelocityFn& velocity) {
    const size_t n = 2 * static_cast<size_t>(params_.tailSamples);
    for (size_t g = 0; g < ghosts_.size();) {
        Ghost& d = ghosts_[g];
        d.life -= dt;
        if (d.life <= 0.0) {
            const size_t last = ghosts_.size() - 1;
            ghosts_[g] = ghosts_[last];
            ghosts_.pop_back();
            std::copy(ghostTail_.begin() + last * n, ghostTail_.begin() + (last + 1) * n,
                      ghostTail_.begin() + g * n);
            ghostTail_.resize(last * n);
            continue;
        }
        double u0, v0, u1, v1;
        velocity(d.x, d.y, u0, v0);
        velocity(d.x + 0.5 * dt * u0, d.y + 0.5 * dt * v0, u1, v1);
        d.x += dt * u1;
        d.y += dt * v1;
        d.age += dt;
        keepInside(d.x, d.y);
        ++g;
    }
}

void Tracers::advance(double dt, const VelocityFn& velocity) {
    advanceGhosts(dt, velocity);
    const double kick = jitter_ && params_.diffusivity > 0.0 ?
        std::sqrt(2.0 * params_.diffusivity * dt) : 0.0;
    for (size_t i = 0; i < x_.size();) {
        double u0, v0, u1, v1;
        velocity(x_[i], y_[i], u0, v0);
        double xm = x_[i] + 0.5 * dt * u0, ym = y_[i] + 0.5 * dt * v0;
        velocity(xm, ym, u1, v1);
        double x = x_[i] + dt * u1, y = y_[i] + dt * v1;
        if (kick > 0.0) { x += kick * rng_.normal(); y += kick * rng_.normal(); }
        age_[i] += dt;
        if (age_[i] > params_.maxLifetime) { removeAt(i); continue; }
        if (y < 0.0 || (x > rightWall_ && y < rightOutletTop_)) {   // out of the outlet
            addGhost(i, x, y);
            removeAt(i);
            continue;
        }
        keepInside(x, y);
        if (hasView_ && (x > viewRight_ || y > viewBottom_)) {   // out of view
            addGhost(i, x, y);
            removeAt(i);
            continue;
        }
        x_[i] = x;
        y_[i] = y;
        ++i;
    }
    tailClock_ += dt;
    clock_ += dt;
    recordTails();
}

const std::vector<float>& Tracers::tailAges() {
    const int k = params_.tailSamples;
    tailAges_.resize(std::max(k, 0));
    if (tailTimes_.empty() && k > 0) tailTimes_.assign(k, 0.0);
    for (int s = 0; s < k; ++s) {   // newest first
        const int slot = ((tailHead_ - s) % k + k) % k;
        tailAges_[s] = static_cast<float>(clock_ - tailTimes_[slot]);
    }
    return tailAges_;
}

const std::vector<float>& Tracers::packed() {
    const int k = params_.tailSamples;
    const size_t stride = static_cast<size_t>(packedStride());
    packed_.resize(stride * static_cast<size_t>(packedCount()));
    auto write = [&](size_t at, double x, double y, float tone, double age,
                     const std::vector<float>& tails, size_t t) {
        float* out = &packed_[stride * at];
        out[0] = static_cast<float>(x);
        out[1] = static_cast<float>(y);
        out[2] = tone;
        out[3] = static_cast<float>(age);
        for (int s = 0; s < k; ++s) {   // newest first
            const int slot = ((tailHead_ - s) % k + k) % k;
            out[4 + 2 * s] = tails[(t * k + slot) * 2];
            out[5 + 2 * s] = tails[(t * k + slot) * 2 + 1];
        }
    };
    for (size_t i = 0; i < x_.size(); ++i) write(i, x_[i], y_[i], tone_[i], age_[i], tail_, i);
    for (size_t g = 0; g < ghosts_.size(); ++g) {
        const Ghost& d = ghosts_[g];
        write(x_.size() + g, d.x, d.y, d.tone, d.age, ghostTail_, g);
    }
    return packed_;
}

}  // namespace lbm
