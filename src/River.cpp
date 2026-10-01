#include "lbm/River.hpp"

#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace lbm {

RiverSimulation::RiverSimulation(std::vector<BankSegment> bank, double width, double height,
                                 const RiverParams& params)
    : params_(params), bank_(std::move(bank)) {
    build(width, height);
}

RiverSimulation::RiverSimulation(std::vector<BankSegment> left, std::vector<BankSegment> right,
                                 double width, double height, const RiverParams& params)
    : params_(params), bank_(std::move(left)),
      right_(std::make_unique<BankSpline>(std::move(right))) {
    build(width, height);
}

void RiverSimulation::build(double width, double height) {
    const bool side = params_.side.enabled;
    {
        const ObstacleParams& box = params_.box;
        if (!(box.width >= 0.0) || !(box.height >= 0.0) ||
            (box.width > 0.0 && box.height > 0.0 &&
             !(box.x >= 0.0 && box.y >= 0.0 && box.x + box.width <= width &&
               box.y + box.height <= height)))
            throw std::invalid_argument("RiverSimulation: the box must lie inside the page");
    }
    {
        const CornerJetParams& jet = params_.cornerJet;
        if (!(jet.width >= 0.0) || !(jet.speed > 0.0) || !(jet.angle > -80.0 && jet.angle < 80.0) ||
            !(jet.edge >= 0.0))
            throw std::invalid_argument(
                "RiverSimulation: corner jet width, edge >= 0, speed > 0, angle in (-80, 80)");
    }
    if (params_.converge.enabled &&
        (!(params_.converge.y < height) ||
         !(params_.converge.maxAngle > 0.0 && params_.converge.maxAngle < 89.0)))
        throw std::invalid_argument(
            "RiverSimulation: the inflow focus must be above the inlet, max angle in (0, 89)");
    if (side) {
        if (right_) throw std::invalid_argument("RiverSimulation: the side inflow needs one bank");
        if (!(params_.side.top >= 0.0 && params_.side.top < std::min(params_.side.bottom, height)) ||
            !(params_.side.margin >= 0.0))
            throw std::invalid_argument("RiverSimulation: the side inflow band must be inside the page");
        // The domain's right edge sits `margin` past the page's, out of sight.
        width += params_.side.margin;
    }
    if (params_.cell > 0.0)
        params_.rows = std::max(8, static_cast<int>(std::lround(height / params_.cell)));
    if (params_.rows < 8) throw std::invalid_argument("RiverSimulation: too few rows");

    // ---- units ----
    RiverUnits& u = units_;
    u.width = width;
    u.height = height;
    u.h = height / params_.rows;
    u.ny = params_.rows;
    const int baseRows = params_.nestedBaseRows ? params_.nestedBaseRows : params_.rows;
    if (baseRows < 8 || params_.rows % baseRows != 0)
        throw std::invalid_argument("RiverSimulation: rows must be a multiple of nestedBaseRows");
    const double baseH = height / baseRows;
    if (right_) {
        // Two banks: the lattice spans both with a margin; no fluid reaches its sides.
        const int baseNx = static_cast<int>(
            std::ceil((right_->maxX() - bank_.minX() + 4.0 * baseH) / baseH));
        u.nx = baseNx * (params_.rows / baseRows);
        u.x0 = bank_.minX() - 2.0 * baseH;
        u.inletWidth = right_->x(height) - bank_.x(height);
        if (!(u.inletWidth > 0.0)) throw std::invalid_argument("RiverSimulation: banks cross at the inlet");
    } else {
        // The lattice's right side is the straight wall at x = width.
        const int baseNx = static_cast<int>(std::ceil((width - bank_.minX() + 2.0 * baseH) / baseH));
        u.nx = baseNx * (params_.rows / baseRows);
        u.x0 = width - baseNx * baseH;
        u.inletWidth = width - bank_.x(height);
    }
    u.y0 = 0.0;
    const bool topOutlet = right_ || (side && params_.side.topOutlet);
    if (topOutlet) {
        u.ny += params_.outletCells;
        u.y0 = -params_.outletCells * u.h;
    }
    u.stepsPerSecond = params_.inletSpeed / (params_.latticeInletSpeed * u.h);
    u.dt = 1.0 / u.stepsPerSecond;
    u.nu = params_.latticeInletSpeed * (u.inletWidth / u.h) / params_.reynolds;
    u.tau = tauFromViscosity(u.nu);
    u.velocityScale = params_.inletSpeed / params_.latticeInletSpeed;

    // ---- geometry: fluid right of the bank; outlet on top (YMin), inlet at the bottom ----
    DomainSpec spec;
    spec.nx = u.nx;
    spec.ny = u.ny;
    const double x0 = u.x0, y0 = u.y0, h = u.h;
    const BankSpline& b = bank_;
    const BankSpline* r = right_.get();
    if (r)
        spec.phi = [x0, y0, h, &b, r](double X, double Y) {
            // Above the page (y < 0) both banks continue along their end tangents.
            const double x = x0 + X * h, y = y0 + Y * h;
            return std::min(x - b.x(y), r->x(y) - x);
        };
    else if (side) {
        const double lid = params_.side.outletRight;
        spec.phi = [x0, y0, h, &b, lid](double X, double Y) {
            // Above the page (y < 0) the bank continues along its end tangent, and the
            // lid (if any) runs parallel to it, so the outlet channel keeps its width.
            // (A vertical lid met the leaning bank and sealed the channel.)
            const double x = x0 + X * h, y = y0 + Y * h;
            const double fromBank = x - b.x(y);
            if (y >= 0.0) return fromBank;
            return std::min(fromBank, lid + (b.x(y) - b.x(0.0)) - x);
        };
    }
    else
        spec.phi = [x0, h, &b](double X, double Y) { return (x0 + X * h) - b.x(Y * h); };
    const ObstacleParams box = params_.box;
    if (box.width > 0.0 && box.height > 0.0) {
        // The solid box: the signed distance to it (negative inside), met with the rest.
        auto outside = spec.phi;
        spec.phi = [outside, box, x0, y0, h](double X, double Y) {
            const double x = x0 + X * h, y = y0 + Y * h;
            const double dx = std::max({box.x - x, 0.0, x - (box.x + box.width)});
            const double dy = std::max({box.y - y, 0.0, y - (box.y + box.height)});
            const double d = dx > 0.0 || dy > 0.0 ? std::hypot(dx, dy) :
                -std::min({x - box.x, box.x + box.width - x, y - box.y, box.y + box.height - y});
            return std::min(outside(X, Y), d);
        };
    }
    spec.xMin = SideKind::Wall;
    spec.xMax = SideKind::Wall;
    spec.yMin = SideKind::PressureOutlet;
    spec.yMax = SideKind::MovingWall;
    if (side) {
        // The right edge: inflow below `top`, outflow above.
        const double top = params_.side.top, bottom = params_.side.bottom;
        const double outletEnd = params_.side.outletBottom;
        const SideKind above = params_.side.outlet ? SideKind::PressureOutlet : SideKind::Wall;
        spec.xMax = SideKind::MovingWall;
        spec.sideKindAt = [y0, h, top, bottom, outletEnd, above, topOutlet](Side s, double,
                                                                            double Y) {
            const double y = y0 + Y * h;
            switch (s) {
                case XMax:
                    if (y < top) return y < outletEnd ? above : SideKind::Wall;
                    return y > bottom ? SideKind::Wall : SideKind::MovingWall;
                case YMin: return topOutlet ? SideKind::PressureOutlet : SideKind::Wall;
                case YMax: return SideKind::MovingWall;
                default: return SideKind::Wall;
            }
        };
    }
    domain_ = std::make_unique<Domain>(spec);

    CollisionConfig cfg = params_.collision;
    cfg.tau = u.tau;
    cfg.kT = params_.thermalKT;
    solver_ = std::make_unique<Solver>(*domain_, cfg, 1.0, hashKey(params_.seed, 4, 0, 0));
    solver_->setThreadCount(params_.threads);
    if (right_ || side) solver_->setOutletBackflowClamp(true);

    // Exact wall distance (lattice units) for closures that want it.
    std::vector<double> wallDistance(domain_->fluidCount());
    for (int n = 0; n < domain_->fluidCount(); ++n) {
        const double x = x0 + (domain_->k(n) + 0.5) * h, y = y0 + (domain_->j(n) + 0.5) * h;
        wallDistance[n] = std::min(r ? r->distance(x, y) : width - x, bank_.distance(x, y)) / h;
    }
    solver_->setWallDistance(std::move(wallDistance));

    // ---- inlet ----
    const auto& links = domain_->links();
    for (size_t li = 0; li < links.size(); ++li) {
        if (links[li].kind != LinkKind::MovingWall) continue;
        const BoundaryLink& L = links[li];
        if (L.side == XMax) {   // the side inflow (designed river only)
            sideLinks_.push_back(static_cast<int>(li));
            sideLinkY_.push_back(y0 + (domain_->j(L.node) + 0.5 + 0.5 * cy[L.dir]) * h);
            continue;
        }
        inletLinks_.push_back(static_cast<int>(li));
        inletLinkX_.push_back(x0 + (domain_->k(L.node) + 0.5 + 0.5 * cx[L.dir]) * h);
    }
    std::vector<double> nodeS, nodeWidth;
    for (int k = 0; k < u.nx; ++k) {
        const int n = domain_->nodeAt(k, u.ny - 1);
        if (n < 0) continue;
        inletNodes_.push_back(n);
        nodeS.push_back(inflowS(x0 + (k + 0.5) * h));
        nodeWidth.push_back(h);
    }
    inflow_ = std::make_unique<InflowModel>(params_.inflow, params_.seed, nodeS, nodeWidth);
    if (topOutlet)
        for (int k = 0; k < u.nx; ++k)
            if (domain_->nodeAt(k, 0) >= 0) outletNodes_.push_back(domain_->nodeAt(k, 0));
    if (topOutlet && outletNodes_.empty())
        throw std::invalid_argument("RiverSimulation: the outlet is closed (no water in the top row)");
    if (side) {
        // The right column: nodes by the side inlet (with the inflow's own random
        // profile, from a separate stream) and by the side outlet above it.
        std::vector<double> sideNodeS, sideWidth;
        for (int j = 0; j < u.ny; ++j) {
            const int n = domain_->nodeAt(u.nx - 1, j);
            if (n < 0) continue;
            const double y = y0 + (j + 0.5) * h;
            if (y < params_.side.top) {
                if (params_.side.outlet && y < params_.side.outletBottom)
                    sideOutletNodes_.push_back(n);
            } else if (y <= params_.side.bottom) {
                sideInletNodes_.push_back(n);
                sideNodeS.push_back(sideS(y));
                sideWidth.push_back(h);
            }
        }
        if (sideNodeS.empty()) throw std::invalid_argument("RiverSimulation: the side inflow is empty");
        if (!topOutlet && sideOutletNodes_.empty())
            throw std::invalid_argument(
                "RiverSimulation: with the top closed, the right edge needs an outlet");
        sideInflow_ = std::make_unique<InflowModel>(params_.inflow, hashKey(params_.seed, 9, 0, 0),
                                                    sideNodeS, sideWidth);
    }

    // ---- stirring (channel-fitted coordinates) ----
    const double H = height, W = width;
    stirring_ = std::make_unique<Stirring>(
        params_.stirring, params_.seed, *domain_,
        [x0, y0, h, &b, r, W, H](double X, double Y, double& xi, double& eta) {
            const double x = x0 + X * h, y = y0 + Y * h;
            const double xb = b.x(y);
            xi = (x - xb) / ((r ? r->x(y) : W) - xb);
            eta = y / H;
        });
    setupForces();
    setupDepth();

    tracers_ = std::make_unique<Tracers>(params_.tracers, params_.seed, bank_, width, height);
    if (box.width > 0.0 && box.height > 0.0)
        tracers_->setObstacle(box.x, box.y, box.x + box.width, box.y + box.height);
    if (r) tracers_->setRightBank(r);
    if (side && params_.side.outlet)
        tracers_->setRightOutlet(std::min(params_.side.top, params_.side.outletBottom));
    gridU_.assign(static_cast<size_t>(u.nx) * u.ny, 0.0);
    gridV_ = gridU_;
}

// The body forces: stirring alone keeps its original path; a side push is precomputed
// and, with stirring on too, summed into one field whenever stirring changes.
void RiverSimulation::setupForces() {
    const PushParams& push = params_.push;
    if (push.strength <= 0.0) {
        if (params_.stirring.enabled) solver_->setForceField(&stirring_->fx(), &stirring_->fy());
        return;
    }
    const double dx = push.aimX - push.x, dy = push.aimY - push.y;
    const double length = std::hypot(dx, dy);
    if (!(push.radius > 0.0) || !(length > 0.0))
        throw std::invalid_argument("RiverSimulation: the push needs a radius and an aim point");
    if (!(push.pulse >= 0.0 && push.pulse <= 1.0) || !(push.period > 0.0))
        throw std::invalid_argument("RiverSimulation: push pulse must be in [0, 1], period > 0");
    // Lattice force at the centre: `strength` inlet speeds of velocity per second.
    const double f0 = push.strength * params_.latticeInletSpeed / units_.stepsPerSecond;
    const int count = domain_->fluidCount();
    pushFx_.assign(count, 0.0);
    pushFy_.assign(count, 0.0);
    for (int n = 0; n < count; ++n) {
        const double x = units_.x0 + (domain_->k(n) + 0.5) * units_.h;
        const double y = units_.y0 + (domain_->j(n) + 0.5) * units_.h;
        const double r2 = (x - push.x) * (x - push.x) + (y - push.y) * (y - push.y);
        const double w = f0 * std::exp(-r2 / (2.0 * push.radius * push.radius));
        pushFx_[n] = w * dx / length;
        pushFy_[n] = w * dy / length;   // lattice rows grow downward, like page y
    }
    if (params_.stirring.enabled || push.pulse > 0.0) {
        forceFx_ = pushFx_;
        forceFy_ = pushFy_;
        solver_->setForceField(&forceFx_, &forceFy_);
    } else {
        solver_->setForceField(&pushFx_, &pushFy_);
    }
}

namespace {

// 0 below -width/2, 1 above width/2, smooth in between (a step for width 0).
double smoothStep(double d, double width) {
    if (!(width > 0.0)) return d > 0.0 ? 1.0 : 0.0;
    const double t = std::clamp(d / width + 0.5, 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

}  // namespace

RiverSimulation::DepthWeights RiverSimulation::depthWeights(double x, double y) const {
    const DepthParams& d = params_.depth;
    const double right = smoothStep(x - d.sillX, d.stepWidth);         // right of the tip
    const double above = smoothStep(d.sillY - y, d.stepWidth);         // above the sill line
    const double upper = smoothStep(d.curveBottom - y, d.stepWidth);   // above layer 2's edge
    DepthWeights w{};
    w.layer[2] = right * above;
    w.layer[1] = (1.0 - right) * upper;
    w.layer[0] = 1.0 - w.layer[1] - w.layer[2];
    if (d.sillWidth > 0.0) {
        // The sill: a band along y = sillY from the right edge to a rounded tip at sillX.
        const double half = 0.5 * d.sillWidth;
        const double along = std::max(0.0, d.sillX + half - x);
        const double distance = std::hypot(along, y - d.sillY);
        w.sill = 1.0 - smoothStep(distance - half, d.stepWidth);
    }
    return w;
}

double RiverSimulation::depthAt(double x, double y) const {
    const DepthParams& d = params_.depth;
    const DepthWeights w = depthWeights(x, y);
    const double layers = w.layer[0] * d.depth1 + w.layer[1] * d.depth2 + w.layer[2] * d.depth3;
    return (1.0 - w.sill) * layers + w.sill * d.sillDepth;
}

int RiverSimulation::layerAt(double x, double y) const {
    const DepthWeights w = depthWeights(x, y);
    if (w.sill >= 0.5) return 0;
    return static_cast<int>(std::max_element(w.layer, w.layer + 3) - w.layer) + 1;
}

// The bottom drag of the hidden bed, per water node: bedDrag / depth^2, in lattice units
// per step (capped at 1, which already stops the water within a few steps).
void RiverSimulation::setupDepth() {
    const DepthParams& d = params_.depth;
    if (!d.enabled) return;
    if (!(d.bedDrag >= 0.0) || !(d.depth1 > 0.0) || !(d.depth2 > 0.0) || !(d.depth3 > 0.0) ||
        !(d.sillDepth > 0.0) || !(d.sillWidth >= 0.0) || !(d.stepWidth >= 0.0))
        throw std::invalid_argument(
            "RiverSimulation: depths must be > 0; bed drag, sill width and step width >= 0");
    const int count = domain_->fluidCount();
    drag_.assign(count, 0.0);
    for (int n = 0; n < count; ++n) {
        const double x = units_.x0 + (domain_->k(n) + 0.5) * units_.h;
        const double y = units_.y0 + (domain_->j(n) + 0.5) * units_.h;
        const double depth = depthAt(x, y);
        drag_[n] = std::min(1.0, d.bedDrag / (depth * depth) * units_.dt);
    }
    solver_->setDragField(&drag_);
}

// The combined field: the push, scaled by its pulse, plus stirring (if on).
void RiverSimulation::updateForces() {
    const PushParams& push = params_.push;
    const double g = push.pulse > 0.0 ?
        1.0 + push.pulse * std::sin(2.0 * pi * time_ / push.period) : 1.0;
    const bool stir = params_.stirring.enabled;
    for (size_t n = 0; n < forceFx_.size(); ++n) {
        forceFx_[n] = g * pushFx_[n] + (stir ? stirring_->fx()[n] : 0.0);
        forceFy_[n] = g * pushFy_[n] + (stir ? stirring_->fy()[n] : 0.0);
    }
}

double RiverSimulation::inflowS(double x) const {
    const double xb = bank_.x(units_.height);
    return (x - xb) / (rightX(units_.height) - xb);
}

double RiverSimulation::meanInletSpeed(double t) const {
    const double T = params_.rampTime;
    const double ramp = (T > 0.0 && t < T) ? 0.5 * (1.0 - std::cos(pi * t / T)) : 1.0;
    return params_.inletSpeed * ramp;
}

void RiverSimulation::bottomInflow(double x, double ubar, double& across, double& up) const {
    const double s = inflowS(x);
    across = inflow_->transverse(s, ubar);
    up = inflow_->streamwise(s, ubar);
    const ConvergeParams& focus = params_.converge;
    const CornerJetParams& jet = params_.cornerJet;
    if (!focus.enabled && !(jet.width > 0.0)) return;
    double tilt = 0.0;   // sideways over up
    if (focus.enabled) {
        const double maxTilt = std::tan(focus.maxAngle * pi / 180.0);
        tilt = std::clamp((focus.x - x) / (units_.height - focus.y), -maxTilt, maxTilt);
    }
    if (jet.width > 0.0) {
        const double w = 1.0 - smoothStep(x - (bank_.x(units_.height) + jet.width), jet.edge);
        up *= 1.0 + w * (jet.speed / params_.inletSpeed - 1.0);
        tilt += w * (std::tan(jet.angle * pi / 180.0) - tilt);
    }
    across += up * tilt;
}

void RiverSimulation::updateInlet() {
    const double ubar = meanInletSpeed(time_);
    const double scale = 1.0 / units_.velocityScale;
    const ConvergeParams& focus = params_.converge;
    const double maxTilt = std::tan(focus.maxAngle * pi / 180.0);
    for (size_t i = 0; i < inletLinks_.size(); ++i) {
        double across, up;
        bottomInflow(inletLinkX_[i], ubar, across, up);
        // Inflow moves toward decreasing y (up the page).
        solver_->setWallVelocity(inletLinks_[i], across * scale, -up * scale);
    }
    if (sideInflow_) {
        // Leftward (toward decreasing x), tilted up (toward decreasing y) by the angle,
        // or toward the focus.
        const double sideBar = params_.side.speed * ubar / params_.inletSpeed;   // same ramp
        const double tilt = std::tan(params_.side.angle * pi / 180.0);
        const double edge = units_.x0 + units_.nx * units_.h;
        for (size_t i = 0; i < sideLinks_.size(); ++i) {
            const double y = sideLinkY_[i];
            const double s = sideS(y);
            const double along = sideInflow_->streamwise(s, sideBar);
            const double lean = focus.enabled && focus.side ?
                std::clamp((y - focus.y) / std::max(edge - focus.x, 1e-9), -maxTilt, maxTilt) : tilt;
            const double up = lean * along + sideInflow_->transverse(s, sideBar);
            solver_->setWallVelocity(sideLinks_[i], -along * scale, -up * scale);
        }
    }
}

void RiverSimulation::stepLattice() {
    const int stride = params_.coupledInflow && params_.nestedBaseRows > 0 ?
        params_.rows / params_.nestedBaseRows : 1;
    if (solver_->time() % stride == 0) {
        inflow_->advance(units_.dt*stride);
        if (sideInflow_) sideInflow_->advance(units_.dt*stride);
    }
    updateInlet();
    if (params_.stirring.enabled) {
        // Force that changes the velocity by ~one inlet speed over one correlation time.
        const double forceScale =
            params_.latticeInletSpeed / (params_.stirring.correlationTime * units_.stepsPerSecond);
        if (solver_->time() % stride == 0) {
            stirring_->advance(units_.dt*stride, forceScale);
        }
    }
    if (!forceFx_.empty()) updateForces();
    solver_->step();
    time_ += units_.dt;
    if (!outletNodes_.empty() || !sideOutletNodes_.empty()) {
        // Two-bank rivers: an absorbing outlet. The big lake and the outlet channel form
        // a slow sloshing mode (the lattice fluid is slightly compressible), which a
        // fixed outlet pressure reflects and the start-up ramp can pump until the flow
        // breaks down. Setting the outlet density from the excess outflow, with the
        // plane-wave impedance (delta rho = rho delta u / c_s), lets that energy leave.
        // In steady flow the outflow equals the inflow and the density returns to 1.
        const auto& rho = solver_->rho();
        const auto& uy = solver_->uy();
        const auto& ux = solver_->ux();
        double out = 0.0, in = 0.0;
        for (int n : outletNodes_) out -= rho[n] * uy[n];
        for (int n : inletNodes_) in -= rho[n] * uy[n];
        for (int n : sideOutletNodes_) out += rho[n] * ux[n];   // rightward, out
        for (int n : sideInletNodes_) in -= rho[n] * ux[n];     // leftward, in
        const double excess =
            (out - in) / static_cast<double>(outletNodes_.size() + sideOutletNodes_.size());
        solver_->setOutletDensity(solver_->rho0() + excess / std::sqrt(cs2));
    }
}

void RiverSimulation::advanceLatticeSteps(int steps) {
    if (steps < 0) throw std::invalid_argument("negative research step count");
    for (int i=0;i<steps;++i) stepLattice();
    updateVelocityGrid();
}

void RiverSimulation::warmUp() {
    while (time_ < params_.warmupTime) stepLattice();
    finishWarmUp();
}

void RiverSimulation::exportFlow(float* out) const {
    const size_t cells = static_cast<size_t>(units_.nx) * units_.ny;
    std::fill(out, out + cells, std::numeric_limits<float>::quiet_NaN());
    std::fill(out + cells, out + 3 * cells, 0.0f);
    for (int n = 0; n < domain_->fluidCount(); ++n) {
        const size_t idx = static_cast<size_t>(domain_->j(n)) * units_.nx + domain_->k(n);
        out[idx] = static_cast<float>(solver_->rho()[n]);
        out[cells + idx] = static_cast<float>(solver_->ux()[n]);
        out[2 * cells + idx] = static_cast<float>(solver_->uy()[n]);
    }
}

void RiverSimulation::startFromFlow(const float* in) {
    if (warmedUp_) throw std::logic_error("startFromFlow after the river started");
    const size_t cells = static_cast<size_t>(units_.nx) * units_.ny;
    for (int n = 0; n < domain_->fluidCount(); ++n) {
        const size_t idx = static_cast<size_t>(domain_->j(n)) * units_.nx + domain_->k(n);
        if (std::isnan(in[idx]))
            solver_->initializeNode(n, solver_->rho0(), 0.0, 0.0);
        else
            solver_->initializeNode(n, in[idx], in[cells + idx], in[2 * cells + idx]);
    }
    time_ = std::max(params_.warmupTime, params_.rampTime);
    finishWarmUp();
}

void RiverSimulation::finishWarmUp() {
    warmedUp_ = true;
    updateVelocityGrid();
    const FillParams& fill = params_.fill;
    if (fill.radius > 0.0) {
        tracers_->fillDisc(fill.x, fill.y, fill.radius,
                           [this](double x, double y) { return fillInside(x, y); });
        // Cells about two dot spacings across: small enough to find a bare patch.
        fillCell_ = 2.0 * params_.tracers.spacing;
        fillX0_ = fill.x - fill.radius;
        fillY0_ = fill.y - fill.radius;
        fillNx_ = fillNy_ = static_cast<int>(std::ceil(2.0 * fill.radius / fillCell_));
        fillTargets_ = tracers_->countCells(fillX0_, fillY0_, fillCell_, fillNx_, fillNy_);
    }
}

// Water that the fill may place dots in: between the banks, on the page, in view.
bool RiverSimulation::fillInside(double x, double y) const {
    const bool water = y > 0.0 && y < units_.height && x > bank_.x(y) && x < rightX(y);
    const bool visible = !hasView_ || (x <= viewRight_ && y <= viewBottom_);
    return water && visible;
}

void RiverSimulation::advance(double seconds) {
    if (!warmedUp_) warmUp();
    seconds = std::min(seconds, 0.1);   // a stalled tab must not trigger a huge catch-up
    carry_ += seconds;
    const long steps = static_cast<long>(carry_ / units_.dt);
    carry_ -= steps * units_.dt;
    for (long s = 0; s < steps; ++s) stepLattice();
    updateVelocityGrid();
    tracers_->advance(seconds, [this](double x, double y, double& u, double& v) {
        velocityAt(x, y, u, v);
    });
    emitTracers(seconds);
    const FillParams& fill = params_.fill;
    if (!fillTargets_.empty() && time_ < params_.warmupTime + fill.hold) {
        // Keep the disc's total at its count after the fill, placing the replacements in
        // the cells furthest below their own count after the fill (new dots fade in).
        // Capping the total matters: the swirl moves dots between cells, and topping up
        // every short cell without it would let the swirl pile up dense clumps.
        const std::vector<int> counts =
            tracers_->countCells(fillX0_, fillY0_, fillCell_, fillNx_, fillNy_);
        int missing = 0;
        std::vector<std::pair<int, size_t>> short_;   // (deficit, cell)
        for (size_t c = 0; c < counts.size(); ++c) {
            missing += fillTargets_[c] - counts[c];
            if (counts[c] < fillTargets_[c]) short_.push_back({fillTargets_[c] - counts[c], c});
        }
        std::sort(short_.begin(), short_.end(), std::greater<>());
        const auto inside = [this](double x, double y) { return fillInside(x, y); };
        for (const auto& [deficit, c] : short_) {
            if (missing <= 0) break;
            const int cx = static_cast<int>(c % fillNx_), cy = static_cast<int>(c / fillNx_);
            missing -= tracers_->fillBox(fillX0_ + cx * fillCell_, fillY0_ + cy * fillCell_,
                                         fillCell_, fillCell_, inside, std::min(deficit, missing));
        }
    }
}

void RiverSimulation::updateVelocityGrid() {
    std::fill(gridU_.begin(), gridU_.end(), 0.0);
    std::fill(gridV_.begin(), gridV_.end(), 0.0);
    const double scale = units_.velocityScale;
    for (int n = 0; n < domain_->fluidCount(); ++n) {
        const size_t idx = static_cast<size_t>(domain_->j(n)) * units_.nx + domain_->k(n);
        gridU_[idx] = solver_->ux()[n] * scale;
        gridV_[idx] = solver_->uy()[n] * scale;
    }
}

void RiverSimulation::velocityAt(double x, double y, double& u, double& v) const {
    const RiverUnits& U = units_;
    const double gx = (x - U.x0) / U.h - 0.5, gy = (y - U.y0) / U.h - 0.5;
    const int k0 = static_cast<int>(std::floor(gx)), j0 = static_cast<int>(std::floor(gy));
    const double tx = gx - k0, ty = gy - j0;
    const double ubar = meanInletSpeed(time_);

    auto node = [&](int k, int j, double& nu, double& nv) {
        if (j >= U.ny) {   // below the inlet: the prescribed inflow
            double up;
            bottomInflow(U.x0 + (k + 0.5) * U.h, ubar, nu, up);
            nv = -up;
            return;
        }
        if (j < 0) j = 0;   // above the outlet: zero gradient
        double sign = 1.0;
        if (k >= U.nx) { k = 2 * U.nx - 1 - k; sign = -1.0; }   // mirror: zero at the wall
        if (k < 0) { nu = nv = 0.0; return; }
        const size_t idx = static_cast<size_t>(j) * U.nx + k;
        nu = sign * gridU_[idx];
        nv = sign * gridV_[idx];
    };
    double u00, v00, u10, v10, u01, v01, u11, v11;
    node(k0, j0, u00, v00);
    node(k0 + 1, j0, u10, v10);
    node(k0, j0 + 1, u01, v01);
    node(k0 + 1, j0 + 1, u11, v11);
    u = (1 - ty) * ((1 - tx) * u00 + tx * u10) + ty * ((1 - tx) * u01 + tx * u11);
    v = (1 - ty) * ((1 - tx) * v00 + tx * v10) + ty * ((1 - tx) * v01 + tx * v11);
}

void RiverSimulation::setView(double right, double bottom) {
    hasView_ = true;
    viewRight_ = right;
    viewBottom_ = bottom;
    tracers_->setView(right, bottom);
    int visible = 0;
    for (int n = 0; n < domain_->fluidCount(); ++n) {
        const double x = units_.x0 + (domain_->k(n) + 0.5) * units_.h;
        const double y = units_.y0 + (domain_->j(n) + 0.5) * units_.h;
        if (x <= right && y <= bottom) ++visible;
    }
    viewArea_ = visible * units_.h * units_.h;
}

// Dots enter where water flows into the view: up through its bottom edge and left
// through its right edge, in proportion to the inward speed.
void RiverSimulation::emitIntoView(double dt, double ubar) {
    const RiverUnits& U = units_;
    const double h = U.h;
    std::vector<EntrySegment> segments;
    const double yb = std::min(viewBottom_, U.height);
    if (viewBottom_ >= U.height) {
        // The view reaches the inlet: use the prescribed inflow, as without a view.
        for (int n : inletNodes_) {
            const double x = U.x0 + (domain_->k(n) + 0.5) * h;
            if (x > viewRight_) continue;
            double across, up;
            bottomInflow(x, ubar, across, up);
            segments.push_back({x, U.height, h, up, true});
        }
    } else {
        const double xr = std::min(viewRight_, rightX(yb));
        for (double x = bank_.x(yb) + 0.5 * h; x < xr; x += h) {
            double vx, vy;
            velocityAt(x, yb, vx, vy);
            if (vy < 0.0) segments.push_back({x, yb, h, -vy, true});
        }
    }
    for (double y = 0.5 * h; y < yb; y += h) {
        if (viewRight_ <= bank_.x(y) || viewRight_ >= rightX(y)) continue;
        // The designed river: beyond the view's right edge, new water only comes from
        // the side inflow band. Above or below it, water crossing into view is the
        // river's own recirculation, and brings no dots of its own.
        if (sideInflow_ && (y < params_.side.top || y > sideBottom())) continue;
        double vx, vy;
        velocityAt(viewRight_, y, vx, vy);
        if (vx < 0.0) segments.push_back({viewRight_, y, h, -vx, false});
    }
    tracers_->emitAlong(segments, dt);
}

void RiverSimulation::emitTracers(double dt) {
    const double ubar = meanInletSpeed(time_);
    if (hasView_) {
        emitIntoView(dt, ubar);
        return;
    }
    if (sideInflow_) {
        // Without a view, dots also enter with the side inflow at the right edge.
        std::vector<EntrySegment> segments;
        const double sideBar = params_.side.speed * ubar / params_.inletSpeed;
        for (int n : sideInletNodes_) {
            const double y = units_.y0 + (domain_->j(n) + 0.5) * units_.h;
            segments.push_back({units_.width, y, units_.h, sideInflow_->streamwise(sideS(y), sideBar), false});
        }
        tracers_->emitAlong(segments, dt);
    }
    std::vector<InletColumn> columns;
    columns.reserve(inletNodes_.size());
    for (int n : inletNodes_) {
        const double x = units_.x0 + (domain_->k(n) + 0.5) * units_.h;
        double across, up;
        bottomInflow(x, ubar, across, up);
        columns.push_back({x, units_.h, up});
    }
    tracers_->emit(columns, dt);
}

double RiverSimulation::area() const {
    return domain_->fluidCount() * units_.h * units_.h;
}

RiverDiagnostics RiverSimulation::diagnostics() const {
    RiverDiagnostics d{};
    d.time = time_;
    d.steps = solver_->time();
    d.targetFlux = meanInletSpeed(time_) * units_.inletWidth;
    if (params_.cornerJet.width > 0.0)   // the jet's extra flux over its (blended) width
        d.targetFlux += (params_.cornerJet.speed / params_.inletSpeed - 1.0) *
                        meanInletSpeed(time_) * params_.cornerJet.width;
    const auto& rho = solver_->rho();
    const auto& ux = solver_->ux();
    const auto& uy = solver_->uy();
    const double scale = units_.velocityScale * units_.h;
    d.finite = true;
    for (int n = 0; n < domain_->fluidCount(); ++n) {
        const double speed = std::sqrt(ux[n] * ux[n] + uy[n] * uy[n]);
        if (!std::isfinite(rho[n]) || !std::isfinite(speed)) d.finite = false;
        d.maxMach = std::max(d.maxMach, speed / std::sqrt(cs2));
        d.maxDensityDeviation = std::max(d.maxDensityDeviation, std::abs(rho[n] - 1.0));
        // Upward (toward y = 0) volume flux through the first and last rows.
        if (domain_->j(n) == units_.ny - 1) d.inflow += -rho[n] * uy[n] * scale;
        if (domain_->j(n) == 0) d.outflow += -rho[n] * uy[n] * scale;
    }
    if (sideInflow_) {
        // The designed river's right edge: side inflow below `top`, outflow above.
        for (int n : sideInletNodes_) d.inflow += -rho[n] * ux[n] * scale;
        for (int n : sideOutletNodes_) d.outflow += rho[n] * ux[n] * scale;
        d.targetFlux += params_.side.speed * meanInletSpeed(time_) / params_.inletSpeed *
                        (sideBottom() - params_.side.top);
    }
    d.maxTauEffective = solver_->maxTauEffective();
    d.tracers = tracers_->count();
    d.expectedTracers = viewArea() / (params_.tracers.spacing * params_.tracers.spacing);
    return d;
}

}  // namespace lbm
