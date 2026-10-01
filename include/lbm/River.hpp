#pragma once

// The homepage river: a lattice Boltzmann channel between the spline bank (left) and
// either the page's right edge (a straight wall, the default) or a second spline bank
// (the lake river), with an inlet along the bottom and an outlet along the top, plus
// the dots. Everything is in the page's fixed reference frame
// (x right, y down); resizing or zooming the page never touches this object.
//
// Units. With lattice spacing h = height / rows page units, a chosen mean inlet speed
// U (page units/s) and its lattice value u:
//   dt = h u / U seconds per lattice step
//   nu_lattice = u (W_in / h) / Re,  tau = 3 nu_lattice + 1/2
// In 2D the flux is the same through every cross-section, so Re = Q / nu holds
// everywhere along the river.

#include "lbm/BankSpline.hpp"
#include "lbm/Collision.hpp"
#include "lbm/Domain.hpp"
#include "lbm/Solver.hpp"
#include "lbm/Stochastic.hpp"
#include "lbm/Tracers.hpp"

#include <cstdint>
#include <memory>
#include <vector>

namespace lbm {

// A push that drives the water toward an aim point, to force eddies where visitors
// look. It is a body force with a Gaussian profile; no water is added. The defaults
// sit at the mouth of the pocket behind the left bank's upper bend and push into it:
// the water swirls into the pocket, down along the bank and back out to the outlet
// (the author's choice; placed just past the view's right edge, it only made an
// eddy that the current carried out of the top-right corner).
struct PushParams {
    double strength = 0.0;             // velocity gained per second at the centre, in
                                       // inlet speeds; 0 = off
    double x = 1000.0, y = 110.0;      // centre, page units
    double radius = 80.0;              // Gaussian radius, page units
    double aimX = 880.0, aimY = 320.0; // it pushes from the centre toward this point
    // Pulsing: strength x (1 + pulse sin(2 pi t / period)). A steady push only makes a
    // steady current; each surge of a pulsing one sheds a pair of eddies.
    double pulse = 0.0;                // 0 = steady, 1 = swings between 0 and 2x
    double period = 4.0;               // seconds
};

// Dots placed at once when the dots start (after the warm-up), inside a soft-edged
// disc, so a swirl is populated from the first frame instead of waiting for the dots
// that enter at the view's edges to reach it. Only water inside the view is filled.
struct FillParams {
    double radius = 0.0;               // page units; 0 = off
    double x = 1010.0, y = 210.0;      // centre, page units (the upper bend's pocket)
    // For this long after the dots start, the disc is kept populated: it is split into
    // small cells, and a cell that has lost dots since the fill is topped back up (new
    // dots fade in). A nearly closed swirl rarely takes in dots from the main current,
    // so without this it slowly empties.
    double hold = 0.0;                 // seconds; 0 = fill once
};

// The designed river: instead of a wall, the right edge takes water in below `top`
// (flowing left, tilted up by `angle`) and lets it out above. With the bottom inlet
// this makes two flows collide, and both leave through the top and the upper right.
// The domain extends `margin` past the page's right edge and continues the left bank
// above the page for an outlet channel (`outletCells` rows), both out of sight.
struct SideInflowParams {
    bool enabled = false;
    double top = 500.0;       // page y: inflow below, outflow above
    double bottom = 1e9;      // page y where the side inflow ends (a wall below it);
                              // past the page bottom = all the way down
    double speed = 30.0;      // mean leftward speed of the side inflow, page units/s
    double angle = 20.0;      // degrees the side inflow tilts upward from leftward
    double margin = 120.0;    // page units the domain extends past the page width
    // Above `top`, the right edge is an outlet (true) or a wall (false). A wall makes
    // the side flow cross the river to the top outlet instead of short-circuiting out.
    bool outlet = true;
    // With the outlet on: page y where it ends, with a wall from there down to `top`, so
    // the outlet can sit above a sill while the inflow band stays low. Default: no wall.
    double outletBottom = 1e9;
    // false: the page top is a wall (no outlet channel above it), so all water leaves
    // through the right edge's outlet and crosses the upper right on the way.
    bool topOutlet = true;
    // An invisible lid above the page: the top outlet ends at this page x, so water
    // right of it must travel left under the lid before it can leave. Default: none.
    double outletRight = 1e9;
};

// A hidden bed depth, the author's idea for steering the flow without visible
// geometry. The solver is 2D, so depth acts only as bottom drag: the rate
// bedDrag / depth^2 (per second) slows shallow water most, so water favours deep
// paths. (Real shallow water also speeds up over shallows; that is not modelled.)
// Three layers, blended over stepWidth: the inflow (1), the upper curve (2) and the
// outflow (3). Right of sillX, layer 3 lies above sillY and layer 1 below it; left of
// sillX, layer 2 lies above curveBottom and layer 1 below it. Where layers 1 and 3
// meet, a shallow sill runs from the right edge left to its tip at (sillX, sillY), so
// water cannot skip layer 2: it passes under the sill, up past the curve, and back
// over the sill to the outflow.
struct DepthParams {
    bool enabled = false;
    // At 0.5/s the shallow inflow layer lost its eddies; 0.15/s keeps them.
    double bedDrag = 0.15;             // drag rate at depth 1, per second
    double depth1 = 1.0, depth2 = 2.0, depth3 = 3.0;   // relative depths of the layers
    // Measured on the page river: 60 wide at depth 0.05 let 38% of the flow across;
    // 100 wide at depth 0.02 lets about 4% across.
    double sillDepth = 0.02;           // relative depth on the sill
    double sillX = 1250.0, sillY = 420.0;   // the sill's tip, page units
    double sillWidth = 100.0;          // page units across the sill; 0 = no sill
    double curveBottom = 650.0;        // page y of layer 2's lower edge, left of sillX
    double stepWidth = 40.0;           // page units over which the depth changes
};

// The inflow aims at a focus point, like a lake draining into a river's mouth (the
// author's drawing: the bottom flow runs up along the bank near it and up-left near the
// right edge, and the currents gather and rise). Each point of the bottom inlet, and
// with `side` each point of the side band, tilts toward the focus, up to maxAngle.
// Only the direction changes; the flux through the edge does not.
struct ConvergeParams {
    bool enabled = false;
    double x = 1450.0, y = -300.0;     // the focus, page units (y < 0: above the page)
    double maxAngle = 45.0;            // largest tilt from the edge's normal, degrees
    bool side = false;                 // aim the side inflow at the focus instead of its angle
};

// A faster stretch of the bottom inflow next to the left bank, angled away from it:
// a second, smaller side jet at the bottom-left corner (the author's request). Its
// shear against the slower water beside it rolls up into eddies, and a small swirl
// forms between it and the bank. It keeps the inflow's wandering lanes (it scales
// their speed), so its eddies never repeat. Only the inflow changes.
struct CornerJetParams {
    double width = 0.0;        // page units of the bottom inlet from the bank; 0 = off
    double speed = 80.0;       // its mean speed, page units/s (the rest keeps inletSpeed)
    double angle = 45.0;       // degrees to the right of straight up, away from the bank
    double edge = 30.0;        // page units over which it blends into the rest
};

// A solid box in the water: the page's control box (the author's "Flow" toggle) floats
// on it, so the water parts around it and sheds eddies. Page units; a zero width or
// height means none.
struct ObstacleParams {
    double x = 0.0, y = 0.0;           // top-left corner
    double width = 0.0, height = 0.0;
};

struct RiverParams {
    int rows = 112;                    // lattice rows over the page height
    double cell = 0.0;                 // > 0: page units per lattice cell (sets rows)
    // Two-bank rivers only: the banks continue along their end tangents for this many
    // lattice rows above the page, so eddies shed near the page top die out before the
    // outlet. (Recirculation reaching the pressure outlet made the flow blow up.)
    int outletCells = 24;
    int nestedBaseRows = 0;            // >0: integer refinement of this exact grid
    int threads = 1;                   // native only; WASM stays single-threaded
    bool coupledInflow = false;        // update OU amplitudes at nested-base cadence
    double inletSpeed = 30.0;          // mean inflow speed, page units per second
    double reynolds = 1000.0;          // Re = Q / nu
    double latticeInletSpeed = 0.015;  // mean inflow in lattice units (sets the Mach number)
    CollisionConfig collision = defaultCollision();   // tau is derived from reynolds
    // The ramp must be long compared with the channel's sound period (~3 s here), or it
    // leaves the outlet flow sloshing by several percent for minutes.
    double rampTime = 6.0;             // seconds to ramp the inflow up from rest
    double warmupTime = 7.0;           // hidden seconds before dots appear
    std::uint64_t seed = 1;
    InflowNoiseParams inflow;
    StirringParams stirring;
    PushParams push;
    FillParams fill;
    SideInflowParams side;
    DepthParams depth;
    ConvergeParams converge;
    CornerJetParams cornerJet;
    ObstacleParams box;
    double thermalKT = 0.0;            // lattice units; 0 = off (research only)
    TracerParams tracers;

    static CollisionConfig defaultCollision() {
        CollisionConfig c;
        c.model = CollisionModel::MRT;
        c.closure = ClosureKind::Smagorinsky;
        c.smagorinsky = 0.12;
        // Near tau = 1/2 the magic-parameter rate for the odd ghosts falls to ~0.07, and
        // those barely-damped modes drew grid-scale stripes behind the bend. A fixed
        // rate removes them (Poiseuille exactness is not needed at this resolution).
        c.oddGhostRate = 1.2;
        return c;
    }
};

struct RiverUnits {
    double width, height;   // page reference frame
    double h;               // page units per lattice spacing
    double x0;              // page x of the lattice's left edge
    double y0;              // page y of the lattice's top edge (< 0 with an outlet section)
    int nx, ny;
    double dt;              // seconds per lattice step
    double stepsPerSecond;
    double inletWidth;      // page units
    double nu, tau;         // lattice units
    double velocityScale;   // page units per second, per lattice velocity unit
};

struct RiverDiagnostics {
    double time;            // simulated seconds
    long long steps;
    double inflow, outflow; // page units^2 / s through the rows next to inlet and outlet
    double targetFlux;      // inlet speed x inlet width, once ramped
    double maxMach;         // max |u| / c_s
    double maxDensityDeviation;   // max |rho - rho0| / rho0
    double maxTauEffective;
    bool finite;
    int tracers;
    double expectedTracers; // density x river area
};

class RiverSimulation {
public:
    RiverSimulation(std::vector<BankSegment> bank, double width, double height,
                    const RiverParams& params);
    // The lake river: a right bank x = g(y) replaces the straight wall at x = width.
    RiverSimulation(std::vector<BankSegment> left, std::vector<BankSegment> right,
                    double width, double height, const RiverParams& params);

    // Run the hidden warm-up (called automatically by the first advance()).
    void warmUp();
    // A saved flow: density and velocity (lattice units) on the full grid, one plane
    // each, row-major (row j, column k at j * nx + k); NaN density marks a cell that is
    // not water. exportFlow() writes 3 * nx * ny floats.
    void exportFlow(float* out) const;
    // Start from a saved flow instead of the hidden warm-up: each water node takes its
    // cell's state at equilibrium; one whose cell was not water (the saved flow's own
    // obstacles, e.g. the Flow box elsewhere) starts at rest. The clock moves to the
    // end of the warm-up and the ramp, so the inflow is at full speed, and the dots
    // start as after warmUp(). Before the first advance() only.
    void startFromFlow(const float* in);
    // Advance by a frame of wall-clock time: lattice steps, then the dots.
    void advance(double seconds);
    // Research stepping: exact integer steps, no tracer advancement.
    void advanceLatticeSteps(int steps);
    Solver& researchSolver() { return *solver_; }

    double time() const { return time_; }
    const RiverUnits& units() const { return units_; }
    const RiverParams& params() const { return params_; }
    const Domain& domain() const { return *domain_; }
    const Solver& solver() const { return *solver_; }
    const BankSpline& bank() const { return bank_; }
    Tracers& tracers() { return *tracers_; }
    RiverDiagnostics diagnostics() const;

    // Interpolated flow velocity at a page point (page units per second), no-slip at walls.
    void velocityAt(double x, double y, double& u, double& v) const;
    // Area of the water in page units^2.
    double area() const;

    // Show only the box x <= right, y <= bottom (page units): dots enter along its
    // bottom and right edges where water flows into view, and leave when they exit it.
    // Changing the view never touches the flow. Without a view, dots enter at the inlet.
    void setView(double right, double bottom);
    bool hasView() const { return hasView_; }
    // Water area inside the view (page units^2); the whole area without a view.
    double viewArea() const { return hasView_ ? viewArea_ : area(); }
    const BankSpline* rightBank() const { return right_.get(); }
    // The hidden bed (DepthParams) at a page point: its relative depth, and which part
    // dominates there: 1, 2 or 3 for a layer, 0 for the sill.
    double depthAt(double x, double y) const;
    int layerAt(double x, double y) const;

private:
    void build(double width, double height);
    void setupForces();
    void setupDepth();
    struct DepthWeights { double layer[3], sill; };
    DepthWeights depthWeights(double x, double y) const;
    double rightX(double y) const { return right_ ? right_->x(y) : units_.width; }
    void emitIntoView(double dt, double ubar);
    void stepLattice();
    double meanInletSpeed(double t) const;
    void updateInlet();
    void updateVelocityGrid();
    void emitTracers(double dt);
    double inflowS(double x) const;
    // The prescribed bottom inflow at page x (page units/s): its sideways part (to the
    // right) and its upward speed, with the wandering lanes, the convergence tilt and
    // the corner jet.
    void bottomInflow(double x, double ubar, double& across, double& up) const;

    RiverParams params_;
    BankSpline bank_;
    std::unique_ptr<BankSpline> right_;
    RiverUnits units_{};
    std::unique_ptr<Domain> domain_;
    std::unique_ptr<Solver> solver_;
    std::unique_ptr<InflowModel> inflow_;
    std::unique_ptr<Stirring> stirring_;
    std::unique_ptr<Tracers> tracers_;

    std::vector<int> inletLinks_;          // link indices on the inlet side
    std::vector<double> inletLinkX_;       // page x where each crosses the inlet
    std::vector<int> inletNodes_;          // fluid nodes in the bottom row
    std::vector<int> outletNodes_;         // fluid nodes in the top row (two-bank rivers)
    // The side inflow (designed river): its links and the right column's nodes.
    std::unique_ptr<InflowModel> sideInflow_;
    std::vector<int> sideLinks_;
    std::vector<double> sideLinkY_;        // page y where each crosses the right edge
    std::vector<int> sideInletNodes_, sideOutletNodes_;
    double sideBottom() const { return std::min(params_.side.bottom, units_.height); }
    double sideS(double y) const { return (y - params_.side.top) / (sideBottom() - params_.side.top); }
    std::vector<double> gridU_, gridV_;    // page-unit velocity on the full lattice grid
    std::vector<double> pushFx_, pushFy_;  // the side push (lattice force per water node)
    std::vector<double> forceFx_, forceFy_;  // push (pulsed) + stirring, when combined
    std::vector<double> drag_;             // bottom drag per water node (lattice units)
    void updateForces();
    double time_ = 0.0;
    double carry_ = 0.0;
    bool warmedUp_ = false;
    void finishWarmUp();
    bool fillInside(double x, double y) const;
    // The fill's cells: a square grid over the disc, and each cell's count after the fill.
    double fillX0_ = 0.0, fillY0_ = 0.0, fillCell_ = 0.0;
    int fillNx_ = 0, fillNy_ = 0;
    std::vector<int> fillTargets_;
    bool hasView_ = false;
    double viewRight_ = 0.0, viewBottom_ = 0.0, viewArea_ = 0.0;
};

}  // namespace lbm
