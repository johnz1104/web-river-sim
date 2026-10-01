#pragma once

// Lattice geometry: which nodes hold fluid, where each population streams from,
// and every link that crosses a boundary (with the exact crossing fraction q).
//
// Coordinates are lattice units. Node (k, j) sits at X = k + 1/2, Y = j + 1/2, so the
// box sides X = 0, X = nx, Y = 0, Y = ny lie exactly halfway along boundary links.

#include "lbm/Lattice.hpp"

#include <cstdint>
#include <functional>
#include <vector>

namespace lbm {

enum class SideKind { Periodic, Wall, MovingWall, PressureOutlet };
enum class LinkKind : std::uint8_t { Wall, MovingWall, PressureOutlet };

// Box sides by lattice direction (not "top/bottom": the river's y grows downward).
enum Side : std::int8_t { XMin = 0, XMax = 1, YMin = 2, YMax = 3, Body = 4 };

// A link from fluid node `node` along outgoing direction `dir` that leaves the fluid.
struct BoundaryLink {
    int node;
    std::int8_t dir;
    LinkKind kind;
    Side side;        // which box side it crosses, or Body for the curved geometry
    int body;         // body id for force accounting (0 = not tracked)
    double q;         // fraction of the link inside the fluid, in (0, 1]
};

struct DomainSpec {
    int nx = 0, ny = 0;
    // Signed geometry: phi(X, Y) > 0 is fluid. Leave empty for an all-fluid box.
    // It must be defined (and continuous) slightly beyond the box.
    std::function<double(double, double)> phi;
    // Optional body id of the solid met at a crossing point (for drag/lift).
    std::function<int(double, double)> bodyAt;
    SideKind xMin = SideKind::Wall, xMax = SideKind::Wall;
    SideKind yMin = SideKind::Wall, yMax = SideKind::Wall;
    // Optional: the kind of box side `s` for a link leaving the node at (X, Y), so one
    // side can be part inlet, part outlet. Unset: every link takes the side's kind.
    std::function<SideKind(Side s, double X, double Y)> sideKindAt;
};

class Domain {
public:
    explicit Domain(const DomainSpec& spec);

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int fluidCount() const { return static_cast<int>(nodeK_.size()); }

    // Fluid index of grid node (k, j), or -1 for solid.
    int nodeAt(int k, int j) const { return index_[j * nx_ + k]; }
    int k(int n) const { return nodeK_[n]; }
    int j(int n) const { return nodeJ_[n]; }

    // Streaming source for incoming direction i at node n: a fluid node index (>= 0),
    // or -(link + 1) when the population must come from a boundary rule.
    int pull(int n, int i) const { return pull_[n * Q + i]; }

    const std::vector<BoundaryLink>& links() const { return links_; }
    SideKind side(Side s) const { return sides_[s]; }

private:
    int nx_, ny_;
    SideKind sides_[4];
    std::vector<int> index_;
    std::vector<int> nodeK_, nodeJ_;
    std::vector<int> pull_;
    std::vector<BoundaryLink> links_;
};

}  // namespace lbm
