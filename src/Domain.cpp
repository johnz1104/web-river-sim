#include "lbm/Domain.hpp"

#include <stdexcept>

namespace lbm {

namespace {

// Crossing fraction of phi along the segment p(s) = (X + s dx, Y + s dy), s in [0, 1],
// given phi(p(0)) > 0 >= phi(p(1)). Bisection is exact enough (2^-60) and needs no
// derivative of the geometry.
double crossingFraction(const std::function<double(double, double)>& phi, double X,
                        double Y, double dx, double dy) {
    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (phi(X + mid * dx, Y + mid * dy) > 0.0) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

LinkKind linkKind(SideKind s) {
    switch (s) {
        case SideKind::MovingWall: return LinkKind::MovingWall;
        case SideKind::PressureOutlet: return LinkKind::PressureOutlet;
        default: return LinkKind::Wall;
    }
}

// Where two boundaries are crossed at the same fraction (box corners), walls win so
// corners stay no-slip; then moving walls; then outlets.
int kindPriority(LinkKind k) {
    switch (k) {
        case LinkKind::Wall: return 0;
        case LinkKind::MovingWall: return 1;
        default: return 2;
    }
}

}  // namespace

Domain::Domain(const DomainSpec& spec) : nx_(spec.nx), ny_(spec.ny) {
    if (nx_ <= 0 || ny_ <= 0) throw std::invalid_argument("Domain: empty lattice");
    sides_[XMin] = spec.xMin;
    sides_[XMax] = spec.xMax;
    sides_[YMin] = spec.yMin;
    sides_[YMax] = spec.yMax;
    if ((spec.xMin == SideKind::Periodic) != (spec.xMax == SideKind::Periodic) ||
        (spec.yMin == SideKind::Periodic) != (spec.yMax == SideKind::Periodic))
        throw std::invalid_argument("Domain: periodic sides must come in pairs");

    auto isFluidAt = [&](double X, double Y) { return !spec.phi || spec.phi(X, Y) > 0.0; };

    index_.assign(static_cast<size_t>(nx_) * ny_, -1);
    for (int j = 0; j < ny_; ++j)
        for (int k = 0; k < nx_; ++k)
            if (isFluidAt(k + 0.5, j + 0.5)) {
                index_[j * nx_ + k] = static_cast<int>(nodeK_.size());
                nodeK_.push_back(k);
                nodeJ_.push_back(j);
            }

    const int n = fluidCount();
    pull_.assign(static_cast<size_t>(n) * Q, 0);
    for (int node = 0; node < n; ++node) {
        const int k0 = nodeK_[node], j0 = nodeJ_[node];
        const double X = k0 + 0.5, Y = j0 + 0.5;
        pull_[node * Q + 0] = node;

        for (int o = 1; o < Q; ++o) {   // outgoing direction; it feeds opposite[o] here
            int k1 = k0 + cx[o], j1 = j0 + cy[o];

            // Box sides: wrap periodic ones, record a half-link crossing for the rest.
            bool crossesSide = false;
            LinkKind sideKind = LinkKind::Wall;
            Side sideHit = Body;
            auto considerSide = [&](Side s) {
                const LinkKind lk = linkKind(spec.sideKindAt ? spec.sideKindAt(s, X, Y) : sides_[s]);
                if (!crossesSide || kindPriority(lk) < kindPriority(sideKind)) {
                    sideKind = lk;
                    sideHit = s;
                }
                crossesSide = true;
            };
            if (k1 < 0 || k1 >= nx_) {
                if (sides_[XMin] == SideKind::Periodic) k1 = (k1 + nx_) % nx_;
                else considerSide(k1 < 0 ? XMin : XMax);
            }
            if (j1 < 0 || j1 >= ny_) {
                if (sides_[YMin] == SideKind::Periodic) j1 = (j1 + ny_) % ny_;
                else considerSide(j1 < 0 ? YMin : YMax);
            }

            // Curved geometry: does the link end inside solid?
            bool crossesBody = false;
            double qBody = 1.0;
            if (spec.phi && spec.phi(X + cx[o], Y + cy[o]) <= 0.0) {
                crossesBody = true;
                qBody = crossingFraction(spec.phi, X, Y, cx[o], cy[o]);
            }

            if (!crossesSide && !crossesBody) {
                const int src = index_[j1 * nx_ + k1];
                if (src < 0)
                    throw std::logic_error("Domain: solid neighbour without a crossing");
                pull_[node * Q + opposite[o]] = src;
                continue;
            }

            BoundaryLink link{node, static_cast<std::int8_t>(o), LinkKind::Wall, Body, 0, 0.5};
            if (crossesBody && (!crossesSide || qBody <= 0.5)) {
                link.q = qBody;
                if (spec.bodyAt)
                    link.body = spec.bodyAt(X + qBody * cx[o], Y + qBody * cy[o]);
            } else {
                link.kind = sideKind;
                link.side = sideHit;
            }
            pull_[node * Q + opposite[o]] = -static_cast<int>(links_.size()) - 1;
            links_.push_back(link);
        }
    }
}

}  // namespace lbm
