#include "lbm/Stochastic.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace lbm {

// ---------------------------------------------------------------- inflow

InflowModel::InflowModel(const InflowNoiseParams& params, std::uint64_t seed,
                         std::vector<double> nodeS, std::vector<double> nodeWeight)
    : params_(params), rng_(seed, 1), nodeS_(std::move(nodeS)),
      nodeWeight_(std::move(nodeWeight)) {
    if (nodeS_.size() != nodeWeight_.size() || nodeS_.empty())
        throw std::invalid_argument("InflowModel: need one weight per inlet node");
    for (double wk : nodeWeight_) width_ += wk;

    const int M = params_.enabled ? params_.modes : 0;
    for (int m = 1; m <= M; ++m) {
        OrnsteinUhlenbeck ou;
        ou.sigma = params_.amplitude / m;
        const double f = M > 1 ? (m - 1.0) / (M - 1.0) : 0.0;
        ou.correlationTime = params_.maxCorrelationTime
                             + f * (params_.minCorrelationTime - params_.maxCorrelationTime);
        ou.value = ou.sigma * rng_.normal();   // start in the stationary distribution
        modes_.push_back(ou);
    }
    meander_.sigma = params_.enabled && params_.meander ? params_.meanderAmplitude : 0.0;
    meander_.correlationTime = params_.meanderCorrelationTime;
    meander_.value = meander_.sigma * rng_.normal();
    renormalize();
}

double InflowModel::shape(double s) const {
    if (s <= 0.0 || s >= 1.0) return 0.0;
    const double t = 2.0 * s - 1.0;
    const double t2 = t * t;
    const double g = 1.0 - t2 * t2 * t2;
    double b = 1.0;
    for (size_t m = 0; m < modes_.size(); ++m)
        b += modes_[m].value * std::cos((m + 1.0) * pi * s);
    return g * std::max(b, params_.floor);
}

void InflowModel::renormalize() {
    double sum = 0.0;
    for (size_t k = 0; k < nodeS_.size(); ++k) sum += shape(nodeS_[k]) * nodeWeight_[k];
    area_ = sum / width_;
}

void InflowModel::advance(double dtSeconds) {
    for (auto& ou : modes_) ou.step(dtSeconds, rng_.normal());
    if (meander_.sigma > 0.0) meander_.step(dtSeconds, rng_.normal());
    renormalize();
}

double InflowModel::streamwise(double s, double ubar) const {
    return ubar * shape(s) / area_;
}

double InflowModel::transverse(double s, double ubar) const {
    return meander_.value * streamwise(s, ubar);
}

double InflowModel::discreteFlux(double ubar) const {
    double q = 0.0;
    for (size_t k = 0; k < nodeS_.size(); ++k) q += streamwise(nodeS_[k], ubar) * nodeWeight_[k];
    return q;
}

// ---------------------------------------------------------------- stirring

Stirring::Stirring(const StirringParams& params, std::uint64_t seed, const Domain& domain,
                   const ChannelMap& channel)
    : params_(params), rng_(seed, 2) {
    const int n = domain.fluidCount();
    fx_.assign(n, 0.0);
    fy_.assign(n, 0.0);
    if (!params_.enabled) return;

    // Large-scale (p across, q along) pairs, cycled if more modes are requested.
    static const int pq[][2] = {{1, 2}, {2, 3}, {1, 4}, {2, 5}, {1, 3}, {3, 4}, {1, 5}, {3, 6}};
    for (int k = 0; k < params_.modes; ++k) {
        const int p = pq[k % 8][0], q = pq[k % 8][1];
        const double phase = 2.0 * pi * rng_.uniform();
        auto psi = [&](double X, double Y) {
            double xi, eta;
            channel(X, Y, xi, eta);
            if (params_.periodic)
                return std::sin(2.0*pi*(p*xi+q*eta)+phase);
            if (xi <= 0.0 || xi >= 1.0 || eta <= 0.0 || eta >= 1.0) return 0.0;
            const double m = 4.0 * xi * (1.0 - xi);   // wall fade: zero value and slope
            return m * m * std::sin(pi * p * xi + phase) * std::sin(pi * q * eta);
        };
        std::vector<double> fx(n), fy(n);
        double sumSq = 0.0;
        for (int node = 0; node < n; ++node) {
            const double X = domain.k(node) + 0.5, Y = domain.j(node) + 0.5;
            fx[node] = 0.5 * (psi(X, Y + 1.0) - psi(X, Y - 1.0));
            fy[node] = -0.5 * (psi(X + 1.0, Y) - psi(X - 1.0, Y));
            sumSq += fx[node] * fx[node] + fy[node] * fy[node];
        }
        const double rms = std::sqrt(sumSq / std::max(n, 1));
        if (rms > 0.0)
            for (int node = 0; node < n; ++node) { fx[node] /= rms; fy[node] /= rms; }
        modeFx_.push_back(std::move(fx));
        modeFy_.push_back(std::move(fy));

        OrnsteinUhlenbeck ou;
        ou.sigma = 1.0 / std::sqrt(static_cast<double>(params_.modes));
        ou.correlationTime = params_.correlationTime;
        ou.value = ou.sigma * rng_.normal();
        amplitudes_.push_back(ou);
    }
}

void Stirring::advance(double dtSeconds, double forceScale) {
    if (!params_.enabled) return;
    std::fill(fx_.begin(), fx_.end(), 0.0);
    std::fill(fy_.begin(), fy_.end(), 0.0);
    for (size_t k = 0; k < amplitudes_.size(); ++k) {
        amplitudes_[k].step(dtSeconds, rng_.normal());
        const double b = amplitudes_[k].value * params_.intensity * forceScale;
        const auto& mx = modeFx_[k];
        const auto& my = modeFy_[k];
        for (size_t node = 0; node < fx_.size(); ++node) {
            fx_[node] += b * mx[node];
            fy_[node] += b * my[node];
        }
    }
}

}  // namespace lbm
