#pragma once

// Counter-based random numbers. A value depends only on (seed, key...), never on
// call order, so runs are reproducible with a fixed seed and loops can be
// reordered or parallelised without changing results.

#include <cmath>
#include <cstdint>

namespace lbm {

// SplitMix64 finalizer: a fast, well-mixed 64-bit hash.
inline std::uint64_t mix64(std::uint64_t z) {
    z += 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

inline std::uint64_t hashKey(std::uint64_t seed, std::uint64_t a, std::uint64_t b,
                             std::uint64_t c) {
    return mix64(seed ^ mix64(a ^ mix64(b ^ mix64(c))));
}

// Uniform in (0, 1): never exactly 0, so log() in Box-Muller is safe.
inline double uniformFromBits(std::uint64_t bits) {
    return (static_cast<double>(bits >> 11) + 0.5) * 0x1.0p-53;
}

// Two independent standard normals from one key (Box-Muller).
inline void normalPair(std::uint64_t seed, std::uint64_t a, std::uint64_t b,
                       std::uint64_t c, double& n0, double& n1) {
    const double u1 = uniformFromBits(hashKey(seed, a, b, 2 * c));
    const double u2 = uniformFromBits(hashKey(seed, a, b, 2 * c + 1));
    const double r = std::sqrt(-2.0 * std::log(u1));
    const double theta = 6.283185307179586 * u2;
    n0 = r * std::cos(theta);
    n1 = r * std::sin(theta);
}

// Sequential stream on top of the counter hash, for the few draws per step made by
// the inflow, stirring and tracer modules. Each module owns its own stream id.
class RandomStream {
public:
    RandomStream(std::uint64_t seed = 0, std::uint64_t stream = 0)
        : seed_(seed), stream_(stream) {}

    double uniform() { return uniformFromBits(hashKey(seed_, stream_, counter_++, 0)); }

    double normal() {
        if (hasSpare_) { hasSpare_ = false; return spare_; }
        double n0, n1;
        normalPair(seed_, stream_, counter_++, 1, n0, n1);
        spare_ = n1;
        hasSpare_ = true;
        return n0;
    }

private:
    std::uint64_t seed_, stream_;
    std::uint64_t counter_ = 0;
    double spare_ = 0.0;
    bool hasSpare_ = false;
};

// Ornstein-Uhlenbeck process: zero mean, stationary standard deviation sigma,
// correlation time T. Uses the exact discrete update, valid for any dt:
//   a <- a e^{-dt/T} + sigma sqrt(1 - e^{-2 dt/T}) xi
struct OrnsteinUhlenbeck {
    double sigma = 0.0;
    double correlationTime = 1.0;
    double value = 0.0;

    void step(double dt, double standardNormal) {
        const double decay = std::exp(-dt / correlationTime);
        value = value * decay + sigma * std::sqrt(1.0 - decay * decay) * standardNormal;
    }
};

}  // namespace lbm
