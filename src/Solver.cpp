#include "lbm/Solver.hpp"

#include "lbm/Random.hpp"

#include <algorithm>
#include <stdexcept>
#ifndef __EMSCRIPTEN__
#include <condition_variable>
#include <exception>
#include <mutex>
#include <thread>
#endif

namespace lbm {

struct Solver::WorkerPool {
#ifndef __EMSCRIPTEN__
    Solver& owner;
    std::mutex mutex;
    std::condition_variable start, done;
    std::vector<std::thread> workers;
    std::vector<std::exception_ptr> errors;
    unsigned generation = 0;
    int pending = 0;
    bool stopping = false;
    WorkerPool(Solver& s, int count) : owner(s), errors(count) {
        try {
            for (int id = 0; id < count; ++id) workers.emplace_back([this, id, count] {
                unsigned seen = 0;
                std::unique_lock<std::mutex> lock(mutex);
                for (;;) {
                    start.wait(lock, [&] { return stopping || generation != seen; });
                    if (stopping) return;
                    seen = generation;
                    lock.unlock();
                    try {
                        const int n = owner.domain_.fluidCount();
                        owner.stepRange(n * id / count, n * (id + 1) / count, id);
                    } catch (...) { errors[id] = std::current_exception(); }
                    lock.lock();
                    if (--pending == 0) done.notify_one();
                }
            });
        } catch (...) { stop(); throw; }
    }
    void stop() {
        { std::lock_guard<std::mutex> lock(mutex); stopping = true; }
        start.notify_all();
        for (auto& t : workers) if (t.joinable()) t.join();
    }
    ~WorkerPool() { stop(); }
    void run() {
        std::unique_lock<std::mutex> lock(mutex);
        std::fill(errors.begin(), errors.end(), nullptr);
        pending = static_cast<int>(workers.size());
        ++generation;
        start.notify_all();
        done.wait(lock, [&] { return pending == 0; });
        for (const auto& error : errors) if (error) std::rethrow_exception(error);
    }
#endif
};

Solver::~Solver() = default;

void Solver::setThreadCount(int count) {
    if (count < 1) throw std::invalid_argument("Solver: thread count must be positive");
#ifdef __EMSCRIPTEN__
    count = 1;
#endif
    pool_.reset();
    threadCount_ = count;
    partialMaxTau_.assign(count, 0.0);
#ifndef __EMSCRIPTEN__
    if (count > 1) pool_ = std::make_unique<WorkerPool>(*this, count);
#endif
}

Solver::Solver(const Domain& domain, const CollisionConfig& config, double rho0,
               std::uint64_t seed)
    : domain_(domain), config_(config), rho0_(rho0), seed_(seed), rhoOut_(rho0) {
    if (config_.tau <= 0.5) throw std::invalid_argument("Solver: tau must exceed 1/2");
    const int n = domain_.fluidCount();
    f_.assign(static_cast<size_t>(n) * Q, real(0));
    fNext_ = f_;
    rho_.assign(n, rho0_);
    ux_.assign(n, 0.0);
    uy_.assign(n, 0.0);
    rhoPrev_ = rho_;
    uxPrev_ = ux_;
    uyPrev_ = uy_;
    wallVelocity_.assign(domain_.links().size(), {0.0, 0.0});
    exchange_.resize(domain_.links().size());
    for (int node = 0; node < n; ++node) for (int i = 0; i < Q; ++i) {
        const int src = domain_.pull(node, i);
        if (src >= 0) continue;
        const int li = -src - 1, body = domain_.links()[li].body;
        if (body <= 0) continue;
        bodyLinks_.push_back(li);
        if (static_cast<int>(bodyForce_.size()) <= body) bodyForce_.resize(body + 1);
    }
    setThreadCount(1);
    refreshRates();
}

void Solver::refreshRates() {
    relaxationRates(config_, config_.tau, baseRates_);
    if(config_.stochastic && learnedNoise_.empty()) {
        learnedNoise_.resize(domain_.fluidCount());
        for(int n=0;n<domain_.fluidCount();++n)for(int k=0;k<3;++k)
            normalPair(seed_^0x72ac834619fe02d1ULL,~std::uint64_t(0),n,k,learnedNoise_[n][2*k],learnedNoise_[n][2*k+1]);
    }
    if(!config_.stochastic) learnedNoise_.clear();
    noiseDecay_=config_.noiseCorrelationSteps>0?std::exp(-1/config_.noiseCorrelationSteps):0.;
    noiseInnovation_=std::sqrt(1-noiseDecay_*noiseDecay_);
}

void Solver::initializeNode(int n, double rho, double ux, double uy) {
    for (int i = 0; i < Q; ++i)
        f_[n * Q + i] = static_cast<real>(equilibrium(i, rho, ux, uy) - w[i] * rho0_);
    rho_[n] = rhoPrev_[n] = rho;
    ux_[n] = uxPrev_[n] = ux;
    uy_[n] = uyPrev_[n] = uy;
}

void Solver::initializeUniform(double rho, double ux, double uy) {
    for (int n = 0; n < domain_.fluidCount(); ++n) initializeNode(n, rho, ux, uy);
}

void Solver::initializePopulations(int n, const double* populations) {
    double rho=0., jx=0., jy=0.;
    for (int i=0;i<Q;++i) {
        f_[n*Q+i]=static_cast<real>(populations[i]-w[i]*rho0_);
        rho+=populations[i]; jx+=cx[i]*populations[i]; jy+=cy[i]*populations[i];
    }
    const auto F=force(n);
    rho_[n]=rhoPrev_[n]=rho;
    ux_[n]=uxPrev_[n]=(jx+0.5*F[0])/rho;
    uy_[n]=uyPrev_[n]=(jy+0.5*F[1])/rho;
}

std::array<double, 2> Solver::bodyForce(int body) const {
    if (body < 0 || body >= static_cast<int>(bodyForce_.size())) return {0.0, 0.0};
    return bodyForce_[body];
}

// Population entering node n along direction i through a boundary link.
// Reads only post-collision populations and velocities of the previous step.
double Solver::incoming(int n, int i, const BoundaryLink& link, int linkIndex) const {
    const int o = link.dir;   // outgoing direction; i == opposite[o]
    const double fo = f_[n * Q + o];
    switch (link.kind) {
        case LinkKind::Wall: {
            const double q = link.q;
            if (q < 0.5) {
                // Needs the node one step back along the link (upstream for direction o).
                const int back = domain_.pull(n, o);
                if (back < 0) return fo;   // no second fluid node: plain bounce-back
                return 2.0 * q * fo + (1.0 - 2.0 * q) * f_[back * Q + o];
            }
            return fo / (2.0 * q) + (2.0 * q - 1.0) / (2.0 * q) * f_[n * Q + i];
        }
        case LinkKind::MovingWall: {
            const auto& u = wallVelocity_[linkIndex];
            const double cu = cx[o] * u[0] + cy[o] * u[1];
            return fo - 6.0 * w[o] * rhoPrev_[n] * cu;
        }
        case LinkKind::PressureOutlet: {
            double uwx = uxPrev_[n], uwy = uyPrev_[n];
            const int back = domain_.pull(n, o);
            if (back >= 0) {
                uwx += 0.5 * (uxPrev_[n] - uxPrev_[back]);
                uwy += 0.5 * (uyPrev_[n] - uyPrev_[back]);
            }
            if (outletBackflowClamp_) {
                // Outward normal of the box side this link crosses.
                const double nx = link.side == XMin ? -1.0 : link.side == XMax ? 1.0 : 0.0;
                const double ny = link.side == YMin ? -1.0 : link.side == YMax ? 1.0 : 0.0;
                const double un = uwx * nx + uwy * ny;
                if (un < 0.0) {
                    uwx -= un * nx;
                    uwy -= un * ny;
                }
            }
            const double cu = cx[o] * uwx + cy[o] * uwy;
            const double u2 = uwx * uwx + uwy * uwy;
            // Deviation form of -f_o* + 2 w rho_out [1 + 4.5 cu^2 - 1.5 u^2] - w rho0.
            return -fo + 2.0 * w[o] * (rhoOut_ * (1.0 + 4.5 * cu * cu - 1.5 * u2) - rho0_);
        }
    }
    return fo;
}

void Solver::step() {
    std::swap(rho_, rhoPrev_);
    std::swap(ux_, uxPrev_);
    std::swap(uy_, uyPrev_);
#ifndef __EMSCRIPTEN__
    if (pool_) pool_->run();
    else
#endif
        stepRange(0, domain_.fluidCount(), 0);
    maxTauEff_ = *std::max_element(partialMaxTau_.begin(), partialMaxTau_.end());
    for (auto& fb : bodyForce_) fb = {0.0, 0.0};
    for (int li : bodyLinks_) {
        const auto& link = domain_.links()[li];
        bodyForce_[link.body][0] += cx[link.dir] * exchange_[li];
        bodyForce_[link.body][1] += cy[link.dir] * exchange_[li];
    }
    std::swap(f_, fNext_);
    ++time_;
}

void Solver::stepRange(int begin, int end, int worker) {
    const auto& links = domain_.links();
    double maxTau = 0.0;
    const bool thermal = config_.kT > 0.0;
    const bool hasWallDistance = !wallDistance_.empty();

    for (int n = begin; n < end; ++n) {
        double f[Q];
        for (int i = 0; i < Q; ++i) {
            const int src = domain_.pull(n, i);
            if (src >= 0) {
                f[i] = f_[src * Q + i];
                continue;
            }
            const int li = -src - 1;
            const BoundaryLink& link = links[li];
            f[i] = incoming(n, i, link, li);
            if (link.body > 0) {
                // Momentum exchange: the body receives c_o (f_o* + f_i) per link.
                const int o = link.dir;
                const double exchange = f_[n * Q + o] + f[i] + 2.0 * w[o] * rho0_;
                exchange_[li] = exchange;
            }
        }

        double Fx = uniformFx_, Fy = uniformFy_;
        if (fieldFx_) { Fx += (*fieldFx_)[n]; Fy += (*fieldFy_)[n]; }
        if (dragField_ && (*dragField_)[n] > 0.0) {
            // Implicit drag: with u = (j + F/2) / (rho (1 + k/2)), the drag -k rho u is
            // -c (j + F/2), c = k / (1 + k/2), and the collision's u includes it.
            const double k = (*dragField_)[n];
            double jx = 0.0, jy = 0.0;
            for (int i = 0; i < Q; ++i) { jx += cx[i] * f[i]; jy += cy[i] * f[i]; }
            const double c = k / (1.0 + 0.5 * k);
            Fx -= c * (jx + 0.5 * Fx);
            Fy -= c * (jy + 0.5 * Fy);
        }

        double noise[6];
        if (thermal) {
            for (int p = 0; p < 3; ++p)
                normalPair(seed_, static_cast<std::uint64_t>(time_),
                           static_cast<std::uint64_t>(n), p, noise[2 * p], noise[2 * p + 1]);
        }

        if(config_.stochastic) {
            for(int k=0;k<3;++k) {
                double a,b;normalPair(seed_^0x72ac834619fe02d1ULL,time_,n,k,a,b);
                learnedNoise_[n][2*k]=noiseDecay_*learnedNoise_[n][2*k]+noiseInnovation_*a;
                learnedNoise_[n][2*k+1]=noiseDecay_*learnedNoise_[n][2*k+1]+noiseInnovation_*b;
            }
        }
        const CollisionOutput out =
            collide(f, rho0_, Fx, Fy, config_, baseRates_,
                    hasWallDistance ? wallDistance_[n] : -1.0, thermal ? noise : nullptr,
                    viscousFlux_.empty() ? nullptr : viscousFlux_[n].data(),
                    config_.stochastic ? learnedNoise_[n].data() : nullptr);
        for (int i = 0; i < Q; ++i) fNext_[n * Q + i] = static_cast<real>(f[i]);
        rho_[n] = out.rho;
        ux_[n] = out.ux;
        uy_[n] = out.uy;
        maxTau = std::max(maxTau, out.tauEff);
    }
    partialMaxTau_[worker] = maxTau;
}

}  // namespace lbm
