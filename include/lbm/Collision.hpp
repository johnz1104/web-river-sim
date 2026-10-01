#pragma once

// Collision operators. All three models share one moment-space kernel:
//   m_k* = m_k - s_k (m_k - m_k^eq) + (1 - s_k / 2) F_k + noise_k
// with the rates s_k chosen per model:
//   BGK  every non-conserved moment relaxes at 1/tau
//   TRT  even moments at s+ = 1/tau, odd (ghost) moments at s- from the magic
//        parameter Lambda = (1/s+ - 1/2)(1/s- - 1/2)
//   MRT  shear at 1/tau, odd ghosts from Lambda, bulk and even ghost set freely
// The shear rate sets the kinematic viscosity nu = c_s^2 (tau - 1/2).

#include "lbm/Lattice.hpp"

#include <algorithm>
#include <cmath>
#include <memory>

namespace lbm {

enum class CollisionModel { BGK, TRT, MRT };

// Sub-grid closure: how the local shear relaxation time is chosen.
enum class ClosureKind { Laminar, Smagorinsky, Custom };

// Local, rotation-invariant information a closure may use. The stress is the
// non-equilibrium momentum flux Pi^neq (forcing-corrected), in lattice units.
struct ClosureInputs {
    double rho, ux, uy;
    double Pxx, Pyy, Pxy;
    double wallDistance;   // lattice units; negative when not supplied
};

// Hook for data-driven closures (the planned ML models): return tau_eff >= tau0.
class ClosureModel {
public:
    virtual ~ClosureModel() = default;
    virtual double tauEffective(const ClosureInputs& in, double tau0) const = 0;
};

// Generalized non-conserved collision hook. Conserved modes are inaccessible.
struct CollisionState {
    const double* moments;
    const double* deviation;
    double rho, ux, uy, Pxx, Pyy, Pxy, wallDistance, tau;
};
class CollisionClosure {
public:
    virtual ~CollisionClosure() = default;
    virtual void apply(const CollisionState&, double rates[6], double correction[6]) const = 0;
    virtual void noiseAmplitudes(const CollisionState&, double sigma[6]) const { std::fill(sigma,sigma+6,0.); }
};

struct CollisionConfig {
    CollisionModel model = CollisionModel::MRT;
    double tau = 0.6;              // shear relaxation time (sets nu)
    double magic = 3.0 / 16.0;     // Lambda for the odd moments (TRT, MRT)
    double bulkRate = 1.2;         // MRT: rate of k = 3 (bulk viscosity, damps sound)
    double evenGhostRate = 1.2;    // MRT: rate of k = 8
    double oddGhostRate = 0.0;     // MRT: fixed rate for k = 6, 7; 0 = derive from magic
    ClosureKind closure = ClosureKind::Laminar;
    double smagorinsky = 0.15;     // C_S, with filter width = one lattice spacing
    std::shared_ptr<const ClosureModel> custom;
    std::shared_ptr<const CollisionClosure> learned;
    std::shared_ptr<const CollisionClosure> stochastic;
    double noiseCorrelationSteps = 0.; // zero: white; positive: stationary OU
    double kT = 0.0;               // thermal-noise temperature, lattice units (0 = off)
};

inline double viscosityFromTau(double tau) { return cs2 * (tau - 0.5); }
inline double tauFromViscosity(double nu) { return nu / cs2 + 0.5; }

// Relaxation rates for all nine moments at a given shear relaxation time.
inline void relaxationRates(const CollisionConfig& c, double tau, double s[Q]) {
    const double sPlus = 1.0 / tau;
    const double sMinus = 1.0 / (c.magic / (tau - 0.5) + 0.5);
    s[0] = s[1] = s[2] = 0.0;   // conserved: density and momentum
    switch (c.model) {
        case CollisionModel::BGK:
            for (int k = 3; k < Q; ++k) s[k] = sPlus;
            break;
        case CollisionModel::TRT:
            s[3] = s[4] = s[5] = s[8] = sPlus;
            s[6] = s[7] = sMinus;
            break;
        case CollisionModel::MRT:
            s[3] = c.bulkRate;
            s[4] = s[5] = sPlus;
            s[6] = s[7] = c.oddGhostRate > 0.0 ? c.oddGhostRate : sMinus;
            s[8] = c.evenGhostRate;
            break;
    }
}

// Smagorinsky in closed form (Hou et al. 1996). With S = -Pi^neq / (2 rho c_s^2 tau_eff)
// and nu_eff = nu_0 + (C_S Delta)^2 |S|, |S| = sqrt(2 S:S):
//   tau_eff = 1/2 [tau_0 + sqrt(tau_0^2 + 18 sqrt(2) C_S^2 |Pi^neq| / rho)]
inline double smagorinskyTau(double tau0, double cSmag, double rho, double Pxx,
                             double Pyy, double Pxy) {
    const double piNorm = std::sqrt(Pxx * Pxx + Pyy * Pyy + 2.0 * Pxy * Pxy);
    return 0.5 * (tau0 + std::sqrt(tau0 * tau0
                                   + 18.0 * std::sqrt(2.0) * cSmag * cSmag * piNorm / rho));
}

struct CollisionOutput {
    double rho, ux, uy;   // physical (force-corrected) velocity
    double tauEff;
};

// One collision at one node. f holds deviations f_i - w_i rho0 on entry (pre-collision)
// and on exit (post-collision). F is the body force per unit volume (Guo forcing).
// noise, when non-null, holds six standard normals for the thermal modes k = 3..8.
inline CollisionOutput collide(double f[Q], double rho0, double Fx, double Fy,
                               const CollisionConfig& cfg, const double baseRates[Q],
                               double wallDistance, const double* noise,
                               double* viscousFlux = nullptr, const double* learnedNoise = nullptr) {
    double m[Q];
    toMoments(f, m);   // deviations shift only m_0: m_0 = rho - rho0

    const double rho = rho0 + m[0];
    const double ux = (m[1] + 0.5 * Fx) / rho;
    const double uy = (m[2] + 0.5 * Fy) / rho;
    const EquilibriumMoments eq = equilibriumMoments(rho, ux, uy);

    double s[Q];
    double tauEff = cfg.tau;
    if (cfg.closure == ClosureKind::Laminar) {
        std::copy(baseRates, baseRates + Q, s);
    } else {
        // Non-equilibrium stress from the pre-collision moments, with the Guo forcing
        // correction 1/2 (F_a u_b + u_a F_b).
        const double dBulk = (m[3] - eq.bulk) / 3.0;
        const double dNormal = m[4] - eq.normalDiff;
        const double Pxx = 0.5 * (dBulk + dNormal) + ux * Fx;
        const double Pyy = 0.5 * (dBulk - dNormal) + uy * Fy;
        const double Pxy = (m[5] - eq.shear) + 0.5 * (Fx * uy + ux * Fy);
        if (cfg.closure == ClosureKind::Smagorinsky) {
            tauEff = smagorinskyTau(cfg.tau, cfg.smagorinsky, rho, Pxx, Pyy, Pxy);
        } else if (cfg.custom) {
            tauEff = cfg.custom->tauEffective({rho, ux, uy, Pxx, Pyy, Pxy, wallDistance},
                                              cfg.tau);
        }
        tauEff = std::max(tauEff, 0.5 + 1e-9);
        // Only the rates that depend on tau change; the rest keep their base values.
        std::copy(baseRates, baseRates + Q, s);
        const double sPlus = 1.0 / tauEff;
        if (cfg.model == CollisionModel::BGK) {
            for (int k = 3; k < Q; ++k) s[k] = sPlus;
        } else {
            s[4] = s[5] = sPlus;
            if (cfg.model == CollisionModel::TRT) s[3] = s[8] = sPlus;
            if (cfg.model == CollisionModel::TRT || cfg.oddGhostRate <= 0.0)
                s[6] = s[7] = 1.0 / (cfg.magic / (tauEff - 0.5) + 0.5);
        }
    }

    double correction[6]={};
    if(cfg.learned || (cfg.stochastic && learnedNoise)) {
        double deviation[Q];std::copy(m,m+Q,deviation);
        deviation[3]-=eq.bulk;deviation[4]-=eq.normalDiff;deviation[5]-=eq.shear;
        const double trace=deviation[3]/3.,diff=deviation[4];
        const CollisionState state{m,deviation,rho,ux,uy,.5*(trace+diff)+ux*Fx,
            .5*(trace-diff)+uy*Fy,deviation[5]+.5*(ux*Fy+uy*Fx),wallDistance,cfg.tau};
        if(cfg.learned) cfg.learned->apply(state,s+3,correction);
        if(cfg.stochastic && learnedNoise) {
            double sigma[6];cfg.stochastic->noiseAmplitudes(state,sigma);
            for(int k=0;k<6;++k) correction[k]+=sigma[k]*learnedNoise[k];
        }
        tauEff=1/s[4];
    }

    // Physical viscous momentum flux, pressure excluded. Capture pre-collision
    // moments directly; recovering them from post-collision data fails at s=1.
    if (viscousFlux) {
        const double trace = (1.0 - 0.5*s[3]) *
            ((m[3]-eq.bulk)/3.0 + ux*Fx + uy*Fy);
        const double diff = (1.0 - 0.5*s[4]) *
            (m[4]-eq.normalDiff + ux*Fx - uy*Fy);
        viscousFlux[0] = 0.5*(trace+diff);
        viscousFlux[1] = 0.5*(trace-diff);
        viscousFlux[2] = (1.0 - 0.5*s[5]) *
            (m[5]-eq.shear + 0.5*(ux*Fy+uy*Fx));
    }

    // Guo forcing in moment space: F_hat = (0, Fx, Fy, 6 u.F, 2(ux Fx - uy Fy),
    // ux Fy + uy Fx, 0, 0, 0). Momentum gains the full force (s = 0).
    const double F3 = 6.0 * (ux * Fx + uy * Fy);
    const double F4 = 2.0 * (ux * Fx - uy * Fy);
    const double F5 = ux * Fy + uy * Fx;

    m[1] += Fx;
    m[2] += Fy;
    m[3] += -s[3] * (m[3] - eq.bulk) + (1.0 - 0.5 * s[3]) * F3;
    m[4] += -s[4] * (m[4] - eq.normalDiff) + (1.0 - 0.5 * s[4]) * F4;
    m[5] += -s[5] * (m[5] - eq.shear) + (1.0 - 0.5 * s[5]) * F5;
    m[6] -= s[6] * m[6];
    m[7] -= s[7] * m[7];
    m[8] -= s[8] * m[8];

    if(cfg.learned || (cfg.stochastic && learnedNoise)) {
        for(int k=3;k<Q;++k) m[k]+=correction[k-3];
        if(viscousFlux) {
            viscousFlux[0]+=.25*(correction[0]/3+correction[1]);
            viscousFlux[1]+=.25*(correction[0]/3-correction[1]);
            viscousFlux[2]+=.5*correction[2];
        }
    }

    // Fluctuating LBM (Duenweg, Schiller & Ladd 2007): each non-conserved mode gets
    // noise of variance mu * norm_k * s_k (2 - s_k), with mu = rho kT / c_s^2, which
    // yields equipartition <rho u_a^2> = kT.
    if (noise && cfg.kT > 0.0) {
        const double mu = rho * cfg.kT / cs2;
        for (int k = 3; k < Q; ++k)
            m[k] += std::sqrt(mu * norm[k] * s[k] * (2.0 - s[k])) * noise[k - 3];
    }

    fromMoments(m, f);
    return {rho, ux, uy, tauEff};
}

// Reference BGK written directly on populations (for learning and for checking that
// the moment-space kernel reduces to it exactly). f holds full populations.
inline CollisionOutput collideBGKPopulations(double f[Q], double tau, double Fx,
                                             double Fy) {
    double rho = 0.0, jx = 0.0, jy = 0.0;
    for (int i = 0; i < Q; ++i) { rho += f[i]; jx += f[i] * cx[i]; jy += f[i] * cy[i]; }
    const double ux = (jx + 0.5 * Fx) / rho;
    const double uy = (jy + 0.5 * Fy) / rho;
    for (int i = 0; i < Q; ++i) {
        const double cu = cx[i] * ux + cy[i] * uy;
        const double forcing = w[i] * (3.0 * ((cx[i] - ux) * Fx + (cy[i] - uy) * Fy)
                                       + 9.0 * cu * (cx[i] * Fx + cy[i] * Fy));
        f[i] += -(f[i] - equilibrium(i, rho, ux, uy)) / tau + (1.0 - 0.5 / tau) * forcing;
    }
    return {rho, ux, uy, tau};
}

}  // namespace lbm
