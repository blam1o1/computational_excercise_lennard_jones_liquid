#include "../include/thermostat.h"
#include "../include/units.h"
#include <cmath>
#include <stdexcept>


std::vector<double> suzuki_yoshida_weights(int nsy) {
    if (nsy == 1) return {1.0};
    if (nsy == 3) {
        const double w1 = 1.0 / (2.0 - std::cbrt(2.0));
        return {w1, 1.0 - 2.0 * w1, w1};
    }
    if (nsy == 7) {   // Yoshida (1990), sixth order
        const double w1 = 0.784513610477560;
        const double w2 = 0.235573213359357;
        const double w3 = -1.17767998417887;
        const double w4 = 1.0 - 2.0 * (w1 + w2 + w3);
        return {w1, w2, w3, w4, w3, w2, w1};
    }
    throw std::runtime_error("suzuki_yoshida_weights: nsy must be 1, 3 or 7");
}


NoseHooverChain::NoseHooverChain(int M_, double T_target, int n_dof_, double tau, int nc_, int nsy)
    : M(M_), n_dof(n_dof_), nc(nc_), kT(units::KB * T_target), weights(suzuki_yoshida_weights(nsy)) {
    if (M < 1) throw std::runtime_error("NoseHooverChain: M must be >= 1");
    if (n_dof < 1) throw std::runtime_error("NoseHooverChain: n_dof must be >= 1");
    if (nc < 1) throw std::runtime_error("NoseHooverChain: nc must be >= 1");
    if (T_target <= 0 || tau <= 0) throw std::runtime_error("NoseHooverChain: T and tau must be positive");

    Q.assign(M, kT * tau * tau);   // Q_j = kB T tau^2
    Q[0] *= n_dof;                 // Q_1 = N_f kB T tau^2
    eta.assign(M, 0.0);
    p_eta.assign(M, 0.0);
}

NoseHooverChain::NoseHooverChain(const Params& p, int n_dof_)
    : NoseHooverChain(p.nhc_M, p.T, n_dof_, nhc_tau(p), p.nhc_nc, p.nhc_nsy) {}


// G_1 = 2K - N_f kT ;  G_j = p_eta_{j-1}^2 / Q_{j-1} - kT   (0-based: j = 0 is the first link)
double NoseHooverChain::G(int j, double twoK) const {
    if (j == 0) return twoK - n_dof * kT;
    return p_eta[j - 1] * p_eta[j - 1] / Q[j - 1] - kT;
}

// p_eta_j <- p_eta_j * exp(-(d/4) p_eta_{j+1}/Q_{j+1})
// p_eta_j <- p_eta_j + (d/2) G_j
// p_eta_j <- p_eta_j * exp(-(d/4) p_eta_{j+1}/Q_{j+1})
void NoseHooverChain::friction_kick(int j, double twoK, double d) {
    const double fac = std::exp(-0.25 * d * p_eta[j + 1] / Q[j + 1]);
    p_eta[j] *= fac;
    p_eta[j] += 0.5 * d * G(j, twoK);
    p_eta[j] *= fac;
}


void NoseHooverChain::half_step(System& sys, double dt) {
    // 2K = sum_i m v_i^2 ; tracked through the scalings instead of recomputed
    double twoK = 0.0;
    for (const Vec3& v : sys.velocities) twoK += sys.mass * magnitude_squared(v);

    double scale = 1.0;   // accumulated velocity scaling, applied once at the end
    const int top = M - 1;

    for (int k = 0; k < nc; ++k) {
        for (double w : weights) {
            const double d = w * dt / (2.0 * nc);

            // (a) top of the chain
            p_eta[top] += 0.5 * d * G(top, twoK);
            // (b) down the chain: j = M-1 ... 1 (0-based M-2 ... 0)
            for (int j = M - 2; j >= 0; --j) friction_kick(j, twoK, d);
            // (c) thermostat positions
            for (int j = 0; j < M; ++j) eta[j] += d * p_eta[j] / Q[j];
            // (d) particle velocity scaling
            const double a = std::exp(-d * p_eta[0] / Q[0]);
            scale *= a;
            twoK *= a * a;
            // (e) up the chain: j = 1 ... M-1 (0-based 0 ... M-2)
            for (int j = 0; j <= M - 2; ++j) friction_kick(j, twoK, d);
            // (f) top of the chain
            p_eta[top] += 0.5 * d * G(top, twoK);
        }
    }

    for (Vec3& v : sys.velocities) v *= scale;
}


double NoseHooverChain::energy() const {
    double E = 0.0;
    for (int j = 0; j < M; ++j) E += 0.5 * p_eta[j] * p_eta[j] / Q[j];
    E += n_dof * kT * eta[0];
    for (int j = 1; j < M; ++j) E += kT * eta[j];
    return E;
}