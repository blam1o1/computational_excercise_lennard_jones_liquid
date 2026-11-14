// header file for thermostat.cpp

// Nose-Hoover chain thermostat (Martyna, Klein, Tuckerman 1992) and its half-step propagator

#pragma once
#include <vector>
#include "params.h"
#include "system.h"

// Equations of motion (p_i = m v_i):
//   dp_i/dt     = F_i - (p_eta_1 / Q_1) p_i
//   deta_j/dt   = p_eta_j / Q_j                                  j = 1..M
//   dp_eta_j/dt = G_j - (p_eta_{j+1} / Q_{j+1}) p_eta_j           j = 1..M-1
//   dp_eta_M/dt = G_M
//   G_1 = sum_i p_i^2 / m - N_f kB T,   G_j = p_eta_{j-1}^2 / Q_{j-1} - kB T  (j >= 2)
// Masses: Q_1 = N_f kB T tau^2, Q_j = kB T tau^2.
// Thermostat contribution to the conserved quantity:
//   sum_j p_eta_j^2 / (2 Q_j) + N_f kB T eta_1 + kB T sum_{j>=2} eta_j
struct NoseHooverChain {
    int M = 0;                      // chain length
    int n_dof = 0;                  // N_f = 3N - 3
    int nc = 1;                     // sub-steps per half step
    double kT = 0.0;                // kB * T_target, amu A^2 / ps^2
    std::vector<double> weights;    // Suzuki-Yoshida weights w_nu
    std::vector<double> Q;          // thermostat masses, amu A^2
    std::vector<double> eta;        // thermostat positions (dimensionless), start at 0
    std::vector<double> p_eta;      // thermostat momenta, amu A^2 / ps, start at 0

    NoseHooverChain(int M, double T_target, int n_dof, double tau, int nc, int nsy);
    NoseHooverChain(const Params& p, int n_dof);

    // exp(iL_NHC dt/2): updates eta, p_eta and scales sys.velocities
    void half_step(System& sys, double dt);

    // thermostat part of the conserved quantity H'
    double energy() const;

  private:
    double G(int j, double twoK) const;                 // force on p_eta_j (0-based j)
    void friction_kick(int j, double twoK, double d);   // the 3-line update of steps (b), (e)
};

// Suzuki-Yoshida weights for nsy = 1 (2nd order), 3 (4th order), 7 (6th order)
std::vector<double> suzuki_yoshida_weights(int nsy);