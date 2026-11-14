// header file for observables.cpp

// quantities measured from the state of the system, and their running averages.
// Calculations only: writing them to files lives in output.h.

#pragma once
#include <cstddef>
#include <vector>
#include "vec3.h"

struct System;   // system.h (declared here so including this header stays lightweight)


// ---------------------------------------------------------------- instantaneous quantities

// Number of degrees of freedom once total momentum has been removed: 3N - 3.
// The Nose-Hoover chain must use this same value for its target kinetic energy.
inline int degrees_of_freedom(std::size_t N) {
    return 3 * static_cast<int>(N) - 3;
}

// V_cm = sum_i(m_i * v_i) / sum_i(m_i), all particles sharing the same mass
Vec3 center_of_mass_velocity(const std::vector<Vec3>& velocities, double mass);

// K = 1/2 * sum_i(m_i * |v_i|^2)    [amu A^2 / ps^2]
double kinetic_energy(const std::vector<Vec3>& velocities, double mass);

// T_kin = 2K / ((3N - 3) * kB)    [K]
double temperature(const std::vector<Vec3>& velocities, double mass);


// ---------------------------------------------------------------- running averages

// Energy conservation measure, accumulated every step:
//   Delta E(n) = (1/n) sum_{k=1}^{n} |(E_k - E0) / E0|
struct EnergyDrift {
    double E0 = 0.0;          // reference conserved energy
    long n = 0;               // steps accumulated
    double sum_abs = 0.0;     // sum_k |(E_k - E0) / E0|

    void start(double E_initial);        // set E0 and clear the sum
    void add(double E);                  // add one step
    double deviation(double E) const;    // (E - E0) / |E0|, signed, instantaneous
    double delta_E() const;              // Delta E(n); NaN before the first step
};

// Running mean and standard deviation of a sampled quantity (used for T_kin):
//   <x>_n = (1/n) sum_{j=1}^{n} x_j      (Welford's update, numerically stable)
struct RunningAverage {
    long n = 0;
    double mean_ = 0.0;
    double m2 = 0.0;          // sum of squared deviations from the running mean

    void reset();
    void add(double x);
    double mean() const;      // NaN before the first sample
    double std_dev() const;   // population standard deviation; NaN before the first sample
};


// ---------------------------------------------------------------- structure

// Radial distribution function from a histogram of minimum-image pair distances.
//   bins of width dr on [0, r_max = L/2), bin i centred at r_i = (i + 1/2) dr
//   each sampled frame: for every pair i < j, h[int(r/dr)] += 2
//   g(r_i) = h(i) / (4 pi rho r_i^2 dr N_conf N),   rho = N / L^3,   i < n_bins (r < L/2)
// The histogram extends to the box corner (sqrt(3) L/2) so every pair can be binned without an
// r < L/2 test (that branch is unpredictable and made sampling ~1.7x slower); the bins beyond
// n_bins hold incomplete shells and are never used for g(r).
struct RadialDistribution {
    double dr = 0.0;                    // bin width, A
    double L = 0.0;                     // box length the histogram belongs to, A
    double r_max = 0.0;                 // L/2
    int n_bins = 0;                     // bins used for g(r): floor(r_max / dr)
    long n_conf = 0;                    // frames sampled
    std::vector<long long> histogram;   // h(i), summed over frames (size > n_bins, see above)

    RadialDistribution(double dr, double L);

    void sample(const System& sys);                 // add one frame
    double r(int i) const { return (i + 0.5) * dr; }
    std::vector<double> radii() const;              // all bin centres r_i
    std::vector<double> g(int N) const;             // normalized g(r_i); needs n_conf > 0
};

// Static structure factor from g(r), midpoint rule on the g(r) bins, integrated to r_max:
//   S(q) = 1 + 4 pi rho dr sum_i r_i^2 [g(r_i) - 1] sin(q r_i) / (q r_i)
std::vector<double> structure_factor(const std::vector<double>& r, const std::vector<double>& g,
                                     double dr, double rho, const std::vector<double>& q);

// q_min, q_min + dq, ... up to q_max (inclusive within round-off)
std::vector<double> q_grid(double q_min, double q_max, double dq);