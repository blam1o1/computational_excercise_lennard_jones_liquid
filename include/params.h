// header file for params.cpp
#pragma once
#include <ostream>
#include <string>

// All simulation inputs, read from input.txt. Units: A, amu, K, ps.
// Quantities that follow from these (box length, cutoff in A, ...) are computed by the
// functions below the struct, so they always agree with whatever was read from the file.
struct Params {

    // System
    int N = 864;                    // number of particles, must be 4 n^3 (FCC supercell)
    double density = 0.02;          // number density, A^-3
    double mass = 39.98;            // amu

    // Lennard-Jones potential
    double sigma = 3.405;           // A
    double epsilon_K = 119.8;       // epsilon / kB, K
    double rc_sigma = 2.5;          // cutoff radius, units of sigma
    double lambda_sigma = 0.1;      // healing length of the switching function, units of sigma

    // Ensemble / thermostat
    double T = 100.0;               // target temperature, K
    std::string thermostat = "nhc"; // "nhc" (Nose-Hoover chain) or "none" (plain velocity Verlet)
    int nhc_M = 5;                  // chain length
    double nhc_tau_steps = 100.0;   // thermostat time scale tau = nhc_tau_steps * dt (>= 20)
    int nhc_nc = 1;                 // NHC sub-steps per half step
    int nhc_nsy = 7;                // Suzuki-Yoshida stages: 1, 3 or 7 (7 = sixth order)

    // Integration
    double dt = 0.005;              // time step, ps
    long n_equil = 10000;           // equilibration steps
    long n_prod = 50000;            // production steps

    // Initialization and output
    unsigned seed = 42;
    int thermo_every = 100;         // thermo_equil.dat / thermo_prod.dat line interval, steps
    int traj_every = 1000;          // trajectory.xyz frame interval, steps (whole run)
    int sample_every = 10;          // g(r) sampling interval, steps (production only)
    std::string output_prefix;      // prepended to every output file name (default: none)

    // Observables
    double gr_dr = 0.05;            // g(r) bin width, A
    double sq_qmin = 0.5;           // S(q) grid start, A^-1
    double sq_qmax = 15.0;          // S(q) grid end, A^-1
    double sq_dq = 0.01;            // S(q) grid spacing, A^-1
};

Params read_params(const std::string& filename);
void check_params(const Params& p);
void print_params(const Params& p, std::ostream& os);

// Derived quantities, in internal units
double box_length(const Params& p);        // L = (N / density)^(1/3), A
double cutoff_radius(const Params& p);     // rc = rc_sigma * sigma, A
double healing_length(const Params& p);    // lambda = lambda_sigma * sigma, A
double epsilon_internal(const Params& p);  // epsilon = kB * epsilon_K, amu A^2 / ps^2
double nhc_tau(const Params& p);           // tau = nhc_tau_steps * dt, ps