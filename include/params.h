
#pragma once

#include <string>
#include <ostream>
#include <cmath>

struct Params {

    // System parameters
    int N = 864;
    double mass = 39.98;        // atomic mass units
    double density = 0.02;      // Angstrom^-3

    // Cubic simulation box (Angstrom)
    double L = std::cbrt(N / density);

    // Lennard-Jones potential
    double sigma = 3.405;       // Angstrom
    double epsilon_K = 119.8;   // epsilon / kB, Kelvin
    double rc_sigma = 2.5;      // cutoff in units of sigma

    // Ensemble
    double T = 100.0;           // Kelvin

    // Integration parameters
    double time_step = 0.0005;  // Units must be defined
    int n_steps = 10000;

    // Thermostat
    std::string thermostat = "nhc";

    // Random seed
    unsigned seed = 42;

    // Output parameters
    int thermo_every = 100;
    int traj_every = 1000;
};

Params read_params(const std::string& filename);

void check_params(const Params& params);

void print_params(const Params& params, std::ostream& os);


