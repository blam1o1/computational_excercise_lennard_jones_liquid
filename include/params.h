// header file for params.cpp
#pragma once
#include <string>


// Params struct to store all simulation parameters, defined in input.txt 
struct Params {
    
    // System parameters
    int N = 500;
    double L = 8.5;
    double mass = 39.98; // in atomic mass units
    double density = 0.02 // angstroms^-3


    // Potential 
    double rc = 2.5;

    // Ensemble 
    double T = 1.0;

    // Integration parameters
    double time_step = 0.0005;
    int n_steps = 10000;

    // Thermostat
    std::string thermostat = "nvh";

    // Seed
    unsigned seed = 42;

    // Initialization and Output parameters
    int thermo_every = 100;
    int traj_every = 1000;

};

Params read_params(const std::string& filename);

void check_params(const Params& params);
void print_params(const Params& params, std::ostream& os);

