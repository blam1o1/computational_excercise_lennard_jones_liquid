// header file for force.cpp

// Lennard-Jones pair potential with a switching function, and the force loop

#pragma once
#include "params.h"
#include "system.h"

// u(r)  = 4 epsilon [ (sigma/r)^12 - (sigma/r)^6 ]
// S(r)  = 1                    r < rc - lambda
//       = 1 + R^2 (2R - 3)     rc - lambda <= r <= rc,   R = (r - (rc - lambda)) / lambda
//       = 0                    r > rc
// pair energy  u(r) S(r);  pair force on i from j  F_ij = -[u'(r) S(r) + u(r) S'(r)] r_ij / r
struct LennardJones {
    double sigma = 0.0;      // A
    double epsilon = 0.0;    // amu A^2 / ps^2  (= kB * epsilon_K)
    double rc = 0.0;         // cutoff radius, A
    double lambda = 0.0;     // healing length, A

    // derived
    double sigma2 = 0.0;     // sigma^2
    double rc2 = 0.0;        // rc^2
    double r_switch = 0.0;   // rc - lambda, where switching starts
    double r_switch2 = 0.0;  // (rc - lambda)^2

    LennardJones(double sigma, double epsilon, double rc, double lambda);
    explicit LennardJones(const Params& p);

    double switching(double r) const;             // S(r)
    double switching_derivative(double r) const;  // dS/dr
    double pair_energy(double r) const;           // u(r) S(r)
    double pair_dudr(double r) const;             // d(uS)/dr = u'S + uS'
};

// Fills sys.forces, sys.potential_energy and sys.virial at the current positions.
// All pairs i < j, minimum-image separations; requires rc < L/2.
void compute_forces(System& sys, const LennardJones& lj);