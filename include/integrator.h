// header file for integrator.cpp

// one time step of the factorized propagator
//   exp(iL dt) ~ exp(iL_NHC dt/2) exp(iL_2 dt/2) exp(iL_1 dt) exp(iL_2 dt/2) exp(iL_NHC dt/2)

#pragma once
#include "force.h"
#include "system.h"
#include "thermostat.h"

// NVT step (Nose-Hoover chain). Requires sys.forces to hold F(t) on entry
// (call compute_forces once before the first step); leaves F(t + dt) on exit.
//   1. NHC half step          4. F, U, W <- compute_forces
//   2. v += (dt/2) F / m      5. v += (dt/2) F / m
//   3. r += dt v ; wrap       6. NHC half step
void integrate_step(System& sys, const LennardJones& lj, NoseHooverChain& nhc, double dt);

// NVE step (plain velocity Verlet, steps 2-5 only), for testing
void integrate_step(System& sys, const LennardJones& lj, double dt);

// conserved quantities: H' = K + U + thermostat energy (NVT), H = K + U (NVE)
double conserved_energy(const System& sys, const NoseHooverChain& nhc);
double conserved_energy(const System& sys);