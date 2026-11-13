#pragma once
// header file for init.cpp

// initial conditions: FCC lattice positions and Maxwell-Boltzmann velocities

#pragma once
#include <cstddef>
#include <random>
#include <utility>
#include <vector>
#include "vec3.h"


// ---------------------------------------------------------------- positions

// Box length of a cube holding N particles at number density rho:  L = (N / rho)^(1/3)
double box_length_from_density(int N, double density);

// N particles on a periodic FCC lattice filling a cube of side L, all coordinates in [0, L).
// Requirements (throws std::runtime_error otherwise):
//   - N = 4 * n^3 (a perfect n x n x n supercell of the 4-atom cubic FCC cell)
//   - N / L^3 equals `density` to within relative tolerance `rel_tol`
std::vector<Vec3> init_fcc_positions(int N, double L, double density, double rel_tol = 1e-4);


// ---------------------------------------------------------------- velocities

// One Box-Muller step: two independent standard normal samples (mean 0, sigma 1)
std::pair<double, double> box_muller(std::mt19937_64& rng);

// Subtract the center-of-mass velocity from every particle (removes drift)
void remove_com_velocity(std::vector<Vec3>& velocities, double mass);

// Multiply every velocity by lambda = sqrt(T_target / T_current)
void rescale_to_temperature(std::vector<Vec3>& velocities, double mass, double T_target);

// Maxwell-Boltzmann velocities for the particles at `positions` (only its size is used,
// so the two arrays are guaranteed to match). Steps:
//   1. Box-Muller samples with sigma = 1
//   2. scale by sqrt(kB T / m)
//   3. remove center-of-mass velocity
//   4. rescale so the instantaneous temperature (3N - 3 dof) equals T_target exactly
// mass in amu, T_target in K, returned velocities in A/ps.
std::vector<Vec3> init_velocities(const std::vector<Vec3>& positions, double mass,
                                  double T_target, std::mt19937_64& rng);