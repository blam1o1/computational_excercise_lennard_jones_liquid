// header file for observables.cpp

// quantities measured from the current state of the system. These are used by
// init.cpp now, and will be reused by the thermostat and the output routines.

#pragma once
#include <cstddef>
#include <vector>
#include "vec3.h"

// Number of degrees of freedom once total momentum has been removed: 3N - 3.
// The Nose-Hoover chain must use this same value for its target kinetic energy.
inline int degrees_of_freedom(std::size_t N) {
    return 3 * static_cast<int>(N) - 3;
}

// V_cm = sum_i(m_i * v_i) / sum_i(m_i), all particles sharing the same mass
Vec3 center_of_mass_velocity(const std::vector<Vec3>& velocities, double mass);

// K = 1/2 * sum_i(m_i * |v_i|^2)    [amu A^2 / ps^2]
double kinetic_energy(const std::vector<Vec3>& velocities, double mass);

// T = sum_i(m_i * |v_i|^2) / ((3N - 3) * kB)    [K]
double temperature(const std::vector<Vec3>& velocities, double mass);