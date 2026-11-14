// header file for system.cpp

// the shared state of the simulation, plus periodic boundary helpers

#pragma once
#include <cmath>
#include <random>
#include <vector>
#include "params.h"
#include "vec3.h"

// Everything that describes the particles at the current time. Units: A, amu, ps.
// The thermostat keeps its own state (eta, p_eta) in NoseHooverChain.
struct System {
    // fixed for the whole run
    int N = 0;
    double mass = 0.0;              // amu (shared by all particles)
    double L = 0.0;                 // cubic box length, A
    int n_dof = 0;                  // degrees of freedom, 3N - 3 (total momentum removed)

    // per-particle state
    std::vector<Vec3> positions;    // A, always wrapped into [0, L)
    std::vector<Vec3> velocities;   // A / ps
    std::vector<Vec3> forces;       // amu A / ps^2, forces at the current positions

    // set by compute_forces at the current positions
    double potential_energy = 0.0;  // amu A^2 / ps^2
    double virial = 0.0;            // sum_{i<j} r_ij . F_ij (minimum image), for pressure later

    // time
    long step = 0;
    double time = 0.0;              // ps
};

// FCC positions + Maxwell-Boltzmann velocities at p.T (init.h). Forces are left at zero:
// call compute_forces once before the first integration step.
System make_system(const Params& p, std::mt19937_64& rng);


// ---------------------------------------------------------------- periodic boundaries

// separation to the nearest periodic image: each component mapped into [-L/2, L/2].
// Assumes both particles are inside the box (|d| < L per component), which wrap_positions
// guarantees every step; a single +-L shift then suffices. Written without branches
// (d > L/2 -> d - L, d < -L/2 -> d + L) because in a liquid the branch is unpredictable;
// this is ~3x faster than d - L*round(d/L) and gives identical results.
inline double minimum_image_component(double d, double L) {
    const double half_L = 0.5 * L;
    return d - L * (static_cast<double>(d > half_L) - static_cast<double>(d < -half_L));
}

inline Vec3 minimum_image(Vec3 d, double L) {
    return Vec3{minimum_image_component(d.x, L),
                minimum_image_component(d.y, L),
                minimum_image_component(d.z, L)};
}

// map one coordinate into [0, L): leaving through one face re-enters through the opposite face
inline double wrap_coordinate(double x, double L) {
    x -= L * std::floor(x / L);
    if (x >= L) x -= L;   // guard: a tiny negative x gives x + L, which rounds to exactly L
    return x;
}

// wrap every particle back into the box (velocities are unchanged)
void wrap_positions(System& sys);