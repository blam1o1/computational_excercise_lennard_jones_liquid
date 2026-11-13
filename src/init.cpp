#include "../include/init.h"
#include "../include/observables.h"
#include "../include/units.h"
#include <cmath>
#include <sstream>
#include <stdexcept>


// ---------------------------------------------------------------- positions

double box_length_from_density(int N, double density) {
    if (N <= 0) throw std::runtime_error("box_length_from_density: N must be positive");
    if (density <= 0.0) throw std::runtime_error("box_length_from_density: density must be positive");
    return std::cbrt(static_cast<double>(N) / density);
}


std::vector<Vec3> init_fcc_positions(int N, double L, double density, double rel_tol) {
    if (N <= 0) throw std::runtime_error("init_fcc_positions: N must be positive");
    if (L <= 0.0) throw std::runtime_error("init_fcc_positions: L must be positive");
    if (density <= 0.0) throw std::runtime_error("init_fcc_positions: density must be positive");

    // N, L and density are over-determined (density = N / L^3), so make sure they agree
    const double actual_density = N / (L * L * L);
    if (std::abs(actual_density - density) > rel_tol * density) {
        std::ostringstream msg;
        msg << "init_fcc_positions: N / L^3 = " << actual_density
            << " does not match requested density " << density
            << " (for N = " << N << " use L = " << box_length_from_density(N, density) << ")";
        throw std::runtime_error(msg.str());
    }

    // number of unit cells per side: N = 4 * n_cells^3
    const int n_cells = static_cast<int>(std::lround(std::cbrt(N / 4.0)));
    if (4 * n_cells * n_cells * n_cells != N) {
        std::ostringstream msg;
        msg << "init_fcc_positions: N = " << N << " is not 4*n^3 (e.g. 108, 256, 500, 864, 1372)";
        throw std::runtime_error(msg.str());
    }

    const double a = L / n_cells;  // lattice constant

    // 4-atom basis of the conventional cubic FCC cell, in units of a
    const Vec3 basis[4] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5},
    };

    std::vector<Vec3> positions;
    positions.reserve(N);
    for (int ix = 0; ix < n_cells; ++ix) {
        for (int iy = 0; iy < n_cells; ++iy) {
            for (int iz = 0; iz < n_cells; ++iz) {
                const Vec3 cell_origin{static_cast<double>(ix),
                                       static_cast<double>(iy),
                                       static_cast<double>(iz)};
                for (const Vec3& b : basis) positions.push_back(a * (cell_origin + b));
            }
        }
    }
    return positions;
}


// ---------------------------------------------------------------- velocities

std::pair<double, double> box_muller(std::mt19937_64& rng) {
    std::uniform_real_distribution<double> uniform(0.0, 1.0);

    // u1 must be in (0, 1] so that log(u1) is finite; uniform() is [0, 1)
    double u1;
    do { u1 = uniform(rng); } while (u1 <= 0.0);
    const double u2 = uniform(rng);

    const double radius = std::sqrt(-2.0 * std::log(u1));
    const double theta = 2.0 * units::PI * u2;
    return {radius * std::cos(theta), radius * std::sin(theta)};
}


void remove_com_velocity(std::vector<Vec3>& velocities, double mass) {
    const Vec3 v_cm = center_of_mass_velocity(velocities, mass);
    for (Vec3& v : velocities) v -= v_cm;
}


void rescale_to_temperature(std::vector<Vec3>& velocities, double mass, double T_target) {
    if (T_target <= 0.0) throw std::runtime_error("rescale_to_temperature: T_target must be positive");

    const double T_current = temperature(velocities, mass);
    if (T_current <= 0.0) throw std::runtime_error("rescale_to_temperature: current temperature is zero");

    // T scales as v^2, so (lambda v) has temperature lambda^2 * T_current = T_target
    const double lambda = std::sqrt(T_target / T_current);
    for (Vec3& v : velocities) v *= lambda;
}


std::vector<Vec3> init_velocities(const std::vector<Vec3>& positions, double mass,
                                  double T_target, std::mt19937_64& rng) {
    const std::size_t N = positions.size();
    if (N < 2) throw std::runtime_error("init_velocities: need at least 2 particles");
    if (mass <= 0.0) throw std::runtime_error("init_velocities: mass must be positive");
    if (T_target <= 0.0) throw std::runtime_error("init_velocities: T_target must be positive");

    // 1. 3N standard normal numbers, generated two at a time by Box-Muller
    std::vector<double> z(3 * N);
    for (std::size_t k = 0; k < z.size(); k += 2) {
        const std::pair<double, double> pair = box_muller(rng);
        z[k] = pair.first;
        if (k + 1 < z.size()) z[k + 1] = pair.second;  // drop the spare if 3N is odd
    }

    // 2. scale to the Maxwell-Boltzmann width: each component has sigma_v = sqrt(kB T / m)
    const double sigma_v = std::sqrt(units::KB * T_target / mass);
    std::vector<Vec3> velocities(N);
    for (std::size_t i = 0; i < N; ++i) {
        velocities[i] = sigma_v * Vec3{z[3 * i], z[3 * i + 1], z[3 * i + 2]};
    }

    // 3. remove center-of-mass drift
    remove_com_velocity(velocities, mass);

    // 4. rescale to exactly T_target (scaling keeps V_cm = 0)
    rescale_to_temperature(velocities, mass, T_target);

    return velocities;
}