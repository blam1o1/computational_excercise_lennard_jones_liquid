#include "../include/observables.h"
#include "../include/units.h"
#include <stdexcept>


Vec3 center_of_mass_velocity(const std::vector<Vec3>& velocities, double mass) {
    if (velocities.empty()) throw std::runtime_error("center_of_mass_velocity: no particles");

    Vec3 total_momentum;     // sum_i m_i v_i
    double total_mass = 0.0; // sum_i m_i
    for (const Vec3& v : velocities) {
        total_momentum += mass * v;
        total_mass += mass;
    }
    return total_momentum / total_mass;
}


double kinetic_energy(const std::vector<Vec3>& velocities, double mass) {
    double sum_mv2 = 0.0;
    for (const Vec3& v : velocities) sum_mv2 += mass * magnitude_squared(v);
    return 0.5 * sum_mv2;
}


double temperature(const std::vector<Vec3>& velocities, double mass) {
    const int n_dof = degrees_of_freedom(velocities.size());
    if (n_dof <= 0) throw std::runtime_error("temperature: need at least 2 particles");

    // 2K = sum_i m_i |v_i|^2
    return 2.0 * kinetic_energy(velocities, mass) / (n_dof * units::KB);
}