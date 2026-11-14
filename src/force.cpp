#include "../include/force.h"
#include <cmath>
#include <sstream>
#include <stdexcept>


LennardJones::LennardJones(double sigma_, double epsilon_, double rc_, double lambda_)
    : sigma(sigma_), epsilon(epsilon_), rc(rc_), lambda(lambda_) {
    if (sigma <= 0 || epsilon <= 0 || rc <= 0 || lambda <= 0)
        throw std::runtime_error("LennardJones: sigma, epsilon, rc and lambda must be positive");
    if (lambda >= rc) throw std::runtime_error("LennardJones: lambda must be smaller than rc");
    sigma2 = sigma * sigma;
    rc2 = rc * rc;
    r_switch = rc - lambda;
    r_switch2 = r_switch * r_switch;
}

LennardJones::LennardJones(const Params& p)
    : LennardJones(p.sigma, epsilon_internal(p), cutoff_radius(p), healing_length(p)) {}


double LennardJones::switching(double r) const {
    if (r < r_switch) return 1.0;
    if (r > rc) return 0.0;
    const double R = (r - r_switch) / lambda;
    return 1.0 + R * R * (2.0 * R - 3.0);
}

double LennardJones::switching_derivative(double r) const {
    if (r < r_switch || r > rc) return 0.0;
    const double R = (r - r_switch) / lambda;
    return 6.0 * R * (R - 1.0) / lambda;
}

double LennardJones::pair_energy(double r) const {
    if (r > rc) return 0.0;
    const double s6 = std::pow(sigma / r, 6);
    return 4.0 * epsilon * (s6 * s6 - s6) * switching(r);
}

double LennardJones::pair_dudr(double r) const {
    if (r > rc) return 0.0;
    const double s6 = std::pow(sigma / r, 6);
    const double u = 4.0 * epsilon * (s6 * s6 - s6);
    const double dudr = -24.0 * epsilon / r * (2.0 * s6 * s6 - s6);
    return dudr * switching(r) + u * switching_derivative(r);
}


void compute_forces(System& sys, const LennardJones& lj) {
    const int N = sys.N;
    const double L = sys.L;
    if (lj.rc >= 0.5 * L) {
        std::ostringstream msg;
        msg << "compute_forces: rc = " << lj.rc << " A must be < L/2 = " << 0.5 * L << " A";
        throw std::runtime_error(msg.str());
    }

    std::vector<Vec3>& F = sys.forces;
    const std::vector<Vec3>& pos = sys.positions;
    F.assign(N, Vec3{});
    double U = 0.0, W = 0.0;

    for (int i = 0; i < N - 1; ++i) {
        const Vec3 ri = pos[i];
        Vec3 Fi;                                        // accumulate force on i locally
        for (int j = i + 1; j < N; ++j) {
            const Vec3 d = minimum_image(ri - pos[j], L);   // r_ij to the nearest image of j
            const double r2 = dot(d, d);
            if (r2 >= lj.rc2) continue;

            const double s2 = lj.sigma2 / r2;
            const double s6 = s2 * s2 * s2;
            const double s12 = s6 * s6;
            double u = 4.0 * lj.epsilon * (s12 - s6);
            // f_over_r = -(1/r) d(uS)/dr ; inside the plateau S = 1 so it is -u'/r
            double f_over_r = 24.0 * lj.epsilon * (2.0 * s12 - s6) / r2;

            if (r2 > lj.r_switch2) {                    // switching shell: include u S'
                const double r = std::sqrt(r2);
                const double R = (r - lj.r_switch) / lj.lambda;
                const double S = 1.0 + R * R * (2.0 * R - 3.0);
                const double dS = 6.0 * R * (R - 1.0) / lj.lambda;
                f_over_r = f_over_r * S - u * dS / r;
                u *= S;
            }

            const Vec3 fij = f_over_r * d;              // force on i from j
            Fi += fij;
            F[j] -= fij;                                // Newton's third law
            U += u;
            W += f_over_r * r2;                         // r_ij . F_ij
        }
        F[i] += Fi;
    }

    sys.potential_energy = U;
    sys.virial = W;
}