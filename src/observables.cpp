#include "../include/observables.h"
#include "../include/system.h"
#include "../include/units.h"
#include <cmath>
#include <limits>
#include <stdexcept>

static const double NaN = std::numeric_limits<double>::quiet_NaN();


// ---------------------------------------------------------------- instantaneous quantities

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


// ---------------------------------------------------------------- running averages

void EnergyDrift::start(double E_initial) {
    if (E_initial == 0.0) throw std::runtime_error("EnergyDrift: reference energy E0 must be nonzero");
    E0 = E_initial;
    n = 0;
    sum_abs = 0.0;
}

void EnergyDrift::add(double E) {
    sum_abs += std::abs((E - E0) / E0);
    ++n;
}

double EnergyDrift::deviation(double E) const { return (E - E0) / std::abs(E0); }

double EnergyDrift::delta_E() const { return n > 0 ? sum_abs / n : NaN; }


void RunningAverage::reset() {
    n = 0;
    mean_ = 0.0;
    m2 = 0.0;
}

void RunningAverage::add(double x) {
    ++n;
    const double d = x - mean_;
    mean_ += d / n;
    m2 += d * (x - mean_);
}

double RunningAverage::mean() const { return n > 0 ? mean_ : NaN; }

double RunningAverage::std_dev() const { return n > 0 ? std::sqrt(m2 / n) : NaN; }


// ---------------------------------------------------------------- structure

RadialDistribution::RadialDistribution(double dr_, double L_) : dr(dr_), L(L_), r_max(0.5 * L_) {
    if (dr <= 0.0 || L <= 0.0) throw std::runtime_error("RadialDistribution: dr and L must be positive");
    n_bins = static_cast<int>(r_max / dr);
    if (n_bins < 1) throw std::runtime_error("RadialDistribution: dr must be smaller than L/2");
    // minimum-image distances reach at most sqrt(3) L/2 (box corner); +2 bins of margin
    histogram.assign(static_cast<int>(std::sqrt(3.0) * r_max / dr) + 2, 0);
}


void RadialDistribution::sample(const System& sys) {
    if (std::abs(sys.L - L) > 1e-12 * L)
        throw std::runtime_error("RadialDistribution::sample: box length differs from the histogram's");

    const int N = sys.N;
    const std::vector<Vec3>& pos = sys.positions;
    long long* h = histogram.data();

    for (int i = 0; i < N - 1; ++i) {
        const Vec3 ri = pos[i];
        for (int j = i + 1; j < N; ++j) {
            const Vec3 d = minimum_image(ri - pos[j], L);
            h[static_cast<int>(std::sqrt(dot(d, d)) / dr)] += 2;   // once for each atom of the pair
        }
    }
    ++n_conf;
}


std::vector<double> RadialDistribution::radii() const {
    std::vector<double> r_i(n_bins);
    for (int i = 0; i < n_bins; ++i) r_i[i] = r(i);
    return r_i;
}


std::vector<double> RadialDistribution::g(int N) const {
    if (n_conf == 0) throw std::runtime_error("RadialDistribution::g: no frames sampled");
    const double rho = N / (L * L * L);
    std::vector<double> g_r(n_bins);
    for (int i = 0; i < n_bins; ++i) {
        const double ideal = 4.0 * units::PI * rho * r(i) * r(i) * dr * n_conf * N;
        g_r[i] = histogram[i] / ideal;
    }
    return g_r;
}


std::vector<double> structure_factor(const std::vector<double>& r, const std::vector<double>& g,
                                     double dr, double rho, const std::vector<double>& q) {
    if (r.size() != g.size()) throw std::runtime_error("structure_factor: r and g differ in length");

    std::vector<double> S(q.size());
    for (std::size_t k = 0; k < q.size(); ++k) {
        double integral = 0.0;
        for (std::size_t i = 0; i < r.size(); ++i) {
            const double x = q[k] * r[i];
            const double sinc = (x == 0.0) ? 1.0 : std::sin(x) / x;
            integral += r[i] * r[i] * (g[i] - 1.0) * sinc;
        }
        S[k] = 1.0 + 4.0 * units::PI * rho * dr * integral;
    }
    return S;
}


std::vector<double> q_grid(double q_min, double q_max, double dq) {
    if (dq <= 0.0 || q_max < q_min) throw std::runtime_error("q_grid: need dq > 0 and q_max >= q_min");
    const long n = static_cast<long>(std::floor((q_max - q_min) / dq + 1e-9)) + 1;
    std::vector<double> q(n);
    for (long k = 0; k < n; ++k) q[k] = q_min + k * dq;
    return q;
}