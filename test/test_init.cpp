// Validation of the initial conditions (positions + velocities) before any dynamics.
// Build:  g++ -std=c++17 -O2 -Wall -Wextra -o bin/test_init test/test_init.cpp src/init.cpp src/observables.cpp
// Run:    ./bin/test_init [seed]   (writes init.xyz and speed_hist.dat in the current directory)

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "../include/init.h"
#include "../include/observables.h"
#include "../include/units.h"


// ---------------------------------------------------------------- tiny test harness
static int n_pass = 0, n_fail = 0;

static void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << what << "\n";
    (ok ? n_pass : n_fail)++;
}

template <class F>
static void check_throws(F f, const std::string& what) {
    bool threw = false;
    try { f(); } catch (const std::runtime_error&) { threw = true; }
    check(threw, what);
}

// minimum-image separation in a cubic box (move to system.h when writing the force loop)
static Vec3 minimum_image(Vec3 d, double L) {
    d.x -= L * std::round(d.x / L);
    d.y -= L * std::round(d.y / L);
    d.z -= L * std::round(d.z / L);
    return d;
}

// Maxwell-Boltzmann speed CDF, x = v / sigma_v
static double mb_speed_cdf(double x) {
    return std::erf(x / std::sqrt(2.0)) - std::sqrt(2.0 / units::PI) * x * std::exp(-0.5 * x * x);
}

// Maxwell-Boltzmann speed pdf in 1/(A/ps)
static double mb_speed_pdf(double v, double sigma_v) {
    const double x = v / sigma_v;
    return std::sqrt(2.0 / units::PI) * x * x * std::exp(-0.5 * x * x) / sigma_v;
}


int main(int argc, char* argv[]) {
    // ---------------- inputs
    const int N = 864;               // 6 x 6 x 6 FCC supercell
    const double density = 0.02;     // A^-3
    const double mass = 39.98;       // amu (argon)
    const double T_target = 100.0;   // K
    const double sigma_lj = 3.405;   // A
    const double epsilon_K = 119.8;  // epsilon / kB in K
    const unsigned seed = (argc > 1) ? static_cast<unsigned>(std::stoul(argv[1])) : 42u;

    std::cout << std::setprecision(10);

    // =====================================================================
    std::cout << "\n=== Positions: FCC lattice ===\n";
    const double L = box_length_from_density(N, density);
    const int n_cells = 6;
    const double a = L / n_cells;
    std::cout << "  L = " << L << " A,  a = " << a << " A\n"
              << "  reduced density rho* = rho sigma^3 = " << density * std::pow(sigma_lj, 3) << "\n"
              << "  reduced temperature T* = T / epsilon = " << T_target / epsilon_K << "\n";

    const std::vector<Vec3> r = init_fcc_positions(N, L, density);

    check(static_cast<int>(r.size()) == N, "N positions generated");

    bool inside = true;
    for (const Vec3& p : r)
        inside = inside && p.x >= 0 && p.x < L && p.y >= 0 && p.y < L && p.z >= 0 && p.z < L;
    check(inside, "all coordinates in [0, L)");

    check(std::abs(N / (L * L * L) - density) < 1e-12 * density, "N / L^3 equals requested density");

    // neighbour shells with periodic boundaries. Perfect FCC: every particle has
    // 12 neighbours at a/sqrt(2) and 6 at a, and nothing closer than a/sqrt(2).
    const double d1 = a / std::sqrt(2.0), d2 = a, tol = 1e-8 * a;
    std::vector<int> shell1(N, 0), shell2(N, 0);
    double r_min = L;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            const double d = magnitude(minimum_image(r[i] - r[j], L));
            r_min = std::min(r_min, d);
            if (std::abs(d - d1) < tol) { shell1[i]++; shell1[j]++; }
            if (std::abs(d - d2) < tol) { shell2[i]++; shell2[j]++; }
        }
    }
    const bool fcc1 = std::all_of(shell1.begin(), shell1.end(), [](int c) { return c == 12; });
    const bool fcc2 = std::all_of(shell2.begin(), shell2.end(), [](int c) { return c == 6; });
    check(fcc1, "every particle has 12 nearest neighbours at a/sqrt(2) (periodic)");
    check(fcc2, "every particle has 6 second neighbours at a (periodic)");
    std::cout << "  minimum pair distance = " << r_min << " A = " << r_min / sigma_lj << " sigma\n";
    check(std::abs(r_min - d1) < tol, "minimum pair distance equals a/sqrt(2) (no overlaps/duplicates)");
    check(r_min > std::pow(2.0, 1.0 / 6.0) * sigma_lj, "minimum distance beyond LJ minimum 2^(1/6) sigma");

    check_throws([&] { init_fcc_positions(100, box_length_from_density(100, density), density); },
                 "throws when N is not 4 n^3");
    check_throws([&] { init_fcc_positions(N, 30.0, density); },
                 "throws when N, L and density are inconsistent");

    // =====================================================================
    std::cout << "\n=== Box-Muller: standard normal statistics (1e6 samples) ===\n";
    {
        std::mt19937_64 rng(seed);
        const int n = 1000000;
        std::vector<double> z;
        z.reserve(n);
        while (static_cast<int>(z.size()) < n) {
            const std::pair<double, double> p = box_muller(rng);
            z.push_back(p.first);
            z.push_back(p.second);
        }
        double m1 = 0, m2 = 0, m3 = 0, m4 = 0, in1 = 0, in2 = 0;
        for (double x : z) {
            m1 += x;
            if (std::abs(x) < 1) in1++;
            if (std::abs(x) < 2) in2++;
        }
        m1 /= n;
        for (double x : z) {
            const double d = x - m1;
            m2 += d * d; m3 += d * d * d; m4 += d * d * d * d;
        }
        m2 /= n; m3 /= n; m4 /= n;
        const double skew = m3 / std::pow(m2, 1.5), exkurt = m4 / (m2 * m2) - 3.0;
        std::cout << "  mean = " << m1 << ", var = " << m2 << ", skew = " << skew
                  << ", excess kurtosis = " << exkurt << "\n"
                  << "  P(|z|<1) = " << in1 / n << " (0.6827), P(|z|<2) = " << in2 / n << " (0.9545)\n";
        // tolerances are 5 standard errors
        check(std::abs(m1) < 5.0 / std::sqrt(n), "mean = 0");
        check(std::abs(m2 - 1.0) < 5.0 * std::sqrt(2.0 / n), "variance = 1");
        check(std::abs(skew) < 5.0 * std::sqrt(6.0 / n), "skewness = 0");
        check(std::abs(exkurt) < 5.0 * std::sqrt(24.0 / n), "excess kurtosis = 0");
        check(std::abs(in1 / n - 0.682689) < 5.0 * std::sqrt(0.2173 / n), "68.27% within 1 sigma");
        check(std::abs(in2 / n - 0.954500) < 5.0 * std::sqrt(0.0434 / n), "95.45% within 2 sigma");
    }

    // =====================================================================
    std::cout << "\n=== Velocities: Maxwell-Boltzmann at T = " << T_target << " K ===\n";
    std::mt19937_64 rng(seed);
    const std::vector<Vec3> v = init_velocities(r, mass, T_target, rng);
    const double sigma_v = std::sqrt(units::KB * T_target / mass);
    std::cout << "  sigma_v = sqrt(kB T / m) = " << sigma_v << " A/ps\n";

    check(v.size() == r.size(), "one velocity per particle");

    const Vec3 v_cm = center_of_mass_velocity(v, mass);
    std::cout << "  V_cm = " << v_cm << " A/ps\n";
    check(magnitude(v_cm) < 1e-12 * sigma_v, "center-of-mass velocity is zero (no drift)");

    const double T_now = temperature(v, mass);
    std::cout << "  T (3N-3 dof) = " << T_now << " K,  K = " << kinetic_energy(v, mass) << " amu A^2/ps^2\n";
    check(std::abs(T_now - T_target) < 1e-10 * T_target, "instantaneous temperature equals T_target");

    // equipartition: each Cartesian direction carries (N-1) dof
    double sx = 0, sy = 0, sz = 0, s4 = 0, s2 = 0, mean_speed = 0;
    for (const Vec3& u : v) {
        sx += mass * u.x * u.x; sy += mass * u.y * u.y; sz += mass * u.z * u.z;
        for (double c : {u.x, u.y, u.z}) { s2 += c * c; s4 += c * c * c * c; }
        mean_speed += magnitude(u);
    }
    const double Tx = sx / ((N - 1) * units::KB), Ty = sy / ((N - 1) * units::KB),
                 Tz = sz / ((N - 1) * units::KB);
    std::cout << "  T_x, T_y, T_z = " << Tx << ", " << Ty << ", " << Tz << " K\n";
    const double tol_T = 5.0 * std::sqrt(2.0 / N) * T_target;
    check(std::abs(Tx - T_target) < tol_T && std::abs(Ty - T_target) < tol_T &&
          std::abs(Tz - T_target) < tol_T, "equipartition across x, y, z (within 5 standard errors)");

    const double kurt = (s4 / (3 * N)) / std::pow(s2 / (3 * N), 2);
    std::cout << "  <v^4>/<v^2>^2 = " << kurt << " (Gaussian: 3)\n";
    check(std::abs(kurt - 3.0) < 5.0 * std::sqrt(24.0 / (3 * N)), "velocity components are Gaussian-shaped");

    mean_speed /= N;
    const double mb_mean = std::sqrt(8.0 * units::KB * T_target / (units::PI * mass));
    std::cout << "  <|v|> = " << mean_speed << " A/ps (MB: " << mb_mean << ")\n";
    check(std::abs(mean_speed - mb_mean) < 5.0 * 0.4220 / std::sqrt(N) * mb_mean, "mean speed matches Maxwell-Boltzmann");

    // Kolmogorov-Smirnov test of speeds against the MB speed distribution
    std::vector<double> speeds;
    for (const Vec3& u : v) speeds.push_back(magnitude(u));
    std::sort(speeds.begin(), speeds.end());
    double ks = 0;
    for (int i = 0; i < N; ++i) {
        const double F = mb_speed_cdf(speeds[i] / sigma_v);
        ks = std::max({ks, std::abs(F - static_cast<double>(i) / N), std::abs(F - static_cast<double>(i + 1) / N)});
    }
    const double ks_crit = 1.949 / std::sqrt(N);  // 0.1% significance
    std::cout << "  KS statistic D = " << ks << " (0.1% critical value " << ks_crit << ")\n";
    check(ks < ks_crit, "speed distribution consistent with Maxwell-Boltzmann (KS, 0.1%)");

    // reproducibility
    std::mt19937_64 rng_same(seed), rng_other(seed + 1);
    const std::vector<Vec3> v_same = init_velocities(r, mass, T_target, rng_same);
    const std::vector<Vec3> v_other = init_velocities(r, mass, T_target, rng_other);
    bool identical = true, differs = false;
    for (int i = 0; i < N; ++i) {
        identical = identical && v[i].x == v_same[i].x && v[i].y == v_same[i].y && v[i].z == v_same[i].z;
        differs = differs || v[i].x != v_other[i].x;
    }
    check(identical, "same seed gives identical velocities");
    check(differs && std::abs(temperature(v_other, mass) - T_target) < 1e-10 * T_target,
          "different seed gives different velocities at the same T");

    check_throws([&] { std::vector<Vec3> one(1); init_velocities(one, mass, T_target, rng); },
                 "throws for fewer than 2 particles");
    check_throws([&] { init_velocities(r, -1.0, T_target, rng); }, "throws for non-positive mass");
    check_throws([&] { init_velocities(r, mass, 0.0, rng); }, "throws for non-positive temperature");

    // =====================================================================
    // files for visual inspection
    {
        std::ofstream xyz("init.xyz");  // extended XYZ: opens in OVITO / VMD
        xyz << std::setprecision(12);
        xyz << N << "\n"
            << "Lattice=\"" << L << " 0 0 0 " << L << " 0 0 0 " << L
            << "\" Properties=species:S:1:pos:R:3:velo:R:3\n";
        for (int i = 0; i < N; ++i)
            xyz << "Ar " << r[i].x << " " << r[i].y << " " << r[i].z << " "
                << v[i].x << " " << v[i].y << " " << v[i].z << "\n";

        std::ofstream hist("speed_hist.dat");  // columns: speed, sampled pdf, MB pdf
        const int n_bins = 30;
        const double v_max = 5.0 * sigma_v, dv = v_max / n_bins;
        std::vector<int> counts(n_bins, 0);
        for (double s : speeds) if (s < v_max) counts[static_cast<int>(s / dv)]++;
        hist << "# speed[A/ps]  sampled_pdf  maxwell_boltzmann_pdf\n";
        for (int b = 0; b < n_bins; ++b) {
            const double c = (b + 0.5) * dv;
            hist << c << " " << counts[b] / (N * dv) << " " << mb_speed_pdf(c, sigma_v) << "\n";
        }
    }
    std::cout << "\n  wrote init.xyz and speed_hist.dat\n";

    std::cout << "\n" << n_pass << " passed, " << n_fail << " failed\n";
    return n_fail == 0 ? 0 : 1;
}