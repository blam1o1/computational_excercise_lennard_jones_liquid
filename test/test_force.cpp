// Validation of the switched Lennard-Jones potential, the force loop and periodic boundaries.
// Build: g++ -std=c++17 -O2 -Wall -Wextra -o bin/test_force test/test_force.cpp src/force.cpp src/system.cpp src/params.cpp src/init.cpp src/observables.cpp
// Run:   ./bin/test_force

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "../include/force.h"
#include "../include/init.h"
#include "../include/params.h"
#include "../include/system.h"


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

// a bare System holding the given positions in a box of side L
static System make_box(const std::vector<Vec3>& positions, double L) {
    System sys;
    sys.N = static_cast<int>(positions.size());
    sys.mass = 39.98;
    sys.L = L;
    sys.positions = positions;
    sys.velocities.assign(sys.N, Vec3{});
    sys.forces.assign(sys.N, Vec3{});
    return sys;
}


int main() {
    std::cout << std::setprecision(8);
    const Params p;                         // argon defaults
    const LennardJones lj(p);
    const double sigma = lj.sigma, eps = lj.epsilon, rc = lj.rc, rs = lj.r_switch, lambda = lj.lambda;
    const double L = box_length(p);
    const double force_scale = eps / sigma; // natural force unit

    auto u_bare = [&](double r) { const double s6 = std::pow(sigma / r, 6); return 4 * eps * (s6 * s6 - s6); };
    auto du_bare = [&](double r) { const double s6 = std::pow(sigma / r, 6); return -24 * eps / r * (2 * s6 * s6 - s6); };

    std::cout << "\nrc = " << rc << " A, lambda = " << lambda << " A, switching starts at " << rs
              << " A, L = " << L << " A\n";

    // =====================================================================
    std::cout << "\n=== 1. Switching function S(r) ===\n";
    check(lj.switching(rs - 1e-9) == 1.0 && lj.switching(rs) == 1.0, "S = 1 up to rc - lambda");
    check(lj.switching(rc) == 0.0 && lj.switching(rc + 1e-9) == 0.0, "S = 0 at and beyond rc");
    check(std::abs(lj.switching(rs + 0.5 * lambda) - 0.5) < 1e-14, "S = 1/2 at the middle of the shell");
    // (rc - (rc - lambda)) / lambda is 1 only to round-off, so allow ~1e-13 / lambda
    check(lj.switching_derivative(rs) == 0.0 && std::abs(lj.switching_derivative(rc)) < 1e-12 / lambda,
          "dS/dr = 0 at both ends of the shell (force continuous)");
    bool monotone = true, bounded = true;
    for (int k = 0; k <= 1000; ++k) {
        const double r = rs + lambda * k / 1000.0;
        const double S = lj.switching(r);
        bounded = bounded && S >= 0.0 && S <= 1.0;
        if (k > 0) monotone = monotone && S <= lj.switching(rs + lambda * (k - 1) / 1000.0);
    }
    check(monotone && bounded, "S decreases monotonically from 1 to 0 across the shell");

    // =====================================================================
    std::cout << "\n=== 2. Pair energy u(r) S(r) ===\n";
    check(std::abs(lj.pair_energy(std::pow(2.0, 1.0 / 6.0) * sigma) + eps) < 1e-12 * eps,
          "u(2^(1/6) sigma) = -epsilon (well depth)");
    check(std::abs(lj.pair_energy(sigma)) < 1e-12 * eps, "u(sigma) = 0");
    bool plateau = true;
    for (double r = 0.9 * sigma; r < rs; r += 0.01 * sigma)
        plateau = plateau && std::abs(lj.pair_energy(r) - u_bare(r)) <= 1e-14 * std::abs(u_bare(r));
    check(plateau, "pair energy equals bare LJ for r < rc - lambda");
    check(lj.pair_energy(rc + 1e-6) == 0.0 && lj.pair_energy(2 * rc) == 0.0, "pair energy is 0 beyond rc");

    // =====================================================================
    std::cout << "\n=== 3. Pair force -d(uS)/dr vs finite differences (includes u S') ===\n";
    {
        const double h = 1e-6;
        double worst = 0.0, worst_without_uSprime = 0.0;
        for (double r = 0.9 * sigma; r < rc + 0.1 * sigma; r += 0.0037 * sigma) {
            const double fd = (lj.pair_energy(r + h) - lj.pair_energy(r - h)) / (2 * h);
            const double scale = std::max(std::abs(fd), 1e-3 * force_scale);
            worst = std::max(worst, std::abs(lj.pair_dudr(r) - fd) / scale);
            if (r > rs && r < rc)   // what you'd get if the u S' term were dropped
                worst_without_uSprime = std::max(worst_without_uSprime,
                                                 std::abs(du_bare(r) * lj.switching(r) - fd) / scale);
        }
        std::cout << "  max relative error of d(uS)/dr:            " << worst << "\n"
                  << "  same, with the u S' term left out (shell): " << worst_without_uSprime << "\n";
        check(worst < 1e-6, "d(uS)/dr matches finite differences on [0.9 sigma, rc + 0.1 sigma]");
        check(worst_without_uSprime > 1e3 * worst, "the test is sensitive to the u S' term");
    }

    // =====================================================================
    std::cout << "\n=== 4. Continuity at rc - lambda and rc ===\n";
    {
        const double d = 1e-9;
        check(std::abs(lj.pair_energy(rs + d) - lj.pair_energy(rs - d)) < 1e-6 * eps, "energy continuous at rc - lambda");
        check(std::abs(lj.pair_dudr(rs + d) - lj.pair_dudr(rs - d)) < 1e-6 * std::abs(lj.pair_dudr(rs)),
              "force continuous at rc - lambda");
        check(std::abs(lj.pair_energy(rc - d)) < 1e-12 * eps, "energy -> 0 at rc");
        // just inside rc the force vanishes linearly: |F(rc - delta)| ~ 6 |u(rc)| delta / lambda^2
        const double f6 = std::abs(lj.pair_dudr(rc - 1e-6)), f9 = std::abs(lj.pair_dudr(rc - 1e-9));
        std::cout << "  |F(rc - 1e-6 A)| = " << f6 << ", |F(rc - 1e-9 A)| = " << f9 << "\n";
        check(f9 < 1e-6 * force_scale && std::abs(f6 / f9 - 1000.0) < 10.0, "force -> 0 linearly at rc");
    }

    // =====================================================================
    std::cout << "\n=== 5. Two particles and periodic boundaries ===\n";
    {
        const double r0 = 4.0;   // A, attractive side of the well
        System in_box = make_box({{10, 10, 10}, {10 + r0, 10, 10}}, L);
        compute_forces(in_box, lj);
        check(std::abs(in_box.potential_energy - lj.pair_energy(r0)) < 1e-12 * eps, "U = u(r)S(r) for one pair");
        check(std::abs(in_box.forces[0].x - lj.pair_dudr(r0)) < 1e-12 * force_scale &&
              magnitude(in_box.forces[0] + in_box.forces[1]) < 1e-12 * force_scale,
              "F on particle 1 = -d(uS)/dr along r_12, F_2 = -F_1");

        // same pair split across the x faces: x = 2 and x = L - 2 are 4 A apart through the boundary
        System across = make_box({{2, 10, 10}, {L - 2, 10, 10}}, L);
        compute_forces(across, lj);
        check(std::abs(across.potential_energy - in_box.potential_energy) < 1e-12 * eps,
              "pair across a face has the same energy as the pair inside the box");
        check(std::abs(across.forces[0].x + in_box.forces[0].x) < 1e-12 * force_scale,
              "pair across a face: particle at x=2 is pulled toward its neighbour's image at x=-2");

        // across three faces: (1,1,1) and (L-1,L-1,L-1) are sqrt(12) A apart
        System corner = make_box({{1, 1, 1}, {L - 1, L - 1, L - 1}}, L);
        compute_forces(corner, lj);
        check(std::abs(corner.potential_energy - lj.pair_energy(std::sqrt(12.0))) < 1e-12 * eps,
              "pair across a corner (three faces) uses the nearest image");

        // separation L/2 - 1 > rc, both directly and through the boundary: no interaction
        System far = make_box({{1, 10, 10}, {1 + 0.5 * L - 1, 10, 10}}, L);
        compute_forces(far, lj);
        check(far.potential_energy == 0.0 && magnitude(far.forces[0]) == 0.0, "no interaction beyond rc");
    }

    // =====================================================================
    std::cout << "\n=== 6. Perfect FCC lattice (864 atoms) ===\n";
    const std::vector<Vec3> lattice = init_fcc_positions(p.N, L, p.density);
    {
        System sys = make_box(lattice, L);
        compute_forces(sys, lj);
        double fmax = 0.0;
        for (const Vec3& f : sys.forces) fmax = std::max(fmax, magnitude(f));
        std::cout << "  max |F_i| = " << fmax << " amu A/ps^2 (eps/sigma = " << force_scale << ")\n";
        check(fmax < 1e-9 * force_scale, "all forces vanish on the perfect lattice (inversion symmetry)");

        // independent lattice sum: enumerate FCC lattice vectors directly, no box, no minimum image
        const double a = L / 6.0;
        const int nmax = static_cast<int>(std::ceil(rc / a)) + 1;
        const Vec3 basis[4] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
        double U_per_atom = 0.0;
        for (int i = -nmax; i <= nmax; ++i)
            for (int j = -nmax; j <= nmax; ++j)
                for (int k = -nmax; k <= nmax; ++k)
                    for (const Vec3& b : basis) {
                        const double r = a * magnitude(Vec3{i + b.x, j + b.y, k + b.z});
                        if (r > 0 && r < rc) U_per_atom += 0.5 * lj.pair_energy(r);
                    }
        std::cout << "  U/N = " << sys.potential_energy / p.N / eps << " eps (direct lattice sum: "
                  << U_per_atom / eps << " eps)\n";
        check(std::abs(sys.potential_energy / p.N - U_per_atom) < 1e-10 * eps,
              "lattice energy matches an independent lattice sum (validates minimum image loop)");
    }

    // =====================================================================
    std::cout << "\n=== 7. Disordered configuration: F = -grad U, sum F = 0, virial ===\n";
    {
        std::mt19937_64 rng(7);
        std::uniform_real_distribution<double> shift(-0.3, 0.3);
        std::vector<Vec3> pos = lattice;
        for (Vec3& r : pos) r += Vec3{shift(rng), shift(rng), shift(rng)};
        System sys = make_box(pos, L);
        wrap_positions(sys);
        compute_forces(sys, lj);

        Vec3 total;
        double fmax = 0.0;
        for (const Vec3& f : sys.forces) { total += f; fmax = std::max(fmax, magnitude(f)); }
        std::cout << "  |sum F| = " << magnitude(total) << ", max |F_i| = " << fmax << "\n";
        check(magnitude(total) < 1e-12 * fmax * p.N, "sum of forces vanishes (Newton's third law)");

        // finite-difference gradient for 60 random Cartesian components
        std::uniform_int_distribution<int> pick(0, 3 * p.N - 1);
        const double h = 1e-5;
        double worst = 0.0;
        for (int t = 0; t < 60; ++t) {
            const int k = pick(rng), i = k / 3, c = k % 3;
            auto coord = [&](System& s) -> double& { return c == 0 ? s.positions[i].x : (c == 1 ? s.positions[i].y : s.positions[i].z); };
            System plus = sys, minus = sys;
            coord(plus) += h;
            coord(minus) -= h;
            compute_forces(plus, lj);
            compute_forces(minus, lj);
            const double fd = -(plus.potential_energy - minus.potential_energy) / (2 * h);
            const double analytic = c == 0 ? sys.forces[i].x : (c == 1 ? sys.forces[i].y : sys.forces[i].z);
            worst = std::max(worst, std::abs(fd - analytic) / std::max(std::abs(analytic), 1e-2 * force_scale));
        }
        std::cout << "  max relative error of F vs -dU/dr (60 components): " << worst << "\n";
        check(worst < 1e-5, "forces equal -grad U (finite differences)");

        // virial: scaling all coordinates and L by s gives dU/ds|_{s=1} = -W
        auto U_scaled = [&](double s) {
            std::vector<Vec3> q = sys.positions;
            for (Vec3& r : q) r *= s;
            System scaled = make_box(q, s * L);
            compute_forces(scaled, lj);
            return scaled.potential_energy;
        };
        const double hs = 1e-6;
        const double dUds = (U_scaled(1 + hs) - U_scaled(1 - hs)) / (2 * hs);
        std::cout << "  virial W = " << sys.virial << ", -dU/ds = " << -dUds << "\n";
        check(std::abs(sys.virial + dUds) < 1e-6 * std::abs(sys.virial), "virial W = -dU/ds (minimum image)");
    }

    // =====================================================================
    std::cout << "\n=== 8. Parameter guards ===\n";
    {
        Params q;
        q.N = 100;
        check_throws([&] { check_params(q); }, "rejects N that is not 4 n^3");
        q = Params(); q.N = 32;   // L = 11.7 A, L/2 = 5.85 A < rc = 8.51 A
        check_throws([&] { check_params(q); }, "rejects rc >= L/2");
        check_throws([&] {
            System tiny = make_box(init_fcc_positions(32, box_length(q), q.density), box_length(q));
            compute_forces(tiny, lj);
        }, "compute_forces refuses rc >= L/2");
        q = Params(); q.lambda_sigma = 0.01;
        check_throws([&] { check_params(q); }, "rejects lambda too small to be resolved by dt");
        q = Params(); q.lambda_sigma = 1.5;
        check_throws([&] { check_params(q); }, "rejects switching that starts inside the LJ well");
        q = Params(); q.nhc_tau_steps = 10;
        check_throws([&] { check_params(q); }, "rejects tau < 20 dt");
        q = Params(); q.nhc_nsy = 5;
        check_throws([&] { check_params(q); }, "rejects nsy not in {1, 3, 7}");
        bool defaults_ok = true;
        try { check_params(Params()); } catch (const std::runtime_error&) { defaults_ok = false; }
        check(defaults_ok, "default parameters pass");
    }

    std::cout << "\n" << n_pass << " passed, " << n_fail << " failed\n";
    return n_fail == 0 ? 0 : 1;
}