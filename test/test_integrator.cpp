// Validation of the NHC + velocity-Verlet integrator (Trotter factorization, Suzuki-Yoshida NHC).
// Uses N = 108 (L = 17.54 A > 2 rc) so the finite-difference Jacobian stays cheap.
// Build: g++ -std=c++17 -O2 -Wall -Wextra -o bin/test_integrator test/test_integrator.cpp src/integrator.cpp src/thermostat.cpp src/force.cpp src/system.cpp src/params.cpp src/init.cpp src/observables.cpp
// Run:   ./bin/test_integrator [seed]

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "../include/force.h"
#include "../include/init.h"
#include "../include/integrator.h"
#include "../include/observables.h"
#include "../include/params.h"
#include "../include/system.h"
#include "../include/thermostat.h"
#include "../include/units.h"


// ---------------------------------------------------------------- tiny test harness
static int n_pass = 0, n_fail = 0;

static void check(bool ok, const std::string& what) {
    std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << what << "\n";
    (ok ? n_pass : n_fail)++;
}

// ---------------------------------------------------------------- helpers

static Params small_params() {
    Params p;
    p.N = 108;   // L = 17.54 A, rc = 8.51 A < L/2 = 8.77 A
    return p;
}

// flat phase-space vector: positions, velocities, eta, p_eta  (dimension 6N + 2M)
static std::vector<double> get_state(const System& s, const NoseHooverChain& nhc) {
    std::vector<double> x;
    for (const Vec3& r : s.positions) { x.push_back(r.x); x.push_back(r.y); x.push_back(r.z); }
    for (const Vec3& v : s.velocities) { x.push_back(v.x); x.push_back(v.y); x.push_back(v.z); }
    for (double e : nhc.eta) x.push_back(e);
    for (double pe : nhc.p_eta) x.push_back(pe);
    return x;
}

static void set_state(const std::vector<double>& x, System& s, NoseHooverChain& nhc, const LennardJones& lj) {
    std::size_t k = 0;
    for (Vec3& r : s.positions) { r.x = x[k++]; r.y = x[k++]; r.z = x[k++]; }
    for (Vec3& v : s.velocities) { v.x = x[k++]; v.y = x[k++]; v.z = x[k++]; }
    for (double& e : nhc.eta) e = x[k++];
    for (double& pe : nhc.p_eta) pe = x[k++];
    compute_forces(s, lj);   // F must correspond to the positions
}

// ln|det A| by LU decomposition with partial pivoting
static double log_abs_det(std::vector<double> A, int n, int& sign) {
    sign = 1;
    double logdet = 0.0;
    for (int c = 0; c < n; ++c) {
        int piv = c;
        for (int r = c + 1; r < n; ++r)
            if (std::abs(A[r * n + c]) > std::abs(A[piv * n + c])) piv = r;
        if (piv != c) { for (int k = 0; k < n; ++k) std::swap(A[c * n + k], A[piv * n + k]); sign = -sign; }
        const double d = A[c * n + c];
        if (d < 0) sign = -sign;
        logdet += std::log(std::abs(d));
        for (int r = c + 1; r < n; ++r) {
            const double f = A[r * n + c] / d;
            if (f == 0.0) continue;
            for (int k = c; k < n; ++k) A[r * n + k] -= f * A[c * n + k];
        }
    }
    return logdet;
}


int main(int argc, char* argv[]) {
    std::cout << std::setprecision(6);
    Params p = small_params();
    if (argc > 1) p.seed = static_cast<unsigned>(std::stoul(argv[1]));
    const LennardJones lj(p);
    const double dt = p.dt;

    std::mt19937_64 rng(p.seed);
    const System start = make_system(p, rng);
    const int N = start.N, n_dof = start.n_dof;
    std::cout << "N = " << N << ", L = " << start.L << " A, N_f = " << n_dof << ", dt = " << dt << " ps\n";

    // =====================================================================
    std::cout << "\n=== 1. Time reversibility (n steps forward, flip v and p_eta, n steps back) ===\n";
    {
        // Returns the largest deviation from the starting state after forward + reversed legs.
        // `step` advances (System, NHC) by one time step.
        auto reversal_error = [&](double tau, int n, auto step) {
            System s = start;
            rescale_to_temperature(s.velocities, s.mass, 200.0);   // start hot so the chain is active
            NoseHooverChain nhc(p.nhc_M, p.T, n_dof, tau, p.nhc_nc, p.nhc_nsy);
            compute_forces(s, lj);
            for (int k = 0; k < 100; ++k) integrate_step(s, lj, nhc, dt);
            const System s0 = s;
            const NoseHooverChain n0 = nhc;
            for (int k = 0; k < n; ++k) step(s, nhc);
            for (Vec3& v : s.velocities) v *= -1.0;
            for (double& pe : nhc.p_eta) pe = -pe;
            for (int k = 0; k < n; ++k) step(s, nhc);
            for (Vec3& v : s.velocities) v *= -1.0;
            for (double& pe : nhc.p_eta) pe = -pe;
            double err = 0;
            for (int i = 0; i < N; ++i) {
                err = std::max(err, magnitude(minimum_image(s.positions[i] - s0.positions[i], s.L)));
                err = std::max(err, magnitude(s.velocities[i] - s0.velocities[i]));
            }
            for (int j = 0; j < nhc.M; ++j) {
                err = std::max(err, std::abs(nhc.eta[j] - n0.eta[j]));
                err = std::max(err, std::abs(nhc.p_eta[j] - n0.p_eta[j]) / nhc.Q[j]);
            }
            return err;
        };
        auto proper = [&](System& s, NoseHooverChain& nhc) { integrate_step(s, lj, nhc, dt); };
        // control: same pieces in a NON-symmetric order (full NHC step, then velocity Verlet)
        auto lopsided = [&](System& s, NoseHooverChain& nhc) {
            nhc.half_step(s, 2.0 * dt);
            for (int i = 0; i < N; ++i) s.velocities[i] += 0.5 * dt / s.mass * s.forces[i];
            for (int i = 0; i < N; ++i) s.positions[i] += dt * s.velocities[i];
            wrap_positions(s);
            compute_forces(s, lj);
            for (int i = 0; i < N; ++i) s.velocities[i] += 0.5 * dt / s.mass * s.forces[i];
        };
        // The chain is chaotic at strong coupling, so round-off grows ~10x per 0.25 ps of leg
        // length; legs of 100 steps keep it well below the tolerance. A genuinely
        // irreversible scheme fails by many orders of magnitude (control below).
        const double e_prod = reversal_error(100 * dt, 300, proper);
        const double e_strong = reversal_error(20 * dt, 100, proper);
        const double e_control = reversal_error(20 * dt, 100, lopsided);
        std::cout << "  max deviation, tau = 100 dt, 300 + 300 steps: " << e_prod << "\n"
                  << "  max deviation, tau =  20 dt, 100 + 100 steps: " << e_strong << "\n"
                  << "  control (non-symmetric splitting):            " << e_control << "\n";
        check(e_prod < 1e-10, "trajectory retraces itself, production coupling (round-off level)");
        check(e_strong < 1e-8, "trajectory retraces itself, strong coupling (round-off level)");
        check(e_control > 1e3 * e_strong, "test detects a non-reversible splitting");
    }

    // =====================================================================
    std::cout << "\n=== 2. Generalized Liouville: Jacobian of the 5-step map ===\n";
    {
        System s = start;
        rescale_to_temperature(s.velocities, s.mass, 300.0);
        NoseHooverChain nhc(p.nhc_M, p.T, n_dof, 10 * dt, p.nhc_nc, p.nhc_nsy);
        compute_forces(s, lj);
        for (int k = 0; k < 20; ++k) integrate_step(s, lj, nhc, dt);

        const int n_map = 5, M = nhc.M;
        const std::vector<double> x0 = get_state(s, nhc);
        const int dim = static_cast<int>(x0.size());

        auto flow = [&](const std::vector<double>& x) {
            set_state(x, s, nhc, lj);
            for (int k = 0; k < n_map; ++k) integrate_step(s, lj, nhc, dt);
            return get_state(s, nhc);
        };
        const std::vector<double> xf = flow(x0);

        // per-variable scale (det is unchanged by the similarity transform D^-1 J D)
        std::vector<double> scale(dim, 1.0);
        for (int j = 0; j < M; ++j) scale[6 * N + M + j] = nhc.Q[j];

        // Central differences with h = 1e-6: the force is C1 but its derivative jumps at rc - lambda
        // and rc, so a pair sitting within h of those points spoils one column by O(h); smaller h
        // trades that for round-off. At 1e-6 both stay below ~1e-5 in ln det J.
        std::vector<double> J(static_cast<std::size_t>(dim) * dim);
        for (int k = 0; k < dim; ++k) {
            const double h = 1e-6 * scale[k];
            std::vector<double> xp = x0, xm = x0;
            xp[k] += h;
            xm[k] -= h;
            const std::vector<double> fp = flow(xp), fm = flow(xm);
            for (int i = 0; i < dim; ++i) {
                double diff = fp[i] - fm[i];
                if (i < 3 * N) diff = minimum_image_component(diff, s.L);   // a wrap is a shift by L
                J[static_cast<std::size_t>(i) * dim + k] = diff / (2 * h) * scale[k] / scale[i];
            }
        }
        int sign = 0;
        const double logdet = log_abs_det(J, dim, sign);
        double sum_rest = 0.0;
        for (int j = 1; j < M; ++j) sum_rest += xf[6 * N + j] - x0[6 * N + j];
        const double deta1 = xf[6 * N] - x0[6 * N];
        const double predicted = -(3.0 * N * deta1 + sum_rest);
        const double wrong = -(n_dof * deta1 + sum_rest);
        std::cout << "  dim = " << dim << ", sign(det J) = " << sign << "\n" << std::setprecision(10)
                  << "  ln det J (finite differences)       = " << logdet << "\n"
                  << "  -(3N d_eta1 + sum_{j>=2} d_eta_j)   = " << predicted << "\n"
                  << "  (same with N_f instead of 3N        = " << wrong << ")\n" << std::setprecision(6);
        const double tol = 1e-4;
        check(sign == 1 && std::abs(logdet - predicted) < tol,
              "det J = exp(-(3N d_eta1 + sum d_eta_j)): invariant measure preserved");
        check(std::abs(predicted - wrong) > 100 * tol, "test resolves 3N vs N_f (it is sharp)");
    }

    // =====================================================================
    std::cout << "\n=== 3. Order of the NHC propagator alone (Suzuki-Yoshida) ===\n";
    {
        // NHC flow only (no forces, no drift), 1 ps with outer step h; error vs a fine reference
        auto run = [&](int nsy, int nc) {
            System s = start;
            rescale_to_temperature(s.velocities, s.mass, 250.0);
            const double tau = 0.5;
            NoseHooverChain nhc(5, p.T, n_dof, tau, nc, nsy);
            const double init[5] = {2.0, -1.0, 0.5, -0.2, 0.1};
            for (int j = 0; j < 5; ++j) nhc.p_eta[j] = init[j] * nhc.Q[j] / tau;
            const double h = 0.005;
            for (int k = 0; k < 200; ++k) { nhc.half_step(s, h); nhc.half_step(s, h); }
            std::vector<double> x;
            for (const Vec3& v : s.velocities) { x.push_back(v.x); x.push_back(v.y); x.push_back(v.z); }
            for (int j = 0; j < 5; ++j) { x.push_back(nhc.eta[j]); x.push_back(nhc.p_eta[j] / nhc.Q[j]); }
            return x;
        };
        const std::vector<double> ref = run(7, 64);
        auto err = [&](const std::vector<double>& x) {
            double e = 0;
            for (std::size_t i = 0; i < x.size(); ++i) e = std::max(e, std::abs(x[i] - ref[i]));
            return e;
        };
        for (int nsy : {1, 3, 7}) {
            const double e1 = err(run(nsy, 1)), e2 = err(run(nsy, 2)), e4 = err(run(nsy, 4));
            const double order = 0.5 * (std::log2(e1 / e2) + std::log2(e2 / e4));
            const int expected = nsy == 1 ? 2 : (nsy == 3 ? 4 : 6);
            std::cout << "  nsy = " << nsy << ": error at nc = 1, 2, 4: " << e1 << ", " << e2 << ", " << e4
                      << "  -> ratios " << e1 / e2 << ", " << e2 / e4 << "\n";
            check(std::abs(order - expected) < 0.3,
                  "nsy = " + std::to_string(nsy) + " converges at order " + std::to_string(expected));
        }
    }

    // =====================================================================
    std::cout << "\n=== 4. Global order: H' error vs dt (same initial state, 2 ps) ===\n";
    {
        System s0 = start;
        NoseHooverChain n0(p.nhc_M, p.T, n_dof, 0.5, p.nhc_nc, p.nhc_nsy);   // tau fixed in ps
        compute_forces(s0, lj);
        for (int k = 0; k < 400; ++k) integrate_step(s0, lj, n0, dt);

        std::vector<double> nvt_err, nve_err;
        for (double h : {0.010, 0.005, 0.0025}) {
            System s = s0;
            NoseHooverChain nhc = n0;
            const double H0 = conserved_energy(s, nhc);
            double mean = 0.0;
            const int steps = static_cast<int>(std::lround(2.0 / h));
            for (int k = 0; k < steps; ++k) {
                integrate_step(s, lj, nhc, h);
                mean += std::abs((conserved_energy(s, nhc) - H0) / H0);
            }
            nvt_err.push_back(mean / steps);

            System e = s0;   // same state, thermostat off
            const double E0 = conserved_energy(e);
            mean = 0.0;
            for (int k = 0; k < steps; ++k) {
                integrate_step(e, lj, h);
                mean += std::abs((conserved_energy(e) - E0) / E0);
            }
            nve_err.push_back(mean / steps);
        }
        std::cout << "  NVT <|dH'/H'|> at dt = 10, 5, 2.5 fs: " << nvt_err[0] << ", " << nvt_err[1] << ", "
                  << nvt_err[2] << "  -> ratios " << nvt_err[0] / nvt_err[1] << ", " << nvt_err[1] / nvt_err[2] << "\n"
                  << "  NVE <|dH/H|>  at dt = 10, 5, 2.5 fs: " << nve_err[0] << ", " << nve_err[1] << ", "
                  << nve_err[2] << "  -> ratios " << nve_err[0] / nve_err[1] << ", " << nve_err[1] / nve_err[2] << "\n";
        auto second_order = [](const std::vector<double>& e) {
            return e[0] / e[1] > 3.0 && e[0] / e[1] < 5.0 && e[1] / e[2] > 3.0 && e[1] / e[2] < 5.0;
        };
        check(second_order(nvt_err), "NVT: H' error falls ~4x when dt halves (second order)");
        check(second_order(nve_err), "NVE: H error falls ~4x when dt halves (second order)");
        check(nvt_err[1] < 1e-4, "NVT: H' conserved to < 1e-4 at dt = 5 fs");
    }

    // =====================================================================
    std::cout << "\n=== 5. Canonical sampling, momentum, boundaries (10 ps equilibration + 200 ps) ===\n";
    {
        System s = start;
        NoseHooverChain nhc(p, n_dof);   // input defaults: M = 5, tau = 100 dt, nc = 1, nsy = 7
        compute_forces(s, lj);
        const double H0 = conserved_energy(s, nhc);
        for (int k = 0; k < 2000; ++k) integrate_step(s, lj, nhc, dt);

        const long n_samples = 40000;
        double sumT = 0, sumT2 = 0, max_drift = 0;
        for (long k = 0; k < n_samples; ++k) {
            integrate_step(s, lj, nhc, dt);
            const double T = temperature(s.velocities, s.mass);
            sumT += T;
            sumT2 += T * T;
            max_drift = std::max(max_drift, std::abs((conserved_energy(s, nhc) - H0) / H0));
        }
        const double meanT = sumT / n_samples;
        const double stdT = std::sqrt(sumT2 / n_samples - meanT * meanT);
        const double stdT_canonical = p.T * std::sqrt(2.0 / n_dof);
        std::cout << "  <T> = " << meanT << " K, std(T) = " << stdT << " K (canonical: T sqrt(2/N_f) = "
                  << stdT_canonical << " K)\n"
                  << "  max |H' - H'0| / |H'0| over 210 ps (incl. melting) = " << max_drift << "\n";
        check(std::abs(meanT - p.T) < 2.0, "<T> = target temperature");
        check(std::abs(stdT / stdT_canonical - 1.0) < 0.15, "temperature fluctuations are canonical");
        check(max_drift < 1e-3, "H' conserved over the whole run");

        Vec3 P;
        for (const Vec3& v : s.velocities) P += s.mass * v;
        const double p_scale = s.mass * std::sqrt(units::KB * p.T / s.mass) * std::sqrt(N);
        std::cout << "  |P| = " << magnitude(P) << " amu A/ps (thermal scale " << p_scale << ")\n";
        check(magnitude(P) < 1e-10 * p_scale, "total momentum stays zero");

        bool inside = true;
        for (const Vec3& r : s.positions)
            inside = inside && r.x >= 0 && r.x < s.L && r.y >= 0 && r.y < s.L && r.z >= 0 && r.z < s.L;
        check(inside, "all positions remain in [0, L)");
        check(s.step == 42000 && std::abs(s.time - 42000 * dt) < 1e-9, "step counter and time advance");

        // for contrast: NVE from the same state has smaller (microcanonical) T fluctuations
        System e = s;
        double eT = 0, eT2 = 0;
        for (long k = 0; k < 20000; ++k) {
            integrate_step(e, lj, dt);
            const double T = temperature(e.velocities, e.mass);
            eT += T;
            eT2 += T * T;
        }
        const double eMean = eT / 20000;
        std::cout << "  (NVE from the same state: std(T) = " << std::sqrt(eT2 / 20000 - eMean * eMean)
                  << " K, smaller than canonical as expected)\n";
    }

    // =====================================================================
    std::cout << "\n=== 6. Single Nose-Hoover thermostat (M = 1) code path ===\n";
    {
        System s = start;
        NoseHooverChain nhc(1, p.T, n_dof, 0.5, 1, 7);
        compute_forces(s, lj);
        const double H0 = conserved_energy(s, nhc);
        double worst = 0;
        for (int k = 0; k < 2000; ++k) {
            integrate_step(s, lj, nhc, dt);
            worst = std::max(worst, std::abs((conserved_energy(s, nhc) - H0) / H0));
        }
        std::cout << "  max |dH'/H'| over 10 ps = " << worst << "\n";
        check(worst < 5e-4, "M = 1 conserves H' (max over 10 ps from the lattice, incl. melting)");
    }

    std::cout << "\n" << n_pass << " passed, " << n_fail << " failed\n";
    return n_fail == 0 ? 0 : 1;
}