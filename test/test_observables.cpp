// Validation of the observables (Delta E, <T_kin>, g(r), S(q)) and of the output formats.
// Build: g++ -std=c++17 -O2 -Wall -Wextra -o bin/test_observables test/test_observables.cpp src/output.cpp src/observables.cpp src/integrator.cpp src/thermostat.cpp src/force.cpp src/system.cpp src/params.cpp src/init.cpp
// Run:   ./bin/test_observables [seed]

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include "../include/init.h"
#include "../include/integrator.h"
#include "../include/observables.h"
#include "../include/output.h"
#include "../include/params.h"
#include "../include/system.h"
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

static std::vector<std::string> split_lines(const std::string& text) {
    std::vector<std::string> lines;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line)) lines.push_back(line);
    return lines;
}

static std::vector<std::string> tokens(const std::string& line) {
    std::istringstream in(line);
    std::vector<std::string> t;
    std::string w;
    while (in >> w) t.push_back(w);
    return t;
}

// a bare System holding the given positions in a box of side L
static System make_box(const std::vector<Vec3>& positions, double L) {
    System sys;
    sys.N = static_cast<int>(positions.size());
    sys.mass = 39.98;
    sys.L = L;
    sys.n_dof = degrees_of_freedom(sys.N);
    sys.positions = positions;
    sys.velocities.assign(sys.N, Vec3{});
    sys.forces.assign(sys.N, Vec3{});
    return sys;
}


int main(int argc, char* argv[]) {
    std::cout << std::setprecision(8);
    Params p;
    if (argc > 1) p.seed = static_cast<unsigned>(std::stoul(argv[1]));
    const int N = p.N;
    const double L = box_length(p);
    const double rho = N / (L * L * L);
    const double PI = units::PI;

    // =====================================================================
    std::cout << "\n=== 1. Running accumulators ===\n";
    {
        EnergyDrift drift;
        drift.start(-10.0);
        check(std::isnan(drift.delta_E()), "Delta E is NaN before the first step");
        drift.add(-10.1);
        const double d1 = drift.delta_E();
        drift.add(-9.8);
        const double d2 = drift.delta_E();
        drift.add(-10.0);
        const double d3 = drift.delta_E();
        std::cout << "  E0 = -10, E_k = -10.1, -9.8, -10.0  ->  Delta E = " << d1 << ", " << d2 << ", " << d3 << "\n";
        check(std::abs(d1 - 0.01) < 1e-15 && std::abs(d2 - 0.015) < 1e-15 && std::abs(d3 - 0.01) < 1e-15,
              "Delta E = (1/n) sum |(E_k - E0)/E0|");
        check(std::abs(drift.deviation(-10.1) + 0.01) < 1e-15, "signed deviation (E - E0)/|E0|");
        check_throws([] { EnergyDrift d; d.start(0.0); }, "rejects E0 = 0");

        RunningAverage avg;
        check(std::isnan(avg.mean()), "<x> is NaN before the first sample");
        for (double x : {1.0, 2.0, 3.0, 4.0}) avg.add(x);
        check(std::abs(avg.mean() - 2.5) < 1e-15 && std::abs(avg.std_dev() - std::sqrt(1.25)) < 1e-15,
              "running mean 2.5 and std sqrt(1.25) for 1, 2, 3, 4");
        avg.reset();
        check(avg.n == 0 && std::isnan(avg.mean()), "reset clears the average");
    }

    // =====================================================================
    std::cout << "\n=== 2. g(r) of an ideal gas (uniform random points, 100 frames) ===\n";
    {
        std::mt19937_64 rng(p.seed);
        std::uniform_real_distribution<double> uniform(0.0, L);
        RadialDistribution rdf(p.gr_dr, L);
        std::vector<Vec3> pos(N);
        for (int frame = 0; frame < 100; ++frame) {
            for (Vec3& r : pos) r = Vec3{uniform(rng), uniform(rng), uniform(rng)};
            rdf.sample(make_box(pos, L));
        }
        const std::vector<double> g = rdf.g(N);
        // with rho = N/V the ideal-gas value is (N-1)/N: each atom sees N-1 others
        const double expected = (N - 1.0) / N;
        double sum_h = 0, sum_ideal = 0, worst_z = 0;
        for (int i = 0; i < rdf.n_bins; ++i) {
            if (rdf.r(i) < 1.0) continue;   // too few counts per bin below 1 A
            sum_h += rdf.histogram[i];
            sum_ideal += rdf.histogram[i] / g[i];
            const double sigma = expected * std::sqrt(2.0 / rdf.histogram[i]);   // h counts pairs twice
            worst_z = std::max(worst_z, std::abs(g[i] - expected) / sigma);
        }
        std::cout << "  bins = " << rdf.n_bins << ", overall g = " << sum_h / sum_ideal << " (expected (N-1)/N = "
                  << expected << "), largest deviation = " << worst_z << " standard errors\n";
        check(std::abs(sum_h / sum_ideal - expected) < 2e-3, "normalization: overall g = (N-1)/N");
        check(worst_z < 5.0, "every bin consistent with the ideal-gas value (within 5 standard errors)");
    }

    // =====================================================================
    std::cout << "\n=== 3. g(r) of the perfect FCC lattice: neighbour shells ===\n";
    {
        const std::vector<Vec3> lattice = init_fcc_positions(N, L, p.density);
        RadialDistribution rdf(p.gr_dr, L);
        rdf.sample(make_box(lattice, L));
        const double a = L / 6.0;
        const double shell[4] = {a / std::sqrt(2.0), a, a * std::sqrt(1.5), a * std::sqrt(2.0)};
        const double cumulative_expected[4] = {12, 18, 42, 54};
        bool shells_ok = true;
        for (int s = 0; s < 4; ++s) {
            const double r_cut = 0.5 * (shell[s] + (s < 3 ? shell[s + 1] : shell[s] + 1.0));
            double count = 0;
            for (int i = 0; i < rdf.n_bins && rdf.r(i) < r_cut; ++i) count += rdf.histogram[i];
            count /= N;
            std::cout << "  neighbours within " << r_cut << " A: " << count << " (expected " << cumulative_expected[s] << ")\n";
            shells_ok = shells_ok && std::abs(count - cumulative_expected[s]) < 1e-12;
        }
        check(shells_ok, "FCC shells hold 12, 6, 24, 12 neighbours (minimum image histogram)");

        const int b = static_cast<int>(shell[0] / p.gr_dr);
        const double g_expected = 12.0 / (4.0 * PI * rho * rdf.r(b) * rdf.r(b) * p.gr_dr);
        const std::vector<double> g = rdf.g(N);
        check(std::abs(g[b] - g_expected) < 1e-12 * g_expected,
              "g(r_i) = h(i) / (4 pi rho r_i^2 dr N_conf N) exactly");
        bool empty_below = true;
        for (int i = 0; i < b; ++i) empty_below = empty_below && rdf.histogram[i] == 0;
        check(empty_below, "no counts closer than the nearest-neighbour distance");

        // gr.dat round trip: the n(r) column must give 12 just past the first shell
        const std::string path = "test_gr_tmp.dat";
        write_gr(path, p, rdf, g);
        std::ifstream in(path);
        std::string line;
        double n_first_shell = -1;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ls(line);
            double r, gr, nr;
            ls >> r >> gr >> nr;
            if (r < 0.5 * (shell[0] + shell[1])) n_first_shell = nr;
        }
        std::remove(path.c_str());
        check(std::abs(n_first_shell - 12.0) < 1e-6, "gr.dat n(r) column reads 12 after the first shell");
    }

    // =====================================================================
    std::cout << "\n=== 4. S(q) quadrature against analytic transforms ===\n";
    {
        RadialDistribution grid(p.gr_dr, L);
        const std::vector<double> r = grid.radii();
        const std::vector<double> q = q_grid(p.sq_qmin, p.sq_qmax, p.sq_dq);
        std::cout << "  q grid: " << q.size() << " points, " << q.front() << " ... " << q.back() << " A^-1\n";
        check(q.size() == 1451 && std::abs(q.front() - 0.5) < 1e-12 && std::abs(q.back() - 15.0) < 1e-9,
              "q grid 0.5 ... 15 A^-1 in steps of 0.01");

        // g = 1 everywhere: no structure, S = 1
        const std::vector<double> flat(r.size(), 1.0);
        const std::vector<double> S_flat = structure_factor(r, flat, grid.dr, rho, q);
        double worst_flat = 0;
        for (double s : S_flat) worst_flat = std::max(worst_flat, std::abs(s - 1.0));
        check(worst_flat < 1e-14, "g = 1 gives S = 1 (the self term)");

        // step g(r) = theta(r - d): S(q) = 1 - 4 pi rho [sin(qd) - qd cos(qd)] / q^3
        const double d = 3.4;   // = 68 dr, a bin edge
        std::vector<double> step(r.size());
        for (std::size_t i = 0; i < r.size(); ++i) step[i] = r[i] < d ? 0.0 : 1.0;
        const std::vector<double> S_step = structure_factor(r, step, grid.dr, rho, q);
        double worst = 0;
        for (std::size_t k = 0; k < q.size(); ++k) {
            const double qd = q[k] * d;
            const double exact = 1.0 - 4.0 * PI * rho * (std::sin(qd) - qd * std::cos(qd)) / std::pow(q[k], 3);
            worst = std::max(worst, std::abs(S_step[k] - exact));
        }
        std::cout << "  step g(r): max |S_numeric - S_exact| = " << worst << "\n";
        check(worst < 5e-4, "midpoint quadrature matches the analytic transform of a step g(r)");
    }

    // =====================================================================
    std::cout << "\n=== 5. Output formats ===\n";
    {
        Params tiny = p;
        tiny.N = 4;   // one FCC cell
        std::mt19937_64 rng(1);
        System sys = make_system(tiny, rng);
        sys.step = 7;
        sys.time = 0.035;

        std::ostringstream out;
        write_xyz_frame(out, sys);
        sys.positions[0] = Vec3{1.25, 2.5, 3.75};
        sys.step = 8;
        sys.time = 0.04;
        write_xyz_frame(out, sys);
        const std::vector<std::string> lines = split_lines(out.str());

        bool no_blank = std::none_of(lines.begin(), lines.end(), [](const std::string& l) { return l.empty(); });
        check(lines.size() == 2 * (4 + 2) && no_blank, "two frames = exactly 2 (N + 2) lines, no blank lines");
        check(lines[0] == "4" && lines[6] == "4", "each frame starts with the atom count");
        check(lines[1].rfind("Lattice=\"", 0) == 0 &&
              lines[1].find("Properties=species:S:1:pos:R:3:velo:R:3") != std::string::npos &&
              lines[1].find("Step=7 Time=0.035") != std::string::npos &&
              lines[7].find("Step=8 Time=0.04") != std::string::npos,
              "comment line: Lattice, Properties, Step, Time");
        const std::vector<std::string> atom = tokens(lines[8]);
        check(atom.size() == 7 && atom[0] == "Ar" && std::stod(atom[1]) == 1.25 && std::stod(atom[2]) == 2.5 &&
              std::stod(atom[3]) == 3.75, "atom line: Ar x y z vx vy vz");

        // thermo line: nine numeric columns, nan before the first accumulated step
        EnergyDrift drift;
        RunningAverage T_avg;
        drift.start(-100.0);
        std::ostringstream log;
        write_thermo_line(log, sys, -100.0, drift, T_avg, 1.0);
        const std::vector<std::string> cols = tokens(log.str());
        check(cols.size() == 9 && cols[3] == "nan" && cols[8] == "nan", "thermo line: 9 columns, nan before first step");
    }

    // =====================================================================
    std::cout << "\n=== 6. Observable parameters parse ===\n";
    {
        const std::string path = "test_params_tmp.txt";
        {
            std::ofstream f(path);
            f << "gr_dr 0.1\nsq_qmin 0.25\nsq_qmax 12\nsq_dq 0.02\noutput_prefix run1_\n";
        }
        const Params q = read_params(path);
        std::remove(path.c_str());
        check(q.gr_dr == 0.1 && q.sq_qmin == 0.25 && q.sq_qmax == 12 && q.sq_dq == 0.02 &&
              output_path(q, "gr.dat") == "run1_gr.dat", "gr_dr, sq_qmin/qmax/dq and output_prefix are read");
        Params bad = p;
        bad.sq_qmax = 0.1;
        check_throws([&] { check_params(bad); }, "rejects sq_qmax <= sq_qmin");
        bad = p;
        bad.gr_dr = 20.0;
        check_throws([&] { check_params(bad); }, "rejects gr_dr >= L/2");
    }

    // =====================================================================
    std::cout << "\n=== 7. Liquid: S(q) from g(r) vs direct S(k) = <|sum_j exp(i k.r_j)|^2>/N ===\n";
    {
        std::mt19937_64 rng(p.seed);
        System sys = make_system(p, rng);
        const LennardJones lj(p);
        NoseHooverChain nhc(p, sys.n_dof);
        compute_forces(sys, lj);
        for (int k = 0; k < 3000; ++k) integrate_step(sys, lj, nhc, p.dt);   // 15 ps: melt

        // Wavevectors k = (2 pi / L) n, grouped into windows of |n|^2 so each group holds ~100+
        // independent vectors (only one of each +-k pair: |rho_-k| = |rho_k|). The direct estimate
        // is noisy, especially at the main peak where density fluctuations relax slowly, so the
        // tolerance is its measured statistical error (5 blocks) times 4, plus 0.05 for the
        // truncation of g(r) at L/2 (an absolute error of a few 0.01 in S).
        const int groups[5][2] = {{69, 71}, {116, 122}, {146, 154}, {278, 284}, {496, 502}};
        std::vector<std::vector<Vec3>> kvecs(5);
        for (int a = -23; a <= 23; ++a)
            for (int b = -23; b <= 23; ++b)
                for (int c = -23; c <= 23; ++c) {
                    const bool first_of_pair = a > 0 || (a == 0 && (b > 0 || (b == 0 && c > 0)));
                    const int n2 = a * a + b * b + c * c;
                    for (int s = 0; s < 5; ++s)
                        if (first_of_pair && n2 >= groups[s][0] && n2 <= groups[s][1] && kvecs[s].size() < 150)
                            kvecs[s].push_back((2.0 * PI / L) * Vec3{double(a), double(b), double(c)});
                }

        RadialDistribution rdf(p.gr_dr, L);
        const int n_frames = 200, n_blocks = 5;
        std::vector<std::vector<double>> block_sum(5, std::vector<double>(n_blocks, 0.0));
        for (int frame = 0; frame < n_frames; ++frame) {
            for (int k = 0; k < 20; ++k) integrate_step(sys, lj, nhc, p.dt);
            rdf.sample(sys);
            for (int s = 0; s < 5; ++s) {
                double acc = 0;
                for (const Vec3& k : kvecs[s]) {
                    double re = 0, im = 0;
                    for (const Vec3& r : sys.positions) {
                        const double phase = dot(k, r);
                        re += std::cos(phase);
                        im += std::sin(phase);
                    }
                    acc += (re * re + im * im) / N;
                }
                block_sum[s][frame * n_blocks / n_frames] += acc / kvecs[s].size();
            }
        }
        const std::vector<double> g = rdf.g(N);
        const std::vector<double> r_i = rdf.radii();
        bool all_agree = true;
        for (int s = 0; s < 5; ++s) {
            double mean = 0, var = 0;
            for (double b : block_sum[s]) mean += b / (n_frames / n_blocks) / n_blocks;
            for (double b : block_sum[s]) var += std::pow(b / (n_frames / n_blocks) - mean, 2) / (n_blocks - 1);
            const double se = std::sqrt(var / n_blocks);
            // S(q) from g(r) averaged over the same vectors' |k|
            std::vector<double> q_k;
            for (const Vec3& k : kvecs[s]) q_k.push_back(magnitude(k));
            const std::vector<double> S_k = structure_factor(r_i, g, rdf.dr, rho, q_k);
            double S_from_g = 0;
            for (double v : S_k) S_from_g += v / S_k.size();
            const double tol = 4.0 * se + 0.05;
            const bool ok = std::abs(S_from_g - mean) < tol;
            all_agree = all_agree && ok;
            std::cout << "  |k| = " << std::setprecision(4) << 2 * PI / L * std::sqrt(double(groups[s][0])) << " - "
                      << 2 * PI / L * std::sqrt(double(groups[s][1])) << " A^-1 (" << kvecs[s].size()
                      << " vectors): direct " << std::setprecision(4) << mean << " +- " << se
                      << ", from g(r) " << S_from_g << (ok ? "" : "   <-- outside tolerance") << "\n";
        }
        std::cout << std::setprecision(8);
        check(all_agree, "S(q) from g(r) agrees with the direct S(k) at 1.5 - 4 A^-1 (4 s.e. + 0.05)");

        double tail = 0;
        int n_tail = 0;
        for (int i = 0; i < rdf.n_bins; ++i)
            if (rdf.r(i) > rdf.r_max - 3.0) { tail += g[i]; ++n_tail; }
        tail /= n_tail;
        const int i_peak = static_cast<int>(std::max_element(g.begin(), g.end()) - g.begin());
        std::cout << "  g(r): highest peak at " << rdf.r(i_peak) << " A, g = " << g[i_peak]
                  << "; mean g over the last 3 A = " << tail << "\n";
        check(g[static_cast<int>(0.8 * p.sigma / p.gr_dr)] == 0.0, "g(r) = 0 inside the repulsive core (r < 0.8 sigma)");
        check(std::abs(tail - 1.0) < 0.02, "g(r) -> 1 at large r");
    }

    std::cout << "\n" << n_pass << " passed, " << n_fail << " failed\n";
    return n_fail == 0 ? 0 : 1;
}