// Molecular dynamics of a Lennard-Jones liquid (argon) in the canonical ensemble.
// Usage: ./bin/md [input file]      (default: input.txt)
//
// Output (each name prefixed by output_prefix, if set):
//   trajectory.xyz     frames at t = 0 and every traj_every steps (equilibration + production)
//   thermo_equil.dat   T_kin, <T_kin>, K, U, H', Delta E every thermo_every steps (equilibration)
//   thermo_prod.dat    the same for production (<T_kin> restarts; Delta E runs on from t = 0)
//   gr.dat, sq.dat     g(r) and S(q) from production frames sampled every sample_every steps

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include "../include/force.h"
#include "../include/integrator.h"
#include "../include/observables.h"
#include "../include/output.h"
#include "../include/params.h"
#include "../include/system.h"
#include "../include/thermostat.h"


int main(int argc, char* argv[]) {
    try {
        // ---- parameters
        const std::string input_file = (argc > 1) ? argv[1] : "input.txt";
        const Params p = read_params(input_file);
        check_params(p);
        print_params(p, std::cout);

        // ---- initial state
        std::mt19937_64 rng(p.seed);
        System sys = make_system(p, rng);
        const LennardJones lj(p);
        NoseHooverChain nhc(p, sys.n_dof);
        const bool nvt = (p.thermostat == "nhc");

        compute_forces(sys, lj);   // F(0): the first velocity-Verlet kick needs it

        auto advance = [&]() {
            if (nvt) integrate_step(sys, lj, nhc, p.dt);
            else integrate_step(sys, lj, p.dt);
        };
        auto conserved = [&]() { return nvt ? conserved_energy(sys, nhc) : conserved_energy(sys); };

        const double eps = epsilon_internal(p);
        const double rho = sys.N / (sys.L * sys.L * sys.L);

        // ---- observables
        EnergyDrift drift;                          // Delta E over the whole run
        RunningAverage T_equil, T_prod;             // <T_kin> for each phase
        RadialDistribution rdf(p.gr_dr, sys.L);     // g(r), production only

        std::ofstream traj = open_output(output_path(p, "trajectory.xyz"));
        std::ofstream thermo_equil = open_output(output_path(p, "thermo_equil.dat"));
        std::ofstream thermo_prod = open_output(output_path(p, "thermo_prod.dat"));
        long n_frames = 0;

        // t = 0: reference energy, first log line and first trajectory frame
        double H = conserved();
        drift.start(H);
        write_thermo_header(thermo_equil, p, "equilibration", H);
        write_thermo_line(thermo_equil, sys, H, drift, T_equil, eps);
        write_xyz_frame(traj, sys);
        ++n_frames;

        // bookkeeping after every step, shared by both phases
        auto on_step = [&](RunningAverage& T_avg, std::ofstream& thermo) {
            H = conserved();
            drift.add(H);
            T_avg.add(temperature(sys.velocities, sys.mass));
            if (sys.step % p.traj_every == 0) { write_xyz_frame(traj, sys); ++n_frames; }
            if (sys.step % p.thermo_every == 0) write_thermo_line(thermo, sys, H, drift, T_avg, eps);
        };

        // progress on screen at every 10% of a phase
        auto progress = [&](const char* phase, long s, long n_phase, const RunningAverage& T_avg) {
            const long every = std::max(1L, n_phase / 10);
            if (s % every == 0 || s == n_phase)
                std::printf("  %-6s step %9ld  t = %9.3f ps  T_kin = %8.3f K  <T_kin> = %8.3f K  DeltaE = %.3e\n",
                            phase, sys.step, sys.time, temperature(sys.velocities, sys.mass), T_avg.mean(),
                            drift.delta_E());
        };

        std::printf("# ---- run\n");
        const auto t_start = std::chrono::steady_clock::now();

        // ---- equilibration
        for (long s = 1; s <= p.n_equil; ++s) {
            advance();
            on_step(T_equil, thermo_equil);
            progress("equil", s, p.n_equil, T_equil);
        }

        // ---- production
        write_thermo_header(thermo_prod, p, "production", drift.E0);
        for (long s = 1; s <= p.n_prod; ++s) {
            advance();
            on_step(T_prod, thermo_prod);
            if (s % p.sample_every == 0) rdf.sample(sys);
            progress("prod", s, p.n_prod, T_prod);
        }

        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();

        // ---- structural results
        int i_peak = -1;
        std::size_t k_peak = 0;
        std::vector<double> g, q, S;
        if (rdf.n_conf > 0) {
            g = rdf.g(sys.N);
            write_gr(output_path(p, "gr.dat"), p, rdf, g);
            q = q_grid(p.sq_qmin, p.sq_qmax, p.sq_dq);
            S = structure_factor(rdf.radii(), g, rdf.dr, rho, q);
            write_sq(output_path(p, "sq.dat"), p, rdf, q, S);
            i_peak = static_cast<int>(std::max_element(g.begin(), g.end()) - g.begin());
            k_peak = static_cast<std::size_t>(std::max_element(S.begin(), S.end()) - S.begin());
        }

        // ---- summary
        const long total_steps = p.n_equil + p.n_prod;
        std::printf("# ---- summary\n");
        std::printf("#   %ld steps in %.2f s (%.3f ms/step)\n", total_steps, seconds,
                    total_steps > 0 ? 1e3 * seconds / total_steps : 0.0);
        std::printf("#   Delta E over the whole run (%ld steps): %.4e\n", drift.n, drift.delta_E());
        std::printf("#   equilibration <T_kin> = %.3f K (std %.3f K)\n", T_equil.mean(), T_equil.std_dev());
        std::printf("#   production    <T_kin> = %.3f K (std %.3f K; canonical T sqrt(2/N_f) = %.3f K)\n",
                    T_prod.mean(), T_prod.std_dev(), p.T * std::sqrt(2.0 / sys.n_dof));
        if (rdf.n_conf > 0) {
            std::printf("#   g(r): N_conf = %ld, highest peak at r = %.3f A (%.3f sigma), g = %.3f\n", rdf.n_conf,
                        rdf.r(i_peak), rdf.r(i_peak) / p.sigma, g[i_peak]);
            std::printf("#   S(q): highest peak at q = %.3f A^-1, S = %.3f\n", q[k_peak], S[k_peak]);
        } else {
            std::printf("#   g(r), S(q): no production frames sampled, gr.dat and sq.dat not written\n");
        }
        std::string files = output_path(p, "trajectory.xyz") + " (" + std::to_string(n_frames) + " frames), " +
                            output_path(p, "thermo_equil.dat") + ", " + output_path(p, "thermo_prod.dat");
        if (rdf.n_conf > 0) files += ", " + output_path(p, "gr.dat") + ", " + output_path(p, "sq.dat");
        std::printf("#   wrote %s\n", files.c_str());

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}