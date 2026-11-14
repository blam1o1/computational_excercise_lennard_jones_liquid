// Molecular dynamics of a Lennard-Jones liquid (argon) in the canonical ensemble.
// Usage: ./bin/md [input file]      (default: input.txt)

#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include "../include/force.h"
#include "../include/integrator.h"
#include "../include/observables.h"
#include "../include/params.h"
#include "../include/system.h"
#include "../include/thermostat.h"


// Minimal progress log until the observables module is designed:
// energies per particle in units of epsilon, temperature in K,
// relative change of the conserved quantity since the start of the current phase.
static void log_thermo(const char* phase, const System& sys, double E_conserved, double E0, double eps) {
    const double K = kinetic_energy(sys.velocities, sys.mass);
    std::printf("%-6s %9ld %10.3f %9.3f %11.5f %11.5f %12.6f %12.3e\n", phase, sys.step, sys.time,
                temperature(sys.velocities, sys.mass), K / (sys.N * eps),
                sys.potential_energy / (sys.N * eps), E_conserved / (sys.N * eps),
                (E_conserved - E0) / std::abs(E0));
}


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
        std::printf("# %-4s %9s %10s %9s %11s %11s %12s %12s\n", "phase", "step", "time[ps]", "T[K]",
                    "K/(N eps)", "U/(N eps)", "H'/(N eps)", "dH'/|H'0|");

        const auto t_start = std::chrono::steady_clock::now();

        // ---- equilibration
        double E0 = conserved();
        log_thermo("equil", sys, E0, E0, eps);
        for (long s = 1; s <= p.n_equil; ++s) {
            advance();
            if (s % p.thermo_every == 0) log_thermo("equil", sys, conserved(), E0, eps);
        }

        // ---- production
        E0 = conserved();
        // obs.begin_production(sys, nhc);              <- observables hook (later)
        for (long s = 1; s <= p.n_prod; ++s) {
            advance();
            // obs.on_step(sys, nhc);                   <- H', delta E, K, <K> (later)
            // if (s % p.sample_every == 0) obs.sample_structure(sys);   <- g(r), S(q) (later)
            // if (s % p.traj_every == 0) write_xyz(sys);                <- trajectory (later)
            if (s % p.thermo_every == 0) log_thermo("prod", sys, conserved(), E0, eps);
        }

        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
        const long total_steps = p.n_equil + p.n_prod;
        std::printf("# %ld steps in %.2f s (%.3f ms/step)\n", total_steps, seconds,
                    total_steps > 0 ? 1e3 * seconds / total_steps : 0.0);

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}