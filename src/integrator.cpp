#include "../include/integrator.h"
#include "../include/observables.h"


// steps 2-5: exp(iL_2 dt/2) exp(iL_1 dt) exp(iL_2 dt/2)
static void velocity_verlet(System& sys, const LennardJones& lj, double dt) {
    const double half_dt_over_m = 0.5 * dt / sys.mass;
    const int N = sys.N;

    // 2. half kick with F(t)
    for (int i = 0; i < N; ++i) sys.velocities[i] += half_dt_over_m * sys.forces[i];

    // 3. drift, then periodic boundaries
    for (int i = 0; i < N; ++i) sys.positions[i] += dt * sys.velocities[i];
    wrap_positions(sys);

    // 4. forces at the new positions
    compute_forces(sys, lj);

    // 5. half kick with F(t + dt)
    for (int i = 0; i < N; ++i) sys.velocities[i] += half_dt_over_m * sys.forces[i];
}


void integrate_step(System& sys, const LennardJones& lj, NoseHooverChain& nhc, double dt) {
    nhc.half_step(sys, dt);          // 1.
    velocity_verlet(sys, lj, dt);    // 2-5.
    nhc.half_step(sys, dt);          // 6.
    sys.step += 1;
    sys.time += dt;
}


void integrate_step(System& sys, const LennardJones& lj, double dt) {
    velocity_verlet(sys, lj, dt);
    sys.step += 1;
    sys.time += dt;
}


double conserved_energy(const System& sys, const NoseHooverChain& nhc) {
    return conserved_energy(sys) + nhc.energy();
}

double conserved_energy(const System& sys) {
    return kinetic_energy(sys.velocities, sys.mass) + sys.potential_energy;
}