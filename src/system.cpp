#include "../include/system.h"
#include "../include/init.h"
#include "../include/observables.h"


System make_system(const Params& p, std::mt19937_64& rng) {
    System sys;
    sys.N = p.N;
    sys.mass = p.mass;
    sys.L = box_length(p);
    sys.n_dof = degrees_of_freedom(p.N);

    sys.positions = init_fcc_positions(p.N, sys.L, p.density);
    sys.velocities = init_velocities(sys.positions, p.mass, p.T, rng);
    sys.forces.assign(p.N, Vec3{});

    sys.potential_energy = 0.0;
    sys.virial = 0.0;
    sys.step = 0;
    sys.time = 0.0;
    return sys;
}


void wrap_positions(System& sys) {
    const double L = sys.L;
    for (Vec3& r : sys.positions) {
        r.x = wrap_coordinate(r.x, L);
        r.y = wrap_coordinate(r.y, L);
        r.z = wrap_coordinate(r.z, L);
    }
}