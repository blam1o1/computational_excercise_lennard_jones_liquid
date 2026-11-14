#include "../include/output.h"
#include <cstdio>
#include <iomanip>
#include <stdexcept>


std::string output_path(const Params& p, const std::string& name) { return p.output_prefix + name; }


std::ofstream open_output(const std::string& path) {
    std::ofstream file(path);
    if (!file) throw std::runtime_error("Cannot open output file: " + path);
    return file;
}


// ---------------------------------------------------------------- trajectory

void write_xyz_frame(std::ostream& os, const System& sys) {
    const std::streamsize old_precision = os.precision(12);
    const double L = sys.L;
    os << sys.N << '\n'
       << "Lattice=\"" << L << " 0 0 0 " << L << " 0 0 0 " << L
       << "\" Properties=species:S:1:pos:R:3:velo:R:3"
       << " Step=" << sys.step << " Time=" << sys.time << '\n';
    for (int i = 0; i < sys.N; ++i) {
        const Vec3& r = sys.positions[i];
        const Vec3& v = sys.velocities[i];
        os << "Ar " << r.x << ' ' << r.y << ' ' << r.z << ' ' << v.x << ' ' << v.y << ' ' << v.z << '\n';
    }
    os.precision(old_precision);
}


// ---------------------------------------------------------------- thermodynamic log

void write_thermo_header(std::ostream& os, const Params& p, const std::string& phase, double E0) {
    const double eps = epsilon_internal(p);
    os << "# " << phase << " phase: N = " << p.N << ", L = " << box_length(p) << " A, T_target = " << p.T
       << " K, dt = " << p.dt << " ps, line every " << p.thermo_every << " steps\n"
       << "# energies per particle in units of epsilon (epsilon/kB = " << p.epsilon_K << " K)\n"
       << "# H' = K + U + thermostat energy; E0 = H'(t = 0) = " << std::setprecision(10)
       << E0 / (p.N * eps) << " eps per particle\n"
       << "# DeltaE = (1/n) sum_{k=1..n} |(H'_k - E0)/E0|, accumulated every step since t = 0 (whole run)\n"
       << "# <T_kin> = running mean of T_kin = 2K/((3N-3) kB) over every step of this phase\n"
       << "#      step   time[ps]   T_kin[K] <T_kin>[K]    K/N[eps]    U/N[eps]    H'/N[eps]"
          "  (H'-E0)/|E0|        DeltaE\n";
}


void write_thermo_line(std::ostream& os, const System& sys, double H, const EnergyDrift& drift,
                       const RunningAverage& T_avg, double epsilon) {
    const double per_particle = 1.0 / (sys.N * epsilon);
    char line[256];
    std::snprintf(line, sizeof(line), "%11ld %10.4f %10.4f %10.4f %11.6f %11.6f %12.7f %13.5e %13.5e\n",
                  sys.step, sys.time, temperature(sys.velocities, sys.mass), T_avg.mean(),
                  kinetic_energy(sys.velocities, sys.mass) * per_particle,
                  sys.potential_energy * per_particle, H * per_particle, drift.deviation(H), drift.delta_E());
    os << line;
}


// ---------------------------------------------------------------- structure

void write_gr(const std::string& path, const Params& p, const RadialDistribution& rdf,
              const std::vector<double>& g) {
    std::ofstream file = open_output(path);
    file << "# radial distribution function, production phase\n"
         << "# N = " << p.N << ", L = " << rdf.L << " A, rho = " << p.N / (rdf.L * rdf.L * rdf.L)
         << " A^-3, dr = " << rdf.dr << " A, r_max = L/2 = " << rdf.r_max << " A, N_conf = " << rdf.n_conf
         << " (every " << p.sample_every << " steps)\n"
         << "# g(r_i) = h(i) / (4 pi rho r_i^2 dr N_conf N),  r_i = (i + 1/2) dr\n"
         << "# n(r) = running coordination number: neighbours within r_i + dr/2\n"
         << "#     r[A]          g(r)          n(r)\n";
    double cumulative = 0.0;
    char line[128];
    for (int i = 0; i < rdf.n_bins; ++i) {
        cumulative += static_cast<double>(rdf.histogram[i]) / (static_cast<double>(rdf.n_conf) * p.N);
        std::snprintf(line, sizeof(line), "%10.4f %13.6f %13.6f\n", rdf.r(i), g[i], cumulative);
        file << line;
    }
}


void write_sq(const std::string& path, const Params& p, const RadialDistribution& rdf,
              const std::vector<double>& q, const std::vector<double>& S) {
    std::ofstream file = open_output(path);
    file << "# static structure factor from g(r), production phase\n"
         << "# S(q) = 1 + 4 pi rho dr sum_i r_i^2 [g(r_i) - 1] sin(q r_i)/(q r_i),  r_i < r_max = L/2 = "
         << rdf.r_max << " A\n"
         << "# rho = " << p.N / (rdf.L * rdf.L * rdf.L) << " A^-3, N_conf = " << rdf.n_conf << "\n"
         << "# note: g(r) is known only to L/2, so values below ~1 A^-1 carry errors of a few 0.01\n"
         << "#   q[A^-1]          S(q)\n";
    char line[128];
    for (std::size_t k = 0; k < q.size(); ++k) {
        std::snprintf(line, sizeof(line), "%10.4f %13.6f\n", q[k], S[k]);
        file << line;
    }
}