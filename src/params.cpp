#include "../include/params.h"
#include "../include/units.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>


// ---------------------------------------------------------------- derived quantities

double box_length(const Params& p) { return std::cbrt(p.N / p.density); }
double cutoff_radius(const Params& p) { return p.rc_sigma * p.sigma; }
double healing_length(const Params& p) { return p.lambda_sigma * p.sigma; }
double epsilon_internal(const Params& p) { return units::KB * p.epsilon_K; }
double nhc_tau(const Params& p) { return p.nhc_tau_steps * p.dt; }


// ---------------------------------------------------------------- reading

// read parameters from input.txt: one "key value" pair per line, '#' starts a comment
Params read_params(const std::string& filename) {

    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open parameter file: " + filename);

    Params p;
    std::string line;

    while (std::getline(file, line)) {

        // delete comments ('#')
        std::size_t pound_sign = line.find('#');
        if (pound_sign != std::string::npos) line.erase(pound_sign);

        std::istringstream param(line);
        std::string key;
        if (!(param >> key)) continue;   // blank line

        if      (key == "N")             param >> p.N;
        else if (key == "density")       param >> p.density;
        else if (key == "mass")          param >> p.mass;
        else if (key == "sigma")         param >> p.sigma;
        else if (key == "epsilon_K")     param >> p.epsilon_K;
        else if (key == "rc_sigma")      param >> p.rc_sigma;
        else if (key == "lambda_sigma")  param >> p.lambda_sigma;
        else if (key == "T")             param >> p.T;
        else if (key == "thermostat")    param >> p.thermostat;
        else if (key == "nhc_M")         param >> p.nhc_M;
        else if (key == "nhc_tau_steps") param >> p.nhc_tau_steps;
        else if (key == "nhc_nc")        param >> p.nhc_nc;
        else if (key == "nhc_nsy")       param >> p.nhc_nsy;
        else if (key == "dt")            param >> p.dt;
        else if (key == "n_equil")       param >> p.n_equil;
        else if (key == "n_prod")        param >> p.n_prod;
        else if (key == "seed")          param >> p.seed;
        else if (key == "thermo_every")  param >> p.thermo_every;
        else if (key == "traj_every")    param >> p.traj_every;
        else if (key == "sample_every")  param >> p.sample_every;
        else if (key == "output_prefix") param >> p.output_prefix;
        else if (key == "gr_dr")         param >> p.gr_dr;
        else if (key == "sq_qmin")       param >> p.sq_qmin;
        else if (key == "sq_qmax")       param >> p.sq_qmax;
        else if (key == "sq_dq")         param >> p.sq_dq;
        else if (key == "L")
            throw std::runtime_error("'L' is not an input: it is derived as (N / density)^(1/3)");
        else throw std::runtime_error("Unknown parameter: " + key);

        // value must parse as the type declared in Params, with nothing left over
        if (param.fail())
            throw std::runtime_error("Failed to parse parameter (bad or missing value): " + key);
        std::string extra;
        if (param >> extra)
            throw std::runtime_error("Unexpected extra value '" + extra + "' for parameter: " + key);
    }

    return p;
}


// ---------------------------------------------------------------- validation

static void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error("Invalid parameters: " + message);
}

void check_params(const Params& p) {
    // system
    require(p.N > 0, "N must be positive");
    const int n_cells = static_cast<int>(std::lround(std::cbrt(p.N / 4.0)));
    require(4 * n_cells * n_cells * n_cells == p.N, "N must be 4 n^3 (108, 256, 500, 864, 1372, ...)");
    require(p.density > 0, "density must be positive");
    require(p.mass > 0, "mass must be positive");

    // potential
    require(p.sigma > 0, "sigma must be positive");
    require(p.epsilon_K > 0, "epsilon_K must be positive");
    require(p.rc_sigma > 0, "rc_sigma must be positive");
    require(p.lambda_sigma > 0, "lambda_sigma must be positive");

    const double L = box_length(p), rc = cutoff_radius(p), lambda = healing_length(p);
    std::ostringstream msg;

    msg << "rc = " << rc << " A must be < L/2 = " << 0.5 * L << " A (minimum image convention)";
    require(rc < 0.5 * L, msg.str());

    msg.str("");
    msg << "switching must start beyond the LJ minimum: rc - lambda = " << rc - lambda
        << " A must be > 2^(1/6) sigma = " << std::pow(2.0, 1.0 / 6.0) * p.sigma << " A";
    require(rc - lambda > std::pow(2.0, 1.0 / 6.0) * p.sigma, msg.str());

    // ensemble / thermostat
    require(p.T > 0, "T must be positive");
    require(p.thermostat == "nhc" || p.thermostat == "none", "thermostat must be 'nhc' or 'none'");
    require(p.nhc_M >= 1, "nhc_M must be >= 1");
    require(p.nhc_tau_steps >= 20, "nhc_tau_steps must be >= 20 (tau >= 20 dt)");
    require(p.nhc_nc >= 1, "nhc_nc must be >= 1");
    require(p.nhc_nsy == 1 || p.nhc_nsy == 3 || p.nhc_nsy == 7, "nhc_nsy must be 1, 3 or 7");

    // integration
    require(p.dt > 0, "dt must be positive");
    require(p.n_equil >= 0, "n_equil must be >= 0");
    require(p.n_prod >= 0, "n_prod must be >= 0");

    // switching shell must be resolved: a pair at the typical relative speed sqrt(2 kB T / m)
    // takes at least 10 time steps to cross it
    const double lambda_min = 10.0 * std::sqrt(2.0 * units::KB * p.T / p.mass) * p.dt;
    msg.str("");
    msg << "healing length lambda = " << lambda << " A is too small for dt = " << p.dt
        << " ps; need lambda >= " << lambda_min << " A (" << lambda_min / p.sigma << " sigma)";
    require(lambda >= lambda_min, msg.str());

    // output
    require(p.thermo_every > 0, "thermo_every must be positive");
    require(p.traj_every > 0, "traj_every must be positive");
    require(p.sample_every > 0, "sample_every must be positive");

    // observables
    require(p.gr_dr > 0, "gr_dr must be positive");
    require(p.gr_dr < 0.5 * L, "gr_dr must be smaller than L/2");
    require(p.sq_qmin > 0, "sq_qmin must be positive");
    require(p.sq_qmax > p.sq_qmin, "sq_qmax must be larger than sq_qmin");
    require(p.sq_dq > 0, "sq_dq must be positive");
}


// ---------------------------------------------------------------- printing

void print_params(const Params& p, std::ostream& os) {
    os << "# ---- parameters\n"
       << "#   N = " << p.N << ", density = " << p.density << " A^-3, mass = " << p.mass << " amu\n"
       << "#   L = " << box_length(p) << " A (derived)\n"
       << "#   sigma = " << p.sigma << " A, epsilon/kB = " << p.epsilon_K << " K ("
       << epsilon_internal(p) << " amu A^2/ps^2)\n"
       << "#   rc = " << p.rc_sigma << " sigma = " << cutoff_radius(p) << " A, lambda = "
       << p.lambda_sigma << " sigma = " << healing_length(p) << " A\n"
       << "#   T = " << p.T << " K, thermostat = " << p.thermostat << "\n";
    if (p.thermostat == "nhc")
        os << "#   NHC: M = " << p.nhc_M << ", tau = " << p.nhc_tau_steps << " dt = " << nhc_tau(p)
           << " ps, nc = " << p.nhc_nc << ", nsy = " << p.nhc_nsy << "\n";
    os << "#   dt = " << p.dt << " ps, n_equil = " << p.n_equil << " (" << p.n_equil * p.dt
       << " ps), n_prod = " << p.n_prod << " (" << p.n_prod * p.dt << " ps)\n"
       << "#   seed = " << p.seed << ", thermo_every = " << p.thermo_every
       << ", traj_every = " << p.traj_every << ", sample_every = " << p.sample_every << "\n"
       << "#   g(r): dr = " << p.gr_dr << " A, r_max = L/2 = " << 0.5 * box_length(p) << " A; S(q): q = "
       << p.sq_qmin << " ... " << p.sq_qmax << " A^-1, dq = " << p.sq_dq << " A^-1\n";
    if (!p.output_prefix.empty()) os << "#   output_prefix = " << p.output_prefix << "\n";
}