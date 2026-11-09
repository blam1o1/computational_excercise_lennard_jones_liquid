#include "../include/params.h"
#include <fstream>
#include <sstream>
#include <stdexcept>


// read parameters from input.txt
Params read_params(const std::string& filename) {

    // open input.txt
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open parameter file:" + filename);

    Params p;
    std::string line;

    while (std::getline(file, line)) {

        // delete comments ('#') as we read input.txt
        std::size_t pound_sign = line.find('#');
        if (pound_sign != std::string::npos) line.erase(pound_sign);

        // parse input.txt and update parameter values
        std::istringstream param(line);
        std::string key;

        if (!(param >> key)) continue;
        if (key == "N") param >> p.N;
        else if (key == "L") param >> p.L;
        else if (key == "rc") param >> p.rc;
        else if (key == "T") param >> p.T;
        else if (key == "time_step") param >> p.time_step;
        else if (key == "n_steps") param >> p.n_steps;
        else if (key == "thermostat") param >> p.thermostat;
        else if (key == "seed") param >> p.seed;
        else if (key == "thermo_every") param >> p.thermo_every;
        else if (key == "traj_every") param >> p.traj_every;
        else throw std::runtime_error("Unknown parameter: " + key);


        // check if parameter parsing succeeded according to type defined in Params struct
        if (param.fail()) {
            throw std::runtime_error("Failed to parse parameter (bad or missing value): " + key);
        }
    }

    return p;
}

// check if parameters are valid, like we cant have negative particles etc...
void check_params(const Params& p) {
    if (p.N <= 0) throw std::runtime_error("Number of particles must be positive");
    if (p.L <= 0) throw std::runtime_error("Box length must be positive");
    if (p.rc <= 0) throw std::runtime_error("'rc' must be positive");
    if (p.T <= 0) throw std::runtime_error("Temperature must be positive");
    if (p.time_step <= 0) throw std::runtime_error("'time_step' must be positive");
    if (p.n_steps <= 0) throw std::runtime_error("'n_steps' must be positive");
    if (p.thermostat != "nhc" && p.thermostat != "nhc") throw std::runtime_error("Invalid thermostat");
    if (p.thermo_every <= 0) throw std::runtime_error("'thermo_every' must be positive");
    if (p.traj_every <= 0) throw std::runtime_error("'traj_every' must be positive");    
}

// print function 
void print_params(const Params& params, std::ostream& os) {
    os << "Number of particles: " << params.N << std::endl;
    os << "Box length: " << params.L << std::endl;
    os << "Cutoff radius: " << params.rc << std::endl;
    os << "Temperature: " << params.T << std::endl;
    os << "Time step: " << params.time_step << std::endl;
    os << "Number of steps: " << params.n_steps << std::endl;
    os << "Thermostat: " << params.thermostat << std::endl;
    os << "Thermo every: " << params.thermo_every << std::endl;
    os << "Traj every: " << params.traj_every << std::endl;
}