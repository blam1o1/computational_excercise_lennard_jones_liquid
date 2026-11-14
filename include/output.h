// header file for output.cpp

// writing simulation output: trajectory frames, thermodynamic log lines, g(r) and S(q) files

#pragma once
#include <fstream>
#include <ostream>
#include <string>
#include <vector>
#include "observables.h"
#include "params.h"
#include "system.h"

// p.output_prefix + name
std::string output_path(const Params& p, const std::string& name);

// open a file for writing (truncating any old one); throws if it cannot be opened
std::ofstream open_output(const std::string& path);

// One extended-XYZ frame:
//   line 1      atom count
//   line 2      Lattice="L 0 0 0 0 L 0 0 0 L" Properties=species:S:1:pos:R:3:velo:R:3 Step=<n> Time=<ps>
//   N lines     Ar x y z vx vy vz      (A, A/ps)
// Frames written back to back form a valid multi-frame .xyz file (no blank lines between frames).
void write_xyz_frame(std::ostream& os, const System& sys);

// Thermodynamic log, one line per call:
//   step  time[ps]  T_kin[K]  <T_kin>[K]  K/N  U/N  H'/N  (H'-E0)/|E0|  DeltaE
// energies per particle in units of epsilon. Cumulative columns print "nan" before their first sample.
void write_thermo_header(std::ostream& os, const Params& p, const std::string& phase, double E0);
void write_thermo_line(std::ostream& os, const System& sys, double H, const EnergyDrift& drift,
                       const RunningAverage& T_avg, double epsilon);

// Final structural results
//   gr.dat:  r[A]  g(r)  n(r) = number of neighbours within r + dr/2
//   sq.dat:  q[A^-1]  S(q)
void write_gr(const std::string& path, const Params& p, const RadialDistribution& rdf,
              const std::vector<double>& g);
void write_sq(const std::string& path, const Params& p, const RadialDistribution& rdf,
              const std::vector<double>& q, const std::vector<double>& S);