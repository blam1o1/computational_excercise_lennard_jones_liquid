// header file for physical constants and the internal unit system
#pragma once

// Internal unit system (chosen so the user-facing inputs need no conversion):
//   length       : Angstrom (A)
//   mass         : atomic mass unit (amu)
//   time         : picosecond (ps)
//   temperature  : Kelvin (K)
//   velocity     : A / ps          (1 A/ps = 100 m/s)
//   energy       : amu A^2 / ps^2  (1 amu A^2/ps^2 = 1.66053906660e-23 J)
//
// Lennard-Jones epsilon given as epsilon/kB in K converts to internal energy
// via  epsilon_internal = KB * epsilon_K.

namespace units {

// Boltzmann constant in amu * A^2 * ps^-2 * K^-1
//   = 1.380649e-23 J/K / 1.66053906660e-27 kg * (1e20 A^2/m^2) * (1e-24 s^2/ps^2)
constexpr double KB = 0.8314462621026538;

constexpr double PI = 3.14159265358979323846;

}  // namespace units