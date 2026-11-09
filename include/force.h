// header file for force.cpp

// abstract of lennard jones force potential and force calculation

#pragma once
#include <cmath>

// Lennard Jones force potential and force calculation
struct LennardJones {
    double sigma = 1.0;
    double epsilon = 1.0;
    double rc = 2.5;
};

// Lennard Jones force calculation