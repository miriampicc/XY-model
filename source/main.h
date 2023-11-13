#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include "rng.h"
#include<cstdio>

namespace paths_dir{
    inline std::string DIR_IN;
    inline std::string DIR_OUT;
}


// Define simulation parameters

extern size_t L;            // Number of spins in each dimension
extern size_t n_steps;   // Number of simulation steps
extern double T;        // Temperature
extern bool restart;    // where does the program have to take the initial file

void mc_step(std::vector<double> & spins, double T, int &acc, double thetabox);
double local_energy(const std::vector<double> & spins);
void myhelp(int argd, char** argu);
void mainloop (std::vector<double> & spins, double T, int n_steps, size_t N, std::string directory_write);
