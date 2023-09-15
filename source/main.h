#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <cstdlib>
#include <ctime>
#include<fstream>

// Define simulation parameters

extern size_t L;            // Number of spins in each dimension
extern size_t n_steps;   // Number of simulation steps
extern double T;        // Temperature
extern bool restart;    // where does the program have to take the initial file

void mc_step(std::vector<double> & spins, double T, int &acc, double thetabox);
double local_energy(const std::vector<double> & spins);
int vortex (std::vector<double> spins);
double J_param (std::vector<double> cc, double T);


