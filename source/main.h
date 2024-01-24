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
#include "initialization.h"

namespace paths_dir{
    inline std::string DIR_IN;
    inline std::string DIR_OUT;
}


// Define simulation parameters

extern size_t L;            // Number of spins in each dimension
extern size_t n_steps;   // Number of simulation steps
extern double T;        // Temperature
extern bool restart;    // where does the program have to take the initial file

void myhelp(int argd, char** argu);
void mainloop (std::vector <Node> &Site, double T, struct MC_parameters &MC, size_t N, struct H_parameters &Hp, std::string directory_write, int NSTART);
