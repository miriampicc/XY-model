#pragma once

#include<cstdio>
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
#include "parallel_temp.h"
#include "montecarlo.h"
#include "o2.h"
#include "initialization.h"
#include "robust_filesystem.h"
#include <mpi.h>



namespace paths_dir{
    inline std::string DIR_IN;
    inline std::string DIR_OUT;
}


// Define simulation parameters

extern size_t L;            // Number of spins in each dimension
extern double T;        // Temperature
extern bool restart;    // where does the program have to take the initial file

void myhelp(int argd, char** argu);
void mainloop(std::vector <Node> &Site, double T_temp, struct MC_parameters &MC, int &my_ind, double &my_beta, struct PT_parameters PTp, struct PTroot_parameters PTrood, size_t N, struct H_parameters &Hp, std::string directory_write, int NSTART);
