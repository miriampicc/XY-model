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
void mainloop(std::vector <Node> &Site, struct MC_parameters &MC, int &my_ind, double &my_beta, struct PT_parameters PTp, struct PTroot_parameters PTrood, size_t N, struct H_parameters &Hp, std::string directory_write, int NSTART);
void update_file_path (const std::string& base_dir, std::ofstream& file_energy, std::ofstream& file_mag1, std::ofstream& file_kin, std::ofstream& file_josph, std::ofstream& file_B, std::ofstream& file_hel1, std::ofstream& file_hel2, std::ofstream& file_trsb, std::ofstream& file_ds, std::ofstream& file_rank);