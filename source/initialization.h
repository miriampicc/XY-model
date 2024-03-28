//
// Created by Mirimi on 16/01/24.
//
#pragma once

#include "robust_filesystem.h"
#include "rng.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include<array>
#include "o2.h"


struct H_parameters {
    double K ;
    double J1;
    double J2;
    double e;
};

struct MC_parameters {
    int n_steps; //Number of Monte Carlo steps we want to perform
    int transient; // estimate of the thermalization time
    double tau;
    long double theta_box;
};


struct Node {
    //complex order parameter of the NC and SC components
   std::array <O2, 2> Psi {};

};

void initialize_Lattice (std::vector <Node> &Site, const fs::path & directory_read, const fs::path & directory_write, int restart, struct H_parameters &Hp);