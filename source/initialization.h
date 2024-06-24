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
    double e;
    double beta_high;
    double beta_low;
    double b1;
    double b2;
};

struct MC_parameters {
    int n_steps; //Number of Monte Carlo steps we want to perform
    int transient; // estimate of the thermalization time
    double tau;
    double theta_box;
    double theta_box_A;
    double theta_box_density;
};


struct Node {
    //complex order parameter of the NC and SC components
   std::array <O2, 2> Psi {};
   //Fluctuating vector potential A, with 2 spatial dimensions
   std::array <double, 2> A{};
};

void initialize_Lattice (std::vector <Node> &Site, const fs::path & directory_param_beta, int restart, struct H_parameters &Hp);