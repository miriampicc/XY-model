//
// Created by Mirimi on 11/04/24.
//

#ifndef CMT_PARALLEL_TEMP_H
#define CMT_PARALLEL_TEMP_H

#include "parallel_temp.h"
#include <mpi.h>
#include <fstream>
#include "initialization.h"

struct PT_parameters{
    /*Parallel Tempering parameters*/
    int np;
    int rank;
    int root=0;
};

struct PTroot_parameters{
    /*Arrays root needs to handle the swaps*/
    std::vector <double> beta;
    std::vector <double> All_Energies;
    std::vector <int> ind_to_rank;
    std::vector <int> rank_to_ind;
};

void initialize_PTarrays(struct PT_parameters &PTp, struct PTroot_parameters &PTroot, struct H_parameters &Hp);
void parallel_temp(double &my_E , double &my_beta,  int &my_ind, struct PT_parameters &PTp, struct PTroot_parameters &PTroot);

#endif //CMT_PARALLEL_TEMP_H
