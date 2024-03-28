//
// Created by Mirimi on 17/01/24.
//

#pragma once

#include "main.h"
#include "initialization.h"
#include "rng.h"
#include <iostream>
#include <cstring>
#include "o2.h"
#include "measures.h"

void metropolis(std::vector<Node> &Site, struct MC_parameters &MC, struct H_parameters &Hp, double T);
double local_energy (std::array<O2, 2> &Psi, int i, H_parameters &Hp, const std::vector<Node> &Site);