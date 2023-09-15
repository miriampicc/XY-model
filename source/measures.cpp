#include "measures.h"

// Define energy function
void energy(const std::vector<double>& spins, struct Measures &mis ) {
    double sum_cosines = 0.0;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            sum_cosines += cos(spins[i+j*L] - spins[((i + 1) % L)+j*L]) +
                           cos(spins[i+j*L] - spins[i+((j + 1)% L)*L]);
        }
    }
    mis.E= -sum_cosines;
}

// Function to calculate the total magnetization of the lattice

double magnetization (const std::vector<double>& spins, struct Measures &mis, int N ) {
    double M_x = 0.0, M_y = 0.0;
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            M_x += cos(spins[i+j*L]);
            M_y += sin(spins[i+j*L]);
        }
    }
    M_x /= (N);
    M_y /= (N);

    mis.M = M_x * M_x + M_y * M_y;
}


//Function to calculate the superfluid stiffness (or better the spin stiffness)
double helicity_modulus (const std::vector<double>& spins, struct Measures &mis, int N ){

    double sum_sines = 0.0;

    //calcolato solo lungo x, per ora

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            int ip =i+1;
            if( ip==L){
                ip =0;
            }
            sum_sines += sin(spins[i+j*L] - spins[ip+j*L]) ;
            mis.Jd += cos(spins[i+j*L] - spins[ip+j*L]) ;
        }
    }
    mis.Ic = sum_sines / N;
    mis.Jd /= N;
}
