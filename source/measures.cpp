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

    sum_cosines = 0.0;
    std::cout<< sum_cosines << std::endl;
    mis.E= -sum_cosines;
}

// Function to calculate the total magnetization of the lattice

void magnetization (const std::vector<double>& spins, struct Measures &mis, int N ) {
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
void helicity_modulus (const std::vector<double>& spins, struct Measures &mis, int N ){

    double sum_sines = 0.0;
    double sum_cos = 0.0;

    //calcolato solo lungo x, per ora

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {

            //sum_sines += sin(spins[i+j*L] - spins[((i + 1) % L)+j*L]); //+
                         //sin(-spins[i+j*L] + spins[i+((j + 1)% L)*L]); //In this way I am taking into account periodic boundary cond
            //sum_cos += cos(spins[i+j*L] - spins[((i + 1) % L)+j*L]); // +   //NB: I HAVE TO THINK ABOUT A BETTER WAY TO DO IT!!
                      //cos(spins[i+j*L] - spins[i+((j + 1)% L)*L]);
        }
    }

    sum_sines = 0.0;
    sum_cos = 0.0;

    mis.Ic = sum_sines / (N) ;
    mis.Jd = sum_cos / (N);

}

void vortex (const std::vector<double>& spins, struct Measures &mis, int N ) {

    double sq; //this should represent the sum of the angle of every square
    int i, j;
    int n_plus, n_minus, n_tot;
    double phi_1, phi_1_fin,  phi_2, phi_2_fin, phi_3, phi_3_fin, phi_4, phi_4_fin;


    n_plus = 0;
    n_minus = 0;
    n_tot = 0;

    for ( i=0; i<L; i++ ) {

        for ( j=0; j<L; j++ ) {

            sq = 0 ;
            phi_1 = 0;
            phi_2 = 0;
            phi_3 = 0;
            phi_4 = 0;

            if ( j == (L-1) & i != (L-1)){

                phi_1 = spins[(i+1)+j*L] - spins[i+j*L] ;
                phi_1_fin = wrapToPi(phi_1);

                phi_2 = spins[(i+1)+(j+1-L)*L] - spins[(i+1)+j*L] ;
                phi_2_fin = wrapToPi(phi_2);

                phi_3 = spins[i+(j+1-L)*L] - spins[(i+1)+(j+1-L)*L] ;
                phi_3_fin = wrapToPi(phi_3);

                phi_4 = spins[i+j*L] - spins[i+(j+1-L)*L] ;
                phi_4_fin = wrapToPi(phi_4);

            } else if ( j != (L-1) & i == (L-1) ) {

                phi_1 = spins[(i+1-L)+j*L] - spins[i+j*L] ;
                phi_1_fin = wrapToPi(phi_1);

                phi_2 = spins[(i+1-L)+(j+1)*L] - spins[(i+1-L)+j*L] ;
                phi_2_fin = wrapToPi(phi_2);

                phi_3 = spins[i+(j+1)*L] - spins[(i+1-L)+(j+1)*L] ;
                phi_3_fin = wrapToPi(phi_3);

                phi_4 = spins[i+j*L] - spins[i+(j+1)*L] ;
                phi_4_fin = wrapToPi(phi_4);


            } else if ( j == (L-1) & i == (L-1) ) {

                phi_1 = spins[(i+1-L)+j*L] - spins[i+j*L] ;
                phi_1_fin = wrapToPi(phi_1);

                phi_2 = spins[(i+1-L)+(j+1-L)*L] - spins[(i+1-L)+j*L] ;
                phi_2_fin = wrapToPi(phi_2);

                phi_3 = spins[i+(j+1-L)*L] - spins[(i+1-L)+(j+1-L)*L] ;
                phi_3_fin = wrapToPi(phi_3);

                phi_4 = spins[i+j*L] - spins[i+(j+1-L)*L] ;
                phi_4_fin = wrapToPi(phi_4);


            } else {

                phi_1 = spins[(i+1)+j*L] - spins[i+j*L] ;

                // Wrap phi_1 to the range [-π, π]
                phi_1_fin = wrapToPi(phi_1);


                phi_2 = spins[(i+1)+(j+1)*L] - spins[(i+1)+j*L] ;
                phi_2_fin = wrapToPi(phi_2);


                phi_3 = spins[i+(j+1)*L] - spins[(i+1)+(j+1)*L] ;
                phi_3_fin = wrapToPi(phi_3);

                phi_4 = spins[i+j*L] - spins[i+(j+1)*L] ;
                phi_4_fin = wrapToPi(phi_4);

            }

            sq = phi_1_fin + phi_2_fin + phi_3_fin + phi_4_fin ;

            sq /= 2*M_PI;

            if ( sq == 1 ) {
                n_plus ++;
            } else if ( sq == -1) {
                n_minus++;
            }

        }
    }

    mis.n_vort = n_plus;
    mis.n_antivort = n_minus;

}



double wrapToPi(double angle) {

    if (angle >= M_PI) {
        angle -= 2*M_PI;
    }   else if (angle < -M_PI) {
        angle += 2*M_PI;
    }
    return angle;
}
