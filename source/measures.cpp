#include "measures.h"
#include <iostream>
#include <cmath>
#include <fstream>


// Define energy function
void energy(struct Measures &mis, struct H_parameters &Hp, std::vector<Node> &Site ) {

    double first_layer = 0.0;
    double second_layer = 0.0;
    double interaction = 0.0;

    for (size_t i = 0; i < L; i++) {
        for (size_t j = 0; j < L; j++) {
            first_layer += Site[i+j*L].Psi[0].r * Site[((i + 1) % L)+j*L].Psi[0].r *(cos(Site[i+j*L].Psi[0].t - Site[((i + 1) % L)+j*L].Psi[0].t))
                            +Site[i+j*L].Psi[0].r*Site[i+((j + 1)% L)*L].Psi[0].r*(cos(Site[i+j*L].Psi[0].t - Site[i+((j + 1)% L)*L].Psi[0].t));
            second_layer += Site[i+j*L].Psi[1].r* Site[((i + 1) % L)+j*L].Psi[1].r *(cos(Site[i+j*L].Psi[1].t - Site[((i + 1) % L)+j*L].Psi[1].t)
                            + Site[i+j*L].Psi[1].r * Site[i+((j + 1)% L)*L].Psi[1].r *cos(Site[i+j*L].Psi[1].t - Site[i+((j + 1)% L)*L].Psi[1].t));
            interaction +=  Site[i+j*L].Psi[1].r * Site[i+j*L].Psi[0].r * cos(2*(Site[i+j*L].Psi[1].t-Site[i+j*L].Psi[0].t));
        }
    }
    mis.E = - first_layer - second_layer + Hp.K * interaction;
}

// Function to calculate the total magnetization of the lattice

void single_magnetization (std::vector<Node> &Site, struct Measures &mis, int N ) {

    double Mx[2]={0}, My[2]={0};

    for (auto & s: Site){
        for (int alpha=0; alpha<2; alpha++){
            Mx[alpha] += cos (s.Psi[alpha].t);
            My[alpha] += sin (s.Psi[alpha].t);
        }
    }

    for(int alpha=0; alpha<2; alpha++){
        Mx[alpha] /=N;
        My[alpha] /=N;
    }

    for(int alpha=0; alpha<2; alpha++) {
        mis.m_phase[alpha] = sqrt( Mx[alpha]*Mx[alpha] + My[alpha]*My[alpha]);
    }
}

void trsb_magnetization(struct Measures &mis, const std::vector<Node> &Site, int N) {
    //The Ising parameter m(x,y)=+/-1 indicates the chirality between the two phases.

    long double phi_shifted = 0.;
    for (size_t iy = 0; iy < L; iy++) {
        for (size_t ix = 0; ix < L; ix++) {
            size_t i = ix + L * (iy);
            phi_shifted = Site[i].Psi[1].t - Site[i].Psi[0].t;
            while (phi_shifted >= M_PI) {
                phi_shifted -= 2 * M_PI;
            }
            while (phi_shifted < -M_PI) {
                phi_shifted += 2 * M_PI;
            }
            if (phi_shifted > 0) {
                mis.trsb_m += 1;
            } else if (phi_shifted < 0) {
                mis.trsb_m += (-1);
            }
        }
    }
}




void helicity_modulus (struct H_parameters &Hp, const std::vector<Node> &Site, struct Measures &mis, int N ){

    double sum_sines[2] = {0};
    double sum_cos[2] = {0};


    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {

            for (int alpha = 0; alpha < 2; ++alpha) {
                sum_sines[alpha] += Site[i+j*L].Psi[alpha].r * Site[((i + 1) % L)+j*L].Psi[alpha].r * sin(Site[i+j*L].Psi[alpha].t - Site[((i + 1) % L)+j*L].Psi[alpha].t);
                sum_cos[alpha] += Site[i+j*L].Psi[alpha].r * Site[((i + 1) % L)+j*L].Psi[alpha].r * cos(Site[i+j*L].Psi[alpha].t - Site[((i + 1) % L)+j*L].Psi[alpha].t);
            }
        }
    }
    sum_sines[0] /= N;
    sum_cos[0] /= N;

    sum_sines[1] /= N;
    sum_cos[1] /= N;

    for (int alpha = 0; alpha < 2; ++alpha) {

        mis.Ic[alpha] = sum_sines[alpha];
        mis.Jd[alpha] = sum_cos[alpha];
    }
}

void vortex (const std::vector<Node> &Site, struct Measures &mis, int N ) {

    double sq; //this should represent the sum of the angle of every square
    double phi_1, phi_1_fin,  phi_2, phi_2_fin, phi_3, phi_3_fin, phi_4, phi_4_fin;


    int n_plus [2] = {0};
    int n_minus [2] = {0};
    int n_tot [2] = {0};

    for (int i=0; i<L; i++ ) {

        for (int j=0; j<L; j++ ) {

            for (int alpha = 0; alpha < 2; ++alpha) {

                sq = 0 ;
                phi_1 = 0;
                phi_2 = 0;
                phi_3 = 0;
                phi_4 = 0;

                if ( j == (L-1) & i != (L-1)){

                    phi_1 = Site[(i+1)+j*L].Psi[alpha].t - Site[i+j*L].Psi[alpha].t ;
                    phi_1_fin = wrapToPi(phi_1);

                    phi_2 = Site[(i+1)+(j+1-L)*L].Psi[alpha].t - Site[(i+1)+j*L].Psi[alpha].t ;
                    phi_2_fin = wrapToPi(phi_2);

                    phi_3 = Site[i+(j+1-L)*L].Psi[alpha].t - Site[(i+1)+(j+1-L)*L].Psi[alpha].t ;
                    phi_3_fin = wrapToPi(phi_3);

                    phi_4 = Site[i+j*L].Psi[alpha].t - Site[i+(j+1-L)*L].Psi[alpha].t ;
                    phi_4_fin = wrapToPi(phi_4);

                } else if ( j != (L-1) & i == (L-1) ) {

                    phi_1 = Site[(i+1-L)+j*L].Psi[alpha].t - Site[i+j*L].Psi[alpha].t ;
                    phi_1_fin = wrapToPi(phi_1);

                    phi_2 = Site[(i+1-L)+(j+1)*L].Psi[alpha].t - Site[(i+1-L)+j*L].Psi[alpha].t ;
                    phi_2_fin = wrapToPi(phi_2);

                    phi_3 = Site[i+(j+1)*L].Psi[alpha].t - Site[(i+1-L)+(j+1)*L].Psi[alpha].t ;
                    phi_3_fin = wrapToPi(phi_3);

                    phi_4 = Site[i+j*L].Psi[alpha].t - Site[i+(j+1)*L].Psi[alpha].t ;
                    phi_4_fin = wrapToPi(phi_4);


                } else if ( j == (L-1) & i == (L-1) ) {

                    phi_1 = Site[(i+1-L)+j*L].Psi[alpha].t - Site[i+j*L].Psi[alpha].t ;
                    phi_1_fin = wrapToPi(phi_1);

                    phi_2 = Site[(i+1-L)+(j+1-L)*L].Psi[alpha].t - Site[(i+1-L)+j*L].Psi[alpha].t ;
                    phi_2_fin = wrapToPi(phi_2);

                    phi_3 = Site[i+(j+1-L)*L].Psi[alpha].t - Site[(i+1-L)+(j+1-L)*L].Psi[alpha].t ;
                    phi_3_fin = wrapToPi(phi_3);

                    phi_4 = Site[i+j*L].Psi[alpha].t - Site[i+(j+1-L)*L].Psi[alpha].t;
                    phi_4_fin = wrapToPi(phi_4);


                } else {

                    phi_1 = Site[(i+1)+j*L].Psi[alpha].t - Site[i+j*L].Psi[alpha].t ;

                    // Wrap phi_1 to the range [-π, π]
                    phi_1_fin = wrapToPi(phi_1);


                    phi_2 = Site[(i+1)+(j+1)*L].Psi[alpha].t - Site[(i+1)+j*L].Psi[alpha].t ;
                    phi_2_fin = wrapToPi(phi_2);


                    phi_3 = Site[i+(j+1)*L].Psi[alpha].t - Site[(i+1)+(j+1)*L].Psi[alpha].t ;
                    phi_3_fin = wrapToPi(phi_3);

                    phi_4 = Site[i+j*L].Psi[alpha].t - Site[i+(j+1)*L].Psi[alpha].t ;
                    phi_4_fin = wrapToPi(phi_4);

                }

                sq = phi_1_fin + phi_2_fin + phi_3_fin + phi_4_fin ;

                sq /= 2*M_PI;

                if ( sq == 1 ) {
                    n_plus[alpha] ++;
                } else if ( sq == -1) {
                    n_minus[alpha] ++;
                }

            }
        }
    }

    for (int alpha = 0; alpha < 2; ++alpha) {

        mis.vortices[alpha] = n_plus[alpha];
        mis.antivortices[alpha] = n_minus[alpha];
    }
}



double wrapToPi(double angle) {

    if (angle >= M_PI) {
        angle -= 2*M_PI;
    }   else if (angle < -M_PI) {
        angle += 2*M_PI;
    }
    return angle;
}


void save_lattice(const std::vector<Node> &Site, const fs::path &directory_write, const std::string &configuration) {
    fs::path outputFile = directory_write / configuration;

    // Open the file for writing in text mode
    FILE *fPsi = fopen(outputFile.c_str(), "w");

    // Check if the file is opened successfully
    if ((fPsi != nullptr)) {

        for(auto & s : Site){
            // Write the Psi field data of the node to the file in text format
            fprintf(fPsi, "%.8lf %.8lf %.8lf %.8lf\n", s.Psi[0].t, s.Psi[0].r, s.Psi[1].t, s.Psi[1].r);
        }
        fclose(fPsi);
    } else {
        // Handle the case where the file couldn't be opened
        std::cerr << "Error opening file for writing: " << outputFile << std::endl;
    }
}
