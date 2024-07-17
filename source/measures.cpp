#include "measures.h"
#include <iostream>
#include <cmath>
#include <fstream>


// Define energy function
void energy(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site ) {

    double interaction = 0.0, dens_fluct = 0.0;
    double gauge_phase , h_Kinetic=0.;
    double A_plaq, A_2=0.;
    size_t ip, nn_ip;

    for (size_t iy = 0; iy < Ly; iy++) {
        for (size_t ix = 0; ix < Ly; ix++) {

            size_t i=ix + Lx * ( iy );

            for (size_t alpha=0; alpha<2; alpha++){

                for (size_t vec2 = 0; vec2 < 2; vec2++) {
                    if (vec2 == 0) {
                        ip = (ix == Lx - 1 ? 0 : ix + 1);
                        nn_ip = ip + Lx * (iy);
                    } else if (vec2 == 1) {
                        ip = (iy == Ly - 1 ? 0 : iy + 1);
                        nn_ip = ix + Lx * ip;
                    }
                    gauge_phase = Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t + Hp.e * Site[i].A[vec2];
                    h_Kinetic -= (Site[i].Psi[alpha].r * Site[nn_ip].Psi[alpha].r) * cos(gauge_phase);
                }
            }

            //interaction +=  2 * Hp.b2 *(Site[i].Psi[1].r * Site[i].Psi[0].r) * (Site[i].Psi[1].r * Site[i].Psi[0].r) *(cos(2*(Site[i].Psi[0].t- Site[i].Psi[1].t))-1.);
            //dens_fluct -=  ((Site[i].Psi[0].r * Site[i].Psi[0].r) + (Site[i].Psi[1].r * Site[i].Psi[1].r)) * ( 1 - (Hp.b1 + Hp.b2) * ((Site[i].Psi[0].r * Site[i].Psi[0].r) + (Site[i].Psi[1].r * Site[i].Psi[1].r)) ) ;
            interaction +=  2 * Hp.K *(Site[i].Psi[1].r * Site[i].Psi[0].r) * (Site[i].Psi[1].r * Site[i].Psi[0].r) *(cos(2*(Site[i].Psi[0].t- Site[i].Psi[1].t))-1.);
            dens_fluct +=  Hp.a * ((Site[i].Psi[0].r * Site[i].Psi[0].r) + (Site[i].Psi[1].r * Site[i].Psi[1].r)) * ( - 1 + 0.5 * ((Site[i].Psi[0].r * Site[i].Psi[0].r) + (Site[i].Psi[1].r * Site[i].Psi[1].r)) ) ;


            if (Hp.e != 0) {
                for(size_t vec1=0; vec1<2; vec1++){
                    for (size_t vec2 = vec1+1; vec2 < 2; vec2++) {
                        //F_{alpha,vec}= A_alpha(r_i) + A_vec(ri+alpha) - A_alpha(r_i+vec) - A_vec(ri)
                        A_plaq = (Site[i].A[vec1] + Site[nn(i, vec1, 1)].A[vec2] - Site[nn(i, vec2, 1)].A[vec1] -
                               Site[i].A[vec2]);

                        A_2 +=  0.5 * (A_plaq * A_plaq);
                    }
                }

            }
        }
    }
    //mis.E =  h_Kinetic + interaction + A_2;
    mis.E_kinetic = h_Kinetic;
    mis.E_josephson = interaction;
    mis.E_B = A_2;
    mis.density_fluct = dens_fluct;

    mis.E =(mis.E_kinetic  + mis.E_josephson + mis.E_B + mis.density_fluct );
}

// Function to calculate the total magnetization of the lattice

void single_magnetization (std::vector<Node> &Site, struct Measures &mis, size_t N ) {

    double Mx[2]={0}, My[2]={0};

    for (auto & s: Site){
        for (int alpha=0; alpha<2; alpha++){
            Mx[alpha] += cos (s.Psi[alpha].t);
            My[alpha] += sin (s.Psi[alpha].t);
        }
    }

    for(int alpha=0; alpha<2; alpha++){
        Mx[alpha] /= static_cast<double>(N);
        My[alpha] /= static_cast<double>(N);
    }

    for(int alpha=0; alpha<2; alpha++) {
        mis.m_phase[alpha] = sqrt( Mx[alpha]*Mx[alpha] + My[alpha]*My[alpha]);
    }
}

void trsb_magnetization(struct Measures &mis, const std::vector<Node> &Site) {
    //The Ising parameter m(x,y)=+/-1 indicates the chirality between the two phases.

    long double phi_shifted = 0.;

    for (size_t iy = 0; iy < Ly; iy++) {
        for (size_t ix = 0; ix < Lx; ix++) {
            size_t i = ix + Lx * (iy);
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




void helicity_modulus (struct H_parameters &Hp, const std::vector<Node> &Site, struct Measures &mis, size_t N ){

    double J_alpha=0., DJ_alpha_Dd=0.;
    //double sum_sines[2] = {0};
    //double sum_cos[2] = {0};
    double gauge_phase;
    int vec=0; //in this case we are calculating the helicity modulus along the x direction


    for (int iy = 0; iy < Ly; iy++) {
        for (int ix = 0; ix < Lx; ix++) {

            size_t i = ix + Lx * iy;
            size_t ip = (ix == Lx-1 ? 0: ix+1);
            size_t nn_ip = ip + Lx * (iy);

            for (int alpha = 0; alpha < 2; ++alpha) {

                //gauge_phase = Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t + Hp.e * Site[i].A[vec];
                //sum_sines[alpha] += Site[i].Psi[alpha].r * Site[nn_ip].Psi[alpha].r * sin(gauge_phase);
                //sum_cos[alpha] += Site[i].Psi[alpha].r * Site[nn_ip].Psi[alpha].r * cos(gauge_phase);

                gauge_phase = Site[nn_ip].Psi[alpha].t - Site[i].Psi[alpha].t + Hp.e*Site[i].A[vec];
                J_alpha = (Site[i].Psi[alpha].r * Site[nn_ip].Psi[alpha].r) * sin(gauge_phase);
                DJ_alpha_Dd = (Site[i].Psi[alpha].r * Site[nn_ip].Psi[alpha].r) * cos(gauge_phase);
                mis.Ic[alpha] += ( J_alpha )/ N;
                mis.Jd[alpha]+= ( DJ_alpha_Dd ) / N;
            }
        }
    }
}

void dual_stiffness (struct Measures &mis, const std::vector<Node> &Site) {

    long double qx_min = (2 * M_PI) / static_cast<long double>(Lx);
    long double normalization_const = 1./((2*M_PI)*(2*M_PI)*static_cast<long double>(Lx)*static_cast<long double>(Lx));
    long double Im_dual=0., Re_dual=0.;
    long double Dx_Ay, Dy_Ax;
    int ix, iy;

    for (iy=0; iy<Ly; iy++){
        for (ix=0; ix<Lx; ix++){

            size_t i=ix +Lx*(iy);

            Dx_Ay = Site[nn(i, 0, 1)].A[1] - Site[i].A[1];
            Dy_Ax = Site[nn(i, 1, 1)].A[0] - Site[i].A[0];

            Im_dual += (Dx_Ay - Dy_Ax) * sin ((double)(qx_min * ix));
            Re_dual += (Dx_Ay - Dy_Ax) * cos ((double)(qx_min * ix));

        }
    }

    mis.dual_stiff_Z = normalization_const * ((Im_dual * Im_dual) + (Re_dual * Re_dual));

}
/*

void vortex (const std::vector<Node> &Site, struct Measures &mis) {

    double sq; //this should represent the sum of the angle of every square
    double phi_1, phi_1_fin,  phi_2, phi_2_fin, phi_3, phi_3_fin, phi_4, phi_4_fin;


    int n_plus [2] = {0};
    int n_minus [2] = {0};

    for (int i=0; i<L; i++ ) {

        for (int j=0; j<L; j++ ) {

            for (int alpha = 0; alpha < 2; ++alpha) {

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
*/



double wrapToPi(double angle) {

    if (angle >= M_PI) {
        angle -= 2*M_PI;
    }   else if (angle < -M_PI) {
        angle += 2*M_PI;
    }
    return angle;
}


void save_lattice(const std::vector<Node> &Site, const fs::path &directory_write, const std::string &configuration, struct H_parameters &Hp) {


    // Check if the directory exists, and create it if it doesn't
    if (!fs::exists(directory_write))
        fs::create_directories(directory_write);

    std::string sPsi;
    sPsi= std::string("Lattice_")+ configuration + std::string(".txt");
    fs::path psi_file = directory_write / sPsi;

    // Open the file for writing in text mode
    FILE *fPsi = fopen(psi_file.c_str(), "w");

    // Check if the file is opened successfully
    if ((fPsi != nullptr)) {

        for(auto & s : Site){
            // Write the Psi field data of the node to the file in text format
            fprintf(fPsi, "%.8lf %.8lf %.8lf %.8lf\n", s.Psi[0].t, s.Psi[0].r, s.Psi[1].t, s.Psi[1].r);
        }
        fclose(fPsi);
    } else {
        // Handle the case where the file couldn't be opened
        std::cerr << "Error opening file for writing: " << psi_file << std::endl;
    }

    if (Hp.e != 0) {

        std::string sA;
        sA= std::string("A_")+ configuration + std::string(".txt");
        fs::path a_file = directory_write / sA;

        FILE *fA = fopen(a_file.c_str(), "w");

        if ((fA != nullptr)) {

            for(auto & s: Site){
                fprintf(fA, "%.8lf %.8lf\n", s.A[0], s.A[1]);
            }
            fclose(fA);
        } else {
            // Handle the case where the file couldn't be opened
            std::cerr << "Error opening file for writing: " << a_file << std::endl;
        }
    }
}


size_t nn (size_t i, size_t coord, int dir ){
    size_t ix = i % Lx;
    size_t iy = (i/Ly) % Ly;

    if(coord==0){
        int ix_new = static_cast<int>(ix)+ (dir == 0 ? 0 : (dir > 0 ? 1 : -1));
        if(ix_new==Lx) { ix_new=0;}
        if(ix_new < 0){ ix_new=static_cast<int>(Lx-1);}
        int iy_new=static_cast<int>(iy);
        return (static_cast<size_t>(ix_new) + Lx * (static_cast<size_t>(iy_new)));

    }
    if(coord==1){
        int iy_new= static_cast<int>(iy) + (dir == 0 ? 0 : (dir > 0 ? 1 : -1));
        if(iy_new==static_cast<int>(Ly)){ iy_new=0;}
        if(iy_new<0){ iy_new=static_cast<int>(Ly-1);}
        int ix_new=static_cast<int>(ix);
        return (static_cast<size_t>(ix_new) + Ly * (static_cast<size_t>(iy_new)));
    }
    return 1;
}
