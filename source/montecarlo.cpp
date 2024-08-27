//
// Created by Mirimi on 17/01/24.
//
#include "montecarlo.h"
#include "main.h"
#include "rng.h"
#include "initialization.h"

void metropolis(std::vector<Node> &Site, struct MC_parameters &MC, struct H_parameters &Hp,  double my_beta) {

    double rand, d_A, d_theta, d_X, X;
    double acc_rate = 0.5, acc_theta = 0., acc_A = 0., acc_density = 0.;
    std::array<O2, 2> NewPsi{};
    std::array<O2, 2> OldPsi{};
    double OldA, NewA;
    double newE, oldE, deltaE;
    double target_acc_rate = 0.3;
    size_t N = Lx * Ly;


    for (size_t iy = 0; iy < Ly; iy++) {
        for (size_t ix = 0; ix < Lx; ix++) {

            size_t i = ix + Lx * (iy);
            /*************PSI UPDATE: first component density update and then second component, no total density contraint *********/

            for (int alpha = 0; alpha < 2; alpha++) {

                OldPsi[0] = Site[i].Psi[0];
                OldPsi[1] = Site[i].Psi[1];
                NewPsi[0] = Site[i].Psi[0];
                NewPsi[1] = Site[i].Psi[1];

                //In this way we have the square

                X = Site[i].Psi[alpha].r * Site[i].Psi[alpha].r;
                //d_X = rn::uniform_real_box(-0.25 * X, 0.25 * X);
                d_X = rn::uniform_real_box(-MC.theta_box_density * X, MC.theta_box_density * X);
                //std::cout << "  d X  = " << d_X << std::endl;

                //d_X = rn::uniform_real_box(-MC.theta_box_density * X, MC.theta_box_density * X);
                //NewPsi[alpha].r = sqrt(fabs(X + d_X)); // Ensure positive densities

                //std::cout << "  d X  = " << d_X << std::endl;

                //std::cout << "  Theta box density   = " << MC.theta_box_density  << std::endl;


                NewPsi[alpha].r = sqrt(X + d_X);

                if (std::isnan(NewPsi[alpha].r )) {

                    std::cout << "  Errore, NewPsi[alpha].r   = " << NewPsi[alpha].r  << std::endl;
                    exit(1);
                }

                 //std::cout << "  Update Psi = " << NewPsi[alpha].r << "  Old Psi = "<< OldPsi[alpha].r << std::endl;


                oldE = local_energy(OldPsi, i, Hp, Site);
                newE = local_energy(NewPsi, i, Hp, Site);

                //std::cout << "  Old E  = " << oldE  << "   New E = " << newE  << std::endl;

                deltaE = (newE - oldE);
                //std::cout << "Delta E = " << deltaE << std::endl;
                if (deltaE < 0) {
                    Site[i].Psi[alpha] = NewPsi[alpha];
                    acc_density++;
                    //std::cout << "YES  " << "Delta E = " << deltaE  << std::endl;
                } else {
                    rand = rn::uniform_real_box(0, 1);
                    if (rand < exp(-my_beta * deltaE)) {
                        Site[i].Psi[alpha] = NewPsi[alpha];
                        acc_density++;
                        //std::cout << "YES  " << "Delta E = " << deltaE << std::endl;
                    }
                }
            }

                //And, remember, you are interested in the total density fluctuations, thus  NewPsi[0].r^2+ NewPsi[1].r^2

                /*************PSI UPDATE: phase update **********/

                for ( int alpha = 0; alpha < 2; alpha++) {
                    OldPsi[0] = Site[i].Psi[0];
                    OldPsi[1] = Site[i].Psi[1];
                    NewPsi[0] = Site[i].Psi[0];
                    NewPsi[1] = Site[i].Psi[1];
                    d_theta = rn::uniform_real_box(-MC.theta_box, MC.theta_box);
                    NewPsi[alpha].t = fmod(OldPsi[alpha].t + d_theta, 2 * M_PI);
                    //NewPsi[alpha].r = OldPsi[alpha].r;  ////da togliere

                    oldE = local_energy(OldPsi, i, Hp, Site);
                    newE = local_energy(NewPsi, i, Hp, Site);
                    deltaE = (newE - oldE);
                    if (deltaE < 0) {
                        Site[i].Psi[alpha] = NewPsi[alpha];
                        acc_theta++;
                    } else {
                        rand = rn::uniform_real_box(0, 1);

                        if (rand < exp(-my_beta * deltaE)) {
                            Site[i].Psi[alpha] = NewPsi[alpha];
                            acc_theta++;
                        }
                    }
                }
            }
        }


        /***********VECTOR POTENTIAL*******/
        if (Hp.e != 0) {   //then we also have to update the vector potential

            for (int iy = 0; iy < Ly; iy++) {
                for (int ix = 0; ix < Lx; ix++) {

                    size_t i = ix + Lx * (iy);
                    for (int alpha = 0; alpha < 2; alpha++) {

                        OldA = Site[i].A[alpha];
                        d_A = rn::uniform_real_box(-MC.theta_box_A, MC.theta_box_A);
                        NewA = OldA + d_A;

                        oldE = local_energy_A(OldA, i, alpha, Hp, Site);
                        newE = local_energy_A(NewA, i, alpha, Hp, Site);

                        deltaE = (newE - oldE);
                        if (deltaE < 0.) {
                            Site[i].A[alpha] = NewA;
                            acc_A++;
                            //std::cout<<"Updating A form "<< OldA <<" to "<< NewA << " on site "<< i<<std::endl;
                        } else {
                            rand = rn::uniform_real_box(0, 1);
                            if (rand < exp(-my_beta * deltaE)) {
                                Site[i].A[alpha] = NewA;
                                acc_A++;
                                //std::cout<<"Updating A form "<< OldA <<" to "<< NewA << " on site "<< i<<std::endl;
                            }
                        }
                    }
                }
            }
        }

        acc_theta = (double) acc_theta / static_cast<double>(2 * N);
        acc_A = (double) acc_A / static_cast<double>(2 * N);
        acc_density = (double) acc_density / static_cast<double>(2 * N);

        //std::cout << "Acc = " << acc_theta << std::endl;


        MC.theta_box = MC.theta_box * ((0.5 * acc_theta / acc_rate) + 0.5);
        MC.theta_box_A = MC.theta_box_A * ((0.5 * acc_A / acc_rate) + 0.5);
        //MC.theta_box_density = MC.theta_box_density * ((0.5 * acc_density / acc_rate) + 0.5);

        //std::cout<<" Update theta_box_desity  "<< MC.theta_box_density << std::endl;

        // Ensure the step sizes do not become too small
        //MC.theta_box_density = std::max(MC.theta_box_density, min_theta_box_density);

    }

    double local_energy(std::array<O2, 2> &Psi, size_t i, H_parameters &Hp, const std::vector<Node> &Site) {

        double h_Kinetic = 0., h_Josephson = 0., tot_energy, dens_fluct = 0.;
        double gauge_phase1, gauge_phase2;
        size_t ix, iy;
        size_t nn_ip, nn_im;

        ix = i % Lx;
        iy = i / Ly;

        size_t ip = (ix == Lx - 1 ? 0 : ix + 1);
        size_t ipx = ip + Lx * (iy);                    //Relevant
        size_t jp = (iy == Ly - 1 ? 0 : iy + 1);
        size_t ipy = ix + (Lx * jp);                    //Relevant
        size_t imx = (ix == 0 ? Lx - 1 : ix - 1) + Lx * (iy); //Relevant
        size_t imy = ix + Lx * ((iy == 0 ? Ly - 1 : iy - 1)); //Relevant

        for (int alpha = 0; alpha < 2; alpha++) {

            for (int vec = 0; vec < 2; vec++) {
                if (vec == 0) {
                    nn_ip = ipx;
                    nn_im = imx;
                } else if (vec == 1) {
                    nn_ip = ipy;
                    nn_im = imy;
                }
                //NB: vec is useful to consider the potential vector in both the components
                gauge_phase1 = Site[nn_ip].Psi[alpha].t - Psi[alpha].t + Hp.e * Site[i].A[vec];
                gauge_phase2 = Psi[alpha].t - Site[nn_im].Psi[alpha].t + Hp.e * Site[nn_im].A[vec];
                h_Kinetic += 0.5 * (Psi[alpha].r * Psi[alpha].r + Site[nn_ip].Psi[alpha].r * Site[nn_ip].Psi[alpha].r - 2 * (Psi[alpha].r * Site[nn_ip].Psi[alpha].r) * cos(gauge_phase1));
                h_Kinetic += 0.5 * (Psi[alpha].r * Psi[alpha].r + Site[nn_im].Psi[alpha].r * Site[nn_im].Psi[alpha].r - 2 * (Psi[alpha].r * Site[nn_im].Psi[alpha].r) * cos(gauge_phase2));
            }
        }

        //h_Josephson +=  2 * Hp.b2 * (Psi[0].r * Psi[1].r) * (Psi[0].r * Psi[1].r) * (cos(2*(Psi[0].t -Psi[1].t)) - 1. );
        //dens_fluct -=  ((Psi[0].r * Psi[0].r) + (Psi[1].r * Psi[1].r)) * ( 1 - (Hp.b1 + Hp.b2) * ((Psi[0].r * Psi[0].r) + (Psi[1].r * Psi[1].r)) ) ;

        //h_Josephson += 2 * Hp.K * (Psi[0].r * Psi[1].r) * (Psi[0].r * Psi[1].r) * (cos(2 * (Psi[0].t - Psi[1].t)) - 1.);
        //dens_fluct += Hp.a * ((Psi[0].r * Psi[0].r) + (Psi[1].r * Psi[1].r)) * (-1 + 0.5 * ((Psi[0].r * Psi[0].r) + (Psi[1].r * Psi[1].r)));

        h_Josephson += Hp.a * Hp.K * (Psi[0].r * Psi[1].r) * (Psi[0].r * Psi[1].r) * (cos(2 * (Psi[0].t - Psi[1].t)) - 1.);
        dens_fluct += Hp.a * ((Psi[0].r * Psi[0].r) + (Psi[1].r * Psi[1].r)) * (-1 + 0.5 * ((Psi[0].r * Psi[0].r) + (Psi[1].r * Psi[1].r)));


        //double diff = (Psi[0].t - Psi[1].t);
        //std::cout << "difference phase = " << diff  << std::endl;

        //std::cout << "Dens fluct = " << dens_fluct << "  Josephson = " << h_Josephson << std::endl;

        tot_energy = h_Kinetic + h_Josephson + dens_fluct;

        return tot_energy;
    }

    double local_energy_A(double A, size_t i, int alpha, H_parameters &Hp, const std::vector<Node> &Site) {

        double h_Kinetic = 0., tot_energy;
        double gauge_phase1, A_plaq, A_2 = 0.;
        size_t ix, iy;
        size_t nn_ip, nn_ipl, nn_iml;

        ix = i % Lx;
        iy = i / Ly;

        size_t ip = (ix == Lx - 1 ? 0 : ix + 1);
        size_t ipx = ip + Lx * (iy);                    //Relevant
        size_t jp = (iy == Ly - 1 ? 0 : iy + 1);
        size_t ipy = ix + (Lx * jp);                    //Relevant
        size_t imx = (ix == 0 ? Lx - 1 : ix - 1) + Lx * (iy); //Relevant
        size_t imy = ix + Lx * ((iy == 0 ? Ly - 1 : iy - 1)); //Relevant

        if (alpha == 0) {
            nn_ip = ipx;
        } else if (alpha == 1) {
            nn_ip = ipy;
        }

        for (int vec = 0; vec < 2; vec++) {
            gauge_phase1 = Site[nn_ip].Psi[vec].t - Site[i].Psi[vec].t + Hp.e * A;
            h_Kinetic -= 0.5 * (Site[i].Psi[vec].r * Site[i].Psi[vec].r + Site[nn_ip].Psi[vec].r * Site[nn_ip].Psi[vec].r - 2 * (Site[i].Psi[vec].r * Site[nn_ip].Psi[vec].r) * cos(gauge_phase1));
        }

        //considering now all the plaquettes that involve the vectors

        for (int l = 0; l < 2; l++) {
            if (l == 0) {
                nn_ipl = ipx;
            } else if (l == 1) {
                nn_ipl = ipy;
            }
            if (l != alpha) {
                A_plaq = A + Site[nn_ip].A[l] - Site[nn_ipl].A[alpha] - Site[i].A[l];
                A_2 += A_plaq * A_plaq;
            }
        }

        for (int l = 0; l < 2; l++) {
            if (l == 0) {
                nn_iml = imx;
            } else if (l == 1) {
                nn_iml = imy;
            }
            if (l != alpha) {
                A_plaq = Site[nn_iml].A[alpha] + Site[nn(nn_iml, alpha, 1)].A[l] - A - Site[nn_iml].A[l];
                A_2 += A_plaq * A_plaq;
            }
        }

        tot_energy = h_Kinetic + 0.5 * A_2;

        return tot_energy;

    }
