//
// Created by Mirimi on 17/01/24.
//
#include "montecarlo.h"
#include "main.h"
#include "rng.h"

void metropolis(std::vector<Node> &Site, struct MC_parameters &MC, struct H_parameters &Hp,  double T){

    double l, d_theta, rand;
    double acc_rate=0.5, acc_theta=0.;
    std::array<O2, 2> NewPsi;
    std::array<O2, 2> OldPsi;
    double newE, oldE, deltaE;
    int N= L*L;
    int k,i,j;

    for (int iy = 0; iy < L; iy++) {
        for (int ix = 0; ix < L; ix++) {
            /*choose randomly a site of the lattice*/
            i = rn::uniform_integer_box(0, N-1);

            /*************PSI UPDATE: density update with total density contraint **********/

            OldPsi[0] = Site[i].Psi[0];
            OldPsi[1] = Site[i].Psi[1];
            NewPsi[0] = Site[i].Psi[0];
            NewPsi[1] = Site[i].Psi[1];

            l = rn::uniform_real_box(0, 1);
            NewPsi[0].r = sqrt(l);
            NewPsi[1].r = sqrt(1-l);

            oldE = local_energy(OldPsi, i, Hp, Site);
            newE = local_energy(NewPsi, i, Hp, Site);
            deltaE = (newE - oldE);

            if (deltaE < 0) {
                Site[i].Psi[0] = NewPsi[0];
                Site[i].Psi[1] = NewPsi[1];
            } else {
                rand = rn::uniform_real_box(0, 1);
                if (rand < exp(-1/T * deltaE)) {

                    Site[i].Psi[0] = NewPsi[0];
                    Site[i].Psi[1] = NewPsi[1];
                }
            }

            /*************PSI UPDATE: phase update **********/

            for (int alpha = 0; alpha < 2; alpha++) {
                OldPsi[0] = Site[i].Psi[0];
                OldPsi[1] = Site[i].Psi[1];
                NewPsi[0] = Site[i].Psi[0];
                NewPsi[1] = Site[i].Psi[1];
                d_theta = rn::uniform_real_box(-MC.theta_box, MC.theta_box);
                NewPsi[alpha].t = fmod(OldPsi[alpha].t + d_theta, 2*M_PI);
                NewPsi[alpha].r = OldPsi[alpha].r;

                oldE = local_energy(OldPsi, i, Hp, Site);
                newE = local_energy(NewPsi, i, Hp, Site);
                deltaE = (newE - oldE);
                if (deltaE < 0) {
                    Site[i].Psi[alpha] = NewPsi[alpha];
                    acc_theta++;
                } else {
                    rand = rn::uniform_real_box(0, 1);

                    if (rand < exp(-1/T * deltaE)) {
                        Site[i].Psi[alpha] = NewPsi[alpha];
                        acc_theta++;
                    }
                }
            }
        }
    }

    acc_theta=(double) acc_theta/(2*N);
    MC.theta_box= MC.theta_box*((0.5*acc_theta/acc_rate)+0.5);

}

double local_energy(std::array<O2, 2> &Psi, int i, H_parameters &Hp, const std::vector<Node> &Site) {

    double h_Kinetic=0., h_Josephson=0., tot_energy;
    double gauge_phase1, gauge_phase2;
    int ix, iy;
    int nn_ip, nn_im;

    ix = i % L;
    iy = i / L;

    int ip=(ix == L-1 ? 0: ix+1);
    int ipx= ip+L*(iy);                    //Relevant
    int jp=(iy == L-1 ? 0: iy+1);
    int ipy= ix+(L*jp);                    //Relevant
    int imx= (ix == 0 ? L-1: ix-1)+L*(iy); //Relevant
    int imy= ix+L*((iy == 0 ? L-1: iy-1)); //Relevant

    for(int alpha=0; alpha<2; alpha ++) {

        for (int vec = 0; vec < 2; vec++) {
            if (vec == 0) {
                 nn_ip = ipx;
                 nn_im = imx;
            } else if (vec == 1) {
                 nn_ip = ipy;
                 nn_im = imy;
            }
            gauge_phase1 = Site[nn_ip].Psi[alpha].t - Psi[alpha].t ;
            gauge_phase2 = Psi[alpha].t - Site[nn_im].Psi[alpha].t ;
            h_Kinetic -=  (Psi[alpha].r * Site[nn_ip].Psi[alpha].r) * cos(gauge_phase1);
            h_Kinetic -=  (Psi[alpha].r * Site[nn_im].Psi[alpha].r) * cos(gauge_phase2);
        }
    }

    h_Josephson +=  Hp.K * (Psi[0].r * Psi[1].r) * (Psi[0].r * Psi[1].r) * (cos(2*(Psi[0].t -Psi[1].t)));

    tot_energy=  h_Kinetic + h_Josephson;

    return tot_energy;
}
