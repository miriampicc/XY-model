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
            //i = k % L;
            //j = k / L;

            /*************PSI UPDATE: phase update **********/

            for (int alpha = 0; alpha < 2; alpha++) {
                OldPsi[0] = Site[i].Psi[0];
                OldPsi[1] = Site[i].Psi[1];
                NewPsi[0] = Site[i].Psi[0];
                NewPsi[1] = Site[i].Psi[1];
                d_theta = rn::uniform_real_box(-MC.theta_box, MC.theta_box);
                NewPsi[alpha].t = fmod(OldPsi[alpha].t + d_theta, 2*M_PI);

                oldE = local_energy(OldPsi, i, Hp, Site, alpha);
                newE = local_energy(NewPsi, i, Hp, Site, alpha);
                deltaE = (newE - oldE);
                if (deltaE < 0) {
                    Site[i].Psi[alpha] = NewPsi[alpha];
                    acc_theta++;
                } else {
                    rand = rn::uniform_real_box(0, 1);
                    //Boltzmann weight: exp(-\beta \Delta E) E= h³ \sum_i E(i)
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

double local_energy(std::array<O2, 2> &Psi, int i, H_parameters &Hp, const std::vector<Node> &Site, int alpha) {


    double cos_in_plane, interaction, tot_energy;
    int ix, iy;

    ix = i % L;
    iy = i / L;

    int ip=(ix == L-1 ? 0: ix+1);
    int ipx= ip+L*(iy);                    //Relevant
    int jp=(iy == L-1 ? 0: iy+1);
    int ipy= ix+(L*jp);                    //Relevant
    int imx= (ix == 0 ? L-1: ix-1)+L*(iy); //Relevant
    int imy= ix+L*((iy == 0 ? L-1: iy-1)); //Relevant


    if (alpha == 0) {
        cos_in_plane = - Hp.J1 * (cos (Psi[alpha].t-Site[ipx].Psi[alpha].t) + cos (Psi[alpha].t-Site[ipy].Psi[alpha].t) +
                                  cos (Site[imy].Psi[alpha].t-Psi[alpha].t) + cos (Site[imx].Psi[alpha].t-Psi[alpha].t) );
        interaction = + Hp.K * (cos (2*(Psi[alpha].t-Psi[alpha+1].t)));

        tot_energy = cos_in_plane + interaction;
    }
    else{
        cos_in_plane = - Hp.J2 * (cos (Psi[alpha].t-Site[ipx].Psi[alpha].t) + cos (Psi[alpha].t-Site[ipy].Psi[alpha].t) +
                                  cos (Site[imy].Psi[alpha].t-Psi[alpha].t) + cos (Site[imx].Psi[alpha].t-Psi[alpha].t) );
        interaction = + Hp.K * (cos (2*(Psi[alpha].t-Psi[alpha-1].t)));

        tot_energy = cos_in_plane + interaction;
    }

    return tot_energy;
}
