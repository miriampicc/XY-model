//
// Created by Mirimi on 16/01/24.
//
#include "initialization.h"
#include "measures.h"
#include "robust_filesystem.h"


void initialize_Lattice ( std::vector <Node> &Site, const fs::path & directory_read, const fs::path & directory_write, int restart, struct H_parameters &Hp) {


    fs::path inputFile1 = directory_read / "Lattice_dir_fin.txt";

    if (restart == 1) {

        if((fs::exists(inputFile1))){

            FILE *fPsi1= nullptr;

            if((fPsi1=fopen(inputFile1.c_str(), "r"))) {
                for(auto & s : Site){

                    double angle1, angle2, modulo1, modulo2;

                    if (fscanf(fPsi1, "%lf %lf %lf %lf", &angle1, &modulo1, &angle2, &modulo2) !=4) {
                        // Handle read error for Psi[0].t
                        std::cerr << "Error reading Psi[0].t from fPsi." << std::endl;
                        break;
                    }
                    s.Psi[0].t = angle1;
                    s.Psi[0].r = modulo1;
                    s.Psi[1].t = angle2;
                    s.Psi[1].r = modulo2;
                }

                fclose(fPsi1);
            }
        }

    } else {

        for(auto & s : Site){
            s.Psi[0].t = rn::uniform_real_box(0, 2*M_PI);
            s.Psi[1].t = rn::uniform_real_box(0, 2*M_PI);
            s.Psi[0].r = sqrt(rn::uniform_real_box(0, 1));
            s.Psi[1].r= (1. - s.Psi[0].r*s.Psi[0].r);
        }
    }

    save_lattice(Site, directory_write, std::string("Lattice_dir_init.txt"));

}