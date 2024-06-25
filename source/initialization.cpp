//
// Created by Mirimi on 16/01/24.
//
#include "initialization.h"
#include "measures.h"
#include "robust_filesystem.h"


void initialize_Lattice ( std::vector <Node> &Site, const fs::path & directory_param_beta, int restart_1, struct H_parameters &Hp) {

    size_t pos = Lx + 2;
    size_t a = 2, b = Lx / 4 ;
    fs::path inputFile1 = directory_param_beta / "Lattice_dir_fin.txt";
    fs::path inputFile_A = directory_param_beta / "A_dir_fin.txt";
    fs::path inputFile_Defects = directory_param_beta / "Defect_dir_fin.txt";


    if (restart_1 == 1) {

        if ((fs::exists(inputFile1))) {

            FILE *fPsi1 = nullptr;

            if ((fPsi1 = fopen(inputFile1.c_str(), "r"))) {
                for (auto &s: Site) {

                    double angle1, angle2, modulo1, modulo2;

                    if (fscanf(fPsi1, "%lf %lf %lf %lf", &angle1, &modulo1, &angle2, &modulo2) != 4) {
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

        for (auto &s: Site) {
            s.Psi[0].t = rn::uniform_real_box(0, 2 * M_PI);
            s.Psi[1].t = rn::uniform_real_box(0, 2 * M_PI);
            s.Psi[0].r = sqrt(rn::uniform_real_box(0, 1));
            s.Psi[1].r = (1. - s.Psi[0].r * s.Psi[0].r);
        }
    }

    if (Hp.e != 0) {     //in case we have a superconductor thus we will also have the vector potential file!

        if(restart_1 == 1) {

            if ((fs::exists(inputFile_A))) {

                FILE *fA = nullptr;

                if ((fA = fopen(inputFile_A.c_str(), "r"))) {
                    for (auto &s: Site) {

                        double vec_pot_A1, vec_pot_A2;

                        if (fscanf(fA, "%lf %lf", &vec_pot_A1, &vec_pot_A2) != 2) {
                            // Handle read error for Psi[0].t
                            std::cerr << "Error reading Psi[0].t from fPsi." << std::endl;
                            break;
                        }
                        s.A[0] = vec_pot_A1;
                        s.A[1] = vec_pot_A2;
                    }

                    fclose(fA);
                }
            }

        } else {

            for (auto &s: Site) {
                s.A[0] = 0.0;
                s.A[1] = 0.0;
            }
        }
    }

    if (restart == 1) {

        if ((fs::exists(inputFile_Defects))) {

            FILE *fDef = nullptr;

            if ((fDef = fopen(inputFile_Defects.c_str(), "r"))) {
                for (auto &s: Site) {

                    double defect_1, defect_2;

                    if (fscanf(fDef, "%lf %lf", &defect_1, &defect_2) != 2) {
                        // Handle read error for Psi[0].t
                        std::cerr << "Error reading Psi[0].t from fPsi." << std::endl;
                        break;
                    }
                    s.defect[0] = defect_1;
                    s.defect[1] = defect_2;
                }

                fclose(fDef);
            }
        }
    }
    else{

        for (auto &s:Site){
            s.defect[0] = 0.0;
            s.defect[1] = 0.0;
        }

        int ix = pos % Lx;
        int iy = pos / Ly;
        Site [ix + iy * Ly].defect[0] = 1;

        int iy_r = iy + a;
        Site [ix + iy_r * Ly].defect[1] = 1;

        int ix_d = ix + b ;
        Site [ix_d + iy * Ly].defect[1] = 1;

        Site [ix_d + iy_r * Ly].defect[0] = 1;

    }


}