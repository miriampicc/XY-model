//
// Created by Mirimi on 16/01/24.
//
#include "initialization.h"
#include "measures.h"
#include "robust_filesystem.h"


void initialize_Lattice ( std::vector <Node> &Site, const fs::path & directory_read, const fs::path & directory_write, int restart, struct H_parameters &Hp) {


    fs::path inputFile1 = directory_read / "Lattice_dir_fin.txt";
    //fs::path outputFile1 = directory_write / "Lattice_dir_init.txt";


    if (restart == 1) {

        if((fs::exists(inputFile1))){

            FILE *fPsi1= nullptr;
            //FILE *fOutput1 = fopen(outputFile1.c_str(), "a");

            if((fPsi1=fopen(inputFile1.c_str(), "r"))) {
                for(auto & s : Site){

                    double value1, value2;

                    if (fscanf(fPsi1, "%lf %lf", &value1, &value2) !=2) {
                        // Handle read error for Psi[0].t
                        std::cerr << "Error reading Psi[0].t from fPsi." << std::endl;
                        break;
                    }
                    s.Psi[0].t = value1;
                    s.Psi[1].t = value2;

                    /*// Now, write the modified values to outputFile1 and outputFile2

                    if (fOutput1 != nullptr) {
                        fprintf(fOutput1, "%.8lf %.8lf\n", s.Psi[0].t, s.Psi[1].t);
                    } else {
                        std::cerr << "Error opening output files for writing." << std::endl;
                        break;
                    }*/
                }

                fclose(fPsi1);
                //fclose(fOutput1);
            }
        }

    } else {

        //FILE *fSpin1Output = fopen(outputFile1.c_str(), "a");  // "a" appends to the file


        for(auto & s : Site){

            s.Psi[0].t = rn::uniform_real_box(0, 2*M_PI);
            s.Psi[1].t = rn::uniform_real_box(0, 2*M_PI);

            // Write Psi[0].t to Spin_1_dir_init.txt
            /* if (fSpin1Output != nullptr) {
                fprintf(fSpin1Output, "%.8lf %.8lf\n", s.Psi[0].t, s.Psi[1].t);
            } else {
                std::cerr << "Error opening " << outputFile1 << " for writing." << std::endl;
            }  */
        }
        //fclose(fSpin1Output);
    }

    save_lattice(Site, directory_write, std::string("Lattice_dir_init.txt"));


}