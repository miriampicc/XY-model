#ifndef XY_MODEL_MEASURES_H
#define XY_MODEL_MEASURES_H

#include "main.h"

struct Measures {
    double E=0; //Energy
    double E_kinetic=0;
    double E_josephson=0;
    double E_B=0;
    double density_fluct=0;
    double m_phase[2] ={0}; //x-component of the single phase. We are defining a sort of vector with two components
    double trsb_m=0;
    double Ic[2]= {0};
    double Jd[2]= {0};
    int vortices[2]= {0};
    int antivortices[2] = {0};
    long double dual_stiff_Z =0.; //dual stiffness along Z
    int my_rank=0;

    void reset(){
        *this = Measures();
    }
};

void single_magnetization (std::vector<Node> &Site, struct Measures &mis, size_t N );
void energy(struct Measures &mis, struct H_parameters &Hp, const std::vector<Node> &Site );
void trsb_magnetization(struct Measures &mis, const std::vector<Node> &Site);
void helicity_modulus (struct H_parameters &Hp, const std::vector<Node> &Site, struct Measures &mis, size_t N );
void dual_stiffness (struct Measures &mis, const std::vector<Node> &Site);
//void vortex (const std::vector<Node> &Site, struct Measures &mis);
void save_lattice(const std::vector<Node> &Site, const fs::path &directory_write, const std::string &configuration, struct H_parameters &Hp);
double wrapToPi(double angle);
size_t nn (size_t i, size_t coord, int dir );


#endif //XY_MODEL_MEASURES_H
