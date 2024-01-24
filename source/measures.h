#ifndef XY_MODEL_MEASURES_H
#define XY_MODEL_MEASURES_H

#include "main.h"

struct Measures {
    double E=0; //Energy
    double m_phase[2] ={0}; //x-component of the single phase. We are defining a sort of vector with two components
    double trsb_m=0;
    double Ic[2]= {0};
    double Jd[2]= {0};
    int vortices[2]= {0};
    int antivortices[2] = {0};

    void reset(){
        *this = Measures();
    }
};

void single_magnetization (std::vector<Node> &Site, struct Measures &mis, int N );
void energy(struct Measures &mis, struct H_parameters &Hp, std::vector<Node> &Site );
void trsb_magnetization(struct Measures &mis, const std::vector<Node> &Site, int N);
void helicity_modulus (struct H_parameters &Hp, const std::vector<Node> &Site, struct Measures &mis, int N );
void vortex (const std::vector<Node> &Site, struct Measures &mis, int N );
void save_lattice(const std::vector<Node> &Site, const fs::path &directory_write, const std::string &configuration);
double wrapToPi(double angle);


#endif //XY_MODEL_MEASURES_H
