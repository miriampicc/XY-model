#ifndef XY_MODEL_MEASURES_H
#define XY_MODEL_MEASURES_H

#include "main.h"

struct Measures {
    double E=0;
    double M=0;
    double Ic=0;
    double Jd=0;
    int n_vort=0;
    int n_antivort=0;

    void reset(){
        *this = Measures();
    }
};

void magnetization(const std::vector<double>& spins, struct Measures &mis, int N );
void energy(const std::vector<double>& spins, struct Measures &mis );
void helicity_modulus (const std::vector<double>& spins, struct Measures &mis, int N );
void vortex (const std::vector<double>& spins, struct Measures &mis, int N );
double wrapToPi(double angle);


#endif //XY_MODEL_MEASURES_H
