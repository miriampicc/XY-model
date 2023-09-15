#ifndef XY_MODEL_MEASURES_H
#define XY_MODEL_MEASURES_H

#include "main.h"

struct Measures {
    double E=0;
    double M=0;
    double Ic=0;
    double Jd=0;

    void reset(){
        *this = Measures();
    }
};

double magnetization(const std::vector<double>& spins, struct Measures &mis, int N );
void energy(const std::vector<double>& spins, struct Measures &mis );
double helicity_modulus (const std::vector<double>& spins, struct Measures &mis, int N );


#endif //XY_MODEL_MEASURES_H
