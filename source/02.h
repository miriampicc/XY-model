//
// Created by Mirimi on 17/01/24.
//

#ifndef CMT_02_H
#define CMT_02_H

#include <math.h>

//#define C_TWO_PI (6.2831853071795864769252867665590058L)
//#define C_PI (3.1415926535897932384626433832795029L)

typedef double _FPTYPE;

struct O2 {
    _FPTYPE x; /* X */
    _FPTYPE y; /* Y */
    _FPTYPE t; /* angle */
    _FPTYPE r; /* modulus */
};

#define O2prod(__a,__b) ((__a).x*(__b).x+(__a).y*(__b).y)
#define O2vprod(__a,__b) ((__a).x*(__b).y-(__a).y*(__b).x) //\a\\b\sin(theta_b -theta_a)
#define O2norm2(__a) O2prod((__a),(__a))
#define O2norm(__a) (sqrt(O2norm2((__a))))



#endif //CMT_02_H
