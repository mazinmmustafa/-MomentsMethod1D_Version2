#ifndef MM_ENGINE_H
#define MM_ENGINE_H

// Includes
#include "myLib.h"
#include "Matrix.h"
#include "Shape.h"

// Definitions
typedef struct integrandArgs integrandArgs;
struct integrandArgs{
    Basis basis_m, basis_n;    
    Vector Lm_p, Lm_m, Ln_p, Ln_m;
    double a;
};

// Functions
void settings(int warningsIn, int NmaxIn);
void MM_deltaGap(Shape *myShape, double tol, Matrix *Zmn, Matrix *Vm, int isWarnings, int Nmax);
complex double I1(double L, double a, double tol);
complex double I2(double L, double a, double tol);
complex double I3(double L, double a, double tol);
//
complex double phi_mn_pp(void *args, double tol);
complex double phi_mn_pm(void *args, double tol);
complex double phi_mn_mp(void *args, double tol);
complex double phi_mn_mm(void *args, double tol);
complex double psi_mn_pp(void *args, double tol);
complex double psi_mn_pm(void *args, double tol);
complex double psi_mn_mp(void *args, double tol);
complex double psi_mn_mm(void *args, double tol);

#endif
