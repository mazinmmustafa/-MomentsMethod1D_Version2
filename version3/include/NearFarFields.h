#ifndef NEARFARFIELDS_H
#define NEARFARFIELDS_H

// Includes
#include "myLib.h"

// Definitions
typedef struct FieldComponents FieldComponents;
struct FieldComponents{
    complex double x, y, z;
};

typedef struct FarFieldComponents FarFieldComponents;
struct FarFieldComponents{
    complex double theta, phi;
};

// Functions
void NearFields(Shape *myShape, Matrix *In, double x, double y, double z, FieldComponents *E, 
    double tol, int warnings, int Nmax);
void FarFields(Shape *myShape, Matrix *In, double theta, double phi, FarFieldComponents *E);
double Prad(Shape *myShape, Matrix *In, double tol, int warnings, int Nmax);

#endif
