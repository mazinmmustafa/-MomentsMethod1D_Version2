#ifndef SHAPE_H
#define SHAPE_H

// Includes
#include "myLib.h"
#include "Vector.h"

// Definitions
typedef struct Basis Basis;
struct Basis{
    Vector rm, rn, rp;
    Vector Lm, Lp;
    double lm, lp;
    int isPort, portNumber;
};

typedef struct Port Port;
struct Port{
    int portNumber;
    complex double Vin, ZL;
};

typedef struct Shape Shape;
struct Shape{
    double lambda;
    int NBasis, NPorts;
    Basis *basisList;
    Port *portsList;
    double a;
};
#define DefaultShape {1.0, 0, 0, NULL, NULL, 0.0}

typedef struct YagiElements YagiElements;
struct YagiElements{
    int NElements;
    double *length;
    double *xLocation;
    int *isExcitaiton;
};

// Functions
void getShape(Shape *myShape, int NPorts, double lambda, double a);
void setPort(Shape *myShape, Port newPort);
void deleteShape(Shape *myShape);
void logShape(Shape *myShape);
void createVerticalDipoleDeltaGapCenter(double L, double lc);
void createCircularLoopDeltaGapCenter(double r, double lc);
void createTransmissionLineDeltaGapCenter(double L, double H, double lc);
void createYagiAntenna(YagiElements *elements, double segments_per_lambda);

#endif