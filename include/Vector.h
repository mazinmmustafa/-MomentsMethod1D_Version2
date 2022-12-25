#ifndef VECTOR_H
#define VECTOR_H

// Includes
#include "myLib.h"

// Definitions
typedef struct Vector Vector;
struct Vector{
    double x, y, z;
};

// Functions
double magVector(Vector A);
double dotVector(Vector A, Vector B);
Vector crossVector(Vector A, Vector B);
Vector unitVector(Vector A);
Vector addVector(Vector A, Vector B);
Vector subVector(Vector A, Vector B);
Vector scaleVector(Vector A, double a);
int isEqualVector(Vector A, Vector B);

#endif