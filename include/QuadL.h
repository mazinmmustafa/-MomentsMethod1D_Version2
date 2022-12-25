#ifndef QUADL_H
#define QUADL_H

// Includes 
#include "myLib.h"

// Functions
void QuadL(int N, double *x, double *w);
void printQuadLRule(int Nmax);
complex double QuadL_1D(complex double func(complex double, void*), 
	void *args, double a, double b, double tol, int warnings, int Nmax);
complex double QuadL_2D(complex double func(complex double, complex double, void*), 
	void *args, double a2, double b2, double a1, double b1, double tol, int warnings, int Nmax);

#endif