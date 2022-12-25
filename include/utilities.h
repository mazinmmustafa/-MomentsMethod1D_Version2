#ifndef UTILITIES_H
#define UTILITIES_H

// Include
#include "myLib.h"

// Definitions
typedef struct Timer Timer;

#ifdef _WIN32
struct Timer{
    time_t start, stop;
    double elapsed;
    int isSet;
};
#endif

#ifdef linux
struct Timer{
    struct timespec start, stop;
    double elapsed;
    int isSet;
};
#endif

// Functions
void showComplex(complex double z);
void progressBar(int n, int N, char *message);
void setTimer(Timer *T);
void unsetTimer(Timer *T);
double deg2rad(double x);
double rad2deg(double x);
void setRandomSeed();
int randInt(int a, int b);
double randDouble(double a, double b);

#endif