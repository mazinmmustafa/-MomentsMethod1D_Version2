//
#include "myLib.h"

double magVector(Vector A){
    return sqrt(A.x*A.x+A.y*A.y+A.z*A.z);
}

double dotVector(Vector A, Vector B){
    return A.x*B.x+A.y*B.y+A.z*B.z;
}

Vector crossVector(Vector A, Vector B){
    Vector C={A.y*B.z-A.z*B.y,
              A.z*B.x-A.x*B.z,
              A.x*B.y-A.y*B.x};
    return C;
}

Vector unitVector(Vector A){
    double a=magVector(A);
    Vector n={A.x/a, A.y/a, A.z/a};
    return n;
}

Vector addVector(Vector A, Vector B){
    Vector C={A.x+B.x, A.y+B.y, A.z+B.z};
    return C;
}

Vector subVector(Vector A, Vector B){
    Vector C={A.x-B.x, A.y-B.y, A.z-B.z};
    return C;
}

Vector scaleVector(Vector A, double a){
    Vector B={a*A.x, a*A.y, a*A.z};
    return B;
}

int isEqualVector(Vector A, Vector B){
    return ((A.x==B.x)&&(A.y==B.y)&&(A.z==B.z));
}