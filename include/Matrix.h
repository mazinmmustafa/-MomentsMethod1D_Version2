#ifndef MATRIX_H
#define MATRIX_H

// Includes
#include "myLib.h"

// Definitions
typedef struct Matrix Matrix;
struct Matrix{
    int rows, cols;
    complex double **data;
    int isAllocated;
};
#define DefaultMatrix {0, 0, NULL, 0}

// Functions
void allocateMatrix(Matrix *A, int rows, int cols);
void deallocateMatrix(Matrix *A);
void copyMatrix(const Matrix *A, Matrix *B);
void showMatrix(const Matrix *A);
void zerosMatrix(Matrix *A);
void onesMatrix(Matrix *A);
void eyeMatrix(Matrix *A);
void addMatrix(const Matrix *A, const Matrix *B, Matrix *C);
void subMatrix(const Matrix *A, const Matrix *B, Matrix *C);
void multMatrix(const Matrix *A, const Matrix *B, Matrix *C);
void LUP_Decompose(Matrix *A, int *P);
void LUP_Solve(const Matrix *A, const int *P, const Matrix *b, Matrix *x);
void LUP_Invert(Matrix *A, const int *P, Matrix *IA);
complex double LUP_Determinant(const Matrix *A, const int *P);

#endif