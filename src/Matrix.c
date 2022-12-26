//
#include "Matrix.h"

void allocateMatrix(Matrix *A, int rows, int cols){
    assert(A->isAllocated==0);
    assert(rows>0&&cols>0);
    A->rows = rows;
    A->cols = cols;
    A->data = (complex double**)malloc(rows*sizeof(complex double*));
    for (int i=0; i<rows; i++){
        A->data[i] = (complex double*)malloc(cols*sizeof(complex double));
    }
    A->isAllocated = 1;
    zerosMatrix(A);
}

void deallocateMatrix(Matrix *A){
    assert(A->isAllocated==1);
    for (int i=0; i<A->rows; i++){
        free(A->data[i]);
        A->data[i] = NULL;
    }
    free(A->data);
    A->data = NULL;
}

void copyMatrix(const Matrix *A, Matrix *B){
    assert(A->rows>0&&A->cols>0);
    assert(B->rows>0&&B->cols>0);
    assert(A->isAllocated==1);
    assert(B->isAllocated==1);
    assert(A->rows&&B->rows);
    assert(A->cols&&B->cols);
    for (int i=0; i<A->rows; i++){
        for (int j=0; j<A->cols; j++){
            B->data[i][j] = A->data[i][j];
        }
    }
}

void showMatrix(const Matrix *A){
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    for (int i=0; i<A->rows; i++){
        for (int j=0; j<A->cols; j++){
            printf("(%9.2E, %9.2E) ", creal(A->data[i][j]), cimag(A->data[i][j]));
        }
        printf("\n");
    }
}

void zerosMatrix(Matrix *A){
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    for (int i=0; i<A->rows; i++){
        for (int j=0; j<A->cols; j++){
            A->data[i][j] = 0.0E0;
        }
    }
}

void onesMatrix(Matrix *A){
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    for (int i=0; i<A->rows; i++){
        for (int j=0; j<A->cols; j++){
            A->data[i][j] = 1.0E0;
        }
    }
}

void eyeMatrix(Matrix *A){
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    for (int i=0; i<A->rows; i++){
        for (int j=0; j<A->cols; j++){
            A->data[i][j] = i==j ? 1.0E0 : 0.0E0;
        }
    }
}

void addMatrix(const Matrix *A, const Matrix *B, Matrix *C){
    // C=A+B
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    assert(B->rows>0&&B->cols>0);
    assert(B->isAllocated==1);
    assert(C->rows>0&&C->cols>0);
    assert(C->isAllocated==1);
    assert((A->rows==B->rows)&&(A->rows==C->rows));
    assert((A->cols==B->cols)&&(A->cols==C->cols));
    for (int i=0; i<A->rows; i++){
        for (int j=0; j<A->cols; j++){
            C->data[i][j] = A->data[i][j]+B->data[i][j];
        }
    }
}

void subMatrix(const Matrix *A, const Matrix *B, Matrix *C){
    // C=A-B
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    assert(B->rows>0&&B->cols>0);
    assert(B->isAllocated==1);
    assert(C->rows>0&&C->cols>0);
    assert(C->isAllocated==1);
    assert((A->rows==B->rows)&&(A->rows==C->rows));
    assert((A->cols==B->cols)&&(A->cols==C->cols));
    for (int i=0; i<A->rows; i++){
        for (int j=0; j<A->cols; j++){
            C->data[i][j] = A->data[i][j]-B->data[i][j];
        }
    }
}

void multMatrix(const Matrix *A, const Matrix *B, Matrix *C){
    // C=A*B
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    assert(B->rows>0&&B->cols>0);
    assert(B->isAllocated==1);
    assert(C->rows>0&&C->cols>0);
    assert(C->isAllocated==1);
    assert((A->rows==C->rows)&&(B->cols==C->cols)&&(A->cols==B->rows));
    complex double sum;
    for (int i=0; i<A->rows; i++){
        for (int j=0; j<B->cols; j++){
            sum = 0.0E0;
            for (int k=0; k<A->cols; k++){
                sum+=A->data[i][k]*B->data[k][j];
            }
            C->data[i][j] = sum;
        }
    }
}

void LUP_Decompose(Matrix *A, int *P){
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    assert(A->rows==A->cols);
    assert(P!=NULL);
    int N=A->rows;
    int imax, itemp; 
    complex double *ptr;
    double maxA, absA;
    for (int i=0; i<=N; i++){
        P[i]=i;
    }
    for (int i=0; i<N; i++){
        maxA=0.0E0;
        imax=i;
        for (int k=i; k<N; k++){
            if ((absA=cabs(A->data[k][i]))>maxA){ 
                maxA = absA;
                imax = k;
            }
        }
        if (maxA==0.0){
            printf("Warning: Matrix Is Degenerate!\n");
            break;
        }
        if (imax!=i) {
            itemp = P[i];
            P[i] = P[imax];
            P[imax] = itemp;
            ptr = A->data[i];
            A->data[i] = A->data[imax];
            A->data[imax] = ptr;
            P[N]++;
        }
        for (int j=i+1; j<N; j++){
            A->data[j][i]/=A->data[i][i];
            for (int k=i+1; k<N; k++){
                A->data[j][k]-=A->data[j][i]*A->data[i][k];
            }    
        }
    }
}

void LUP_Solve(const Matrix *A, const int *P, const Matrix *b, Matrix *x){
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    assert(b->rows>0&&b->cols>0);
    assert(b->isAllocated==1);
    assert(x->rows>0&&x->cols>0);
    assert(x->isAllocated==1);
    assert(A->rows==A->cols);
    assert(b->rows==A->rows);
    assert(x->rows==x->rows);
    assert(b->cols==1);
    assert(x->cols==1);
    assert(P!=NULL);
    int N=A->rows;
    for (int i=0; i<N; i++){
        x->data[i][0] = b->data[P[i]][0];
        for (int k=0; k<i; k++){
            x->data[i][0]-=A->data[i][k]*x->data[k][0];
        }
    }
    for (int i=N-1; i>=0; i--){
        for (int k=i+1; k<N; k++){
            x->data[i][0]-=A->data[i][k]*x->data[k][0];
        }
        x->data[i][0]/=A->data[i][i];
    }
}

void LUP_Invert(Matrix *A, const int *P, Matrix *IA){
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    assert(IA->rows>0&&IA->cols>0);
    assert(IA->isAllocated==1);
    assert(A->rows==A->cols);
    assert(IA->rows==IA->cols);
    assert(A->rows==IA->rows);
    assert(A->cols==IA->cols);
    assert(P!=NULL);
    int N=A->rows;
    for (int j=0; j<N; j++) {
        for (int i=0; i<N; i++) {
            IA->data[i][j]= P[i] == j ? 1.0 : 0.0;
            for (int k=0; k<i; k++){
                IA->data[i][j]-=A->data[i][k]*IA->data[k][j];
            }
        }
        for (int i=N-1; i>=0; i--) {
            for (int k=i+1; k<N; k++){
                IA->data[i][j]-=A->data[i][k]*IA->data[k][j];
            }
            IA->data[i][j]/=A->data[i][i];
        }
    }
}

complex double LUP_Determinant(const Matrix *A, const int *P){
    assert(A->rows>0&&A->cols>0);
    assert(A->isAllocated==1);
    assert(A->rows==A->cols);
    assert(P!=NULL);
    int N=A->rows;
    complex double det=A->data[0][0];
    for (int i=1; i<N; i++){
        det*=A->data[i][i];
    }
    return (P[N]-N)%2 == 0 ? det : -det;
}
