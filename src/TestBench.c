//
#include "TestBench.h"
#include "Matrix.h"
#include "Vector.h"
#include "Shape.h"
#include "MM_Engine.h"
#include "NearFarFields.h"
#include "QuadL.h"

void TestConstants(){

    double I=pi;
    printf("%21.14E\n", I);
    printf("%21.14E\n", c0);
    printf("%21.14E\n", mu0);
    printf("%21.14E\n", eps0);
    printf("%21.14E\n", eta0);

}

void TestMatrix(){
    
    // Test Basics
    int N=3;
    Matrix A=DefaultMatrix;
    Matrix b=DefaultMatrix;
    Matrix x=DefaultMatrix;
    Matrix IA=DefaultMatrix;
    allocateMatrix(&A, N, N);
    allocateMatrix(&b, N, 1);
    allocateMatrix(&x, N, 1);
    allocateMatrix(&IA, N, N);

    A.data[0][0] = 0.0; A.data[0][1] = 1.0; A.data[0][2] = 1.0;
    A.data[1][0] = 1.0; A.data[1][1] = 1.0; A.data[1][2] = 1.0; 
    A.data[2][0] = 2.0; A.data[2][1] = 1.0; A.data[2][2] = 0.0; 

    b.data[0][0] = 3.0;
    b.data[1][0] = 1.0;
    b.data[2][0] = 2.0;

    int *P=(int*)malloc((N+1)*sizeof(int));

    LUP_Decompose(&A, P);
    showMatrix(&A);
    for (int i=0; i<=N; i++){
        printf("%d\n", P[i]);
    }

    LUP_Solve(&A, P, &b, &x);
    showMatrix(&x);

    LUP_Invert(&A, P, &IA);
    showMatrix(&IA);

    double complex detA=LUP_Determinant(&A, P);
    printf("(%21.14E, %21.14E)\n", creal(detA), cimag(detA));

    free(P);
    P = NULL;
    deallocateMatrix(&A);
    deallocateMatrix(&b);
    deallocateMatrix(&x);
    deallocateMatrix(&IA);
    
}

void TestUtilities(){

    complex double j=csqrt(-1.0);
    complex double z=1.2-j*2.4;
    showComplex(z);

    int N=1000;
    Timer T;
    setTimer(&T);
    for (int i=0; i<N; i++){
        usleep(6000);
        progressBar(i, N, "Computing...");
    }
    unsetTimer(&T);

    time_t t;
    srand((unsigned)time(&t));

    setRandomSeed();

    int Ns=10;
    int a=4, b=-1;
    for (int i=0; i<Ns; i++){
        printf("%d, ", randInt(a, b));
    }
    printf("\n");

    Ns = 10000;
    double a_=-1.0, b_=+1.0;
    FILE *file=fopen("Data/TestBench/Random/data.dat", "w");
    assert(file!=NULL);
    for (int i=0; i<Ns; i++){
        fprintf(file, "%21.14E %21.14E\n", randDouble(a_, b_), randDouble(a_, b_));
    }
    fclose(file);
    
}


typedef struct func1_Args func1_Args;
struct func1_Args{
    double a;
};

typedef struct func2_Args func2_Args;
struct func2_Args{
    double a, b;
};

complex double func1(complex double x, void *Args){
    func1_Args *myArgs=(func1_Args*)Args;
    double a=myArgs->a;
    return ccos(a*x);
}

complex double func2(complex double x, complex double y, void *Args){
    func2_Args *myArgs=(func2_Args*)Args;
    double a=myArgs->a;
    double b=myArgs->b;
    return ccos(a*x+cexp(-b*y));
}

void TestQuadL(){

    double a1=0.1;
    double b1=0.6;
    double a2=0.3;
    double b2=0.4;

    double tol=1.0E-14;
    double a=2.0;
    double b=4.7;

    complex double ans;

    func1_Args myArgs1={a};
    ans = QuadL_1D(func1, &myArgs1, a1, b1, tol, 1, 2048);
    showComplex(ans);

    func2_Args myArgs2={a, b};
    ans = QuadL_2D(func2, &myArgs2, a2, b2, a1, b1, tol, 1, 2048);
    showComplex(ans);

}

void TestSingularIntegrands(){
    double tol=1.0E-4;
    double L=0.01;
    double a=1.0E-4;
    complex double ans;
    ans = I1(L, a, tol); showComplex(ans);
    ans = I2(L, a, tol); showComplex(ans);
    ans = I3(L, a, tol); showComplex(ans);
}

void TestBasicSolution(){
    int warnings=1;
    int Nmax=32;
    double tol=1.0E-3;
    double lambda=1.0;
    double L=0.1*lambda;
    double a=5.0E-3*lambda;
    double segments_per_lambda=5.0;
    double dL=L/segments_per_lambda;
    // Create Shape
    createVerticalDipoleDeltaGapCenter(L, dL);
    Shape myShape=DefaultShape;
    getShape(&myShape, 1, lambda, a);
    Port newPort={1, 1.0, 0.0};
    setPort(&myShape, newPort);
    logShape(&myShape);
    int N=myShape.NBasis;
    Matrix Zmn=DefaultMatrix;
    Matrix Vm=DefaultMatrix;
    Matrix In=DefaultMatrix;
    allocateMatrix(&Zmn, N, N);
    allocateMatrix(&Vm, N, 1);
    allocateMatrix(&In, N, 1);
    MM_deltaGap(&myShape, tol, &Zmn, &Vm, warnings, Nmax);
    FILE *file=fopen("Data/TestBench/Singular/Zmn.dat", "w");
    printf("N = %d\n", myShape.NBasis);
    for (int m=0; m<N; m++){
        for (int n=0; n<N; n++){
			printf("Term (%3d, %3d) = %21.14E < %21.14E\n", m, n, cabs(Zmn.data[m][n]), rad2deg(carg(Zmn.data[m][n])));
            fprintf(file, "Term (%3d, %3d) = (%21.14E, %21.14E)\n", m, n, creal(Zmn.data[m][n]), cimag(Zmn.data[m][n]));
        }
    }
    fclose(file);
    int *P=(int*)malloc((N+1)*sizeof(int));
    LUP_Decompose(&Zmn, P);
    LUP_Solve(&Zmn, P, &Vm, &In);
    free(P);
    P = NULL;
	for (int n=0; n<N; n++){
		printf("%21.14E < %21.14E\n", cabs(In.data[n][0]), rad2deg(carg(In.data[n][0])));
	}
    complex double Zin;
    for (int n=0; n<N; n++){
        for (int i=0; i<myShape.NPorts; i++){
            if ((myShape.basisList[n].portNumber==myShape.portsList[i].portNumber)&&
                (myShape.basisList[n].isPort)){
                Zin = (myShape.portsList[i].Vin/In.data[n][0])-myShape.portsList[i].ZL;
                showComplex(Zin);
            }
        }
    }
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
}

void TestNearTerms(){
    settings(1, 256);
    double tol=1.0E-4;
    double L=0.0001;
    double theta=deg2rad(90.0);
    double a=1.0E-4;
    Vector rm_m={0.0, 0.0, 0.0};
    Vector rm_n={L, 0.0, 0.0};
    Vector rm_p={2.0*L, 0.0, 0.0};
    Vector Lm_m=subVector(rm_n, rm_m);
    Vector Lm_p=subVector(rm_p, rm_n);
    Basis basis_m={rm_m, rm_n, rm_p, Lm_m, Lm_p, magVector(Lm_m), magVector(Lm_p), 0, 0};
    Vector rn_m={2.0*L*cos(theta), 2.0*L*sin(theta), 0.0};
    Vector rn_n={L*cos(theta), L*sin(theta), 0.0};
    Vector rn_p={0.0, 0.0, 0.0};
    Vector Ln_m=subVector(rn_n, rn_m);
    Vector Ln_p=subVector(rn_p, rn_n);
    Basis basis_n={rn_m, rn_n, rn_p, Ln_m, Ln_p, magVector(Ln_m), magVector(Ln_p), 0, 0};
    integrandArgs args={basis_m, basis_n, Lm_p, Lm_m, Ln_p, Ln_m, a};
    complex double ans;

    ans = phi_mn_pp(&args, tol);
    showComplex(ans);

    ans = phi_mn_pm(&args, tol);
    showComplex(ans);
    
    ans = phi_mn_mp(&args, tol);
    showComplex(ans);

    ans = phi_mn_mm(&args, tol);
    showComplex(ans);

    ans = psi_mn_pp(&args, tol);
    showComplex(ans);

    ans = psi_mn_pm(&args, tol);
    showComplex(ans);
    
    ans = psi_mn_mp(&args, tol);
    showComplex(ans);

    ans = psi_mn_mm(&args, tol);
    showComplex(ans);

}


void TestShapeProcessing(){
    double lambda=1.0;
    double L=0.5*lambda;
    double a=1.0E-3*lambda;
    double segments_per_lambda=11.0;
    double dL=L/segments_per_lambda;
    // Create Shape
    createVerticalDipoleDeltaGapCenter(L, dL);
    Shape myShape=DefaultShape;
    getShape(&myShape, 1, lambda, a);
    Port newPort={1, 1.0, 0.0};
    setPort(&myShape, newPort);

    logShape(&myShape);

    deleteShape(&myShape);
}

void TestVerticalDipoleDeltaGapCenter(){
    int warnings=1;
    int Nmax=256;
    int Ns=201;
    double tol=1.0E-4;
    double lambda=1.0;
    double L_min=1.0E-3*lambda;
    double L_max=2.0*lambda;
    double segments_per_lambda=31.0;
    double dL=(L_max-L_min)/(Ns-1.0);
    double lc;
    double L, a;
    //
    Timer T;
    setTimer(&T);
    FILE *file=fopen("Data/TestBench/DipoleAdmittance/Yin.dat", "w");
    assert(file!=NULL);
    for (int i=0; i<Ns; i++){
        L = L_min+i*dL;
        a = L/(2.0*74.2);
        lc = L/segments_per_lambda;
        // Create Shape
        if (i<Ns-1){progressBar(i, Ns, "Creating Shape ...");}
        createVerticalDipoleDeltaGapCenter(L, lc);
        Shape myShape=DefaultShape;
        getShape(&myShape, 1, lambda, a);
        Port newPort={1, 1.0, 0.0};
        setPort(&myShape, newPort);
        // Allocate Matricies
        int N=myShape.NBasis;
        Matrix Zmn=DefaultMatrix;
        Matrix Vm=DefaultMatrix;
        Matrix In=DefaultMatrix;
        allocateMatrix(&Zmn, N, N);
        allocateMatrix(&Vm, N, 1);
        allocateMatrix(&In, N, 1);
        // Solve MM
        if (i<Ns-1){progressBar(i, Ns, "Builing Matricies ...");}
        MM_deltaGap(&myShape, tol, &Zmn, &Vm, warnings, Nmax);
        if (i<Ns-1){progressBar(i, Ns, "Solving Linear System ...");}
        int *P=(int*)malloc((N+1)*sizeof(int));
        LUP_Decompose(&Zmn, P);
        LUP_Solve(&Zmn, P, &Vm, &In);
        free(P);
        P = NULL;
        // Get Results
        if (i<Ns-1){progressBar(i, Ns, "Writing Results ...");}
        complex double Zin;
        for (int n=0; n<N; n++){
            for (int i=0; i<myShape.NPorts; i++){
                if ((myShape.basisList[n].portNumber==myShape.portsList[i].portNumber)&&
                    (myShape.basisList[n].isPort)){
                    Zin = (myShape.portsList[i].Vin/In.data[n][0])-myShape.portsList[i].ZL;
                    fprintf(file, "%21.14E %21.14E %21.14E\n", L, creal(1.0/Zin), cimag(1.0/Zin));
                }
            }
        }
        // Free The Memory
        progressBar(i, Ns, "Free The Memory ...");
        deallocateMatrix(&Zmn);
        deallocateMatrix(&Vm);
        deallocateMatrix(&In);
        deleteShape(&myShape);
    }
    fclose(file);
    unsetTimer(&T);
}

void TestCircularLoopDeltaGapCenter(){
    int warnings=1;
    int Nmax=256;
    int Ns=601;
    double tol=1.0E-4;
    double lambda=1.0;
    double a=1.0E-4;
    double S_min=1.0E-3*lambda;
    double S_max=2.5*lambda;
    double segments_per_lambda=31.0;
    double dS=(S_max-S_min)/(Ns-1.0);
    double lc;
    double S, r;
    //
    Timer T;
    setTimer(&T);
    FILE *file=fopen("Data/TestBench/CircularLoopAdmittance/Zin.dat", "w");
    assert(file!=NULL);
    for (int i=0; i<Ns; i++){
        S = S_min+i*dS;
        r = S/(2.0*pi);
        lc = S/segments_per_lambda;
        // Create Shape
        if (i<Ns-1){progressBar(i, Ns, "Creating Shape ...");}
        createCircularLoopDeltaGapCenter(r, lc);
        Shape myShape=DefaultShape;
        getShape(&myShape, 1, lambda, a);
        Port newPort={1, 1.0, 0.0};
        setPort(&myShape, newPort);
        // Allocate Matricies
        int N=myShape.NBasis;
        Matrix Zmn=DefaultMatrix;
        Matrix Vm=DefaultMatrix;
        Matrix In=DefaultMatrix;
        allocateMatrix(&Zmn, N, N);
        allocateMatrix(&Vm, N, 1);
        allocateMatrix(&In, N, 1);
        // Solve MM
        if (i<Ns-1){progressBar(i, Ns, "Builing Matricies ...");}
        MM_deltaGap(&myShape, tol, &Zmn, &Vm, warnings, Nmax);
        if (i<Ns-1){progressBar(i, Ns, "Solving Linear System ...");}
        int *P=(int*)malloc((N+1)*sizeof(int));
        LUP_Decompose(&Zmn, P);
        LUP_Solve(&Zmn, P, &Vm, &In);
        free(P);
        P = NULL;
        // Get Results
        if (i<Ns-1){progressBar(i, Ns, "Writing Results ...");}
        complex double Zin;
        for (int n=0; n<N; n++){
            for (int i=0; i<myShape.NPorts; i++){
                if ((myShape.basisList[n].portNumber==myShape.portsList[i].portNumber)&&
                    (myShape.basisList[n].isPort)){
                    Zin = (myShape.portsList[i].Vin/In.data[n][0])-myShape.portsList[i].ZL;
                    fprintf(file, "%21.14E %21.14E %21.14E\n", S, creal(Zin), cimag(Zin));
                }
            }
        }
        // Free The Memory
        progressBar(i, Ns, "Free The Memory ...");
        deallocateMatrix(&Zmn);
        deallocateMatrix(&Vm);
        deallocateMatrix(&In);
        deleteShape(&myShape);
    }
    fclose(file);
    unsetTimer(&T);
}