//
#include "programs.h"
#include "MM_Engine.h"
#include "NearFarFields.h"

void VerticalDipoleNearFields(){
    int warnings=1;
    int Nmax=256;
    int Ns=100;
    double tol=1.0E-4;
    double lambda=1.0;
    double L=0.47*lambda;
    double a=5.0E-3*lambda;
    double segments_per_lambda=31.0;
    double lc=L/segments_per_lambda;

    double x_min=-6.0*lambda;
    double x_max=+6.0*lambda;
    double z_min=-6.0*lambda;
    double z_max=+6.0*lambda;
    double y=0.0;

    double dx, dz, x, z;
    dx = (x_max-x_min)/(Ns-1.0);
    dz = (z_max-z_min)/(Ns-1.0);

    // Create Shape
    createVerticalDipoleDeltaGapCenter(L, lc);
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
    int *P=(int*)malloc((N+1)*sizeof(int));
    LUP_Decompose(&Zmn, P);
    LUP_Solve(&Zmn, P, &Vm, &In);
    free(P);
    P = NULL;

    FILE *fileX=fopen("Data/TestBench/NearField/DataX.dat", "w");
    FILE *fileY=fopen("Data/TestBench/NearField/DataY.dat", "w");
    FILE *fileZ=fopen("Data/TestBench/NearField/DataZ.dat", "w");
    FILE *fileDataXr=fopen("Data/TestBench/NearField/DataXr.dat", "w");
    FILE *fileDataXi=fopen("Data/TestBench/NearField/DataXi.dat", "w");
    FILE *fileDataYr=fopen("Data/TestBench/NearField/DataYr.dat", "w");
    FILE *fileDataYi=fopen("Data/TestBench/NearField/DataYi.dat", "w");
    FILE *fileDataZr=fopen("Data/TestBench/NearField/DataZr.dat", "w");
    FILE *fileDataZi=fopen("Data/TestBench/NearField/DataZi.dat", "w");
    assert(fileX!=NULL);
    assert(fileY!=NULL);
    assert(fileZ!=NULL);
    assert(fileDataXr!=NULL);
    assert(fileDataXi!=NULL);
    assert(fileDataYr!=NULL);
    assert(fileDataYi!=NULL);
    assert(fileDataZr!=NULL);
    assert(fileDataZi!=NULL);
    FieldComponents E;
    Timer T;
    setTimer(&T);
    for (int i=0; i<Ns; i++){
        x = x_min+i*dx;
        for (int j=0; j<Ns; j++){
            z = z_min+j*dz;
            NearFields(&myShape, &In, x, y, z, &E, tol, warnings, Nmax);
            fprintf(fileX, "%21.14E ", x);
            fprintf(fileY, "%21.14E ", y);
            fprintf(fileZ, "%21.14E ", z);
            fprintf(fileDataXr, "%21.14E ", creal(E.x));
            fprintf(fileDataXi, "%21.14E ", cimag(E.x));
            fprintf(fileDataYr, "%21.14E ", creal(E.y));
            fprintf(fileDataYi, "%21.14E ", cimag(E.y));
            fprintf(fileDataZr, "%21.14E ", creal(E.z));
            fprintf(fileDataZi, "%21.14E ", cimag(E.z));
        }
        fprintf(fileX, "\n");
        fprintf(fileY, "\n");
        fprintf(fileZ, "\n");
        fprintf(fileDataXr, "\n");
        fprintf(fileDataXi, "\n");
        fprintf(fileDataYr, "\n");
        fprintf(fileDataYi, "\n");
        fprintf(fileDataZr, "\n");
        fprintf(fileDataZi, "\n");
        progressBar(i, Ns, "Computing Near Fields ...");
    }
    unsetTimer(&T);
    fclose(fileX);
    fclose(fileY);
    fclose(fileZ);
    fclose(fileDataXr);
    fclose(fileDataXi);
    fclose(fileDataYr);
    fclose(fileDataYi);
    fclose(fileDataZr);
    fclose(fileDataZi);
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
}

void TransmissionLineNearFields(){
    int warnings=1;
    int Nmax=256;
    int Nsx=4000;
    int Nsz=2000;
    double tol=1.0E-4;
    double GHz=1.0E9;
    double mm=1.0E-3;
    double cm=1.0E-2;

    double freq=2.0*GHz;
    double L=120.0*cm;
    double H=10.0*cm;
    double a=0.8*mm;
    double Z0=50.0;
    double segments_per_lambda=31.0;

    double lambda=c0/freq;
    double lc=L/segments_per_lambda;

    double x_min=-200.0*cm;
    double x_max=+200.0*cm;
    double z_min=-0.0*cm;
    double z_max=+200.0*cm;
    double y=0.0;

    double dx, dz, x, z;
    dx = (x_max-x_min)/(Nsx-1.0);
    dz = (z_max-z_min)/(Nsz-1.0);

    // Create Shape
    createTransmissionLineDeltaGapCenter(L, H, lc);
    Shape myShape=DefaultShape;
    getShape(&myShape, 2, lambda, a);
    Port newPort1={1, 1.0, Z0};
    Port newPort2={2, 0.0, Z0};
    setPort(&myShape, newPort1);
    setPort(&myShape, newPort2);
    logShape(&myShape);
    int N=myShape.NBasis;
    Matrix Zmn=DefaultMatrix;
    Matrix Vm=DefaultMatrix;
    Matrix In=DefaultMatrix;
    allocateMatrix(&Zmn, N, N);
    allocateMatrix(&Vm, N, 1);
    allocateMatrix(&In, N, 1);
    MM_deltaGap(&myShape, tol, &Zmn, &Vm, warnings, Nmax);
    int *P=(int*)malloc((N+1)*sizeof(int));
    LUP_Decompose(&Zmn, P);
    LUP_Solve(&Zmn, P, &Vm, &In);
    free(P);
    P = NULL;

    FILE *fileX=fopen("Data/TestBench/NearField/DataX.dat", "w");
    FILE *fileY=fopen("Data/TestBench/NearField/DataY.dat", "w");
    FILE *fileZ=fopen("Data/TestBench/NearField/DataZ.dat", "w");
    FILE *fileDataXr=fopen("Data/TestBench/NearField/DataXr.dat", "w");
    FILE *fileDataXi=fopen("Data/TestBench/NearField/DataXi.dat", "w");
    FILE *fileDataYr=fopen("Data/TestBench/NearField/DataYr.dat", "w");
    FILE *fileDataYi=fopen("Data/TestBench/NearField/DataYi.dat", "w");
    FILE *fileDataZr=fopen("Data/TestBench/NearField/DataZr.dat", "w");
    FILE *fileDataZi=fopen("Data/TestBench/NearField/DataZi.dat", "w");
    assert(fileX!=NULL);
    assert(fileY!=NULL);
    assert(fileZ!=NULL);
    assert(fileDataXr!=NULL);
    assert(fileDataXi!=NULL);
    assert(fileDataYr!=NULL);
    assert(fileDataYi!=NULL);
    assert(fileDataZr!=NULL);
    assert(fileDataZi!=NULL);
    Timer T;
    setTimer(&T);
    FieldComponents E;
    for (int i=0; i<Nsx; i++){
        x = x_min+i*dx;
        for (int j=0; j<Nsz; j++){
            z = z_min+j*dz;
            NearFields(&myShape, &In, x, y, z, &E, tol, warnings, Nmax);
            fprintf(fileX, "%21.14E ", x);
            fprintf(fileY, "%21.14E ", y);
            fprintf(fileZ, "%21.14E ", z);
            fprintf(fileDataXr, "%21.14E ", creal(E.x));
            fprintf(fileDataXi, "%21.14E ", cimag(E.x));
            fprintf(fileDataYr, "%21.14E ", creal(E.y));
            fprintf(fileDataYi, "%21.14E ", cimag(E.y));
            fprintf(fileDataZr, "%21.14E ", creal(E.z));
            fprintf(fileDataZi, "%21.14E ", cimag(E.z));
        }
        fprintf(fileX, "\n");
        fprintf(fileY, "\n");
        fprintf(fileZ, "\n");
        fprintf(fileDataXr, "\n");
        fprintf(fileDataXi, "\n");
        fprintf(fileDataYr, "\n");
        fprintf(fileDataYi, "\n");
        fprintf(fileDataZr, "\n");
        fprintf(fileDataZi, "\n");
        progressBar(i, Nsx, "Computing Near Fields ...");
    }
    unsetTimer(&T);
    fclose(fileX);
    fclose(fileY);
    fclose(fileZ);
    fclose(fileDataXr);
    fclose(fileDataXi);
    fclose(fileDataYr);
    fclose(fileDataYi);
    fclose(fileDataZr);
    fclose(fileDataZi);
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
}

void VerticalDipoleNearFieldsCut(){
    int warnings=1;
    int Nmax=256;
    int Ns=1001;
    double tol=1.0E-4;
    double lambda=1.0;
    double L=0.47*lambda;
    double a=5.0E-3*lambda;
    double segments_per_lambda=31.0;
    double lc=L/segments_per_lambda;

    double x=+4.0*lambda;
    double y=+4.0*lambda;
    double z_min=-10.0*lambda;
    double z_max=+10.0*lambda;
    
    double dz, z;
    dz = (z_max-z_min)/(Ns-1.0);

    // Create Shape
    createVerticalDipoleDeltaGapCenter(L, lc);
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
    int *P=(int*)malloc((N+1)*sizeof(int));
    LUP_Decompose(&Zmn, P);
    LUP_Solve(&Zmn, P, &Vm, &In);
    free(P);
    P = NULL;

    FILE *fileX=fopen("Data/TestBench/NearFieldCut/DataX.dat", "w");
    FILE *fileY=fopen("Data/TestBench/NearFieldCut/DataY.dat", "w");
    FILE *fileZ=fopen("Data/TestBench/NearFieldCut/DataZ.dat", "w");
    FILE *fileDataXr=fopen("Data/TestBench/NearFieldCut/DataXr.dat", "w");
    FILE *fileDataXi=fopen("Data/TestBench/NearFieldCut/DataXi.dat", "w");
    FILE *fileDataYr=fopen("Data/TestBench/NearFieldCut/DataYr.dat", "w");
    FILE *fileDataYi=fopen("Data/TestBench/NearFieldCut/DataYi.dat", "w");
    FILE *fileDataZr=fopen("Data/TestBench/NearFieldCut/DataZr.dat", "w");
    FILE *fileDataZi=fopen("Data/TestBench/NearFieldCut/DataZi.dat", "w");
    assert(fileX!=NULL);
    assert(fileY!=NULL);
    assert(fileZ!=NULL);
    assert(fileDataXr!=NULL);
    assert(fileDataXi!=NULL);
    assert(fileDataYr!=NULL);
    assert(fileDataYi!=NULL);
    assert(fileDataZr!=NULL);
    assert(fileDataZi!=NULL);
    FieldComponents E;
    Timer T;
    setTimer(&T);
    
    for (int i=0; i<Ns; i++){
        z = z_min+i*dz;
        NearFields(&myShape, &In, x, y, z, &E, tol, warnings, Nmax);
        fprintf(fileX, "%21.14E \n", x);
        fprintf(fileY, "%21.14E \n", y);
        fprintf(fileZ, "%21.14E \n", z);
        fprintf(fileDataXr, "%21.14E \n", creal(E.x));
        fprintf(fileDataXi, "%21.14E \n", cimag(E.x));
        fprintf(fileDataYr, "%21.14E \n", creal(E.y));
        fprintf(fileDataYi, "%21.14E \n", cimag(E.y));
        fprintf(fileDataZr, "%21.14E \n", creal(E.z));
        fprintf(fileDataZi, "%21.14E \n", cimag(E.z));
        progressBar(i, Ns, "Computing Near Fields ...");
    }
       
    unsetTimer(&T);
    fclose(fileX);
    fclose(fileY);
    fclose(fileZ);
    fclose(fileDataXr);
    fclose(fileDataXi);
    fclose(fileDataYr);
    fclose(fileDataYi);
    fclose(fileDataZr);
    fclose(fileDataZi);
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
}

void CircularLoopNearFieldsCut(){
    int warnings=1;
    int Nmax=256;
    int Ns=1001;
    double tol=1.0E-4;
    double lambda=1.0;
    double r=0.2*lambda;
    double a=5.0E-3*lambda;
    double segments_per_lambda=31.0;
    double lc=(2.0*pi*r)/segments_per_lambda;

    double x=+4.0*lambda;
    double y=+4.0*lambda;
    double z_min=-10.0*lambda;
    double z_max=+10.0*lambda;
    
    double dz, z;
    dz = (z_max-z_min)/(Ns-1.0);

    // Create Shape
    createCircularLoopDeltaGapCenter(r, lc);
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
    int *P=(int*)malloc((N+1)*sizeof(int));
    LUP_Decompose(&Zmn, P);
    LUP_Solve(&Zmn, P, &Vm, &In);
    free(P);
    P = NULL;

    FILE *fileX=fopen("Data/TestBench/NearFieldCut/DataX.dat", "w");
    FILE *fileY=fopen("Data/TestBench/NearFieldCut/DataY.dat", "w");
    FILE *fileZ=fopen("Data/TestBench/NearFieldCut/DataZ.dat", "w");
    FILE *fileDataXr=fopen("Data/TestBench/NearFieldCut/DataXr.dat", "w");
    FILE *fileDataXi=fopen("Data/TestBench/NearFieldCut/DataXi.dat", "w");
    FILE *fileDataYr=fopen("Data/TestBench/NearFieldCut/DataYr.dat", "w");
    FILE *fileDataYi=fopen("Data/TestBench/NearFieldCut/DataYi.dat", "w");
    FILE *fileDataZr=fopen("Data/TestBench/NearFieldCut/DataZr.dat", "w");
    FILE *fileDataZi=fopen("Data/TestBench/NearFieldCut/DataZi.dat", "w");
    assert(fileX!=NULL);
    assert(fileY!=NULL);
    assert(fileZ!=NULL);
    assert(fileDataXr!=NULL);
    assert(fileDataXi!=NULL);
    assert(fileDataYr!=NULL);
    assert(fileDataYi!=NULL);
    assert(fileDataZr!=NULL);
    assert(fileDataZi!=NULL);
    FieldComponents E;
    Timer T;
    setTimer(&T);
    
    for (int i=0; i<Ns; i++){
        z = z_min+i*dz;
        NearFields(&myShape, &In, x, y, z, &E, tol, warnings, Nmax);
        fprintf(fileX, "%21.14E \n", x);
        fprintf(fileY, "%21.14E \n", y);
        fprintf(fileZ, "%21.14E \n", z);
        fprintf(fileDataXr, "%21.14E \n", creal(E.x));
        fprintf(fileDataXi, "%21.14E \n", cimag(E.x));
        fprintf(fileDataYr, "%21.14E \n", creal(E.y));
        fprintf(fileDataYi, "%21.14E \n", cimag(E.y));
        fprintf(fileDataZr, "%21.14E \n", creal(E.z));
        fprintf(fileDataZi, "%21.14E \n", cimag(E.z));
        progressBar(i, Ns, "Computing Near Fields ...");
    }
       
    unsetTimer(&T);
    fclose(fileX);
    fclose(fileY);
    fclose(fileZ);
    fclose(fileDataXr);
    fclose(fileDataXi);
    fclose(fileDataYr);
    fclose(fileDataYi);
    fclose(fileDataZr);
    fclose(fileDataZi);
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
}

void TransmissionLineNearFieldsCut(){
    int warnings=1;
    int Nmax=256;
    int Ns=1001;
    double tol=1.0E-4;
    double GHz=1.0E9;
    double mm=1.0E-3;
    double cm=1.0E-2;

    double freq=1.0*GHz;
    double L=30.0*cm;
    double H=1.0*cm;
    double a=0.8*mm;
    double Z0=50.0;
	complex double Vin=sqrt(8*Z0);
    double segments_per_lambda=31.0;

    double lambda=c0/freq;
    double lc=L/segments_per_lambda;

    double x=15.0*cm;
    double y=10.0*cm;
    double z_min=-30.0*cm;
    double z_max=+30.0*cm;

    double dz, z;
    dz = (z_max-z_min)/(Ns-1.0);

    // Create Shape
    createTransmissionLineDeltaGapCenter(L, H, lc);
    Shape myShape=DefaultShape;
    getShape(&myShape, 1, lambda, a);
    Port newPort1={1, Vin, Z0};
    setPort(&myShape, newPort1);
    logShape(&myShape);
    int N=myShape.NBasis;
    Matrix Zmn=DefaultMatrix;
    Matrix Vm=DefaultMatrix;
    Matrix In=DefaultMatrix;
    allocateMatrix(&Zmn, N, N);
    allocateMatrix(&Vm, N, 1);
    allocateMatrix(&In, N, 1);
    MM_deltaGap(&myShape, tol, &Zmn, &Vm, warnings, Nmax);
    int *P=(int*)malloc((N+1)*sizeof(int));
    LUP_Decompose(&Zmn, P);
    LUP_Solve(&Zmn, P, &Vm, &In);
    free(P);
    P = NULL;

    FILE *fileX=fopen("Data/TestBench/NearFieldCut/DataX.dat", "w");
    FILE *fileY=fopen("Data/TestBench/NearFieldCut/DataY.dat", "w");
    FILE *fileZ=fopen("Data/TestBench/NearFieldCut/DataZ.dat", "w");
    FILE *fileDataXr=fopen("Data/TestBench/NearFieldCut/DataXr.dat", "w");
    FILE *fileDataXi=fopen("Data/TestBench/NearFieldCut/DataXi.dat", "w");
    FILE *fileDataYr=fopen("Data/TestBench/NearFieldCut/DataYr.dat", "w");
    FILE *fileDataYi=fopen("Data/TestBench/NearFieldCut/DataYi.dat", "w");
    FILE *fileDataZr=fopen("Data/TestBench/NearFieldCut/DataZr.dat", "w");
    FILE *fileDataZi=fopen("Data/TestBench/NearFieldCut/DataZi.dat", "w");
    assert(fileX!=NULL);
    assert(fileY!=NULL);
    assert(fileZ!=NULL);
    assert(fileDataXr!=NULL);
    assert(fileDataXi!=NULL);
    assert(fileDataYr!=NULL);
    assert(fileDataYi!=NULL);
    assert(fileDataZr!=NULL);
    assert(fileDataZi!=NULL);
    FieldComponents E;
    Timer T;
    setTimer(&T);
    for (int i=0; i<Ns; i++){
        z = z_min+i*dz;
        NearFields(&myShape, &In, x, y, z, &E, tol, warnings, Nmax);
        fprintf(fileX, "%21.14E \n", x);
        fprintf(fileY, "%21.14E \n", y);
        fprintf(fileZ, "%21.14E \n", z);
        fprintf(fileDataXr, "%21.14E \n", creal(E.x));
        fprintf(fileDataXi, "%21.14E \n", cimag(E.x));
        fprintf(fileDataYr, "%21.14E \n", creal(E.y));
        fprintf(fileDataYi, "%21.14E \n", cimag(E.y));
        fprintf(fileDataZr, "%21.14E \n", creal(E.z));
        fprintf(fileDataZi, "%21.14E \n", cimag(E.z));
        progressBar(i, Ns, "Computing Near Fields ...");
    }
    unsetTimer(&T);
    fclose(fileX);
    fclose(fileY);
    fclose(fileZ);
    fclose(fileDataXr);
    fclose(fileDataXi);
    fclose(fileDataYr);
    fclose(fileDataYi);
    fclose(fileDataZr);
    fclose(fileDataZi);
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
}

void TransmissionLineNearInputImpedance(){
    int warnings=1;
    int Nmax=256;
    int Ns=401;
    double tol=1.0E-4;
    double MHz=1.0E6;
    double kHz=1.0E3;
    double mm=1.0E-3;
    double cm=1.0E-2;

    double freq_min=100.0*kHz;
    double freq_max=1000.0*MHz;
    double L=120.0*cm;
    double H=10.0*cm;
    double a=0.8*mm;
    double Z0=50.0;
    complex double ZL=50.0;

    // double freq_min=100.0*kHz;
    // double freq_max=100.0*MHz;
    // double L=250.0*cm;
    // double H=0.5*3.0*cm;
    // double a=1.0*mm;
    // double Z0=50.0;
    // complex double ZL=200.0;

    double segments_per_lambda=31.0;
    double lambda, lc=L/segments_per_lambda;

    double freq, dfreq;
    // dfreq = (log10(freq_max)-log10(freq_min))/(Ns-1.0);
    dfreq = (freq_max-freq_min)/(Ns-1.0);

    FILE *file=fopen("Data/TestBench/LoadedTL/Data.dat", "w");
    assert(file!=NULL);
    Timer T;
    setTimer(&T);
    for (int i=0; i<Ns; i++){

        // freq = pow(10.0, log10(freq_min)+i*dfreq);
        freq = freq_min+i*dfreq;
        lambda = c0/freq;

        // Create Shape
        createTransmissionLineDeltaGapCenter(L, H, lc);
        Shape myShape=DefaultShape;
        getShape(&myShape, 2, lambda, a);
        Port newPort1={1, 1.0, Z0};
        Port newPort2={2, 0.0, ZL};
        setPort(&myShape, newPort1);
        setPort(&myShape, newPort2);
        logShape(&myShape);
        int N=myShape.NBasis;
        Matrix Zmn=DefaultMatrix;
        Matrix Vm=DefaultMatrix;
        Matrix In=DefaultMatrix;
        allocateMatrix(&Zmn, N, N);
        allocateMatrix(&Vm, N, 1);
        allocateMatrix(&In, N, 1);
        MM_deltaGap(&myShape, tol, &Zmn, &Vm, warnings, Nmax);
        int *P=(int*)malloc((N+1)*sizeof(int));
        LUP_Decompose(&Zmn, P);
        LUP_Solve(&Zmn, P, &Vm, &In);
        free(P);
        P = NULL;

        complex double Zin;
        fprintf(file, "%21.14E ", freq);
        for (int n=0; n<N; n++){
            for (int i=0; i<myShape.NPorts; i++){
                if ((myShape.basisList[n].portNumber==myShape.portsList[i].portNumber)&&
                    (myShape.basisList[n].isPort)){
                    Zin = (myShape.portsList[i].Vin/In.data[n][0])-myShape.portsList[i].ZL;
                    fprintf(file, "%21.14E %21.14E ", creal(Zin), cimag(Zin));
                    fprintf(file, "%21.14E %21.14E ", creal(In.data[n][0]*myShape.portsList[i].ZL), 
                            cimag(In.data[n][0]*myShape.portsList[i].ZL));
                }
            }
        }
        fprintf(file, "\n");
        
        deallocateMatrix(&Zmn);
        deallocateMatrix(&Vm);
        deallocateMatrix(&In);
        deleteShape(&myShape);

        progressBar(i, Ns, "Computing Input Impedance ...");
    }
    unsetTimer(&T);
    fclose(file);
}

void VerticalDipoleDeltaGapCenterRadiationResistance(){
    int warnings=1;
    int Nmax=256;
    int Ns=1001;
    double tol=1.0E-4;
    double lambda=1.0;
    double L_min=0.01*lambda;
    double L_max=3.0*lambda;
    double segments_per_lambda=31.0;
    double dL=(L_max-L_min)/(Ns-1.0);
    double lc;
    double L, a;
    a = (1.0E-3)*lambda;
    //
    Timer T;
    setTimer(&T);
    FILE *file=fopen("Data/TestBench/RadiationResistance/Data.dat", "w");
    assert(file!=NULL);
    for (int i=0; i<Ns; i++){
        L = L_min+i*dL;
        lc = L/segments_per_lambda;
        // Create Shape
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
        MM_deltaGap(&myShape, tol, &Zmn, &Vm, warnings, Nmax);
        int *P=(int*)malloc((N+1)*sizeof(int));
        LUP_Decompose(&Zmn, P);
        LUP_Solve(&Zmn, P, &Vm, &In);
        free(P);
        P = NULL;
        // Get Results
        double Pradiated=Prad(&myShape, &In, tol, warnings, Nmax);
        double Rrad;
        complex double Zin;
        fprintf(file, "%21.14E ", L);
        double I_max=0.0;
        for (int n=0; n<N; n++){
            if (cabs(In.data[n][0])>I_max){
                I_max = cabs(In.data[n][0]);
            }
        }
        Rrad = 2.0*Pradiated/(I_max*I_max);
        fprintf(file, "%21.14E ", Rrad);
        for (int n=0; n<N; n++){
            for (int i=0; i<myShape.NPorts; i++){
                if ((myShape.basisList[n].portNumber==myShape.portsList[i].portNumber)&&
                    (myShape.basisList[n].isPort)){
                    Rrad = 2.0*Pradiated/(cabs(In.data[n][0])*cabs(In.data[n][0]));
                    Zin = (myShape.portsList[i].Vin/In.data[n][0])-myShape.portsList[i].ZL;
                    fprintf(file, "%21.14E %21.14E ", Rrad, creal(Zin));
                }
            }
        }
        fprintf(file, "\n");
        // Free The Memory
        progressBar(i, Ns, "Computing Far Fields ...");
        deallocateMatrix(&Zmn);
        deallocateMatrix(&Vm);
        deallocateMatrix(&In);
        deleteShape(&myShape);
    }
    fclose(file);
    unsetTimer(&T);
}

void YagiAntennaRadiationResistance(){
    int warnings=1;
    int Nmax=256;
    int Ns=101;

    double MHz=1.0E6;
    double mm=1.0E-3;
    double tol=1.0E-4;

    double freq_min=144*MHz;
    double freq_max=148*MHz;
    double segments_per_lambda=31.0;

    double a=0.5*6.35*mm;

    int NElements=9;
    double length[]={
        1038*mm,
        955*mm,
        956*mm,
        932*mm,
        916*mm,
        906*mm,
        897*mm,
        891*mm,
        887*mm
    };
    double xLocation[]={
        0*mm,
        312*mm,
        447*mm,
        699*mm,
        1050*mm,
        1482*mm,
        1986*mm,
        2553*mm,
        3168*mm
    };
    int isExcitaiton[]={
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        0
    };
    YagiElements elements={NElements, length, xLocation, isExcitaiton};
    double lambda, freq;
    double dfreq=(freq_max-freq_min)/(Ns-1.0);
    //
    Timer T;
    setTimer(&T);
    FILE *file=fopen("Data/TestBench/Yagi/Data.dat", "w");
    assert(file!=NULL);
    for (int i=0; i<Ns; i++){
        freq = freq_min+i*dfreq;
        lambda = c0/freq;
        // Create Shape
        createYagiAntenna(&elements, segments_per_lambda);
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
        MM_deltaGap(&myShape, tol, &Zmn, &Vm, warnings, Nmax);
        int *P=(int*)malloc((N+1)*sizeof(int));
        LUP_Decompose(&Zmn, P);
        LUP_Solve(&Zmn, P, &Vm, &In);
        free(P);
        P = NULL;
        // Get Results
        double Pradiated=Prad(&myShape, &In, tol, warnings, Nmax);
        double Rrad;
        complex double Zin;
        fprintf(file, "%21.14E ", freq);
        for (int n=0; n<N; n++){
            for (int i=0; i<myShape.NPorts; i++){
                if ((myShape.basisList[n].portNumber==myShape.portsList[i].portNumber)&&
                    (myShape.basisList[n].isPort)){
                    Rrad = 2.0*Pradiated/(cabs(In.data[n][0])*cabs(In.data[n][0]));
                    Zin = (myShape.portsList[i].Vin/In.data[n][0])-myShape.portsList[i].ZL;
                    fprintf(file, "%21.14E %21.14E %21.14E ", Rrad, creal(Zin), cimag(Zin));
                }
            }
        }
        fprintf(file, "\n");
        // Free The Memory
        progressBar(i, Ns, "Computing Far Fields ...");
        deallocateMatrix(&Zmn);
        deallocateMatrix(&Vm);
        deallocateMatrix(&In);
        deleteShape(&myShape);
    }
    fclose(file);
    unsetTimer(&T);
}

void YagiAntennaRadiationPattern(){
    int warnings=1;
    int Nmax=256;
    int Ns=1001;

    double MHz=1.0E6;
    double mm=1.0E-3;
    double tol=1.0E-4;

    double freq=146*MHz;
    double lambda=c0/freq;
    double segments_per_lambda=31.0;

    double a=0.5*6.35*mm;

    int NElements=9;
    double length[]={
        1038*mm,
        955*mm,
        956*mm,
        932*mm,
        916*mm,
        906*mm,
        897*mm,
        891*mm,
        887*mm
    };
    double xLocation[]={
        0*mm,
        312*mm,
        447*mm,
        699*mm,
        1050*mm,
        1482*mm,
        1986*mm,
        2553*mm,
        3168*mm
    };
    int isExcitaiton[]={
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        0
    };
    YagiElements elements={NElements, length, xLocation, isExcitaiton};
    
    double theta, phi;
    double dtheta=2.0*pi/(Ns-1.0);
    double dphi=2.0*pi/(Ns-1.0);
    //
    Timer T;
    setTimer(&T);
    // Create Shape
    createYagiAntenna(&elements, segments_per_lambda);
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
    progressBar(0, Ns, "Solving MM System ...");
    MM_deltaGap(&myShape, tol, &Zmn, &Vm, warnings, Nmax);
    int *P=(int*)malloc((N+1)*sizeof(int));
    LUP_Decompose(&Zmn, P);
    LUP_Solve(&Zmn, P, &Vm, &In);
    free(P);
    P = NULL;
    FILE *file=fopen("Data/TestBench/Yagi/DataFarFieldTheta.dat", "w");
    assert(file!=NULL);
    FarFieldComponents E;
    for (int i=0; i<Ns; i++){
        phi = 0.0;
        theta = -pi+i*dtheta;
        FarFields(&myShape, &In, theta, phi, &E);
        fprintf(file, "%21.14E %21.14E\n", theta, cabs(E.theta));
        progressBar(i, Ns, "Computing Far Fields Theta Cut ...");
    }
    fclose(file);
    file=fopen("Data/TestBench/Yagi/DataFarFieldPhi.dat", "w");
    assert(file!=NULL);
    for (int i=0; i<Ns; i++){
        phi = i*dphi;
        theta = pi/2.0;
        FarFields(&myShape, &In, theta, phi, &E);
        fprintf(file, "%21.14E %21.14E\n", phi, cabs(E.theta));
        progressBar(i, Ns, "Computing Far Fields Theta Cut ...");
    }
    fclose(file);
    // Free The Memory
    deallocateMatrix(&Zmn);
    deallocateMatrix(&Vm);
    deallocateMatrix(&In);
    deleteShape(&myShape);
    
    unsetTimer(&T);
}