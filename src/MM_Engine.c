//
#include "MM_Engine.h"
#include "QuadL.h"

// 
int warnings=0;
int Nmax=32;

void settings(int warningsIn, int NmaxIn){
    warnings = warningsIn;
    Nmax = NmaxIn;
}

double Rmn_pp(double alpha, double alpha_, Basis basis_m, Basis basis_n, double a){
    Vector R=basis_m.rp;
    R = subVector(R, basis_n.rp);
    R = subVector(R, scaleVector(basis_m.Lp, alpha));
    R = addVector(R, scaleVector(basis_n.Lp, alpha_));
    double magR=magVector(R);
    return sqrt(magR*magR+a*a);
}

double Rmn_pm(double alpha, double alpha_, Basis basis_m, Basis basis_n, double a){
    Vector R=basis_m.rp;
    R = subVector(R, basis_n.rm);
    R = subVector(R, scaleVector(basis_m.Lp, alpha));
    R = subVector(R, scaleVector(basis_n.Lm, alpha_));
    double magR=magVector(R);
    return sqrt(magR*magR+a*a);
}

double Rmn_mp(double alpha, double alpha_, Basis basis_m, Basis basis_n, double a){
    Vector R=basis_m.rm;
    R = subVector(R, basis_n.rp);
    R = addVector(R, scaleVector(basis_m.Lm, alpha));
    R = addVector(R, scaleVector(basis_n.Lp, alpha_));
    double magR=magVector(R);
    return sqrt(magR*magR+a*a);
}

double Rmn_mm(double alpha, double alpha_, Basis basis_m, Basis basis_n, double a){
    Vector R=basis_m.rm;
    R = subVector(R, basis_n.rm);
    R = addVector(R, scaleVector(basis_m.Lm, alpha));
    R = subVector(R, scaleVector(basis_n.Lm, alpha_));
    double magR=magVector(R);
    return sqrt(magR*magR+a*a);
}

// Green's Function

complex double gmn_pp(double alpha, double alpha_, Basis basis_m, Basis basis_n, double a){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    double R=Rmn_pp(alpha, alpha_, basis_m, basis_n, a);
    return cexp(-j*k*R)/(4.0*pi*R);
}

complex double gmn_pm(double alpha, double alpha_, Basis basis_m, Basis basis_n, double a){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    double R=Rmn_pm(alpha, alpha_, basis_m, basis_n, a);
    return cexp(-j*k*R)/(4.0*pi*R);
}

complex double gmn_mp(double alpha, double alpha_, Basis basis_m, Basis basis_n, double a){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    double R=Rmn_mp(alpha, alpha_, basis_m, basis_n, a);
    return cexp(-j*k*R)/(4.0*pi*R);
}

complex double gmn_mm(double alpha, double alpha_, Basis basis_m, Basis basis_n, double a){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    double R=Rmn_mm(alpha, alpha_, basis_m, basis_n, a);
    return cexp(-j*k*R)/(4.0*pi*R);
}

// Terms ++

complex double phi_mn_integrand_pp(complex double alpha, complex double alpha_, void *args){
    integrandArgs *myArgs=(integrandArgs*)args;
    Basis basis_m=myArgs->basis_m;
    Basis basis_n=myArgs->basis_n;
    double a=myArgs->a;
    return gmn_pp(alpha, alpha_, basis_m, basis_n, a);
}

complex double phi_mn_pp(void *args, double tol){
    return QuadL_2D(phi_mn_integrand_pp, args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax);
}

complex double psi_mn_integrand_pp(complex double alpha, complex double alpha_, void *args){
    integrandArgs *myArgs=(integrandArgs*)args;
    Basis basis_m=myArgs->basis_m;
    Basis basis_n=myArgs->basis_n;
    Vector Lm_p=myArgs->Lm_p;
    Vector Ln_p=myArgs->Ln_p;
    double a=myArgs->a;
    return dotVector(Lm_p, Ln_p)*alpha*alpha_*gmn_pp(alpha, alpha_, basis_m, basis_n, a);
}

complex double psi_mn_pp(void *args, double tol){
    return QuadL_2D(psi_mn_integrand_pp, args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax);
}

// Terms +-

complex double phi_mn_integrand_pm(complex double alpha, complex double alpha_, void *args){
    integrandArgs *myArgs=(integrandArgs*)args;
    Basis basis_m=myArgs->basis_m;
    Basis basis_n=myArgs->basis_n;
    double a=myArgs->a;
    return gmn_pm(alpha, alpha_, basis_m, basis_n, a);
}

complex double phi_mn_pm(void *args, double tol){
    return QuadL_2D(phi_mn_integrand_pm, args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax);
}

complex double psi_mn_integrand_pm(complex double alpha, complex double alpha_, void *args){
    integrandArgs *myArgs=(integrandArgs*)args;
    Basis basis_m=myArgs->basis_m;
    Basis basis_n=myArgs->basis_n;
    Vector Lm_p=myArgs->Lm_p;
    Vector Ln_m=myArgs->Ln_m;
    double a=myArgs->a;
    return dotVector(Lm_p, Ln_m)*alpha*alpha_*gmn_pm(alpha, alpha_, basis_m, basis_n, a);
}

complex double psi_mn_pm(void *args, double tol){
    return QuadL_2D(psi_mn_integrand_pm, args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax);
}

// Terms -+

complex double phi_mn_integrand_mp(complex double alpha, complex double alpha_, void *args){
    integrandArgs *myArgs=(integrandArgs*)args;
    Basis basis_m=myArgs->basis_m;
    Basis basis_n=myArgs->basis_n;
    double a=myArgs->a;
    return gmn_mp(alpha, alpha_, basis_m, basis_n, a);
}

complex double phi_mn_mp(void *args, double tol){
    return QuadL_2D(phi_mn_integrand_mp, args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax);
}

complex double psi_mn_integrand_mp(complex double alpha, complex double alpha_, void *args){
    integrandArgs *myArgs=(integrandArgs*)args;
    Basis basis_m=myArgs->basis_m;
    Basis basis_n=myArgs->basis_n;
    Vector Lm_m=myArgs->Lm_m;
    Vector Ln_p=myArgs->Ln_p;
    double a=myArgs->a;
    return dotVector(Lm_m, Ln_p)*alpha*alpha_*gmn_mp(alpha, alpha_, basis_m, basis_n, a);
}

complex double psi_mn_mp(void *args, double tol){
    return QuadL_2D(psi_mn_integrand_mp, args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax);
}

// Terms --

complex double phi_mn_integrand_mm(complex double alpha, complex double alpha_, void *args){
    integrandArgs *myArgs=(integrandArgs*)args;
    Basis basis_m=myArgs->basis_m;
    Basis basis_n=myArgs->basis_n;
    double a=myArgs->a;
    return gmn_mm(alpha, alpha_, basis_m, basis_n, a);
}

complex double phi_mn_mm(void *args, double tol){
    return QuadL_2D(phi_mn_integrand_mm, args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax);
}

complex double psi_mn_integrand_mm(complex double alpha, complex double alpha_, void *args){
    integrandArgs *myArgs=(integrandArgs*)args;
    Basis basis_m=myArgs->basis_m;
    Basis basis_n=myArgs->basis_n;
    Vector Lm_m=myArgs->Lm_m;
    Vector Ln_m=myArgs->Ln_m;
    double a=myArgs->a;
    return dotVector(Lm_m, Ln_m)*alpha*alpha_*gmn_mm(alpha, alpha_, basis_m, basis_n, a);
}

complex double psi_mn_mm(void *args, double tol){
    return QuadL_2D(psi_mn_integrand_mm, args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax);
}

//

complex double sinc(complex double z){
    return cabs(z)<1.0E-8 ? 1.0 : csin(z)/z;
}

typedef struct SingularArgs SingularArgs;
struct SingularArgs{
    double L, a;
};

complex double I1_integrand(complex double alpha, complex double alpha_, void *args){
    SingularArgs *myArgs=(SingularArgs*)args;
    double L=myArgs->L;
    double a=myArgs->a;
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    double R=sqrt(L*L*(alpha-alpha_)*(alpha-alpha_)+a*a);
    return (-j*k*L*L/(4.0*pi))*alpha*alpha_*cexp(-j*k*R/2.0)*sinc(k*R/2.0);
}

complex double I2_integrand(complex double alpha, void *args){
    SingularArgs *myArgs=(SingularArgs*)args;
    double L=myArgs->L;
    double a=myArgs->a;
    double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a)+(1.0-alpha)*L;
    double B=sqrt(alpha*alpha*L*L+a*a)-alpha*L;
    return (L/(4.0*pi))*alpha*alpha*log(A/B);
}

complex double I3_integrand(complex double alpha, void *args){
    SingularArgs *myArgs=(SingularArgs*)args;
    double L=myArgs->L;
    double a=myArgs->a;
    double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a);
    double B=sqrt(alpha*alpha*L*L+a*a);
    return (1.0/(4.0*pi))*alpha*(A-B);
}

complex double I4_integrand(complex double alpha, complex double alpha_, void *args){
    SingularArgs *myArgs=(SingularArgs*)args;
    double L=myArgs->L;
    double a=myArgs->a;
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    double R=sqrt(L*L*(alpha-alpha_)*(alpha-alpha_)+a*a);
    return (-j*k/(4.0*pi))*cexp(-j*k*R/2.0)*sinc(k*R/2.0);
}

complex double I5_integrand(complex double alpha, void *args){
    SingularArgs *myArgs=(SingularArgs*)args;
    double L=myArgs->L;
    double a=myArgs->a;
    double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a)+(1.0-alpha)*L;
    double B=sqrt(alpha*alpha*L*L+a*a)-alpha*L;
    return (1.0/(4.0*pi*L))*log(A/B);
}

complex double I6_integrand(complex double alpha, complex double alpha_, void *args){
    SingularArgs *myArgs=(SingularArgs*)args;
    double L=myArgs->L;
    double a=myArgs->a;
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    double R=sqrt(L*L*(1-alpha-alpha_)*(1-alpha-alpha_)+a*a);
    return (-j*k*L*L/(4.0*pi))*alpha*alpha_*cexp(-j*k*R/2.0)*sinc(k*R/2.0);
}

complex double I7_integrand(complex double alpha, void *args){
    SingularArgs *myArgs=(SingularArgs*)args;
    double L=myArgs->L;
    double a=myArgs->a;
    double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a)+(1.0-alpha)*L;
    double B=sqrt(alpha*alpha*L*L+a*a)-alpha*L;
    return (L/(4.0*pi))*alpha*(1.0-alpha)*log(A/B);
}

complex double I8_integrand(complex double alpha, void *args){
    SingularArgs *myArgs=(SingularArgs*)args;
    double L=myArgs->L;
    double a=myArgs->a;
    double A=sqrt((1.0-alpha)*(1.0-alpha)*L*L+a*a);
    double B=sqrt(alpha*alpha*L*L+a*a);
    return (-1.0/(4.0*pi))*alpha*(A-B);
}

complex double I1(double L, double a, double tol){
    SingularArgs args={L, a};
    return QuadL_2D(I1_integrand, &args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax)+
        QuadL_1D(I2_integrand, &args, 0.0, 1.0, tol, warnings, Nmax)+
        QuadL_1D(I3_integrand, &args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double I2(double L, double a, double tol){
    SingularArgs args={L, a};
    return QuadL_2D(I4_integrand, &args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax)+
        QuadL_1D(I5_integrand, &args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double I3(double L, double a, double tol){
    SingularArgs args={L, a};
    return QuadL_2D(I6_integrand, &args, 0.0, 1.0, 0.0, 1.0, tol, warnings, Nmax)+
        QuadL_1D(I7_integrand, &args, 0.0, 1.0, tol, warnings, Nmax)+
        QuadL_1D(I8_integrand, &args, 0.0, 1.0, tol, warnings, Nmax);
}

void MM_deltaGap(Shape *myShape, double tol, Matrix *Zmn, Matrix *Vm, int isWarnings, int NmaxIn){
    warnings = isWarnings;
    Nmax = NmaxIn;
    double k=2.0*pi;
    complex double j=csqrt(-1.0);
    int NBasis=myShape->NBasis;
    int NPorts=myShape->NPorts;
    complex double A, B;
    int flag;
    for (int m=0; m<NBasis; m++){
        for (int i=0; i<NPorts; i++){
            if ((myShape->basisList[m].portNumber==myShape->portsList[i].portNumber)&&
                (myShape->basisList[m].isPort)){
                Vm->data[m][0] = myShape->portsList[i].Vin;
                Zmn->data[m][m] = myShape->portsList[i].ZL;
            }
        }
        for (int n=0; n<NBasis; n++){
            flag = 0;
            integrandArgs args={myShape->basisList[m], myShape->basisList[n], 
                myShape->basisList[m].Lp, myShape->basisList[m].Lm,
                myShape->basisList[n].Lp, myShape->basisList[n].Lm,
                myShape->a};
            if (isEqualVector(myShape->basisList[m].rm, myShape->basisList[n].rm)&&
                isEqualVector(myShape->basisList[m].rn, myShape->basisList[n].rn)&&
                isEqualVector(myShape->basisList[m].rp, myShape->basisList[n].rp)){
                // printf("Scenario I\n");
                A = I1(myShape->basisList[m].lp, myShape->a, tol)+
                    psi_mn_pm(&args, tol)+
                    psi_mn_mp(&args, tol)+
                    I1(myShape->basisList[m].lm, myShape->a, tol);
                B = I2(myShape->basisList[m].lp, myShape->a, tol)+
                    -phi_mn_pm(&args, tol)+
                    -phi_mn_mp(&args, tol)+
                    I2(myShape->basisList[m].lm, myShape->a, tol);
                flag++;
            }else
            if (isEqualVector(myShape->basisList[m].rn, myShape->basisList[n].rm)&&
                isEqualVector(myShape->basisList[m].rp, myShape->basisList[n].rn)){
                // printf("Scenario II\n");
                A = psi_mn_pp(&args, tol)+
                    I3(myShape->basisList[m].lp, myShape->a, tol)+
                    psi_mn_mp(&args, tol)+
                    psi_mn_mm(&args, tol);
                B = phi_mn_pp(&args, tol)+
                    -I2(myShape->basisList[m].lp, myShape->a, tol)+
                    -phi_mn_mp(&args, tol)+
                    phi_mn_mm(&args, tol);
                flag++;
            }else 
            if (isEqualVector(myShape->basisList[m].rn, myShape->basisList[n].rp)&&
                isEqualVector(myShape->basisList[m].rm, myShape->basisList[n].rn)){
                // printf("Scenario III\n");
                A = psi_mn_pp(&args, tol)+
                    psi_mn_pm(&args, tol)+
                    I3(myShape->basisList[m].lm, myShape->a, tol)+
                    psi_mn_mm(&args, tol);
                B = phi_mn_pp(&args, tol)+
                    -phi_mn_pm(&args, tol)+
                    -I2(myShape->basisList[m].lm, myShape->a, tol)+
                    phi_mn_mm(&args, tol);
                flag++;
            }else
            if (isEqualVector(myShape->basisList[m].rn, myShape->basisList[n].rm)&&
                isEqualVector(myShape->basisList[m].rm, myShape->basisList[n].rn)){
                // printf("Scenario IV\n");
                A = psi_mn_pp(&args, tol)+
                    psi_mn_pm(&args, tol)+
                    psi_mn_mp(&args, tol)+
                    -I3(myShape->basisList[m].lm, myShape->a, tol);
                B = phi_mn_pp(&args, tol)+
                    -phi_mn_pm(&args, tol)+
                    -phi_mn_mp(&args, tol)+
                    +I2(myShape->basisList[m].lm, myShape->a, tol);
                flag++;
            }else
            if (isEqualVector(myShape->basisList[m].rn, myShape->basisList[n].rp)&&
                isEqualVector(myShape->basisList[m].rp, myShape->basisList[n].rn)){
                // printf("Scenario V\n");
                A = -I3(myShape->basisList[m].lp, myShape->a, tol)+
                    psi_mn_pm(&args, tol)+
                    psi_mn_mp(&args, tol)+
                    psi_mn_mm(&args, tol);
                B = I2(myShape->basisList[m].lp, myShape->a, tol)+
                    -phi_mn_pm(&args, tol)+
                    -phi_mn_mp(&args, tol)+
                    phi_mn_mm(&args, tol);
                flag++;
            }else{
                // printf("Scenario VI\n");
                A = psi_mn_pp(&args, tol)+
                    psi_mn_pm(&args, tol)+
                    psi_mn_mp(&args, tol)+
                    psi_mn_mm(&args, tol);
                B = phi_mn_pp(&args, tol)+
                    -phi_mn_pm(&args, tol)+
                    -phi_mn_mp(&args, tol)+
                    phi_mn_mm(&args, tol);
                flag++;
            }
            assert(flag==1);
            Zmn->data[m][n]+=j*k*eta0*A-j*(eta0/k)*B;
        }
    }
}
