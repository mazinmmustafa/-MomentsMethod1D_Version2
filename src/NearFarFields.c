//
#include "NearFarFields.h"
#include "QuadL.h"

// 

Vector Rn_m(Vector r, Vector rn_m, Vector Ln_m, double alpha_){
    Vector R=r;
    R = subVector(R, rn_m);
    R = subVector(R, scaleVector(Ln_m, alpha_));
    return R;
}

Vector Rn_p(Vector r, Vector rn_p, Vector Ln_p, double alpha_){
    Vector R=r;
    R = subVector(R, rn_p);
    R = addVector(R, scaleVector(Ln_p, alpha_));
    return R;
}

double magRn(Vector R, double a){
    return sqrt(magVector(R)*magVector(R)+a*a);
}

complex double gn_m(Vector r, Vector rn_m, Vector Ln_m, double alpha_, double a){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    Vector R=Rn_m(r, rn_m, Ln_m, alpha_);
    double magR=magRn(R, a);
    return cexp(-j*k*magR)/(4.0*pi*magR);
}

complex double gn_p(Vector r, Vector rn_p, Vector Ln_p, double alpha_, double a){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    Vector R=Rn_p(r, rn_p, Ln_p, alpha_);
    double magR=magRn(R, a);
    return cexp(-j*k*magR)/(4.0*pi*magR);
}

typedef struct TermsArgs TermsArgs;
struct TermsArgs{
    Vector r, rn_m, rn_p, Ln_m, Ln_p;
    double a;
};

complex double term_1_integrand_x(complex double alpha_, void *args){
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_m=myArgs->rn_m;
    Vector Ln_m=myArgs->Ln_m;
    double a=myArgs->a;
    double factor=Ln_m.x;
    return alpha_*gn_m(r, rn_m, Ln_m, alpha_, a)*factor;
}

complex double term_1_integrand_y(complex double alpha_, void *args){
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_m=myArgs->rn_m;
    Vector Ln_m=myArgs->Ln_m;
    double a=myArgs->a;
    double factor=Ln_m.y;
    return alpha_*gn_m(r, rn_m, Ln_m, alpha_, a)*factor;
}

complex double term_1_integrand_z(complex double alpha_, void *args){
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_m=myArgs->rn_m;
    Vector Ln_m=myArgs->Ln_m;
    double a=myArgs->a;
    double factor=Ln_m.z;
    return alpha_*gn_m(r, rn_m, Ln_m, alpha_, a)*factor;
}

complex double term_2_integrand_x(complex double alpha_, void *args){
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_p=myArgs->rn_p;
    Vector Ln_p=myArgs->Ln_p;
    double a=myArgs->a;
    double factor=Ln_p.x;
    return alpha_*gn_p(r, rn_p, Ln_p, alpha_, a)*factor;
}

complex double term_2_integrand_y(complex double alpha_, void *args){
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_p=myArgs->rn_p;
    Vector Ln_p=myArgs->Ln_p;
    double a=myArgs->a;
    double factor=Ln_p.y;
    return alpha_*gn_p(r, rn_p, Ln_p, alpha_, a)*factor;
}

complex double term_2_integrand_z(complex double alpha_, void *args){
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_p=myArgs->rn_p;
    Vector Ln_p=myArgs->Ln_p;
    double a=myArgs->a;
    double factor=Ln_p.z;
    return alpha_*gn_p(r, rn_p, Ln_p, alpha_, a)*factor;
}

complex double term_3_integrand_x(complex double alpha_, void *args){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_m=myArgs->rn_m;
    Vector Ln_m=myArgs->Ln_m;
    double a=myArgs->a;
    complex double g=gn_m(r, rn_m, Ln_m, alpha_, a);
    Vector R=Rn_m(r, rn_m, Ln_m, alpha_);
    double magR=magRn(R, a);
    double factor = R.x/magR;
    return (+(1.0+j*k*magR)/(magR))*g*factor;
}

complex double term_3_integrand_y(complex double alpha_, void *args){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_m=myArgs->rn_m;
    Vector Ln_m=myArgs->Ln_m;
    double a=myArgs->a;
    complex double g=gn_m(r, rn_m, Ln_m, alpha_, a);
    Vector R=Rn_m(r, rn_m, Ln_m, alpha_);
    double magR=magRn(R, a);
    double factor = R.y/magR;
    return (+(1.0+j*k*magR)/(magR))*g*factor;
}

complex double term_3_integrand_z(complex double alpha_, void *args){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_m=myArgs->rn_m;
    Vector Ln_m=myArgs->Ln_m;
    double a=myArgs->a;
    
    complex double g=gn_m(r, rn_m, Ln_m, alpha_, a);
    Vector R=Rn_m(r, rn_m, Ln_m, alpha_);
    double magR=magRn(R, a);
    double factor = R.z/magR;
    return (+(1.0+j*k*magR)/(magR))*g*factor;
}

complex double term_4_integrand_x(complex double alpha_, void *args){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_p=myArgs->rn_p;
    Vector Ln_p=myArgs->Ln_p;
    double a=myArgs->a;
    complex double g=gn_p(r, rn_p, Ln_p, alpha_, a);
    Vector R=Rn_p(r, rn_p, Ln_p, alpha_);
    double magR=magRn(R, a);
    double factor = R.x/magR;
    return (-(1.0+j*k*magR)/(magR))*g*factor;
}

complex double term_4_integrand_y(complex double alpha_, void *args){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_p=myArgs->rn_p;
    Vector Ln_p=myArgs->Ln_p;
    double a=myArgs->a;
    complex double g=gn_p(r, rn_p, Ln_p, alpha_, a);
    Vector R=Rn_p(r, rn_p, Ln_p, alpha_);
    double magR=magRn(R, a);
    double factor = R.y/magR;
    return (-(1.0+j*k*magR)/(magR))*g*factor;
}

complex double term_4_integrand_z(complex double alpha_, void *args){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    TermsArgs *myArgs=(TermsArgs*)args;
    Vector r=myArgs->r;
    Vector rn_p=myArgs->rn_p;
    Vector Ln_p=myArgs->Ln_p;
    double a=myArgs->a;
    complex double g=gn_p(r, rn_p, Ln_p, alpha_, a);
    Vector R=Rn_p(r, rn_p, Ln_p, alpha_);
    double magR=magRn(R, a);
    double factor = R.z/magR;
    return (-(1.0+j*k*magR)/(magR))*g*factor;
}

complex double term_1_x(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_1_integrand_x, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_1_y(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_1_integrand_y, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_1_z(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_1_integrand_z, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_2_x(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_2_integrand_x, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_2_y(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_2_integrand_y, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_2_z(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_2_integrand_z, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_3_x(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_3_integrand_x, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_3_y(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_3_integrand_y, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_3_z(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_3_integrand_z, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_4_x(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_4_integrand_x, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_4_y(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_4_integrand_y, args, 0.0, 1.0, tol, warnings, Nmax);
}

complex double term_4_z(void *args, double tol, int warnings, int Nmax){
    return QuadL_1D(term_4_integrand_z, args, 0.0, 1.0, tol, warnings, Nmax);
}

void NearFields(Shape *myShape, Matrix *In, double x, double y, double z, FieldComponents *E, 
    double tol, int warnings, int Nmax){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    int N=myShape->NBasis;
    double a=myShape->a;
    double lambda=myShape->lambda;
    Vector r={x/lambda, y/lambda, z/lambda};
    E->x = 0.0;
    E->y = 0.0;
    E->z = 0.0;
    for (int n=0; n<N; n++){
        Vector rn_m=myShape->basisList[n].rm;
        Vector rn_p=myShape->basisList[n].rp;
        Vector Ln_m=myShape->basisList[n].Lm;
        Vector Ln_p=myShape->basisList[n].Lp;
        TermsArgs args={r, rn_m, rn_p, Ln_m, Ln_p, a};
        E->x+=In->data[n][0]*(-j*k*eta0*(term_1_x(&args, tol, warnings, Nmax)+
                                         term_2_x(&args, tol, warnings, Nmax))+
                             +(j*eta0/k)*(term_3_x(&args, tol, warnings, Nmax)+
                                          term_4_x(&args, tol, warnings, Nmax)));
        E->y+=In->data[n][0]*(-j*k*eta0*(term_1_y(&args, tol, warnings, Nmax)+
                                         term_2_y(&args, tol, warnings, Nmax))+
                             +(j*eta0/k)*(term_3_y(&args, tol, warnings, Nmax)+
                                          term_4_y(&args, tol, warnings, Nmax)));
        E->z+=In->data[n][0]*(-j*k*eta0*(term_1_z(&args, tol, warnings, Nmax)+
                                         term_2_z(&args, tol, warnings, Nmax))+
                             +(j*eta0/k)*(term_3_z(&args, tol, warnings, Nmax)+
                                          term_4_z(&args, tol, warnings, Nmax)));
    }
}

complex double kappa_m(double theta, double phi, Vector Ln_m, Vector rn_m){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    Vector r_hat={sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
    double dot=dotVector(Ln_m, r_hat);
    complex double factor=cexp(j*k*dotVector(rn_m, r_hat));
    if (fabs(dot)<1.0E-14){
        return factor/2.0;
    }else{
        return factor*(cexp(+j*k*dot)*(1.0-j*k*dot)-1.0)/(k*k*dot*dot);
    }
}

complex double kappa_p(double theta, double phi, Vector Ln_p, Vector rn_p){
    complex double j=csqrt(-1.0);
    double k=2.0*pi;
    Vector r_hat={sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)};
    double dot=dotVector(Ln_p, r_hat);
    complex double factor=cexp(j*k*dotVector(rn_p, r_hat));
    if (fabs(dot)<1.0E-14){
        return factor/2.0;
    }else{
        return factor*(cexp(-j*k*dot)*(1.0+j*k*dot)-1.0)/(k*k*dot*dot);
    }
}

void FarFields(Shape *myShape, Matrix *In, double theta, double phi, FarFieldComponents *E){
    double k=2.0*pi;
    int N=myShape->NBasis;
    E->theta = 0.0;
    E->phi = 0.0;
    Vector Theta={cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta)};
    Vector Phi={-sin(phi), cos(phi), 0.0};
    for (int n=0; n<N; n++){
        Vector rn_m=myShape->basisList[n].rm;
        Vector rn_p=myShape->basisList[n].rp;
        Vector Ln_m=myShape->basisList[n].Lm;
        Vector Ln_p=myShape->basisList[n].Lp;
        E->theta+=In->data[n][0]*(dotVector(Theta, Ln_m)*kappa_m(theta, phi, Ln_m, rn_m)+
                                    dotVector(Theta, Ln_p)*kappa_p(theta, phi, Ln_p, rn_p));
        E->phi+=In->data[n][0]*(dotVector(Phi, Ln_m)*kappa_m(theta, phi, Ln_m, rn_m)+
                                  dotVector(Phi, Ln_p)*kappa_p(theta, phi, Ln_p, rn_p));
    }
    E->theta*=(k*eta0/(4.0*pi));
    E->phi*=(k*eta0/(4.0*pi));
}

typedef struct TermsArgs2 TermsArgs2;
struct TermsArgs2{
    Shape *myShape;
    Matrix *In;
};

complex double FarFieldIntegrand(complex double theta, complex double phi, void *args){
    TermsArgs2 *myArgs=(TermsArgs2*)args;
    FarFieldComponents E;
    FarFields(myArgs->myShape, myArgs->In, theta, phi, &E);
    double E_theta_mag=cabs(E.theta);
    double E_phi_mag=cabs(E.phi);
    return (E_theta_mag*E_theta_mag+E_phi_mag*E_phi_mag)*sin(theta);
}

double Prad(Shape *myShape, Matrix *In, double tol, int warnings, int Nmax){
    TermsArgs2 args={myShape, In};
    return QuadL_2D(FarFieldIntegrand, &args, 0.0, 2.0*pi, 0.0, pi, tol, warnings, Nmax)/(2.0*eta0);
}