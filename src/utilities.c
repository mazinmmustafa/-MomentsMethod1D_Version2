//
#include "utilities.h"

void showComplex(complex double z){
    printf("(%21.14E, %21.14E)\n", creal(z), cimag(z));
}

#ifdef _WIN32
void progressBar(int n, int N, char *message){
    n++;
    int NBar=25;
    printf("\rProgress:[");
    for (int i=0; i<NBar*n/N; i++){
        printf("#");
    }
    for (int i=NBar*n/N; i<NBar; i++){
        printf("-");
    }
    printf("] %3d%%", 100*n/N);
    if (n==N){
        printf(" Done!\n");
    }else{
        printf(" %s", message);
    }
    fflush(stdout);
}

void setTimer(Timer *T){
    T->start = clock();
    T->isSet = 1;
}

void unsetTimer(Timer *T){
    if (T->isSet==1){
		T->stop = clock();
		T->elapsed = (double)(T->stop-T->start)/CLOCKS_PER_SEC;
        if (T->elapsed<60.0){
            printf("Elapsed time is %4.2f seconds\n", T->elapsed);
        }else
        if (T->elapsed<3600.0){
            printf("Elapsed time is %4.2f minutes\n", T->elapsed/60.0);
        }else{
            printf("Elapsed time is %4.2f hours\n", T->elapsed/3600.0);
        }
        T->isSet = 0;
    }else{
        printf("Warning: Timer Was Not Set!\n");
    }
}
#endif

#ifdef linux
void progressBar(int n, int N, char *message){
    n++;
    int NBar=25;
    printf("\33[2K\rProgress:[");
    for (int i=0; i<NBar*n/N; i++){
        printf("#");
    }
    for (int i=NBar*n/N; i<NBar; i++){
        printf("-");
    }
    printf("] %3d%%", 100*n/N);
    if (n==N){
        printf(" Done!\n");
    }else{
        printf(" %s", message);
    }
    fflush(stdout);
}

void setTimer(Timer *T){
    timespec_get(&(T->start), TIME_UTC);
    T->isSet = 1;
}

void unsetTimer(Timer *T){
    if (T->isSet==1){
        timespec_get(&(T->stop), TIME_UTC);
        T->elapsed = (double)(T->stop.tv_sec-T->start.tv_sec)+
            ((double)(T->stop.tv_nsec-T->start.tv_nsec)/1000000000L);
        if (T->elapsed<60.0){
            printf("Elapsed time is %4.2f seconds\n", T->elapsed);
        }else
        if (T->elapsed<3600.0){
            printf("Elapsed time is %4.2f minutes\n", T->elapsed/60.0);
        }else{
            printf("Elapsed time is %4.2f hours\n", T->elapsed/3600.0);
        }
        T->isSet = 0;
    }else{
        printf("Warning: Timer Was Not Set!\n");
    }
}
#endif

double deg2rad(double x){
    return x*pi/180.0;
}

double rad2deg(double x){
    return x*180.0/pi;
}

void setRandomSeed(){
    time_t t;
    srand(time(&t));
}

int randInt(int a, int b){
    return b>a ? a+rand()%(b-a+1) : 0;
}

double randDouble(double a, double b){
    return b>a ? a+rand()/(RAND_MAX/(b-a)) : 0.0;
}


