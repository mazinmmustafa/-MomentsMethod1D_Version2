//
#include "myLib.h"

int portCounter=0;

void getShape(Shape *myShape, int NPorts, double lambda, double a){
    assert(myShape->basisList==NULL);
    assert(myShape->portsList==NULL);
    myShape->lambda = lambda;
    FILE *file;
    file = fopen("Mesh/basisInfo.dat", "r");
    assert(file!=NULL);
    assert(fscanf(file, "%d", &(myShape->NBasis)));
    fclose(file);
    int NBasis=myShape->NBasis;
    assert(NBasis>0);
    myShape->basisList = (Basis*)malloc(NBasis*sizeof(Basis));
    file = fopen("Mesh/basis.dat", "r");
    assert(file!=NULL);
    double x, y, z;
    int portNumber;
    for (int i=0; i<NBasis; i++){
        assert(fscanf(file, "%lf", &x));
        assert(fscanf(file, "%lf", &y));
        assert(fscanf(file, "%lf", &z));
        Vector rm={x/lambda, y/lambda, z/lambda};
        assert(fscanf(file, "%lf", &x));
        assert(fscanf(file, "%lf", &y));
        assert(fscanf(file, "%lf", &z));
        Vector rn={x/lambda, y/lambda, z/lambda};
        assert(fscanf(file, "%lf", &x));
        assert(fscanf(file, "%lf", &y));
        assert(fscanf(file, "%lf", &z));
        Vector rp={x/lambda, y/lambda, z/lambda};
        assert(fscanf(file, "%d", &portNumber));
        Vector Lm=subVector(rn, rm);
        Vector Lp=subVector(rp, rn);
        double lm=magVector(Lm);
        double lp=magVector(Lp);
        if (portNumber>0){
            Basis newBasis={rm, rn, rp, Lm, Lp, lm, lp, 1, portNumber};
            myShape->basisList[i] = newBasis;
        }else{
            Basis newBasis={rm, rn, rp, Lm, Lp, lm, lp, 0, portNumber};
            myShape->basisList[i] = newBasis;
        }
    }
    fclose(file);   
    myShape->NPorts = NPorts;
    myShape->portsList = (Port*)malloc(NPorts*sizeof(Port));
    myShape->a = a/lambda;
}

void setPort(Shape *myShape, Port newPort){
    assert(myShape->portsList!=NULL);
    assert(portCounter<myShape->NPorts);
	newPort.Vin/=myShape->lambda;
    myShape->portsList[portCounter] = newPort;
    portCounter++;
}

void deleteShape(Shape *myShape){
    assert(myShape->NBasis>0);
    assert(myShape->basisList!=NULL);
    myShape->NBasis = 0;
    free(myShape->basisList);
    myShape->basisList = NULL;
    if (myShape->NPorts>0){
        free(myShape->portsList);
        myShape->portsList  = NULL;
    }
    myShape->NPorts = 0;
    portCounter = 0;
}

void logShape(Shape *myShape){
    assert(myShape->NBasis>0);
    assert(myShape->basisList!=NULL);
    FILE *file=fopen("log/shapeLog.txt", "w");
    fprintf(file, "Number of basis is %d\n", myShape->NBasis);
    for (int i=0; i<myShape->NBasis; i++){
        fprintf(file, "Basis %d: Port %d\n", i+1, myShape->basisList[i].portNumber);
        fprintf(file, "rm: ");
        fprintf(file, "%21.14E %21.14E %21.14E\n", myShape->basisList[i].rm.x, 
            myShape->basisList[i].rm.y, myShape->basisList[i].rm.z);
        fprintf(file, "rn: ");
        fprintf(file, "%21.14E %21.14E %21.14E\n", myShape->basisList[i].rn.x, 
            myShape->basisList[i].rn.y, myShape->basisList[i].rn.z);
        fprintf(file, "rp: ");
        fprintf(file, "%21.14E %21.14E %21.14E\n", myShape->basisList[i].rp.x, 
            myShape->basisList[i].rp.y, myShape->basisList[i].rp.z);
        fprintf(file, "Lm: ");
        fprintf(file, "%21.14E %21.14E %21.14E\n", myShape->basisList[i].Lm.x, 
            myShape->basisList[i].Lm.y, myShape->basisList[i].Lm.z);
        fprintf(file, "Lp: ");
        fprintf(file, "%21.14E %21.14E %21.14E\n", myShape->basisList[i].Lp.x, 
            myShape->basisList[i].Lp.y, myShape->basisList[i].Lp.z);
        fprintf(file, "Lm = %21.14E, Lp = %21.14E \n", myShape->basisList[i].lm,
            myShape->basisList[i].lp);
    }
    fprintf(file, "====================================================\n");
    fprintf(file, "Number of ports is %d\n", myShape->NPorts);
    for (int i=0; i<myShape->NPorts; i++){
        for (int n=0; n<myShape->NBasis; n++){
            if ((myShape->basisList[n].portNumber==myShape->portsList[i].portNumber)&&
                (myShape->basisList[n].isPort)){
                fprintf(file, "Port %d at basis %d\n", myShape->portsList[i].portNumber, n+1);
                fprintf(file, "Vin = (%21.14E, %21.14E)\n", creal(myShape->portsList[i].Vin),
                    cimag(myShape->portsList[i].Vin));
                fprintf(file, "ZL = (%21.14E, %21.14E)\n", creal(myShape->portsList[i].ZL),
                    cimag(myShape->portsList[i].ZL));
            }
        }
    }
    fclose(file);
}

// Shapes Library

void createVerticalDipoleDeltaGapCenter(double L, double lc){
    FILE *file=fopen("Mesh/shape.py", "w");
    assert(file!=NULL);
    //
    fprintf(file, "import meshLib\n");
    fprintf(file, "tol = %21.14E\n", lc/1000.0);
    fprintf(file, "meshLib.setTol(tol)\n");
    fprintf(file, "L = %21.14E\n", L);
    fprintf(file, "lc = %21.14E\n", lc);
    fprintf(file, "\n");
    fprintf(file, "v1 = meshLib.Vertex(+0.0, +0.0, -L/2.0)\n");
    fprintf(file, "v2 = meshLib.Vertex(+0.0, +0.0, +L/2.0)\n");
    fprintf(file, "line1 = meshLib.Line(v1, v2, lc, 1)\n");
    fprintf(file, "\n");
    fprintf(file, "shape = meshLib.Shape()\n");
    fprintf(file, "shape.addLines([line1])\n");
    fprintf(file, "shape.getBasis()\n");
    // 
    fclose(file);
    //
#ifdef linux
    assert(!system("python3 Mesh/shape.py"));
#endif

#ifdef _WIN32
	assert(!system("python Mesh/shape.py"));
#endif
}

void createCircularLoopDeltaGapCenter(double r, double lc){
    FILE *file=fopen("Mesh/shape.py", "w");
    assert(file!=NULL);
    // v1, v2, portIndex, index
    fprintf(file, "import meshLib\n");
    fprintf(file, "import numpy as np\n");
    fprintf(file, "tol = %21.14E\n", lc/1000.0);
    fprintf(file, "Nround = int(-np.log10(tol))\n");
    fprintf(file, "meshLib.setTol(tol)\n");
    fprintf(file, "r = %21.14E\n", r);
    fprintf(file, "lc = %21.14E\n", lc);
    fprintf(file, "N = %d\n", (int)(2.0*pi*r/lc));
    fprintf(file, "dphi = 2.0*np.pi/N\n");
    fprintf(file, "\n");
    fprintf(file, "shape = meshLib.Shape()\n");
    fprintf(file, "segmentsList = []\n");
    fprintf(file, "for i in range(N):\n");
    // fprintf(file, "\tv1 = meshLib.Vertex(np.round(r*np.cos((i+0)*dphi), Nround), np.round(r*np.sin((i+0)*dphi), Nround), 0.0)\n");
    // fprintf(file, "\tv2 = meshLib.Vertex(np.round(r*np.cos((i+1)*dphi), Nround), np.round(r*np.sin((i+1)*dphi), Nround), 0.0)\n");
    fprintf(file, "\tv1 = meshLib.Vertex(r*np.cos((i+0)*dphi), r*np.sin((i+0)*dphi), 0.0)\n");
    fprintf(file, "\tv2 = meshLib.Vertex(r*np.cos((i+1)*dphi), r*np.sin((i+1)*dphi), 0.0)\n");
    // fprintf(file, "\tv1 = meshLib.Vertex(r*np.cos((i+0)*dphi), 0.0, r*np.sin((i+0)*dphi))\n");
    // fprintf(file, "\tv2 = meshLib.Vertex(r*np.cos((i+1)*dphi), 0.0, r*np.sin((i+1)*dphi))\n");
    fprintf(file, "\tif i==0 or i==N-1:\n");
    fprintf(file, "\t\tsegmentsList.append(meshLib.Segment(v1, v2, 1, i))\n");
    fprintf(file, "\telse:\n");
    fprintf(file, "\t\tsegmentsList.append(meshLib.Segment(v1, v2, 0, i))\n");
    fprintf(file, "\tcontinue\n");
    fprintf(file, "\n");
    fprintf(file, "shape.addSegments(segmentsList)\n");
    fprintf(file, "shape.getBasis()\n");
    // 
    fclose(file);
    //
#ifdef linux
    assert(!system("python3 Mesh/shape.py"));
#endif

#ifdef _WIN32
	assert(!system("python Mesh/shape.py"));
#endif
}

void createTransmissionLineDeltaGapCenter(double L, double H, double lc){
    FILE *file=fopen("Mesh/shape.py", "w");
    assert(file!=NULL);
    //
    fprintf(file, "import meshLib\n");
    fprintf(file, "tol = %21.14E\n", lc/1000.0);
    fprintf(file, "meshLib.setTol(tol)\n");
    fprintf(file, "L = %21.14E\n", L);
    fprintf(file, "H = %21.14E\n", H);
    fprintf(file, "lc = %21.14E\n", lc);
    fprintf(file, "\n");
    fprintf(file, "v1 = meshLib.Vertex(-L/2.0, +0.0, -H)\n");
    fprintf(file, "v2 = meshLib.Vertex(-L/2.0, +0.0, +H)\n");
    fprintf(file, "v3 = meshLib.Vertex(+L/2.0, +0.0, +H)\n");
    fprintf(file, "v4 = meshLib.Vertex(+L/2.0, +0.0, -H)\n");
    fprintf(file, "line1 = meshLib.Line(v1, v2, lc, 1)\n");
    // fprintf(file, "line1 = meshLib.Line(v1, v2, lc/10.0, 1)\n");
    fprintf(file, "line2 = meshLib.Line(v2, v3, lc, 0)\n");
    fprintf(file, "line3 = meshLib.Line(v3, v4, lc, 2)\n");
    // fprintf(file, "line3 = meshLib.Line(v3, v4, lc/10.0, 2)\n");
    fprintf(file, "line4 = meshLib.Line(v4, v1, lc, 0)\n");
    fprintf(file, "\n");
    fprintf(file, "shape = meshLib.Shape()\n");
    fprintf(file, "shape.addLines([line1, line2, line3, line4])\n");
    fprintf(file, "shape.getBasis()\n");
    // 
    fclose(file);
    //
#ifdef linux
    assert(!system("python3 Mesh/shape.py"));
#endif

#ifdef _WIN32
	assert(!system("python Mesh/shape.py"));
#endif
}

void createYagiAntenna(YagiElements *elements, double segments_per_lambda){
    FILE *file=fopen("Mesh/shape.py", "w");
    assert(file!=NULL);
    //
    fprintf(file, "import meshLib\n");
    int N=elements->NElements;
    fprintf(file, "lines= []\n");
    double lc_max=0;
    for (int i=0; i<N; i++){
        fprintf(file, "L = %21.14E\n", elements->length[i]);
        fprintf(file, "x = %21.14E\n", elements->xLocation[i]);
        fprintf(file, "lc = %21.14E\n", elements->length[i]/segments_per_lambda);
        lc_max<elements->length[i]/segments_per_lambda ? lc_max=elements->length[i]/segments_per_lambda : lc_max;
        fprintf(file, "\n");
        fprintf(file, "v1 = meshLib.Vertex(x, +0.0, -L/2.0)\n");
        fprintf(file, "v2 = meshLib.Vertex(x, +0.0, +L/2.0)\n");
        if (elements->isExcitaiton[i]==1){
            fprintf(file, "lines.append(meshLib.Line(v1, v2, lc, 1))\n");
        }else{
            fprintf(file, "lines.append(meshLib.Line(v1, v2, lc, 0))\n");
        }
    }
    fprintf(file, "tol = %21.14E\n", lc_max/1000.0);
    fprintf(file, "meshLib.setTol(tol)\n");
    
    fprintf(file, "\n");
    fprintf(file, "shape = meshLib.Shape()\n");
    fprintf(file, "shape.addLines(lines)\n");
    fprintf(file, "shape.getBasis()\n");
    // 
    fclose(file);
    //
#ifdef linux
    assert(!system("python3 Mesh/shape.py"));
#endif

#ifdef _WIN32
	assert(!system("python Mesh/shape.py"));
#endif
}