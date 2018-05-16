#ifndef FVMGUD_H
#define FVMGUD_H 
#include "algebralineal.h"
int ConveccionDifusion1D(int nn, double **mnodos,double v, double k, double u0, double uf, 
    double *Solucion);
int ConveccionTransitoria1D(int nn, double **nnodos,double dt,int tmax, double v, double u0,
    int(g0(int,double*)),double *Solucion);

int inicial(int n, double *u);
#endif 
