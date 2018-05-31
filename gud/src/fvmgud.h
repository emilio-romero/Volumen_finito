#ifndef FVMGUD_H
#define FVMGUD_H 
#include "algebralineal.h"
int ConveccionDifusion1D(int nn, double **mnodos,double v, double k, double u0, double uf, 
    double *Solucion);
int ConveccionTransitoria1D(int nn, double **nnodos,double dt,int tmax, double v, double u0,
    int(g0(int,double*)),double *Solucion);

int ConveccionTransitoriaImp1D(int nn, double **nnodos,double dt,int tmax, double v, double u0,
    int(g0(int,double*)),double *Solucion);

int ConveccionTransitoria2D(int nn, double **nnodos,double dt,int tmax, double vx,double vy, 
    double *x0,double *y0,int(g0(int,double**)),double *Solucion);
int inicial(int n, double *u);
int inicial2(int n, double *u);

int inicio2d(int n, double **u);

int ConveccionNolineal1D(int nn, double **nnodos, double dt, double tmax,
    double u0, double(*f)(double), int(*g0)(int,double**,double*));

double fburger(double u); 
int g0burger(int n, double **nx, double *out);
#endif 
