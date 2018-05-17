#ifndef MALLA_H
#define MALLA_H 
#include "algebralineal.h"

int malla1D(double x0,double xf,int ncv, int **conectividad, double **nodos);
int malla2D(double x0,double xf,double y0,double yf,int nn, 
    int **conectividad, double **nodos);
#endif 
