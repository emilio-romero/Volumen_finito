#include "malla.h"

int malla1D(double x0, double xf,int nn, int **conectividad, double **nodos){
  double dx=(xf-x0)/((double)(nn-1));
  for(int i=0;i<nn;++i){
    nodos[i][0]=x0+(double)i*dx;
    //printf("nodo %d:%lf\n",i,nodos[i][0]);
  }
return(0);}
