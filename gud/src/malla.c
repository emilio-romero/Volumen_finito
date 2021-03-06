#include "malla.h"

int malla1D(double x0, double xf,int nn, int **conectividad, double **nodos){
  double dx=(xf-x0)/((double)(nn-1));
  for(int i=0;i<nn;++i){
    nodos[i][0]=x0+(double)i*dx;
    //printf("nodo %d:%lf\n",i,nodos[i][0]);
  }
return(0);}


int malla2D(double x0,double xf,double y0,double yf,int nn, 
    int **conectividad, double **nodos){
  double dx=(xf-x0)/((double)(nn-1));
  double dy=(yf-y0)/((double)(nn-1));
  for(int i=0;i<nn;++i){
    for(int j=0;j<nn;++j){
    nodos[nn*i+j][0]=x0+(double)i*dx;
    nodos[nn*i+j][1]=y0+(double)j*dy;
    }
  }
  
return(1);}
