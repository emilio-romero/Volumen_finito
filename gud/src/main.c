#include <stdio.h>
#include <stdlib.h> 
#include "fvmgud.h"
#include "malla.h"
int main(int argc, char *argv[]){
  int n=7; 
  if(argc>1) n=atoi(argv[1]); else n=7; 
  double **nodos=crear_matriz(n,1);
  int **mc; 
  malla1D(-1.0,1.0,n,mc,nodos);
  for(int i=0;i<n;++i){
    //printf("nodo %d: %lf\n",i,nodos[i][0]);
  }
  double *sol=crear_vector(n);

  //ConveccionDifusion1D(n,nodos,2.5,0.1,1.0,0.0,sol);
 // for(int i=0;i<n;++i)
 //   printf("%lf %lf\n",nodos[i][0],sol[i]);
  ConveccionTransitoria1D(n,nodos,0.006,5,1.0,0.0,inicial,sol);
  liberar_matriz(nodos,n);
  free(sol);
//printf("Los calculos han terminado\n");
return(0);}
