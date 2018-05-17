#include <stdio.h>
#include <stdlib.h> 
#include "fvmgud.h"
#include "malla.h"
int main(int argc, char *argv[]){
  int n=7; 
  if(argc>1) n=atoi(argv[1]); else n=7; 
  double **nodos1=crear_matriz(n,1);
  double **nodos=crear_matriz(n*n,2);
  int **mc; 
  double *x0=crear_vector(n);
  double *y0=crear_vector(n);
  for(int i=0;i<n;++i){
    //printf("nodo %d: %lf\n",i,nodos[i][0]);
  }
  double *sol=crear_vector(n);
  malla1D(-1.0,1.0,n,mc,nodos1);
  /* //Bloque de difusion conveccion 1D
  ConveccionDifusion1D(n,nodos1,2.5,0.1,1.0,0.0,sol);
  for(int i=0;i<n;++i)
    printf("%lf %lf\n",nodos1[i][0],sol[i]);
  */
  ConveccionTransitoriaImp1D(n,nodos1,0.01,5,1.0,0.0,inicial,sol);
  malla2D(-5,5,-5,5,n,mc,nodos);
 // ConveccionTransitoria2D(n,nodos,0.01,5,2.0,2.0,x0,y0,inicio2d,sol);
  /*double **tu=crear_matriz(n-1,n-1);
  inicio2d(n-1,tu);
  for(int i=0;i<n-1;++i){
    for(int j=0;j<n-1;++j){
      printf("%lf ",tu[i][j]);
    }printf("\n");
  }*/
  liberar_matriz(nodos1,n);
  liberar_matriz(nodos,n*n);
  free(sol);
  free(x0);
  free(y0);
  //printf("Los calculos han terminado\n");
return(0);}
