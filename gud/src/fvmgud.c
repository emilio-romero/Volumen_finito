#include "fvmgud.h"

int ConveccionDifusion1D(int nn, double **nnodos, double v, double k, double u0, double uf,
    double *Solucion){
  int ni=nn-2; 
  double **A=crear_matriz(ni,ni);
  double *F=crear_vector(ni);
  double Pe=v/k, dx=nnodos[1][0]-nnodos[0][0]; 
  double a,b,c,d;
  double *s=crear_vector(ni);
  a=(-1.0/(dx*dx)+3*Pe/(8.0*dx)); 
  b=(3.0*Pe/(8.0*dx)+2.0/(dx*dx));
  c=(-7.0*Pe/(8.0*dx)-1.0/(dx*dx));
  d=Pe/(8.0*dx);
  A[0][0]=b-d; A[0][1]=a; F[0]=-(c+2.0*d)*u0;
  A[1][0]=c; A[1][1]=b; A[1][2]=a; F[1]=-d*u0;
  for(int i=2;i<ni-1;++i){
    A[i][i-2]=d; 
    A[i][i-1]=c; 
    A[i][i]=b; 
    A[i][i+1]=a; 
  }
  A[ni-1][ni-1]=b; A[ni-1][ni-2]=c; A[ni-1][ni-3]=d; 
  F[ni-1]=-a*uf; 
  
  solLU(A,F,ni,ni,s);
  for(int i=0;i<ni;++i){
    Solucion[i+1]=s[i];
  }
  Solucion[0]=u0; 
  Solucion[nn-1]=uf; 
  liberar_matriz(A,ni);
  free(F);
  free(s);
return(1);}


int ConveccionTransitoria1D(int nn, double **nnodos,double dt,int tmax, double v, double u0,
    int(g0(int,double*)),double *Solucion){
  int ni=nn-1; 
  double **A=crear_matriz(ni,ni);
  double *u=crear_vector(ni);
  double *up1=crear_vector(ni);
  double *aux=crear_vector(ni);
  double  dx=nnodos[1][0]-nnodos[0][0];
  double pc=v*dt/(8.0*dx);
  int limt=(int)((double)tmax/dt);
  g0(ni,u); 
  A[0][0]=2.0; A[0][1]=3.0;
  A[1][0]=-7.0; A[1][1]=3.0; A[1][2]=3.0;
  for(int i=2;i<ni-1;++i){
    A[i][i-2]=1.0; 
    A[i][i-1]=-7.0; 
    A[i][i]=3.0; 
    A[i][i+1]=3.0;
  }
  A[ni-1][ni-1]=-2.0; A[ni-1][ni-2]=-6.0; A[ni-11][ni-3]=1.0;
  

  for(int m=0;m<limt;++m){
    //printf("Tiempo: %lf\n",(double)(m+1)*dt);
    u[0]=u[0]+pc*5.0*u0; 
    u[1]=u[1]-pc*u0;
    matriz_vector_mul(A,u,ni,ni,aux);
    vector_escalar(-1.0*pc,aux,ni,aux);
    vector_suma(u,aux,ni,up1);
    printf("%lf %lf\n",nnodos[0][0],u0);
    for(int i=0;i<ni;++i){
      printf("%lf %lf\n",nnodos[i+1][0],up1[i]);
    }
    printf("\n\n");
    vector_copiar(up1,ni,u);
  }
  liberar_matriz(A,ni);
  free(u);
  free(up1);
  free(aux);
return(1);}

int inicial(int n, double *u){
  int nmedios=n/2;
  for(int i=0;i<n;++i){
    if(i<nmedios)
      u[i]=0.0;
    else 
      u[i]=5.0;
  }

return(1);}
