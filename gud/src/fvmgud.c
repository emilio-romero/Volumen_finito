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

int ConveccionTransitoriaImp1D(int nn, double **nnodos,double dt,int tmax, double v, double u0,
    int(g0(int,double*)),double *Solucion){
  int ni=nn-1; 
  double **A=crear_matriz(ni,ni);
  double *u=crear_vector(ni);
  double *up1=crear_vector(ni);
  double *aux=crear_vector(ni);
  double  dx=nnodos[1][0]-nnodos[0][0];
  double pc=8.0*dx/(v*dt);
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
  matriz_diag_suma(pc,A,ni,ni,A);

  for(int m=0;m<limt;++m){
    //printf("Tiempo: %lf\n",(double)(m+1)*dt);
    vector_escalar(pc,u,ni,u);
    u[0]=u[0]+5.0*u0; 
    u[1]=u[1]-u0;
    solLU(A,u,ni,ni,up1);
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


int Conveccion2D(){

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



int ConveccionTransitoria2D(int nn, double **nnodos,double dt,int tmax, double vx,double vy, 
    double *x0,double *y0,int(g0(int,double**)),double *Solucion){
  int ni=nn-1; 
  double **U=crear_matriz(ni,ni);
  double **up1=crear_matriz(ni,ni);
  double *aux=crear_vector(ni);
  double  dx=nnodos[nn][0]-nnodos[0][0];
  double  dy=nnodos[1][1]-nnodos[0][1];
  double alfa=vx*dt/(8.0*dx), beta=vy*dt/(8.0*dy);
  int limt=(int)((double)tmax/dt);
  g0(ni,U); 
 //printf("dx: %lf, alfa: %lf, beta: %lf\n",dx,alfa,beta); 

  for(int m=0;m<limt;++m){
    //printf("Tiempo: %lf\n",(double)(m+1)*dt);
    /*Parte del calculo anterior y parte de cambio en x*/
    for(int i=0;i<ni;++i){
      up1[i][0]=U[i][0]-alfa*(3.0*U[i][1]+2.0*U[i][0]-5.0*x0[i+1]);
      up1[i][1]=U[i][1]-alfa*(-7.0*U[i][0]+3.0*U[i][1]+3.0*U[i][2]+x0[i+1]);
      for(int j=2;j<ni-1;++j){
        up1[i][j]=U[i][j]-alfa*(3.0*U[i][j+1]+3.0*U[i][j]-7.0*U[i][j-1]+U[i][j-2]);
      }
      up1[i][ni-1]=U[i][ni-1]-alfa*(-2.0*U[i][ni-1]-6.0*U[i][ni-2]+U[i][ni-3]);
    }

    /*Parte de cambio en y*/
    for(int j=0;j<ni;++j){
      up1[0][j]=up1[0][j]-beta*(3.0*U[1][j]+2.0*U[0][j]-5.0*y0[j+1]);
      up1[1][j]=up1[1][j]-beta*(-7.0*U[0][j]+3.0*U[1][j]+3.0*U[2][j]+y0[j+1]);
      for(int i=2;i<ni-1;++i){
        up1[i][j]=up1[i][j]-beta*(3.0*U[i+1][j]+3.0*U[i][j]-7.0*U[i-1][j]+U[i-2][j]);
      }
      up1[ni-1][j]=up1[ni-1][j]-beta*(-2.0*U[ni-1][j]-6.0*U[ni-2][j]+U[ni-3][j]);
    }
   
   
    //for(int i=0;i<nn;++i){
    //  printf("%lf %lf %lf\n",nnodos[i][0],nnodos[i][1],x0[i]);
   // }
    for(int i=1;i<ni;++i){  
      for(int j=1;j<ni;++j){
        printf("%lf %lf %lf\n",nnodos[nn*i+1+j][0],nnodos[nn*i+1+j][1],up1[i-1][j-1]);
      }
    }
    printf("\n\n");
    matriz_copiar(up1,ni,ni,U);
  }
  liberar_matriz(U,ni);
  liberar_matriz(up1,ni);
  free(aux);
return(1);}

int inicio2d(int n, double **u){
  double mita=n/3;
  double dmita=2*n/3; 
  for(int i=mita;i<dmita;++i){
    for(int j=mita;j<dmita;++j){
      u[i][j]=15.0;
    }
  }
return(1);}
