#include "algebralineal.h"

/*
 * Crea un vector de tamanio n y llena de ceros y
 * libera vectores (no necesita la cantidad de elementos)
 */
double *crear_vector(int n){
  double *aux=(double*)malloc(n*sizeof(double));
  for(int i=0;i<n;i++) aux[i]=0.0; 
  return(aux);
}

//int liberar_vector();
// Se omitio porque un free() es mas sencillo 
/*
 * Crea una matriz de tamanio nr x nc y llena de ceros y
 * libera matrices se le da la cantidad de filas
 */
double **crear_matriz(int nr, int nc){
  double **aux=(double**)malloc(nr*sizeof(double*));
  for(int i=0;i<nr;i++) aux[i]=(double*)calloc(nc,sizeof(double));
  return(aux);
}

int liberar_matriz(double **x, int nr){
  for(int i=0;i<nr;i++) free(x[i]); 
  free(x);
return(1);}

/*
* Producto punto de dos vectores de tamanio n 
*/

double punto(double *a, double *b, int n){
  double suma=0; 
  for(int i=0;i<n;i++){
    suma+=a[i]*b[i]; 
  }
return(suma);
}

/*
 * Suma de vectores 
*/
int vector_suma(double *x, double *y, int n, double *out){
  for(int i=0;i<n;i++)
    out[i]=x[i]+y[i];
return(1);}

/*
 * Resta de vectores 
*/
int vector_resta(double *x, double *y, int n, double *out){
  for(int i=0;i<n;i++)
    out[i]=x[i]-y[i];
return(1);}

/*
 *  Vector por escalar
*/
int vector_escalar(double a, double *x, int n, double *out){
  for(int i=0;i<n;i++)
    out[i]=a*x[i];
return(1);}

/*
 * Copia un vector en otro 
*/
int vector_copiar(double *original, int n, double *copia){
  for(int i=0;i<n;i++)
    copia[i]=original[i];
return(1);}
/*
 * Multiplicacion de un vector por otro transpuesto
 * x * y^t 
 * devuelve una matriz 
 */
int vector_vector_mul(double *x, double *y, int n, double **out){
  for(int i=0;i<n;++i){
    for(int j=0;j<n;++j){
      out[i][j]=x[i]*y[j];
    }
  }
return(1);}
/* Se hace la suma y mulplicacion 
*  de tal manera que 
* z=a+yx
* con z=out (vector), a (vector), y (escalar), x(vector) 
*/ 
int vector_zapyx(double *a, double y, double *x, int n, double *out){
 for(int i=0;i<n;++i){
   out[i]=a[i]+y*x[i];
 } 
return(1);}

/*
* Producto de una matriz por un vector, devuelve un vector 
* m es para filas de la matriz y n para las columnas de la matriz
* n tambien son las filas del vector
*/
int matriz_vector_mul(double **A, double *b, int m, int n, double *out){
  double suma=0; 
  for(int i=0;i<m;i++){
    suma=0;
    for(int j=0;j<n;j++){
      suma+=A[i][j]*b[j]; 
    }
    out[i]=suma;
  }
return 1;}

/*
* Producto de matrices, devuelve una matriz 
* l y m son las filas y columnas de la primera matriz
* m y n son las filas y columnas de la segunda matriz
*/

int matriz_mul(double **A, double **B, int l, int m, int n, double **out){
  double suma=0; 

  for(int i=0;i<l;i++){
    for(int j=0;j<n;j++){
      suma=0; 
      for(int k=0;k<m;k++){
        suma+=A[i][k]*B[k][j];
      }
      out[i][j]=suma; 
    }
  }

return 1;}
/*
* Suma de matrices
*/
int matriz_suma(double **A, double **B, int nr, int nc, double **out){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      out[i][j]=A[i][j]+B[i][j]; 
    }
  }
return(1);}
/*
* Resta de matrices
*/
int matriz_resta(double **A, double **B, int nr, int nc, double **out){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      out[i][j]=A[i][j]-B[i][j]; 
    }
  }
return(1);}
/*
* Copia una matriz A a otra B
*/
int matriz_copiar(double **original, int nr, int nc, double **copia){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      copia[i][j]=original[i][j];
    }
  }
return(1);}
/*
* Transpone la matriz A a otra B
* nr son las filas de B 
* nc son las columnas de B
*/
int matriz_transponer(double **original, int nr, int nc, double **copia){
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      copia[i][j]=original[j][i];
    }
  }
return(1);}

/*
 * Multiplica todos los elementos de la matriz A, por el escalar a 
 */
int matriz_escalar(double a,double **A, int nr, int nc, double **out){
  for(int i=0;i<nr;++i){
    for(int j=0;j<nc;++j){
      out[i][j]=a*A[i][j];
    }
  }
return(1);}

/*
 * Llena la matriz out de ceros 
 */
int matriz_ceros(int nr, int nc, double **out){
  for(int i=0;i<nr;++i){
    for(int j=0;j<nc;++j){
      out[i][j]=0.0; 
    }
  }
return(1);}
/*
 * Calcula la inversa de la matriz  A de nxn 
 *
 */

int matriz_inversa(double **A, int n, double **out){
  double *auxv=crear_vector(n);
  double *sole=crear_vector(n);
  for(int e=0;e<n;++e){
    auxv[e]=1.0;
    solLU(A,auxv,n,n,sole);
    for(int j=0;j<n;++j){
      out[j][e]=sole[j]; 
    }
    auxv[e]=0.0;
  }
  free(auxv);
  free(sole);
return(1);}

/*
* Llena una matriz out con unos en la diagonal y ceros en las demas partes
*/ 
int matriz_identidad(int nr, int nc, double **out){
  
  for(int i=0;i<nr;++i){ 
    for(int j=0;j<nc;++j){
      out[i][j]=0.0;
    }
    out[i][i]=1.0;
  }
return(1);}

/*
 * Suma (o resta de acuerdo al signo) lamb en la diagonal de A
 *
 */

int matriz_diag_suma(double lamb,double **A, int nr, int nc,double **out){

  for(int i=0;i<nr;++i){
    out[i][i]=A[i][i]+lamb; 
  }

return(1);}

/*
* Calcula la norma 1 de una matriz 
*/
int Norma_1_matriz(double **A, int nr, int nc, double *out){
  double max=0; 
  double aux=0; 
  for(int i=0;i<nr;i++) 
    max+=fabs(A[i][0]); 
  
  for(int j=1;j<nc;j++){
    for(int i=0;i<nr;i++) aux+=fabs(A[i][j]); 
    if(aux>max){
      max=aux; 
    }
    aux=0; 
  }
*out=max; 
}


/*
 * Calcula la norma 1 de un vector := suma(abs(x_i))
 */
double Norma_1_vector(double *x, int n){
  double aux=0; 
  for(int i=0;i<n;i++)
    aux+=fabs(x[i]);
return(aux);} 

/*
 * Calcula la norma 2 de un vector := sqrt(suma(x_i^2))
 */
double Norma_2_vector(double *x, int n){
  double aux=0; 
  for(int i=0;i<n;++i){
    aux=aux+x[i]*x[i];
  }
  aux=sqrt(aux);
return(aux);} 

/*
 * Calcula la norma infinito de un vector := max |x_i|
 */
double Norma_inf_vector(double *x, int n){
  double aux=fabs(x[0]); 
  for(int i=1;i<n;i++){
    if(aux<fabs(x[i]))
      aux=fabs(x[i]);
  } 
return(aux);} 

int es_spd(char *cfile){
  int nr, nc; 
  double **A=readMatrix(cfile, &nr, &nc);
  double **AA=(double**)malloc(nr*sizeof(double*)); 
  double N1;
  double meps=sqrt(DBL_EPSILON);
  for(int i=0;i<nr;i++) AA[i]=(double*)malloc(nc*sizeof(double));
    
  for(int i=0;i<nr;i++){
    for(int j=0;j<nc;j++){
      AA[i][j]=A[i][j]-A[j][i];
    }
  }
  
  Norma_1_matriz(AA,nr,nc,&N1);
  //printf("%g\n",N1);
  if(N1>meps){
    printf("La matriz dada no es simetrica\n");
    printf("Se debe utilizar la matriz (A+A^T)/2\n");
    for(int i=0;i<nr;i++) free(AA[i]); 
    free(AA);
    freeMatrix(A);
    return(0);
  } else{
    printf("La matriz es simetrica, hurra!\n");
    for(int i=0;i<nr;i++) free(AA[i]); 
    free(AA);
    freeMatrix(A);
    return(1);
  }
  /*Liberacion de memoria*/
return 1;}
/* Calcula el numero de condicion \kappa de una matriz A de nxn 
 *
 */
double numero_condicion(double **A, int n){
  double aux; 
  double **Ai=crear_matriz(n,n);
  double nA,nAi; 
  Norma_1_matriz(A,n,n,&nA);
  matriz_inversa(A,n,Ai);
  Norma_1_matriz(Ai,n,n,&nAi); 
  aux=nA*nAi;
return(aux);}



/*
 * Solucionadores basicos: matriz inferior, matriz superior 
 */
int sinferior(double **L, double *b, int n, double *out){
for(int i=0;i<n;i++) out[i]=0; 
double tol=1e-9; 
double suma; 
for(int i=0;i<n;i++){
  if(fabs(L[i][i])<tol){
    L[i][i]=L[i][i]+0.01;
    //printf("\nEl sistema no tiene soluciones, se perturbo el sistema\n");
  }
}

out[0]=b[0]/L[0][0];
for(int k=1;k<n;k++){
  suma=0;
  for(int j=0;j<k;j++){
    suma+=L[k][j]*out[j]; 
  }
  out[k]=(b[k]-suma)/L[k][k];
}

return(1);}

int ssuperior(double **U,double *b, int n, double *out){
for(int i=0;i<n;i++) out[i]=0; 
double tol=10E-10, suma; 
for(int i=0;i<n;i++){
  if(fabs(U[i][i])<tol){
   U[i][i]=U[i][i]+0.01;
   // printf("\nEl sistema no tiene soluciones\n");
  }
}
out[n-1]=b[n-1]/U[n-1][n-1];
for(int i=n-2;i>=0;i--){
  suma=0;
  for(int j=i;j<n;j++){
    suma+=U[i][j]*out[j]; 
  }
  out[i]=(b[i]-suma)/U[i][i];
}

return(1);}



/*
* Factorizacion LU
* esta funcion devuelve las matrices L y U de la factorizacion
*/
int factoLU(double **A, int n, double **L, double **U){
  double tol=1E-10, sumal, sumau; 
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      L[i][j]=0; 
      U[i][j]=0;
    }
    U[i][i]=1.0; 
    L[i][0]=A[i][0]; 
    U[0][i]=A[0][i]/L[0][0];
  }

  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      sumal=0; sumau=0;
      if(i>=j){
        for(int k=0;k<j;k++) sumal+=L[i][k]*U[k][j]; 
        L[i][j]=A[i][j]-sumal; 
      }

      else{
        for(int k=0;k<i;k++) sumau+=L[i][k]*U[k][j]; 
        U[i][j]=(A[i][j]-sumau)/L[i][i]; 
      }
    }
  }

return(1);}

int solLU(double **A,double *b, int nr, int nc,double *out){
  double **L=(double**)malloc(nr*sizeof(double*));
  double **U=(double**)malloc(nr*sizeof(double*));
  double *vaux=(double*)malloc(nr*sizeof(double));
  for(int i=0;i<nr;i++) {
     L[i]=(double*)malloc(nc*sizeof(double));
     U[i]=(double*)malloc(nc*sizeof(double));
  }
  factoLU(A,nr,L,U); //Calculo de L y U
  //Resolver el sistema 
  sinferior(L,b,nr,vaux);
  ssuperior(U,vaux,nr,out);

  for(int i=0;i<nr;i++){ 
    free(L[i]);
    free(U[i]);
  }
  free(L);
  free(U);
  free(vaux);
return(1);}


/*
* Factorizacion Cholesky
* esta funcion devuelve la matriz L de la factorizacion
*/


int Chol(double **A, int n, double **out){
  double restakj=0; 
  double restaikj=0; 
  double meps=sqrt(DBL_EPSILON);
  out[0][0]=sqrt(A[0][0]);
  if(out[0][0]<meps){
    printf("No se ha podido factorizar\n");
    return(0);
  }
  for(int i=1;i<n;i++){
    out[i][0]=(A[i][0])/out[0][0];
  }
  for(int j=1;j<n;j++){
    for(int k=0;k<j;k++){
      restakj+=out[j][k]*out[j][k]; 
    }
    out[j][j]=sqrt(A[j][j]-restakj);
    if(out[j][j]<meps || (A[j][j]-restakj)<0){
      //printf("No se ha podido factorizar\n");
      return(0);
    }
    restakj=0;
    for(int i=j+1;i<n;i++){
      for(int k=0;k<j;k++){
        restaikj+=out[i][k]*out[j][k]; 
      }
      out[i][j]=(A[i][j]-restaikj)/out[j][j];
      restaikj=0; 
    }
  }

return(1);}

int Cholesky(char *cfile, double **out){
  int simetria=es_spd(cfile); 
  int nr, nc; 
  int df, contador=0;
  double error;
  double **A=readMatrix(cfile, &nr, &nc);
  double **AA=(double**)malloc(nr*sizeof(double*)); 
  double **LL=(double**)malloc(nr*sizeof(double*)); 
  for(int i=0;i<nr;i++){ AA[i]=(double*)malloc(nc*sizeof(double));
     LL[i]=(double*)malloc(nc*sizeof(double));}
  double meps=sqrt(DBL_EPSILON);

  if(simetria==0){
    for(int i=0;i<nr;i++){
      for(int j=0;j<nc;j++){
        AA[i][j]=(A[i][j]+A[j][i])/2.0; 
      }
    }
   matriz_copiar(AA,nr,nc,A);
  }

  df=Chol(A,nr,out);
  //printf("%d\n",df);
  while(df==0){
    for(int i=0;i<nr;i++){
      A[i][i]=A[i][i]+0.01; 
    }
    df=Chol(A,nr,out); 
  //printf("%d\n",df);
    contador++; 
  }
  printMatrix(out,nr,nc);
  matriz_transponer(out,nr,nc,AA);
  matriz_mul(out,AA,nr,nc,nr,LL);
  matriz_resta(A,LL,nr,nc,AA);
  Norma_1_matriz(AA,nr,nc,&error);
  printf("La matriz se perturbo %d veces.\n",contador); 
  printf("El error es %g\n",error);
return(1);}


int solLL(double **A,double *b, int nr, int nc,double *out){
  double **L=(double**)malloc(nr*sizeof(double*));
  double **Lt=(double**)malloc(nr*sizeof(double*));
  double *vaux=(double*)malloc(nr*sizeof(double));
  for(int i=0;i<nr;i++) {
     L[i]=(double*)malloc(nc*sizeof(double));
     Lt[i]=(double*)malloc(nc*sizeof(double));
  }
  Chol(A,nr,L); //Calculo de L 
  matriz_transponer(L,nr,nc,Lt); //Calculo de L*
  //Resolver el sistema 
  sinferior(L,b,nr,vaux);
  ssuperior(Lt,vaux,nr,out);

  for(int i=0;i<nr;i++){ 
    free(L[i]);
    free(Lt[i]);
  }
  free(L);
  free(Lt);
  free(vaux);
return(1);}


/*
 * Gradiente conjugado
 *
 */
int GradienteConjugado(double **A, double *b, int nr, int nc,double tol,double *out){
  double *r=crear_vector(nr);
  double *p=crear_vector(nr);
  double *q=crear_vector(nr);
  double err,dk,ak,bk; 
  for(int i=0;i<nr;++i){
    r[i]=b[i]; 
    p[i]=b[i];
  }
  err=sqrt(punto(r,r,nr)/(double)nr);
  for(int k=0;k<nr;++k){
    matriz_vector_mul(A,p,nr,nr,q);
    dk=punto(r,r,nr);
    ak=dk/(punto(p,q,nr));
    for(int i=0;i<nr;++i){
      out[i]=out[i]+ak*p[i];
      r[i]=r[i]-ak*q[i];
    }
    bk=punto(r,r,nr);
    for(int i=0;i<nr;++i){
      p[i]=r[i]+bk*p[i]; 
    }
    err=sqrt(punto(r,r,nr)/(double)nr);
    if(err<tol) break; 
  }
  free(r); free(p); free(q);
return(1);}


/*
* Minimos cuadrados 
* esta funcion devuelve los coeficientes del polinomio
*/

//Correciones por el cambio de LU faltan
double *aproximaPolinomio(int grado, double **data, int npuntos){
double *aux=(double*)malloc((grado+1)*sizeof(double)); 
double *y=(double*)malloc((grado+1)*sizeof(double));
double *g=(double*)malloc(npuntos*sizeof(double));
double **A=(double**)malloc((npuntos)*sizeof(double*));
double **At=(double**)malloc((grado+1)*sizeof(double*));
double **AA=(double**)malloc((grado+1)*sizeof(double*));
  for(int i=0;i<npuntos;i++){
    g[i]=data[i][1];
    A[i]=(double*)malloc((grado+1)*sizeof(double));
    for(int j=0;j<(grado+1);j++) A[i][j]=pow(data[i][0],j);
  }
  for(int i=0;i<(grado+1);i++){
    AA[i]=(double*)malloc((grado+1)*sizeof(double));
    At[i]=(double*)malloc(npuntos*sizeof(double));
    for(int j=0;j<npuntos;j++) At[i][j]=pow(data[j][0],i);
  }
  matriz_mul(At,A,(grado+1),npuntos,(grado+1),AA);
  matriz_vector_mul(At,g,(grado+1),npuntos,y);
  //factLU(AA,y,(grado+1));  
  free(y); free(g); 
  for(int i=0;i<npuntos;i++) free(A[i]);
  for(int i=0;i<=grado;i++) free(At[i]);
  for(int i=0;i<=grado;i++) free(AA[i]);
  free(A); free(At); free(AA);
return(aux);}

/*
 *Miscelanea
 *
 */
double randx(){
  double aux; 
  aux=rand()/((double)RAND_MAX); 
return(aux);}

int reduccionMatriz(double **A, double *a, int nc, int nr, int nnf, double **aux, double *auxv){
  int cont=0; 
  for(int i=0;i<nr;i++){
    if(cont==(nnf-1)) break;
    if(randx()<0.5){
      for(int j=0;j<nc;j++){
        aux[cont][j]=A[i][j];
      }
      auxv[cont]=a[i];
      cont++; 
    }
  }
return(1);}


double minimo2(double a, double b){
 if(a<b) return(a);
 else return(b);
}


double max_diag(double **A, int n){
  double aux=A[0][0]; 
  for(int i=1;i<n;++i){
    if(aux<A[i][i]){
      aux=A[i][i]; 
    }
  }
return(aux);}
