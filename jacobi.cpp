#include "vars.h"

extern int nump;
extern int numiter;

extern int *displ;
extern int *ind;
extern double *val;
extern double *diag;

extern double *ref;

void jacobi(double *x, double *b){

  int numproc;
  int myid;
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  FILE *filef = fopen("/gpfs/alpine/scratch/merth/csc362/jacobi.txt","w");

  double *r = new double[nump];

  double l0b = 0;
  double l1b = 0;
  double l2b = 0;
  for(int n = 0; n < nump; n++){
    double bn = abs(b[n]);
    l2b += pow(bn,2);
    l1b += bn;
    if(bn>l0b)l0b=bn;
  }
  l2b = sqrt(l2b);
  double l0f = 0;
  double l1f = 0;
  double l2f = 0;
  for(int n = 0; n < nump; n++){
    double fn = abs(ref[n]);
    l2f += pow(fn,2);
    l1f += fn;
    if(fn>l0f)l0f=fn;
  }
  l2f = sqrt(l2f);

  double kernel1 = 0;
  double kernel2 = 0;
  double kernel3 = 0;
  double timet;
  MPI_Barrier(MPI_COMM_WORLD);
  double timesol = MPI_Wtime();
  {
    for(int n = 0; n < nump; n++)
      x[n] = 0;
    if(myid==0)fprintf(filef,"%e %e %e %e %e %e\n",1.,1.,1.,1.,1.,1.);
    if(myid==0)printf("%e %e %e %e %e %e\n",1.,1.,1.,1.,1.,1.);
    for(int iter = 1; iter <= numiter; iter++){
      timet = MPI_Wtime();
      #pragma omp parallel for
      for(int m = 0; m < nump; m++){
        double reduce = 0;
        for(int n = displ[m]; n < displ[m+1]; n++)
          if(ind[n]!=m)
            reduce += val[n]*x[ind[n]];
        r[m] = (b[m]-reduce)/diag[m];
      }
      double *temp = x;
      x = r;
      r = temp;
      kernel1 = kernel1 + MPI_Wtime()-timet;
      timet = MPI_Wtime();
      #pragma omp parallel for
      for(int m = 0; m < nump; m++){
        double reduce = 0;
        for(int n = displ[m]; n < displ[m+1]; n++)
          reduce += val[n]*x[ind[n]];
        r[m] = reduce;
      }
      kernel2 = kernel2 + MPI_Wtime()-timet;
      timet = MPI_Wtime();
      //ERROR
      double l0r = 0;
      double l1r = 0;
      double l2r = 0;
      for(int n = 0; n < nump; n++){
        double rn = abs(r[n]-b[n]);
        l2r += pow(rn,2);
        l1r += rn;
        if(rn>l0r)l0r=rn;
      }
      l2r = sqrt(l2r);
      double l0x = 0;
      double l1x = 0;
      double l2x = 0;
      for(int n = 0; n < nump; n++){
        double xn = abs(x[n]-ref[n]);
        l2x += pow(xn,2);
        l1x += xn;
        if(xn>l0x)l0x=xn;
      }
      l2x = sqrt(l2x);
      if(myid==0)fprintf(filef,"%e %e %e %e %e %e\n",l0r/l0b,l1r/l1b,l2r/l2b,l0x/l0f,l1x/l1f,l2x/l2f);
      if(myid==0)printf("%e %e %e %e %e %e\n",l0r/l0b,l1r/l1b,l2r/l2b,l0x/l0f,l1x/l1f,l2x/l2f);
      kernel3 = kernel3 + MPI_Wtime()-timet;
    }
    fclose(filef);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  timesol = MPI_Wtime()-timesol;
  if(myid==0){
    printf("\nJACOBI SOLVER\n");
    printf("KERNEL1 TIME       : %e s\n",kernel1);
    printf("KERNEL2 TIME       : %e s\n",kernel2);
    printf("KERNEL3 TIME       : %e s\n",kernel3);
    printf("SOLUTION TIME      : %e s\n",timesol);
    printf("\n");
  }
}
