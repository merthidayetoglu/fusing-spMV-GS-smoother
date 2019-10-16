#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <limits>
#include <omp.h>

using namespace std;

void setpositions();

void spmv(double*,double*);

void jacobi(double*,double*);
void gaussseidel(double*,double*);
void gradientdescent(double*,double*);
void conjugategradient(double*,double*);

