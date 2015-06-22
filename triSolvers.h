// ##################################################################
//
// triSolvers.h
//
// contains routines to set up the initial mesh and re-organizing
// the strcutures to make is cache friendly
// ##################################################################
#ifndef TRISOLVERS_H
#define TRISOLVERS_H
// ##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// J. Sitaraman
// last updated 06/13
// Penta diagonal solver

namespace TRISOLVER
{

  void pentSolve(double *a,double *b,double *c,double *d,double *e,double *f,int N);

  //general 5*5 inversion
  void matInv5(double *f1, double *f2, double *f3,double *f4,double *f5,int nq);

  //tridiagonal block inversion
  // nq = 5
  void blockInv(double ***a,double ***b, double ***c,double **f, int N,int nq); 

  //tridiagonal block inversion
  // nq = 5
  void blockTridag(double ***a,double ***b, double ***c,double **f, int N,int nq);    

  //tridiagonal block inversion
  // nq = 5
  void blockTridagPeriodic(double ***a,double ***b, double ***c,double **f, int N,int nq);  

  //tridiagonal block inversion
  // nq = 4
  void blockTridag4(double ***a,double ***b, double ***c,double **f, int N,int nq);  
  
  //tridiagonal block inversion
  // nq = 4
  void blockTridagPeriodic4(double ***a,double ***b, double ***c,double **f, int N,int nq);  
  
}
#endif
// ##################################################################
// END OF FILE
// ##################################################################