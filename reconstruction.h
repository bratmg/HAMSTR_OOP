// ##################################################################
//
// reconstruction.h
//
// ##################################################################
#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include "math.h"

namespace RECONSTRUCTION
{

// ##################################################################
//
// MUSCL RECONSTRUCTION
//
// ##################################################################
   extern void MUSCL(double **f,
               double **ql,
               double **qr,
               double **f2,
               int is,
               int ie,
               double th,
               double qt,
               double eps,
               int imax,
               int nq);

// ##################################################################
//
// WENO 5 RECONSTRUCTION
//
// ##################################################################

   extern double WENO5(double a,  //fluxes at each edge
                double b, // left states
                double c, // right states
                double d,   // start index of chain
                double e,   // end index of chain
                double epsw);

// ##################################################################
//
// WENO RECONSTRUCTION (WRAPPER FUNCTION)
//   
// ##################################################################

   extern void WENO(double **f,  //fluxes at each edge
             double **ql, // left states
             double **qr, // right states
             int    is,   // start index of chain
             int    ie,   // end index of chain
             double epsw,  // epsilon in limiter
             int    imax, // chainSize
             int    nq);   // number of variables

// ##################################################################
//
// MUSCL RECONSTRUCTION DERIVATIVE
//
// ##################################################################

   extern void MUSCL_DERIV(double **f,
                    double **ql,
                    double **qr,
                    double **f2,
                    double **dq,
                    int      is,
                    int      ie,
                    double   th,
                    double   qt,
                    double   eps,
                    int      imax,
                    int      nq);

//##########################################################
//
// WENO RECONSTRUCTION DERIVATIVE 
//
//##########################################################
   extern double WENO_5(double a,  //fluxes at each edge
                double b, // left states
                double c, // right states
                double d,   // start index of chain
                double e,   // end index of chain
                double a_d,
                double b_d, // left states
                double c_d, // right states
                double d_d,   // start index of chain
                double e_d,   // end index of chain
                double epsw);
  
// ##################################################################
//
// WENO RECONSTRUCTION DERIVATIVE (WRAPPER)
//
// ##################################################################
  extern void WENO_DERIV(double **f,
                  double **ql,
                  double **qr,
                  double **dq,
                  int      is,
                  int      ie,
                  double   eps,
                  int      imax,
                  int      nq);

};

#endif
// ##################################################################
// END OF FILE
// ##################################################################
