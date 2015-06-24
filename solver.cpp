// ##################################################################
//
// solver.cpp
//
// ##################################################################
#include "solver.h"

using namespace CODEVARS;
// ##################################################################
//
// constuctor
//
// ##################################################################
SOLVER::SOLVER()
{
   //
   // variable allocations
   //
   int i,j,workSize;
   //
   // mb       = mb_in;
   // sb       = sb_in;
   stepType = scheme;
   //
   // allocation of all the solution arrays
   //
   q     = new double [NQ*nCell];
   qt    = new double [NQ*nCell];
   qtt   = new double [NQ*nCell];
   r     = new double [NQ*nCell];
   dq    = new double [NQ*nCell];
   ddq   = new double [NQ*nCell];
   ddqb  = new double [NQ*nCell];
   ddqf  = new double [NQ*nCell];
   dtac  = new double [nCell];   
   r0    = new double [NQ*nCell];
   sigma = new double [nCell];
   D     = new double** [nCell];
   ff    = new FACEMAT [nFace];
   // itag = new int [ncell];

   workSize  =  nMaxChain+5;
   // cindx = new int [workSize];
   // ctype = new int [workSize];
   ql    = new double* [workSize]; 
   qr    = new double* [workSize]; 
   dql   = new double* [workSize]; 
   dqr   = new double* [workSize]; 
   flux  = new double* [workSize]; 
   fv    = new double* [workSize]; 
   df    = new double* [workSize]; 
   flux2 = new double* [workSize]; 
   Q     = new double* [workSize]; 
   F     = new double* [workSize]; 
   A     = new double** [workSize];
   B     = new double** [workSize];
   C     = new double** [workSize];
   //
   for(i=0;i<workSize;i++)
   {
      ql[i]   = new double [NQ];
      qr[i]   = new double [NQ];
      dql[i]  = new double [NQ];
      dqr[i]  = new double [NQ];
      flux[i] = new double [NQ];
      fv[i]   = new double [NQ];
      df[i]   = new double [NQ];
      flux2[i] = new double [NQ];
      F[i]    = new double [NQ];
      Q[i]    = new double [NQ];
      A[i]    = new double* [NQ];
      B[i]    = new double* [NQ];
      C[i]    = new double* [NQ];
      for(j=0;j<NQ;j++)
      {
         A[i][j] = new double[NQ];
         B[i][j] = new double[NQ];
         C[i][j] = new double[NQ];
      }
   }

}

// ##################################################################
// END OF FILE
// ##################################################################
