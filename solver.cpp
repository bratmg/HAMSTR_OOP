// ##################################################################
//
// solver.cpp
//
// ##################################################################
#include "solver.h"
// ##################################################################
// ##################################################################
//
// constuctor
//
// ##################################################################
SOLVER::SOLVER(MESHBLOCK *mb_in, SOLNBLOCK *sb_in)
{
   //
   // variable allocations
   //
   int i,j,workSize;
   //
   mb       = mb_in;
   sb       = sb_in;
   stepType = sb->scheme;
   //
   // allocation of all the solution arrays
   //
   sb->q     = new double [NQ*mb->nCell];
   sb->qt    = new double [NQ*mb->nCell];
   sb->qtt   = new double [NQ*mb->nCell];
   sb->r     = new double [NQ*mb->nCell];
   sb->dq    = new double [NQ*mb->nCell];
   sb->ddq   = new double [NQ*mb->nCell];
   sb->ddqb  = new double [NQ*mb->nCell];
   sb->ddqf  = new double [NQ*mb->nCell];
   sb->dtac  = new double [mb->nCell];   
   sb->r0    = new double [NQ*mb->nCell];
   sb->sigma = new double [mb->nCell];
   sb->D     = new double** [mb->nCell];
   sb->ff    = new FACEMAT [mb->nFace];
   // sb->itag = new int [mb->ncell];

   workSize  =  mb->nMaxChain+5;
   // sb->cindx = new int [workSize];
   // sb->ctype = new int [workSize];
   sb->ql    = new double* [workSize]; 
   sb->qr    = new double* [workSize]; 
   sb->dql   = new double* [workSize]; 
   sb->dqr   = new double* [workSize]; 
   sb->f     = new double* [workSize]; 
   sb->fv    = new double* [workSize]; 
   sb->df    = new double* [workSize]; 
   sb->f2    = new double* [workSize]; 
   sb->Q     = new double* [workSize]; 
   sb->F     = new double* [workSize]; 
   sb->A     = new double** [workSize];
   sb->B     = new double** [workSize];
   sb->C     = new double** [workSize];
   //
   for(i=0;i<workSize;i++)
   {
      sb->ql[i]  = new double [NQ];
      sb->qr[i]  = new double [NQ];
      sb->dql[i] = new double [NQ];
      sb->dqr[i] = new double [NQ];
      sb->f[i]   = new double [NQ];
      sb->fv[i]  = new double [NQ];
      sb->df[i]  = new double [NQ];
      sb->f2[i]  = new double [NQ];
      sb->F[i]   = new double [NQ];
      sb->Q[i]   = new double [NQ];
      sb->A[i]   = new double* [NQ];
      sb->B[i]   = new double* [NQ];
      sb->C[i]   = new double* [NQ];
      for(j=0;j<NQ;j++)
      {
         sb->A[i][j] = new double[NQ];
         sb->B[i][j] = new double[NQ];
         sb->C[i][j] = new double[NQ];
      }
   }

}

// ##################################################################
// END OF FILE
// ##################################################################
