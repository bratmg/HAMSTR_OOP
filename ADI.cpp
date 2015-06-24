// ##################################################################
//
// ADI.cpp
//
// ##################################################################
#include "solver.h"
#include "triSolvers.h"

using namespace CODEVARS;

using TRISOLVER::blockTridag;
using TRISOLVER::blockTridagPeriodic;
using TRISOLVER::blockTridag4;
using TRISOLVER::blockTridagPeriodic4;
// ##################################################################

void SOLVER::ADI(void)
{
   //
   // local variable declaration
   //
   int    i,j,k,l,m,f,n,worksize;
   int    mm1,f1,f2,iface,iPeriodic,chainSize;
   int    node1,node2,node3,node4,leftCell,rightCell,icell;
   int    isweep,ntotal;
   double ds[NQ];
   double lmat[NQ][NQ];
   double rmat[NQ][NQ];
   double dsnorm;
   double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
   double xa,ya,za,xb,yb,zb;
   double r01,r02,r03,a1,a2,a3,b1,b2,b3,c1,c2,c3,pp,dtfac;
   double linearl2rho,linearlinfrho;
// ==================================================================
// Initialization
// ==================================================================

   //
   // one loop per chain to evaluate fluxes
   // on all the faces in the chain
   //
   for(i = 0; i < nCell; i++)
   {
      dtfac = CFL/sigma[i];
      for(m = 0; m < NQ; m++)
      {
         r[NQ*i+m]  *= dtfac;
         r0[NQ*i+m]  = r[NQ*i+m];
         dq[NQ*i+m]  = ZERO;
      }
   }

// ==================================================================
// Perform multiple sweeps and loop through all the chains
// ==================================================================  
   worksize = nMaxChain+5;
   for(isweep = 0; isweep < msweep; isweep++)
   {
      for(i = 0; i < nChain; i++)
      {
         //
         // zero out block matrices
         //
         for(m = 0; m < worksize; m++)
            for(k = 0; k < NQ; k++)
            {
               F[m][k]=0;
               for(j = 0; j < NQ; j++)
                  A[m][k][j] = B[m][k][j] = C[m][k][j] = ZERO;
            }
    
         //
         // collect cells on the loop
         // 
         f1        = faceStartPerChain[i];
         f2        = faceStartPerChain[i+1];
         iPeriodic = (chainConn[f1]==chainConn[f2-1]); // test for periodicity 
        
         m=0;
         chainSize = (f2-iPeriodic-f1);

// ==================================================================  
// Create matrices A,B,C and vector F for inversion for all
// faces in the chain
// ================================================================== 

         for(f = f1; f < f2-iPeriodic; f++)
         {
            iface     = chainConn[f];
            leftCell  = faces[(NFACE+2)*iface + FNODE    ];
            rightCell = faces[(NFACE+2)*iface + FNODE + 2];
            // 
            // construct left and right matrices for this face
            //
            for(j = 0; j < NQ; j++)
               for(k = 0; k < NQ; k++)
               {
                  lmat[j][k] = (ff[iface]).lmat[j][k];
                  rmat[j][k] = (ff[iface]).rmat[j][k];
               }

            //
            // closed loop
            //
            if (iPeriodic==1) 
            {
               mm1 = (m==0) ? chainSize-1 : m-1;//mm1 = mm - 1
                                       // if mm=0, mm1 = nmax (periodic)
               for(j = 0; j < NQ; j++)
               {
                  for(k = 0; k < NQ; k++)
                  {
                     B[m  ][j][k] += (lmat[j][k]);
                     B[mm1][j][k] -= (rmat[j][k]);
                     A[m  ][j][k] += (rmat[j][k]);
                     C[mm1][j][k] -= (lmat[j][k]);
                  }
                  F[m][j] = r[NQ*leftCell+j];
               }
            }
            //
            // open loop
            //
            else
            {
               mm1 = (m-1);
               if (rightCell==-2)
               {

#ifdef Dim3 /* three-dimensional space */
                  a1 = refMtx[iface][0][0]; 
                  b1 = refMtx[iface][0][1];
                  c1 = refMtx[iface][0][2];
                  a2 = refMtx[iface][1][0];
                  b2 = refMtx[iface][1][1];
                  c2 = refMtx[iface][1][2];
                  a3 = refMtx[iface][2][0];
                  b3 = refMtx[iface][2][1];
                  c3 = refMtx[iface][2][2];

               
                  for (j = 0; j < NQ; j++)
                  {
                     r01 = rmat[j][1];
                     r02 = rmat[j][2];
                     r03 = rmat[j][3];

                     rmat[j][1]=a1*r01+b1*r02+c1*r03;
                     rmat[j][2]=a2*r01+b2*r02+c2*r03;
                     rmat[j][3]=a3*r01+b3*r02+c3*r03;
                  }
#else /* two-dimensional space */
                  a1 = refMtx[iface][0][0]; 
                  b1 = refMtx[iface][0][1];
                  a2 = refMtx[iface][1][0];
                  b2 = refMtx[iface][1][1];
               
                  for (j = 0; j < NQ; j++)
                  {
                     r01 = rmat[j][1];
                     r02 = rmat[j][2];

                     rmat[j][1]=a1*r01+b1*r02;
                     rmat[j][2]=a2*r01+b2*r02;
                  }

#endif               
               }
           
               if (rightCell < 0 && iPeriodic==0 && f==f2-1) m--;

               //
               for(j = 0; j < NQ; j++)
               {
                  for(k = 0; k < NQ; k++)
                  {
                     B[m][j][k] += (lmat[j][k]);
                     if (mm1 > -1 && rightCell > -1)
                     {
                       A[m][j][k]   += (rmat[j][k]);
                       B[mm1][j][k] -= (rmat[j][k]);
                       C[mm1][j][k] -= (lmat[j][k]);
                     }
                     else
                     {
                       if (rightCell==-2) B[m][j][k] += (rmat[j][k]);
                     }
                  }
                  F[m][j]=r[NQ*leftCell+j];
               }
            } // iPeriodic     
            m++;
         } //f = f1 ~ f2-1 

// ==================================================================
// Final tweaks before inverting using a block tridiagonal solver
// ==================================================================

         m = 0;
         chainSize = chainSize-(iPeriodic==0);

         for(f = f1; f < f2-1; f++)
         {
            iface    = chainConn[f];    
            leftCell = faces[(NFACE+2)*iface+FNODE];
            dtfac    = CFL/sigma[leftCell];

            for(j = 0;j < NQ; j++)
               for(k=0; k < NQ; k++)
               {
                  A[m][j][k] *= dtfac;
                  B[m][j][k] *= dtfac;
                  C[m][j][k] *= dtfac;
               }

            // add the Identity matrix
            for(j = 0; j < NQ; j++) B[m][j][j] += 1.0;
            m++;
         }

         //
         // invert using appropriate banded block solver
         //
#ifdef Dim3  /* three-dimensional space */
         if (iPeriodic==1) blockTridagPeriodic(A,B,C,F,chainSize,NQ);
         if (iPeriodic==0) blockTridag(A,B,C,F,chainSize,NQ);
#else  /* two-dimensional space */
         if (iPeriodic==1) blockTridagPeriodic4(A,B,C,F,chainSize,NQ);
         if (iPeriodic==0) blockTridag4(A,B,C,F,chainSize,NQ);
#endif
         //
         // reassign values back at the unknown locations
         //
         m = 0;
         for(f = f1; f < f2-1; f++)
         {
            iface    = chainConn[f];
            leftCell = faces[(NFACE+2)*iface+FNODE];

            for(j = 0; j < NQ; j++)
            {
               r[NQ*leftCell+j] = F[m][j];
            }

            m++;
         }      
      } // nchains loop

      ntotal = nCell*NQ;   
      for(i = 0; i < ntotal; i++) dq[i] += r[i];

      computeLinearRHS();

   } // isweep loop

// ==================================================================
// update q
// ==================================================================

   m = 0;
   for(i = 0; i < nCell; i++)
      for(j = 0;j < NQ; j++)
      {
         q[m] += dq[m];
         m++;
      }

}
// ##################################################################
// END OF FILE
// ##################################################################
