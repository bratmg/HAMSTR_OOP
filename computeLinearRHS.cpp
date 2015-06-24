// ##################################################################
//
// computeLinearRHS.cpp
//
// ##################################################################
#include "solver.h"
#include "reconstruction.h"
// ##################################################################

using namespace CODEVARS;

using RECONSTRUCTION::MUSCL_DERIV;
using RECONSTRUCTION::WENO_DERIV;
// ##################################################################
void SOLVER::computeLinearRHS(void)
{
   //
   int    i,j,k,m,f,n;
   int    f1,f2;
   int    is,ie;
   int    iface;
   int    iPeriodic;
   int    chainSize,ntotal;
   int    node1,node2,node3,node4,leftCell,rightCell,icell,iflag,nbase;
   int    iface1,iface2,rightCell1,rightCell2,c1,c2,c3;
   double dbletemp;
   double leftState[NQ];
   double rightState[NQ];
   double dleftState[NQ];
   double drightState[NQ];
   double consVar[NQ];
   double dqvar[NQ];
   double fluxLocal[NQ];
   double specRadius;
   double faceVel=0.;
   double dsnorm;
   double dtfac;
   double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;  
   double xa,ya,za,xb,yb,zb;
   double pp;
   double th,qt,eps;
   double dscheck[2];  
   //
   // add diagonal term to the linear residual
   //
   ntotal = nCell*NQ;
   for(i = 0; i < ntotal; i++) r[i] = (r0[i] - dq[i]);

// ==================================================================
// one loop per chain to evaluate fluxes
// on all the faces in the chain
// ==================================================================

   for(i = 0; i < nChain; i++)
   {
      iflag = 0;
      f1    = faceStartPerChain[i  ];
      f2    = faceStartPerChain[i+1];
      m     = nGhost;
      //
      for(f = f1; f < f2; f++)
      {
         iface       = chainConn[f];
         cindx[m] = faces[(NFACE+2)*iface+FNODE];
         m++;
      }

// ==================================================================
// add buffer cells to the chain and collect the chain indices in
// a contigous array to help with cache
// ==================================================================

      if (chainConn[f1] == chainConn[f2-1])
      {
         //
         // this is a closed chain
         // make it periodic
         //
         iflag         =  0;
         f             =  f1+1;
         iface         =  chainConn[f];
         cindx[m]  =  faces[(NFACE+2)*iface+FNODE];
         m++;
         chainSize     =  m;
         m             =  0;
         for(f = f2-nGhost-1; f < f2-1; f++)
         {
            iface       = chainConn[f];
            cindx[m] = faces[(NFACE+2)*iface+FNODE];
            m++;
         } // f loop
      }
      else
      {
         //
         // this is a open chain
         // -ve index indicates necessity to create
         // ghost cells
         //
         iflag = 1;
         if(test!=1)
         {
            if(order==5) 
            {
               m--;
               cindx[m] = -cindx[m];
               m++;
               cindx[m] = -cindx[m-3];
               m++;
               cindx[m] = -cindx[m-5];

               chainSize = m+1;
               m = 0;
               cindx[m] = -cindx[m+5];
               m = 1;
               cindx[m] = -cindx[m+3];
               m = 2;
               cindx[m] = -cindx[m+1];
            }
            else
            {
               m--;
               cindx[m]=-cindx[m];
               m++;
               cindx[m]=-cindx[m-3];
               chainSize=m+1;
               m=0;
               cindx[m]=-cindx[m+3];
               m=1;
               cindx[m]=-cindx[m+1];
            }


#ifdef useMPI /* if using MPI */ 

#endif

            // periodic bc at only for strand grid 
            if(NDIM==3)
            {
               cout << " Not yet working for strands.\n";
               exit(1);
            }
         }
      } // open loops

// ==================================================================
// Loop through each chain and collect the fluxes
// ==================================================================
      for(j = 0; j <chainSize; j++)
      {
         icell = cindx[j];
         if (icell >=0 ||(icell==0&&j==nGhost)||(icell==0&&j==chainSize-nGhost-1)) 
         {
            m = NQ*icell;
            for(k = 0;k < NQ; k++)
            {
               consVar[k] = q[m];
               dqvar[k]   = dq[m];
               m++;

               flux[j][k]  = consVar[k];
               df[j][k] = dqvar[k];
                  
            }

         }

         if(icell<0||(icell==0&&j==nGhost-1)||(icell==0&&j==chainSize-nGhost))// icell < 0
         {
            //
            // operate on ghost cells
            // based on whether they are on the solid boundary or not
            //
            if (j < nGhost) 
               iface=chainConn[f1  ];
            else 
               iface=chainConn[f2-1];

            rightCell = faces[(NFACE+2)*iface+(FNODE+2)];

            if (rightCell == -2)  /* this is a face on solid wall */
            {
               icell = -icell;
               m     = NQ*icell;
               for(k = 0;k < NQ; k++)
               {
                  consVar[k] = q[m];
                  dqvar[k]   = dq[m];
                  m++;
               }

#ifdef Dim3 /* three-dimensional space */         
               flux[j][0]   =   consVar[0];
               flux[j][1]   =  (consVar[1]*refMtx[iface][0][0]
                             +   consVar[2]*refMtx[iface][0][1]
                             +   consVar[3]*refMtx[iface][0][2]);
               flux[j][2]   =  (consVar[1]*refMtx[iface][1][0]
                             +   consVar[2]*refMtx[iface][1][1]
                             +   consVar[3]*refMtx[iface][1][2]);
               flux[j][3]   =  (consVar[1]*refMtx[iface][2][0]
                             +   consVar[2]*refMtx[iface][2][1]
                             +   consVar[3]*refMtx[iface][2][2]);
               flux[j][4]   =   consVar[4];       

               df[j][0]  =   dqvar[0];
               df[j][1]  =  (dqvar[1]*refMtx[iface][0][0]
                             +   dqvar[2]*refMtx[iface][0][1]
                             +   dqvar[3]*refMtx[iface][0][2]);
               df[j][2]  =  (dqvar[1]*refMtx[iface][1][0]
                             +   dqvar[2]*refMtx[iface][1][1]
                             +   dqvar[3]*refMtx[iface][1][2]);
               df[j][3]  =  (dqvar[1]*refMtx[iface][2][0]
                             +   dqvar[2]*refMtx[iface][2][1]
                             +   dqvar[3]*refMtx[iface][2][2]);

               df[j][4]  =   dqvar[4];         
#else /* two-dimensional space */
               flux[j][0]   =   consVar[0];
               flux[j][1]   =  (consVar[1]*refMtx[iface][0][0]
                             +   consVar[2]*refMtx[iface][0][1]);
               flux[j][2]   =  (consVar[1]*refMtx[iface][1][0]
                             +   consVar[2]*refMtx[iface][1][1]);
               flux[j][3]   =   consVar[3];         

               df[j][0]  =   dqvar[0];
               df[j][1]  =  (dqvar[1]*refMtx[iface][0][0]
                             +   dqvar[2]*refMtx[iface][0][1]);
               df[j][2]  =  (dqvar[1]*refMtx[iface][1][0]
                             +   dqvar[2]*refMtx[iface][1][1]);
               df[j][3]  =   dqvar[3]; 

#endif              
            }
//=========================================================
// add for parallization
//=========================================================
            else if(rightCell == -5)
            {
#ifdef useMPI

#endif           
            } 
            else //this is for far field bc. 
            {
               if(test==1) 
               {
                  printf("Periodic bc has a problem!\n");
                  exit(1);
               } 
#ifdef Dim3 /* three-dimensional space */
               flux[j][0]  =  rinf;
               flux[j][1]  =  rinf*uinf;
               flux[j][2]  =  rinf*vinf;
               flux[j][3]  =  rinf*winf;
               flux[j][4]  =  einf;

               df[j][0] = ZERO;
               df[j][1] = ZERO;
               df[j][2] = ZERO;
               df[j][3] = ZERO;
               df[j][4] = ZERO;
#else /* two-dimensional space */
               flux[j][0]  =  rinf;
               flux[j][1]  =  rinf*uinf;
               flux[j][2]  =  rinf*vinf;
               flux[j][3]  =  einf;

               df[j][0] = ZERO;
               df[j][1] = ZERO;
               df[j][2] = ZERO;
               df[j][3] = ZERO;
#endif               
            }
      
         } // if icell conditions

      } // j loop

// ==================================================================
// Find reconstruction derivatives along each loop
// ==================================================================
      is  = nGhost-1;
      ie  = chainSize-1;
      th  = THIRD;
      qt  = 0.25;
      if (order==1) qt=0.0;
      eps = 1e-10;
      
      if(order==1 || order==3)
         MUSCL_DERIV(flux,dql,dqr,flux2,df,is,ie,th,qt,eps,chainSize,NQ);
    
      if(order==5) 
         WENO_DERIV(flux,dql,dqr,df,is,ie,eps,chainSize,NQ); //5th weno

      n         = is;
      iPeriodic = (chainConn[f1]==chainConn[f2-1]);

      for(f = f1; f < f2-iPeriodic; f++)
      {
         iface     = chainConn[f];
         leftCell  = faces[(NFACE+2)*iface + FNODE    ];
         rightCell = faces[(NFACE+2)*iface + FNODE + 2];

         for(m = 0; m < NQ; m++)
         {
            if (f == f2-iPeriodic-1 && iPeriodic==0) 
            {
               dleftState[m]  = dql[n  ][m];
               drightState[m] = dqr[n+1][m];
            }
            else
            {
               dleftState[m]  = dqr[n+1][m];
               drightState[m] = dql[n  ][m];
            }
         }
         //
         for(j = 0; j < NQ; j++)
         {
            fluxLocal[j] = 0; 
            for(k = 0;k < NQ; k++)
            {
               fluxLocal[j] += (((ff[iface]).lmat[j][k]*dleftState[k])+
                           ((ff[iface]).rmat[j][k]*drightState[k]));
            }
         }
         //
         m     = NQ*leftCell;
         dtfac = CFL/sigma[leftCell];
         for(j = 0; j < NQ; j++)
         {
            r[m] -= (fluxLocal[j]*dtfac);

            m++;
         }
 
         if (rightCell > -1) 
         {
            m=NQ*rightCell;
            dtfac = CFL/sigma[rightCell];
            for(j = 0; j < NQ; j++)
            {
               r[m] += (fluxLocal[j]*dtfac);
               m++;
            }
         } 

         n++;
      } // f loop
   } // i loop (chains)
// ==================================================================
// Compute the residual norms
// ==================================================================

   LInfNorm = L2Norm = ZERO;

   for (i = 0; i < nCell; i++)
   {
      dbletemp = fabs(r[NQ*i]);

      
      if (LInfNorm < dbletemp)
      {
         LInfNorm = dbletemp;
      }

      L2Norm  += r[NQ*i]*r[NQ*i];
   } // i loop
   L2Norm = sqrt(L2Norm)/nCell;


}
// ##################################################################
// END OF FILE
// ##################################################################