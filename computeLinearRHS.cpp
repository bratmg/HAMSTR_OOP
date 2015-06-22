// ##################################################################
//
// computeLinearRHS.cpp
//
// ##################################################################
#include "solver.h"
#include "reconstruction.h"
// ##################################################################

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
   double flux[NQ];
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
   ntotal = mb->nCell*NQ;
   for(i = 0; i < ntotal; i++) sb->r[i] = (sb->r0[i] - sb->dq[i]);

// ==================================================================
// one loop per chain to evaluate fluxes
// on all the faces in the chain
// ==================================================================

   for(i = 0; i < mb->nChain; i++)
   {
      iflag = 0;
      f1    = mb->faceStartPerChain[i  ];
      f2    = mb->faceStartPerChain[i+1];
      m     = sb->nGhost;
      //
      for(f = f1; f < f2; f++)
      {
         iface       = mb->chainConn[f];
         mb->cindx[m] = mb->faces[(NFACE+2)*iface+FNODE];
         m++;
      }

// ==================================================================
// add buffer cells to the chain and collect the chain indices in
// a contigous array to help with cache
// ==================================================================

      if (mb->chainConn[f1] == mb->chainConn[f2-1])
      {
         //
         // this is a closed chain
         // make it periodic
         //
         iflag         =  0;
         f             =  f1+1;
         iface         =  mb->chainConn[f];
         mb->cindx[m]  =  mb->faces[(NFACE+2)*iface+FNODE];
         m++;
         chainSize     =  m;
         m             =  0;
         for(f = f2-sb->nGhost-1; f < f2-1; f++)
         {
            iface       = mb->chainConn[f];
            mb->cindx[m] = mb->faces[(NFACE+2)*iface+FNODE];
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
         if(sb->test!=1)
         {
            if(sb->order==5) 
            {
               m--;
               mb->cindx[m] = -mb->cindx[m];
               m++;
               mb->cindx[m] = -mb->cindx[m-3];
               m++;
               mb->cindx[m] = -mb->cindx[m-5];

               chainSize = m+1;
               m = 0;
               mb->cindx[m] = -mb->cindx[m+5];
               m = 1;
               mb->cindx[m] = -mb->cindx[m+3];
               m = 2;
               mb->cindx[m] = -mb->cindx[m+1];
            }
            else
            {
               m--;
               mb->cindx[m]=-mb->cindx[m];
               m++;
               mb->cindx[m]=-mb->cindx[m-3];
               chainSize=m+1;
               m=0;
               mb->cindx[m]=-mb->cindx[m+3];
               m=1;
               mb->cindx[m]=-mb->cindx[m+1];
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
         icell = mb->cindx[j];
         if (icell >=0 ||(icell==0&&j==sb->nGhost)||(icell==0&&j==chainSize-sb->nGhost-1)) 
         {
            m = NQ*icell;
            for(k = 0;k < NQ; k++)
            {
               consVar[k] = sb->q[m];
               dqvar[k]   = sb->dq[m];
               m++;

               sb->f[j][k]  = consVar[k];
               sb->df[j][k] = dqvar[k];
                  
            }

         }

         if(icell<0||(icell==0&&j==sb->nGhost-1)||(icell==0&&j==chainSize-sb->nGhost))// icell < 0
         {
            //
            // operate on ghost cells
            // based on whether they are on the solid boundary or not
            //
            if (j < sb->nGhost) 
               iface=mb->chainConn[f1  ];
            else 
               iface=mb->chainConn[f2-1];

            rightCell = mb->faces[(NFACE+2)*iface+(FNODE+2)];

            if (rightCell == -2)  /* this is a face on solid wall */
            {
               icell = -icell;
               m     = NQ*icell;
               for(k = 0;k < NQ; k++)
               {
                  consVar[k] = sb->q[m];
                  dqvar[k]   = sb->dq[m];
                  m++;
               }

#ifdef Dim3 /* three-dimensional space */         
               sb->f[j][0]   =   consVar[0];
               sb->f[j][1]   =  (consVar[1]*mb->refMtx[iface][0][0]
                             +   consVar[2]*mb->refMtx[iface][0][1]
                             +   consVar[3]*mb->refMtx[iface][0][2]);
               sb->f[j][2]   =  (consVar[1]*mb->refMtx[iface][1][0]
                             +   consVar[2]*mb->refMtx[iface][1][1]
                             +   consVar[3]*mb->refMtx[iface][1][2]);
               sb->f[j][3]   =  (consVar[1]*mb->refMtx[iface][2][0]
                             +   consVar[2]*mb->refMtx[iface][2][1]
                             +   consVar[3]*mb->refMtx[iface][2][2]);
               sb->f[j][4]   =   consVar[4];       

               sb->df[j][0]  =   dqvar[0];
               sb->df[j][1]  =  (dqvar[1]*mb->refMtx[iface][0][0]
                             +   dqvar[2]*mb->refMtx[iface][0][1]
                             +   dqvar[3]*mb->refMtx[iface][0][2]);
               sb->df[j][2]  =  (dqvar[1]*mb->refMtx[iface][1][0]
                             +   dqvar[2]*mb->refMtx[iface][1][1]
                             +   dqvar[3]*mb->refMtx[iface][1][2]);
               sb->df[j][3]  =  (dqvar[1]*mb->refMtx[iface][2][0]
                             +   dqvar[2]*mb->refMtx[iface][2][1]
                             +   dqvar[3]*mb->refMtx[iface][2][2]);

               sb->df[j][4]  =   dqvar[4];         
#else /* two-dimensional space */
               sb->f[j][0]   =   consVar[0];
               sb->f[j][1]   =  (consVar[1]*mb->refMtx[iface][0][0]
                             +   consVar[2]*mb->refMtx[iface][0][1]);
               sb->f[j][2]   =  (consVar[1]*mb->refMtx[iface][1][0]
                             +   consVar[2]*mb->refMtx[iface][1][1]);
               sb->f[j][3]   =   consVar[3];         

               sb->df[j][0]  =   dqvar[0];
               sb->df[j][1]  =  (dqvar[1]*mb->refMtx[iface][0][0]
                             +   dqvar[2]*mb->refMtx[iface][0][1]);
               sb->df[j][2]  =  (dqvar[1]*mb->refMtx[iface][1][0]
                             +   dqvar[2]*mb->refMtx[iface][1][1]);
               sb->df[j][3]  =   dqvar[3]; 

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
               if(sb->test==1) 
               {
                  printf("Periodic bc has a problem!\n");
                  exit(1);
               } 
#ifdef Dim3 /* three-dimensional space */
               sb->f[j][0]  =  sb->rinf;
               sb->f[j][1]  =  sb->rinf*sb->uinf;
               sb->f[j][2]  =  sb->rinf*sb->vinf;
               sb->f[j][3]  =  sb->rinf*sb->winf;
               sb->f[j][4]  =  sb->einf;

               sb->df[j][0] = ZERO;
               sb->df[j][1] = ZERO;
               sb->df[j][2] = ZERO;
               sb->df[j][3] = ZERO;
               sb->df[j][4] = ZERO;
#else /* two-dimensional space */
               sb->f[j][0]  =  sb->rinf;
               sb->f[j][1]  =  sb->rinf*sb->uinf;
               sb->f[j][2]  =  sb->rinf*sb->vinf;
               sb->f[j][3]  =  sb->einf;

               sb->df[j][0] = ZERO;
               sb->df[j][1] = ZERO;
               sb->df[j][2] = ZERO;
               sb->df[j][3] = ZERO;
#endif               
            }
      
         } // if icell conditions

      } // j loop

// ==================================================================
// Find reconstruction derivatives along each loop
// ==================================================================
      is  = sb->nGhost-1;
      ie  = chainSize-1;
      th  = THIRD;
      qt  = 0.25;
      if (sb->order==1) qt=0.0;
      eps = 1e-10;
      
      if(sb->order==1 || sb->order==3)
         MUSCL_DERIV(sb->f,sb->dql,sb->dqr,sb->f2,sb->df,is,ie,th,qt,eps,chainSize,NQ);
    
      if(sb->order==5) 
         WENO_DERIV(sb->f,sb->dql,sb->dqr,sb->df,is,ie,eps,chainSize,NQ); //5th weno

      n         = is;
      iPeriodic = (mb->chainConn[f1]==mb->chainConn[f2-1]);

      for(f = f1; f < f2-iPeriodic; f++)
      {
         iface     = mb->chainConn[f];
         leftCell  = mb->faces[(NFACE+2)*iface + FNODE    ];
         rightCell = mb->faces[(NFACE+2)*iface + FNODE + 2];

         for(m = 0; m < NQ; m++)
         {
            if (f == f2-iPeriodic-1 && iPeriodic==0) 
            {
               dleftState[m]  = sb->dql[n  ][m];
               drightState[m] = sb->dqr[n+1][m];
            }
            else
            {
               dleftState[m]  = sb->dqr[n+1][m];
               drightState[m] = sb->dql[n  ][m];
            }
         }
         //
         for(j = 0; j < NQ; j++)
         {
            flux[j] = 0; 
            for(k = 0;k < NQ; k++)
            {
               flux[j] += (((sb->ff[iface]).lmat[j][k]*dleftState[k])+
                           ((sb->ff[iface]).rmat[j][k]*drightState[k]));
            }
         }
         //
         m     = NQ*leftCell;
         dtfac = sb->CFL/sb->sigma[leftCell];
         for(j = 0; j < NQ; j++)
         {
            sb->r[m] -= (flux[j]*dtfac);

            m++;
         }
 
         if (rightCell > -1) 
         {
            m=NQ*rightCell;
            dtfac = sb->CFL/sb->sigma[rightCell];
            for(j = 0; j < NQ; j++)
            {
               sb->r[m] += (flux[j]*dtfac);
               m++;
            }
         } 

         n++;
      } // f loop
   } // i loop (chains)
// ==================================================================
// Compute the residual norms
// ==================================================================

   sb->LInfNorm = sb->L2Norm = ZERO;

   for (i = 0; i < mb->nCell; i++)
   {
      dbletemp = fabs(sb->r[NQ*i]);

      
      if (sb->LInfNorm < dbletemp)
      {
         sb->LInfNorm = dbletemp;
      }

      sb->L2Norm  += sb->r[NQ*i]*sb->r[NQ*i];
   } // i loop
   sb->L2Norm = sqrt(sb->L2Norm)/mb->nCell;


}
// ##################################################################
// END OF FILE
// ##################################################################