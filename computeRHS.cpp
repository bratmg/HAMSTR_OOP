// ##################################################################
//
// computeRHS.cpp
//
// ##################################################################
#include "solver.h"
#include "reconstruction.h"
// ##################################################################

using RECONSTRUCTION::MUSCL;
using RECONSTRUCTION::WENO;
// ##################################################################
//
// computeRHS
//
// ##################################################################
void SOLVER::computeRHS(char * option)
{
   // 
   // local variable declarations
   //
   int    i,j,k,f,m,n,temp,iflag;
   int    f1,f2,is,ie,iPeriodic;
   int    iface,chainSize,leftCell,rightCell,icell;
   int    imode = 1;
   double th,qt,eps,ds[NDIM],specRadius,dbletemp;
   double leftState0[NQ],rightState0[NQ];
   double leftState[NQ],rightState[NQ];
   double consVar[NQ],flux[NQ];
   double **lmat,**rmat;
   //
   // allocate
   //
   lmat = new double* [NQ];
   rmat = new double* [NQ];
   for (i = 0; i < NQ; i++)
   {
      lmat[i] = new double [NQ];
      rmat[i] = new double [NQ];
   }
   // 
   // zero out residual and spectral radii
   //
   temp = NQ*mb->nCell;
   for (i = 0; i < temp; i++)      sb->r[i]     = ZERO;
   for (i = 0; i < mb->nCell; i++) sb->sigma[i] = ZERO;      

// ==================================================================
// one loop per chain to evaluate fluxes
// on all the faces in the chain
// ==================================================================

   for (i = 0; i < mb->nChain; i++)
   {
      iflag = 0;

      f1 = mb->faceStartPerChain[i  ];
      f2 = mb->faceStartPerChain[i+1];
      m  = sb->nGhost;

      // loop through the faces of a given loop
      for(f=f1;f<f2;f++)
      {
         iface        = mb->chainConn[f]; //iface = face index
         mb->cindx[m] = mb->faces[(NFACE+2)*iface+FNODE]; //cindx = cell index
         m++;
      }
      

// ==================================================================
// add buffer cells to the chain and collect the chain indices in
// a contigous array to help with cache
// ==================================================================

      if (mb->chainConn[f1] == mb->chainConn[f2-1])
      {
         iflag = 0;
         //
         // this is a closed chain
         // make it periodic
         //
         f            =  f1+1;
         iface        =  mb->chainConn[f];
         mb->cindx[m] =  mb->faces[(NFACE+2)*iface+FNODE];
         m++;
         chainSize    =  m;

         m=0;
         for(f = f2-sb->nGhost-1; f < f2-1; f++)
         {
            iface        =  mb->chainConn[f];
            mb->cindx[m] =  mb->faces[(NFACE+2)*iface+FNODE];
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
         // solid bc 
         if(sb->test != 1)
         {
            if(sb->order == 5) 
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

#ifdef useMPI
      
#endif
            // periodic bc at only for strand grid 
            if(NDIM==3)
            {
               cout << " Not yet working for strands.\n";
               exit(1);
            }
        
      
         } //not for periodic bc condition         

      } // open or closed chains

// ==================================================================
// Loop through each chain and collect the fluxes
// ==================================================================
      for(j = 0; j < chainSize; j++)
      {
         icell = mb->cindx[j];
         //if (icell >=0) 
         if (icell >=0 ||(icell==0&&j==sb->nGhost)||(icell==0&&j==chainSize-sb->nGhost-1)) // not negative
         {
            m = NQ*icell;
            for(k = 0; k < NQ; k++)
            {
               consVar[k]=sb->q[m]; // conservative variables : rho, rho*u, rho*e
               m++;

               sb->f[j][k] = consVar[k];
            } // k loop
         }

         if(icell<0||(icell==0&&j==sb->nGhost-1)||(icell==0&&j==chainSize-sb->nGhost))// icell < 0
         //else   //icell < 0
         {
            //
            // loop over ghost cells
            // based on whether they are on the solid boundary or not
            //
            if (j < sb->nGhost)
               iface = mb->chainConn[f1 ];
            else 
               iface = mb->chainConn[f2-1];

            rightCell = mb->faces[(NFACE+2)*iface+(FNODE+2)];
            //
            // Surface boundary conditions 
            //
            if (rightCell == BC_WALL_INVISC)  /* this is a face on solid wall */ 
            {

               icell = -icell;

               m     = NQ*icell;
               for(k = 0; k < NQ; k++)
               {
                  consVar[k]=sb->q[m];
                  m++;
               }

#ifdef Dim3 /* three-dimensional space */

               sb->f[j][0]  =  consVar[0];
               sb->f[j][1]  =  (consVar[1]*mb->refMtx[iface][0][0]
                            +   consVar[2]*mb->refMtx[iface][0][1]
                            +   consVar[3]*mb->refMtx[iface][0][2]);

               sb->f[j][2]  =  (consVar[1]*mb->refMtx[iface][1][0]
                            +   consVar[2]*mb->refMtx[iface][1][1]
                            +   consVar[3]*mb->refMtx[iface][1][2]);

               sb->f[j][3]  =  (consVar[1]*mb->refMtx[iface][2][0]
                            +   consVar[2]*mb->refMtx[iface][2][1]
                            +   consVar[3]*mb->refMtx[iface][2][2]);

               sb->f[j][4]  =  consVar[4];            

#else /* two-dimensional space */

               sb->f[j][0]  =  consVar[0];
               sb->f[j][1]  =  (consVar[1]*mb->refMtx[iface][0][0]
                            +   consVar[2]*mb->refMtx[iface][0][1]);

               sb->f[j][2]  =  (consVar[1]*mb->refMtx[iface][1][0]
                            +   consVar[2]*mb->refMtx[iface][1][1]);

               sb->f[j][3]  =  consVar[3];

#endif            
            }
// ==================================================================
//
// paralellization
//
// ==================================================================
            else if(rightCell == -5)
            {
#ifdef useMPI
    
#endif
            }
            else // this is for far field b.c, so constant variabls 
            {
               if(sb->test==1) 
               {
                  printf("Periodic bc not yet implemented.\n");
                  exit(1);
               }
               sb->f[j][0]  =  sb->rinf;
               sb->f[j][1]  =  sb->rinf*sb->uinf;
               sb->f[j][2]  =  sb->rinf*sb->vinf;

#ifdef Dim3 /* three-dimensional space */          
               sb->f[j][3]  =  sb->rinf*sb->winf;
               sb->f[j][4]  =  sb->einf;//pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf+s->winf*s->winf);
#else /* two-dimensional space */
               sb->f[j][3]  =  sb->einf;
#endif            
            }
         } //icell

         
      } // j loop (chainsize)

// ==================================================================
// Reconstruction along the loops
// ==================================================================

      is  = sb->nGhost-1;
      ie  = chainSize-1;
      th  = THIRD;
      qt  = 0.25;
      if (sb->order==1) qt=0.0;
      eps = 1e-10;

      if(sb->order==1 || sb->order==3) 
        MUSCL(sb->f,sb->ql,sb->qr,sb->f2,is,ie,th,qt,eps,chainSize,NQ);
      if(sb->order==5) 
        WENO(sb->f,sb->ql,sb->qr,is,ie,eps,chainSize,NQ); 

      n   = is;
      iPeriodic = (mb->chainConn[f1] == mb->chainConn[f2-1]);
      //
      for(f = f1; f < f2-iPeriodic; f++)
      {
         iface     = mb->chainConn[f];
         leftCell  = mb->faces[(NFACE+2)*iface + FNODE    ];
         rightCell = mb->faces[(NFACE+2)*iface + FNODE + 2];

         for(m = 0; m < NQ; m++)
         {
            if (f == f2-iPeriodic-1 && iPeriodic==0) 
            {
               leftState[m]  = sb->ql[n  ][m];
               rightState[m] = sb->qr[n+1][m];
            }
            else
            {
               leftState[m]  = sb->qr[n+1][m];
               rightState[m] = sb->ql[n  ][m];
            }       
            leftState0[m] = sb->q[NQ*leftCell+m]; //g->ql[j][m];              
         
            if (rightCell > -1) 
               rightState0[m] = sb->q[NQ*rightCell+m]; //g->qr[j+1][m];

         } // m loop
        
         if (rightCell==-1 && sb->test!=1) 
         {
            rightState0[0] = sb->rinf;
            rightState0[1] = sb->rinf*sb->uinf;
            rightState0[2] = sb->rinf*sb->vinf;
#ifdef Dim3         
            rightState0[3] = sb->rinf*sb->winf;
            rightState0[4] = sb->einf;
#else
            rightState0[3] = sb->einf;    
#endif
         }
         else if (rightCell==-2 && sb->test!=1) 
         {

#ifdef Dim3         
            rightState0[0]  =  leftState0[0];
            rightState0[1]  =  (leftState0[1]*mb->refMtx[iface][0][0]
                            +   leftState0[2]*mb->refMtx[iface][0][1]
                            +   leftState0[3]*mb->refMtx[iface][0][2]);

            rightState0[2]  =  (leftState0[1]*mb->refMtx[iface][1][0]
                            +   leftState0[2]*mb->refMtx[iface][1][1]
                            +   leftState0[3]*mb->refMtx[iface][1][2]);

            rightState0[3]  =  (leftState0[1]*mb->refMtx[iface][2][0]
                            +   leftState0[2]*mb->refMtx[iface][2][1]
                            +   leftState0[3]*mb->refMtx[iface][2][2]);

            rightState0[4]  =  leftState0[4];             
#else
            rightState0[0]  =  leftState0[0];
            rightState0[1]  =  (leftState0[1]*mb->refMtx[iface][0][0]
                            +   leftState0[2]*mb->refMtx[iface][0][1]);

            rightState0[2]  =  (leftState0[1]*mb->refMtx[iface][1][0]
                            +   leftState0[2]*mb->refMtx[iface][1][1]);

            rightState0[3]  =  leftState0[3];

#endif
         }

// ==================================================================
// Compute Roe's flux and Jacobian at each face
// ==================================================================

#ifndef Dim3      
         ds[0] = mb->normVec[iface][0];
         ds[1] = mb->normVec[iface][1];
         fluxRoe2D(ds,leftState,rightState,flux,&specRadius,sb->gamma);


         if (strcmp(option,"implicit")==0)
         {
            // compute the Roe's flux-based Jacobian
            jacobianRoe2D(ds,leftState,rightState,lmat,rmat,sb->gamma);
            //
            // save the left and right matrix at each face
            for(j = 0; j < NQ; j++)
               for(k = 0; k < NQ; k++)
               {
                  (sb->ff[iface]).lmat[j][k]=lmat[j][k];
                  (sb->ff[iface]).rmat[j][k]=rmat[j][k];   
               }

         }
#else      
         ds[0] = mb->normVec[iface][0];
         ds[1] = mb->normVec[iface][1];
         ds[2] = mb->normVec[iface][2];
         fluxRoe3D(ds,leftState,rightState,flux,&specRadius,sb->gamma);
#endif

// ==================================================================
// Compute the residual and store the spectral radius
// ==================================================================

         // compute residual array for the cells
         m = NQ*leftCell;
         for(j = 0; j < NQ; j++)
         {
            sb->r[m] -= flux[j]; //residual
            m++;
         }

         // accumulate spectral radius for left cell
         sb->sigma[leftCell] += specRadius;

         // if NOT a boundary cell
         if (rightCell > -1) 
         {
            m = NQ*rightCell;
            for(j = 0; j < NQ; j++)
            {
               sb->r[m] += flux[j]; //residual
               m++;
            }
            // spectral radius for right cell
            sb->sigma[rightCell] += specRadius;
         }
         n++;  // n is a running index of the chain faces (is to ie-1) 
      }  // f1-f2
   } // i loop

// ==================================================================
// Compute the residual norms
// ==================================================================

   sb->LInfNorm = sb->L2Norm = ZERO;

   for (i = 0; i < mb->nCell; i++)
   {
      dbletemp = fabs(sb->r[NQ*i]);
      if (sb->LInfNorm < dbletemp)
         sb->LInfNorm = dbletemp;

      sb->L2Norm  += sb->r[NQ*i]*sb->r[NQ*i];
   } // i loop
   sb->L2Norm = sqrt(sb->L2Norm)/mb->nCell;

}

// ##################################################################
// END OF FILE
// ##################################################################
