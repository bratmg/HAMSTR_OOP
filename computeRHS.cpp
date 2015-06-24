// ##################################################################
//
// computeRHS.cpp
//
// ##################################################################
#include "solver.h"
#include "reconstruction.h"
// ##################################################################

using namespace CODEVARS;

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
   double consVar[NQ],fluxLocal[NQ];
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
   temp = NQ*nCell;
   for (i = 0; i < temp; i++)      r[i]     = ZERO;
   for (i = 0; i < nCell; i++) sigma[i] = ZERO;      

// ==================================================================
// one loop per chain to evaluate fluxes
// on all the faces in the chain
// ==================================================================

   for (i = 0; i < nChain; i++)
   {

      iflag = 0;

      f1 = faceStartPerChain[i  ];
      f2 = faceStartPerChain[i+1];
      m  = nGhost;

      // loop through the faces of a given loop
      for(f=f1;f<f2;f++)
      {
         iface        = chainConn[f]; //iface = face index
         cindx[m] = faces[(NFACE+2)*iface+FNODE]; //cindx = cell index
         m++;
      }
      

// ==================================================================
// add buffer cells to the chain and collect the chain indices in
// a contigous array to help with cache
// ==================================================================

      if (chainConn[f1] == chainConn[f2-1])
      {
         iflag = 0;
         //
         // this is a closed chain
         // make it periodic
         //
         f            =  f1+1;
         iface        =  chainConn[f];
         cindx[m] =  faces[(NFACE+2)*iface+FNODE];
         m++;
         chainSize    =  m;

         m=0;
         for(f = f2-nGhost-1; f < f2-1; f++)
         {
            iface        =  chainConn[f];
            cindx[m] =  faces[(NFACE+2)*iface+FNODE];
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
         if(test != 1)
         {
            if(order == 5) 
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
         icell = cindx[j];
         //if (icell >=0) 
         if (icell >=0 ||(icell==0&&j==nGhost)||(icell==0&&j==chainSize-nGhost-1)) // not negative
         {
            m = NQ*icell;
            for(k = 0; k < NQ; k++)
            {
               consVar[k]=q[m]; // conservative variables : rho, rho*u, rho*e
               m++;

               flux[j][k] = consVar[k];
            } // k loop
         }

         if(icell<0||(icell==0&&j==nGhost-1)||(icell==0&&j==chainSize-nGhost))// icell < 0
         //else   //icell < 0
         {
            //
            // loop over ghost cells
            // based on whether they are on the solid boundary or not
            //
            if (j < nGhost)
               iface = chainConn[f1 ];
            else 
               iface = chainConn[f2-1];

            rightCell = faces[(NFACE+2)*iface+(FNODE+2)];
            //
            // Surface boundary conditions 
            //
            if (rightCell == BC_WALL_INVISC)  /* this is a face on solid wall */ 
            {

               icell = -icell;

               m     = NQ*icell;
               for(k = 0; k < NQ; k++)
               {
                  consVar[k]=q[m];
                  m++;
               }

#ifdef Dim3 /* three-dimensional space */

               flux[j][0]  =  consVar[0];
               flux[j][1]  =  (consVar[1]*refMtx[iface][0][0]
                            +   consVar[2]*refMtx[iface][0][1]
                            +   consVar[3]*refMtx[iface][0][2]);

               flux[j][2]  =  (consVar[1]*refMtx[iface][1][0]
                            +   consVar[2]*refMtx[iface][1][1]
                            +   consVar[3]*refMtx[iface][1][2]);

               flux[j][3]  =  (consVar[1]*refMtx[iface][2][0]
                            +   consVar[2]*refMtx[iface][2][1]
                            +   consVar[3]*refMtx[iface][2][2]);

               flux[j][4]  =  consVar[4];            

#else /* two-dimensional space */

               flux[j][0]  =  consVar[0];
               flux[j][1]  =  (consVar[1]*refMtx[iface][0][0]
                            +   consVar[2]*refMtx[iface][0][1]);

               flux[j][2]  =  (consVar[1]*refMtx[iface][1][0]
                            +   consVar[2]*refMtx[iface][1][1]);

               flux[j][3]  =  consVar[3];

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
               if(test==1) 
               {
                  printf("Periodic bc not yet implemented.\n");
                  exit(1);
               }
               flux[j][0]  =  rinf;
               flux[j][1]  =  rinf*uinf;
               flux[j][2]  =  rinf*vinf;

#ifdef Dim3 /* three-dimensional space */          
               flux[j][3]  =  rinf*winf;
               flux[j][4]  =  einf;//pinf/gm1+0.5*rinf*(s->uinf*s->uinf+s->vinf*s->vinf+s->winf*s->winf);
#else /* two-dimensional space */
               flux[j][3]  =  einf;
#endif            
            }
         } //icell

         
      } // j loop (chainsize)

// ==================================================================
// Reconstruction along the loops
// ==================================================================

      is  = nGhost-1;
      ie  = chainSize-1;
      th  = THIRD;
      qt  = 0.25;
      if (order==1) qt=0.0;
      eps = 1e-10;

      if(order==1 || order==3) 
        MUSCL(flux,ql,qr,flux2,is,ie,th,qt,eps,chainSize,NQ);
      if(order==5) 
        WENO(flux,ql,qr,is,ie,eps,chainSize,NQ); 

      n   = is;
      iPeriodic = (chainConn[f1] == chainConn[f2-1]);
      //
      for(f = f1; f < f2-iPeriodic; f++)
      {
         iface     = chainConn[f];
         leftCell  = faces[(NFACE+2)*iface + FNODE    ];
         rightCell = faces[(NFACE+2)*iface + FNODE + 2];

         for(m = 0; m < NQ; m++)
         {
            if (f == f2-iPeriodic-1 && iPeriodic==0) 
            {
               leftState[m]  = ql[n  ][m];
               rightState[m] = qr[n+1][m];
            }
            else
            {
               leftState[m]  = qr[n+1][m];
               rightState[m] = ql[n  ][m];
            }       
            leftState0[m] = q[NQ*leftCell+m]; //g->ql[j][m];              
         
            if (rightCell > -1) 
               rightState0[m] = q[NQ*rightCell+m]; //g->qr[j+1][m];

         } // m loop
        
         if (rightCell==-1 && test!=1) 
         {
            rightState0[0] = rinf;
            rightState0[1] = rinf*uinf;
            rightState0[2] = rinf*vinf;
#ifdef Dim3         
            rightState0[3] = rinf*winf;
            rightState0[4] = einf;
#else
            rightState0[3] = einf;    
#endif
         }
         else if (rightCell==-2 && test!=1) 
         {

#ifdef Dim3         
            rightState0[0]  =  leftState0[0];
            rightState0[1]  =  (leftState0[1]*refMtx[iface][0][0]
                            +   leftState0[2]*refMtx[iface][0][1]
                            +   leftState0[3]*refMtx[iface][0][2]);

            rightState0[2]  =  (leftState0[1]*refMtx[iface][1][0]
                            +   leftState0[2]*refMtx[iface][1][1]
                            +   leftState0[3]*refMtx[iface][1][2]);

            rightState0[3]  =  (leftState0[1]*refMtx[iface][2][0]
                            +   leftState0[2]*refMtx[iface][2][1]
                            +   leftState0[3]*refMtx[iface][2][2]);

            rightState0[4]  =  leftState0[4];             
#else
            rightState0[0]  =  leftState0[0];
            rightState0[1]  =  (leftState0[1]*refMtx[iface][0][0]
                            +   leftState0[2]*refMtx[iface][0][1]);

            rightState0[2]  =  (leftState0[1]*refMtx[iface][1][0]
                            +   leftState0[2]*refMtx[iface][1][1]);

            rightState0[3]  =  leftState0[3];

#endif
         }

// ==================================================================
// Compute Roe's flux and Jacobian at each face
// ==================================================================

#ifndef Dim3      
         ds[0] = normVec[iface][0];
         ds[1] = normVec[iface][1];
         fluxRoe2D(ds,leftState,rightState,fluxLocal,&specRadius,Gamma);


         if (strcmp(option,"implicit")==0)
         {
            // compute the Roe's flux-based Jacobian
            jacobianRoe2D(ds,leftState,rightState,lmat,rmat,Gamma);
            //
            // save the left and right matrix at each face
            for(j = 0; j < NQ; j++)
               for(k = 0; k < NQ; k++)
               {
                  (ff[iface]).lmat[j][k]=lmat[j][k];
                  (ff[iface]).rmat[j][k]=rmat[j][k];   
               }

         }
#else      
         ds[0] = normVec[iface][0];
         ds[1] = normVec[iface][1];
         ds[2] = normVec[iface][2];
         fluxRoe3D(ds,leftState,rightState,fluxLocal,&specRadius,Gamma);
#endif

// ==================================================================
// Compute the residual and store the spectral radius
// ==================================================================

         // compute residual array for the cells
         m = NQ*leftCell;
         for(j = 0; j < NQ; j++)
         {
            r[m] -= fluxLocal[j]; //residual
            m++;
         }

         // accumulate spectral radius for left cell
         sigma[leftCell] += specRadius;

         // if NOT a boundary cell
         if (rightCell > -1) 
         {
            m = NQ*rightCell;
            for(j = 0; j < NQ; j++)
            {
               r[m] += fluxLocal[j]; //residual
               m++;
            }
            // spectral radius for right cell
            sigma[rightCell] += specRadius;
         }
         n++;  // n is a running index of the chain faces (is to ie-1) 
      }  // f1-f2
   } // i loop

// ==================================================================
// Compute the residual norms
// ==================================================================

   LInfNorm = L2Norm = ZERO;

   for (i = 0; i < nCell; i++)
   {
      dbletemp = fabs(r[NQ*i]);
      if (LInfNorm < dbletemp)
         LInfNorm = dbletemp;

      L2Norm  += r[NQ*i]*r[NQ*i];
   } // i loop
   L2Norm = sqrt(L2Norm)/nCell;

}

// ##################################################################
// END OF FILE
// ##################################################################
