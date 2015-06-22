// ###############################################################################
//
// Patch code that calls 3D jacobian routine
//
// ###############################################################################
#include "solver.h"
#define IMODE 1

void SOLVER::jacobianRoe2D(double *nxyz2, 
                           double *ql2, 
                           double *qr2, 
                           double **lmat2, 
                           double **rmat2,
                           double   gam)
{
   //
   // local variable declaration
   //
   int    i,j;
   int    idiv[4];
   double nxyz[3];
   double ql[5],qr[5];
   double area;
   double **lmat,**rmat;
   //
   // allocate
   //
   lmat = new double* [5];
   rmat = new double* [5];
   for (i = 0; i < 5; i++)
   {
      lmat[i] = new double [5];
      rmat[i] = new double [5];
   }
   //
   idiv[0]=0;
   idiv[1]=1;
   idiv[2]=2;
   idiv[3]=4;
   //
   nxyz[0]=nxyz2[0];
   nxyz[1]=nxyz2[1];
   nxyz[2]=0;
   area = sqrt(nxyz2[0]*nxyz2[0]+nxyz2[1]*nxyz2[1]);
   //
   ql[0] = ql2[0];
   ql[1] = ql2[1];
   ql[2] = ql2[2];
   ql[3] = 0;
   ql[4] = ql2[3];
   //
   qr[0] = qr2[0];
   qr[1] = qr2[1];
   qr[2] = qr2[2];
   qr[3] = 0;
   qr[4] = qr2[3];
   //
   jacobianRoe3D(nxyz,ql,qr,lmat,rmat,gam);

   // 
   // reduce and multiply by area
   // 
   for (i = 0; i < 4; i++)
      for (j = 0; j < 4; j++)
      {
         lmat2[i][j]=lmat[idiv[i]][idiv[j]]*area;
         rmat2[i][j]=rmat[idiv[i]][idiv[j]]*area;
      }

}

// ###############################################################################
// #     Computes flux Jacobian matrix for subroutine "flux_roe" that can be
// #               used in Newton solvers and adjoint codes.
// #
// #     nxyz(3)             - Three components of face normal vector.
// #                           These can be dimensional (i.e including face area)
// #                           or non-dimensional.
// #                           Returned flux does not include face area.
// #     ql(5),qr(5)         - Conserved variables (ro, ro*u, ro*v, ro*w, Et)
// #     lmat(5,5),rmat(5,5) - Left and right state flux Jacobian matrices
// #     gam                 - Ratio of specific heats
// #     imode               - 0 = approximate linearization where the eigenvalues
// #                               are treated as constants
// #                           1 = exact linearization
// #
// ###############################################################################

// ################# (C) Copyright Karthik Mani 2006 #############################


void SOLVER::jacobianRoe3D(double  *nxyz,
                           double  *ql,
                           double  *qr,
                           double **lmat,
                           double **rmat,
                           double   gam)
{
   //
   // local variable allocation
   //
   double  gm1;
   double  nxd,nyd,nzd,area,nx,ny,nz;
   double  rol,ul,vl,wl,pl,hl;
   double  ror,ur,vr,wr,pr,hr;
   double  uconl,uconr;
   double  ubar,vbar,wbar,hbar,uconbar,cbar,robar;
   double  dp,dro,du,dv,dw;
   double  eig1,eig2,eig3;

   double  fact,A,B,term1,term2,del1,del2;

   int    k,i,j;
   double  dro_dql[5],dro_dqr[5];
   double  du_dql[5],du_dqr[5];
   double  dv_dql[5],dv_dqr[5];
   double  dw_dql[5],dw_dqr[5];
   double  dp_dql[5],dp_dqr[5];
   double  ducon_dql[5],ducon_dqr[5];
   double  ddel1_dql[5],ddel1_dqr[5];
   double  ddel2_dql[5],ddel2_dqr[5];

   double  dq5_dql[5],dq5_dqr[5];
   double  dh_dql[5],dh_dqr[5];
   double  dfact_dql[5],dfact_dqr[5];
   double  dA_dql[5],dA_dqr[5];
   double  dB_dql[5],dB_dqr[5];
   double  drobar_dql[5],dubar_dql[5],dvbar_dql[5],dwbar_dql[5];
   double  drobar_dqr[5],dubar_dqr[5],dvbar_dqr[5],dwbar_dqr[5];
   double  dhbar_dql[5],duconbar_dql[5],dcbar_dql[5];
   double  dhbar_dqr[5],duconbar_dqr[5],dcbar_dqr[5];

   double  deig1_dql[5],deig2_dql[5],deig3_dql[5];
   double  deig1_dqr[5],deig2_dqr[5],deig3_dqr[5];
   double  dterm1_dql[5],dterm1_dqr[5];
   double  dterm2_dql[5],dterm2_dqr[5];
   double  cl,cr,dc_dql[5],dc_dqr[5];
   double  t1a,t1b,t2a,t2b,t3a,t3b;
   double  eps1,eps2,eps3;

   double  **lmat1,**rmat1,**imat;
   lmat1 = new double* [5];
   rmat1 = new double* [5];
   imat  = new double* [5];
   for (i = 0; i < 5; i++)
   {
      lmat1[i] = new double [5];
      rmat1[i] = new double [5];
      imat[i]  = new double [5];
   }

   //
   // first executable statement
   //
   gm1 = gam - 1.0;

   nxd = nxyz[0];
   nyd = nxyz[1];
   nzd = nxyz[2];

   area = sqrt(nxd*nxd + nyd*nyd + nzd*nzd);

   nx = nxd/area;
   ny = nyd/area;
   nz = nzd/area;
   //
   // back calculate primitive state
   //
   rol = ql[0];
   ul  = ql[1]/ql[0];
   vl  = ql[2]/ql[0];
   wl  = ql[3]/ql[0];
   pl  = gm1*( ql[4] - 0.5 * rol * (ul*ul + vl*vl + wl*wl) );
   hl  = (ql[4] + pl)/rol;
   cl  = sqrt(gam*pl/rol);
     
   ror = qr[0];
   ur  = qr[1]/qr[0];
   vr  = qr[2]/qr[0];
   wr  = qr[3]/qr[0];
   pr  = gm1*( qr[4] - 0.5 * ror * (ur*ur + vr*vr + wr*wr) );
   hr  = (qr[4] + pr)/ror;
   cr  = sqrt(gam*pr/ror);

   //
   // primitive state differences
   //
   dro = ror - rol;
   du  =  ur - ul;
   dv  =  vr - vl;
   dw  =  wr - wl;
   dp  =  pr - pl;
   //
   // face normal velocities
   //
   uconr = ur*nx + vr*ny + wr*nz;
   uconl = ul*nx + vl*ny + wl*nz;

// ==================================================================
// linearization of left and right primitive states
// ==================================================================
  
   //
   // left state
   //

   for (i = 0; i < 5; i++) dro_dql[i] = 0.0;
   dro_dql[0] = 1.0;

   for (i = 0; i < 5; i++) du_dql[i] =  0.0;
   du_dql[0] = -ul /rol;
   du_dql[1] =  1.0/rol;

   for (i = 0; i < 5; i++) dv_dql[i] =  0.0;
   dv_dql[0] = -vl /rol;
   dv_dql[2] =  1.0/rol;

   for (i = 0; i < 5; i++) dw_dql[i] =  0.0;
   dw_dql[0] = -wl /rol;
   dw_dql[3] =  1.0/rol;

   dp_dql[0] =  0.5*gm1*( ul*ul + vl*vl + wl*wl );
   dp_dql[1] = -gm1*ul;
   dp_dql[2] = -gm1*vl;
   dp_dql[3] = -gm1*wl;
   dp_dql[4] =  gm1;

   for (i = 0; i < 5; i++) dq5_dql[i] = 0.0;
   dq5_dql[4] = 1.0;
   
   for (i = 0; i < 5; i++)
   {
      dh_dql[i] = -(ql[4] + pl)*dro_dql[i]/rol/rol 
                + (1.0/rol)*(dq5_dql[i] + dp_dql[i]);
   }

   for (i = 0; i < 5; i++)                
      dc_dql[i] = (0.5*gam/cl)*( (1.0/rol)*dp_dql[i] - (pl/rol/rol)*dro_dql[i] );

   ducon_dql[0] = -uconl/rol;
   ducon_dql[1] =  nx   /rol;
   ducon_dql[2] =  ny   /rol;
   ducon_dql[3] =  nz   /rol;
   ducon_dql[4] =  0.0;


   //
   // right state
   //

   for (i = 0; i < 5; i++) dro_dqr[i] = 0.0;
   dro_dqr[0] = 1.0;

   for (i = 0; i < 5; i++) du_dqr[i] =  0.0;
   du_dqr[0] = -ur /ror;
   du_dqr[1] =  1.0/ror;

   for (i = 0; i < 5; i++) dv_dqr[i] =  0.0;
   dv_dqr[0] = -vr /ror;
   dv_dqr[2] =  1.0/ror;

   for (i = 0; i < 5; i++) dw_dqr[i] =  0.0;
   dw_dqr[0] = -wr /ror;
   dw_dqr[3] =  1.0/ror;

   dp_dqr[0] =  0.5*gm1*( ur*ur + vr*vr + wr*wr);
   dp_dqr[1] = -gm1*ur;
   dp_dqr[2] = -gm1*vr;
   dp_dqr[3] = -gm1*wr;
   dp_dqr[4] =  gm1;

   for (i = 0; i < 5; i++) dq5_dqr[i] = 0.0;
   dq5_dqr[4] = 1.0;

   for (i = 0; i < 5; i++)
   { 
      dh_dqr[i] = -(qr[4] + pr)*dro_dqr[i]/ror/ror 
                + (1.0/ror)*(dq5_dqr[i] + dp_dqr[i]);
   }             

   for (i = 0; i < 5; i++) 
      dc_dqr[i] = (0.5*gam/cr)*( (1.0/ror)*dp_dqr[i] - (pr/ror/ror)*dro_dqr[i] );

   ducon_dqr[0] = -uconr/ror;
   ducon_dqr[1] =  nx   /ror;
   ducon_dqr[2] =  ny   /ror;
   ducon_dqr[3] =  nz   /ror;
   ducon_dqr[4] =  0.0;

// ==================================================================
// Roe average state
// ==================================================================

   fact = sqrt(ror/rol);

   A    = 1.0 /(1.0 + fact);
   B    = fact/(1.0 + fact);

   robar = rol*fact;
   ubar  = ul*A + ur*B;
   vbar  = vl*A + vr*B;
   wbar  = wl*A + wr*B;
   hbar  = hl*A + hr*B;
   cbar = gm1*(hbar - 0.5*(ubar*ubar + vbar*vbar + wbar*wbar));
   cbar = sqrt(cbar);
   uconbar = ubar*nx + vbar*ny + wbar*nz;

// ==================================================================
// Eigenvalues
// ==================================================================

   eig1 = fabs(uconbar);
   eig2 = fabs(uconbar + cbar);
   eig3 = fabs(uconbar - cbar);


// ==================================================================
// Approximate linearization section
// ==================================================================

   if ( IMODE==1 )
   {
      term1 = -eig1 + 0.5*(eig2 + eig3);
      term2 = 0.5*(eig2 - eig3);
      del1  = term1*dp/cbar/cbar + term2*robar*(uconr - uconl)/cbar;
      del2  = term1*(uconr - uconl)*robar + term2*dp/cbar;

      for (i = 0; i < 5; i++)
      {
         ddel1_dql[i] = - term1*dp_dql[i]/cbar/cbar - term2*robar*ducon_dql[i]/cbar;
         ddel1_dqr[i] = + term1*dp_dqr[i]/cbar/cbar + term2*robar*ducon_dqr[i]/cbar;

         ddel2_dql[i] = - term1*ducon_dql[i]*robar - term2*dp_dql[i]/cbar;
         ddel2_dqr[i] = + term1*ducon_dqr[i]*robar + term2*dp_dqr[i]/cbar;
      }
      
   }
   else
   {

// ==================================================================
// Linearization of Roe averaged state
// ==================================================================      

      for (i = 0; i < 5; i++)
      {
         dfact_dql[i] = (0.5/fact)*(-ror/rol/rol)*dro_dql[i];
         dfact_dqr[i] = (0.5/fact)*(1.0/rol)*dro_dqr[i];

         dA_dql[i] = -dfact_dql[i]/(1.0+fact)/(1.0+fact);
         dA_dqr[i] = -dfact_dqr[i]/(1.0+fact)/(1.0+fact);

         dB_dql[i] = dfact_dql[i]/(1.0 + fact)/(1.0 + fact);
         dB_dqr[i] = dfact_dqr[i]/(1.0 + fact)/(1.0 + fact);

         drobar_dql[i] = dro_dql[i]*fact + rol*dfact_dql[i];
         drobar_dqr[i] =                   rol*dfact_dqr[i];

         dubar_dql[i] = du_dql[i]*A + ul*dA_dql[i]               + ur*dB_dql[i];
         dubar_dqr[i] =               ul*dA_dqr[i] + du_dqr[i]*B + ur*dB_dqr[i];

         dvbar_dql[i] = dv_dql[i]*A + vl*dA_dql[i]               + vr*dB_dql[i];
         dvbar_dqr[i] =               vl*dA_dqr[i] + dv_dqr[i]*B + vr*dB_dqr[i];

         dwbar_dql[i] = dw_dql[i]*A + wl*dA_dql[i]               + wr*dB_dql[i];
         dwbar_dqr[i] =               wl*dA_dqr[i] + dw_dqr[i]*B + wr*dB_dqr[i];

         dhbar_dql[i] = dh_dql[i]*A + hl*dA_dql[i]               + hr*dB_dql[i];
         dhbar_dqr[i] =               hl*dA_dqr[i] + dh_dqr[i]*B + hr*dB_dqr[i];

         dcbar_dql[i] = gm1*( dhbar_dql[i] - ubar*dubar_dql[i]     
                                           - vbar*dvbar_dql[i]     
                                           - wbar*dwbar_dql[i] );
         dcbar_dql[i] = dcbar_dql[i]*0.5/cbar;

         dcbar_dqr[i] = gm1*( dhbar_dqr[i] - ubar*dubar_dqr[i]     
                                           - vbar*dvbar_dqr[i]     
                                           - wbar*dwbar_dqr[i] );
         dcbar_dqr[i] = dcbar_dqr[i]*0.5/cbar;

         duconbar_dql[i] = dubar_dql[i]*nx + dvbar_dql[i]*ny + dwbar_dql[i]*nz;
         duconbar_dqr[i] = dubar_dqr[i]*nx + dvbar_dqr[i]*ny + dwbar_dqr[i]*nz;
      } // i loop

// ==================================================================
// Linearization of Eigenvalues
// ==================================================================

      if(uconbar>=0.0)
      {
         for (i = 0; i < 5; i++)
         {
            deig1_dql[i] = duconbar_dql[i];
            deig1_dqr[i] = duconbar_dqr[i];
         }
      }
      else //if(uconbar< 0.0) then
      {
         for (i = 0; i < 5; i++)
         {
            deig1_dql[i] = -duconbar_dql[i];
            deig1_dqr[i] = -duconbar_dqr[i];
         }
      }

      if( (uconbar + cbar) >= 0.0 )
      {
         for (i = 0; i < 5; i++)
         {
            deig2_dql[i] = ( duconbar_dql[i] + dcbar_dql[i] );
            deig2_dqr[i] = ( duconbar_dqr[i] + dcbar_dqr[i] );
         }
      }
      else //if( (uconbar + cbar) < 0.0 ) then
      {
         for (i = 0; i < 5; i++)
         {
            deig2_dql[i] = -( duconbar_dql[i] + dcbar_dql[i] );
            deig2_dqr[i] = -( duconbar_dqr[i] + dcbar_dqr[i] );
         }
      }

      if( (uconbar - cbar) >= 0.0 )
      {
         for (i = 0; i < 5; i++)
         {
            deig3_dql[i] = ( duconbar_dql[i] - dcbar_dql[i] );
            deig3_dqr[i] = ( duconbar_dqr[i] - dcbar_dqr[i] );
         }
      }                
      else //if( (uconbar - cbar) < 0.0 ) then
      {
         for (i = 0; i < 5; i++)
         {
            deig3_dql[i] = -( duconbar_dql[i] - dcbar_dql[i] );
            deig3_dqr[i] = -( duconbar_dqr[i] - dcbar_dqr[i] );
         }
      }

// ==================================================================

// ==================================================================

      term1 = -eig1 + 0.5*(eig2 + eig3);
      term2 = 0.5*(eig2 - eig3);
      del1  = term1*dp/cbar/cbar + term2*robar*(uconr - uconl)/cbar;
      del2  = term1*(uconr - uconl)*robar + term2*dp/cbar;

      for (i = 0; i < 5; i++)
      {
         dterm1_dql[i] = -deig1_dql[i] + 0.5*( deig2_dql[i] + deig3_dql[i] );
         dterm1_dqr[i] = -deig1_dqr[i] + 0.5*( deig2_dqr[i] + deig3_dqr[i] );

         dterm2_dql[i] = 0.5*( deig2_dql[i] - deig3_dql[i] );
         dterm2_dqr[i] = 0.5*( deig2_dqr[i] - deig3_dqr[i] );

         ddel1_dql[i] = dterm1_dql[i]*dp/cbar/cbar 
                      - term1*dp_dql[i]/cbar/cbar 
                      - 2.0*term1*dp*dcbar_dql[i]/cbar/cbar/cbar;

         ddel1_dql[i] = ddel1_dql[i] + dterm2_dql[i]*robar*( uconr-uconl )/cbar 
                                     + term2*drobar_dql[i]*(uconr-uconl)/cbar 
                                     - term2*robar*ducon_dql[i]/cbar 
                                     - dcbar_dql[i]*term2*robar*(uconr-uconl)/cbar/cbar;

         ddel1_dqr[i] = dterm1_dqr[i]*dp/cbar/cbar 
                      + term1*dp_dqr[i]/cbar/cbar 
                      - 2.0*term1*dp*dcbar_dqr[i]/cbar/cbar/cbar;

         ddel1_dqr[i] = ddel1_dqr[i] + dterm2_dqr[i]*robar*( uconr-uconl )/cbar 
                                     + term2*drobar_dqr[i]*(uconr-uconl)/cbar
                                     + term2*robar*ducon_dqr[i]/cbar 
                                     - dcbar_dqr[i]*term2*robar*(uconr-uconl)/cbar/cbar;

         ddel2_dql[i] = dterm1_dql[i]*(uconr-uconl)*robar 
                      - term1*ducon_dql[i]*robar 
                      + term1*(uconr-uconl)*drobar_dql[i];

         ddel2_dql[i] = ddel2_dql[i] + dterm2_dql[i]*dp/cbar 
                                     - term2*dp_dql[i]/cbar 
                                     - dcbar_dql[i]*term2*dp/cbar/cbar;

         ddel2_dqr[i] = dterm1_dqr[i]*(uconr-uconl)*robar 
                      + term1*ducon_dqr[i]*robar 
                      + term1*(uconr-uconl)*drobar_dqr[i];

         ddel2_dqr[i] = ddel2_dqr[i] + dterm2_dqr[i]*dp/cbar 
                                     + term2*dp_dqr[i]/cbar 
                                     - dcbar_dqr[i]*term2*dp/cbar/cbar;
      } // i loop
   } // imode if statement

// ==================================================================
// Roe flux Jacobian
// ==================================================================

   //
   // common linearization terms
   //
   for (i = 0; i < 5; i++)
      for (j = 0; j < 5; j++)
         lmat[i][j] = rmat[i][j] = imat[i][j] = 0.0;
  
   for (i = 0; i < 5; i++) imat[i][i] = 1.0;

   for (i = 0; i < 5; i++)
   {
      for (j = 0; j < 5; j++)
      {
         lmat[i][j] = lmat[i][j] - eig1*imat[i][j];
         rmat[i][j] = rmat[i][j] + eig1*imat[i][j];
      }
   }
   //
   for (i = 0; i < 5; i++)
   {
      lmat[0][i] = lmat[0][i] + ddel1_dql[i];
      rmat[0][i] = rmat[0][i] + ddel1_dqr[i];

      lmat[1][i] = lmat[1][i] + ddel1_dql[i]*ubar + ddel2_dql[i]*nx;
      rmat[1][i] = rmat[1][i] + ddel1_dqr[i]*ubar + ddel2_dqr[i]*nx;

      lmat[2][i] = lmat[2][i] + ddel1_dql[i]*vbar + ddel2_dql[i]*ny;
      rmat[2][i] = rmat[2][i] + ddel1_dqr[i]*vbar + ddel2_dqr[i]*ny;

      lmat[3][i] = lmat[3][i] + ddel1_dql[i]*wbar + ddel2_dql[i]*nz;
      rmat[3][i] = rmat[3][i] + ddel1_dqr[i]*wbar + ddel2_dqr[i]*nz;

      lmat[4][i] = lmat[4][i] + ddel1_dql[i]*hbar + ddel2_dql[i]*uconbar;
      rmat[4][i] = rmat[4][i] + ddel1_dqr[i]*hbar + ddel2_dqr[i]*uconbar;
   }

// // ==================================================================
// // Additional terms for exact linearization
// // ==================================================================

   if (IMODE != 1)
   {   
      for (i = 0; i < 5; i++)
      {
         for (j = 0; j < 5; j++)
         {
            lmat[i][j] = lmat[i][j] + ( qr[i] - ql[i] )* deig1_dql[j];
            rmat[i][j] = rmat[i][j] + ( qr[i] - ql[i] )* deig1_dqr[j];
         } // j loop
      } // i loop

      for (i = 0; i < 5; i++)
      {
         lmat[1][i] = lmat[1][i] + del1*dubar_dql[i];
         rmat[1][i] = rmat[1][i] + del1*dubar_dqr[i];

         lmat[2][i] = lmat[2][i] + del1*dvbar_dql[i];
         rmat[2][i] = rmat[2][i] + del1*dvbar_dqr[i];

         lmat[3][i] = lmat[3][i] + del1*dwbar_dql[i];
         rmat[3][i] = rmat[3][i] + del1*dwbar_dqr[i];

         lmat[4][i] = lmat[4][i] + del1*dhbar_dql[i] + del2*duconbar_dql[i];
         rmat[4][i] = rmat[4][i] + del1*dhbar_dqr[i] + del2*duconbar_dqr[i];
      } // i loop
   }

// ==================================================================
//
// Compute native flux Jacobian
//
// ==================================================================
   //
   // left state
   //

   for (i = 0; i < 5; i++)
   {
      lmat1[0][i] = dro_dql[i]*uconl + rol*ducon_dql[i];

      lmat1[1][i] = dro_dql[i]*uconl*ul + rol*ducon_dql[i]*ul + rol*uconl*du_dql[i];
      lmat1[2][i] = dro_dql[i]*uconl*vl + rol*ducon_dql[i]*vl + rol*uconl*dv_dql[i];
      lmat1[3][i] = dro_dql[i]*uconl*wl + rol*ducon_dql[i]*wl + rol*uconl*dw_dql[i];

      lmat1[2][i] = lmat1[2][i] + ny*dp_dql[i];
      lmat1[1][i] = lmat1[1][i] + nx*dp_dql[i];
      lmat1[3][i] = lmat1[3][i] + nz*dp_dql[i];

      lmat1[4][i] = ( dq5_dql[i] + dp_dql[i] )*uconl + (ql[4] + pl)*ducon_dql[i];
   }


   //
   // right state
   //
   for (i = 0; i < 5; i++)
   {
      rmat1[0][i] = dro_dqr[i]*uconr + ror*ducon_dqr[i];

      rmat1[1][i] = dro_dqr[i]*uconr*ur + ror*ducon_dqr[i]*ur + ror*uconr*du_dqr[i];
      rmat1[2][i] = dro_dqr[i]*uconr*vr + ror*ducon_dqr[i]*vr + ror*uconr*dv_dqr[i];
      rmat1[3][i] = dro_dqr[i]*uconr*wr + ror*ducon_dqr[i]*wr + ror*uconr*dw_dqr[i];

      rmat1[2][i] = rmat1[2][i] + ny*dp_dqr[i];
      rmat1[1][i] = rmat1[1][i] + nx*dp_dqr[i];
      rmat1[3][i] = rmat1[3][i] + nz*dp_dqr[i];

      rmat1[4][i] = ( dq5_dqr[i] + dp_dqr[i] )*uconr + (qr[4] + pr)*ducon_dqr[i];
   }

// ==================================================================
// Final Jacobian
// ==================================================================

   for (i = 0; i < 5; i++)
   {
      for (j = 0; j < 5; j++)
      {
         lmat[i][j] = 0.5*( lmat1[i][j] - lmat[i][j] );
         rmat[i][j] = 0.5*( rmat1[i][j] - rmat[i][j] );
      }
   }

}
// ##################################################################
// END OF FILE
// ##################################################################