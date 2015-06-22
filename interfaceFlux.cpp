// ##############################################################################
// #     Computes total convective flux across a face using Roe's approximate
// #                 Riemann solver given left and right conserved states
// #
// #     nxyz(3)     - Three components of face normal vector. These can be
// #                   dimensional (i.e including face area) or non-dimensional.
// #                   Returned flux does not include face area.
// #     ql(5),qr(5) - Vector of conserved variables (ro, ro*u, ro*v, ro*w, Et)
// #     flux(5)     - Vector of flux across face. Not multiplied by face area.
// #     specRadius  - spectral radius (largest eigen value)
// #     gam         - Ratio of specific heats
// #
// ##############################################################################

// ################# (C) Copyright Karthik Mani 2005 ############################

#include "solver.h"
#include "math.h"
// ##############################################################################
//
// THREE-DIMENSIONAL ROX FLUX
//
// ##############################################################################
void SOLVER::fluxRoe2D(double *nxyz, // normal vector of face
                       double *ql,   // left states
                       double *qr,   // right states
                       double *flux, // 
                       double *specRadius,
                       double  gamma)
{
   //
   // Local variable declaration
   //
   int    i;
   double gm1;
   double nxd,nyd,area,invArea;
   double rol,ul,vl,pl,hl,cl;
   double ror,ur,vr,pr,hr,cr;
   double uconl,uconr;
   double nx,ny,nz;
   double fm[4];
   double robar,srl,srr,denom;
   double ubar,vbar,hbar,uconbar,uvwbar,cbar;
   double dp,dro,du,dv,ducon,duvw,dc;
   double eig1,eig2,eig3;
   double term;
   double lamL,lamR,lamB;
   
   double fact,A,B,term1,term2,del1,del2;
   double eps1,eps2,eps3;
   double eig1L,eig2L,eig3L;
   double eig1R,eig2R,eig3R;
   double fl[4],fr[4];

   gm1 = gamma - 1.0;

   nxd = nxyz[0];
   nyd = nxyz[1];
    
   area    = sqrt(nxd*nxd + nyd*nyd);
   invArea = 1./area;

   nx = nxd*invArea;
   ny = nyd*invArea;
   //
   // back calculate left and right primitive states
   //
   rol = ql[0];
   ul  = ql[1]/ql[0];
   vl  = ql[2]/ql[0];
   pl  = gm1*( ql[3] - 0.5 * rol * (ul*ul + vl*vl) );
   hl  = (ql[3] + pl)/rol;
   cl  = sqrt(gamma*pl/rol);
    
   ror = qr[0];
   ur  = qr[1]/qr[0];
   vr  = qr[2]/qr[0];
   pr  = gm1*( qr[3] - 0.5 * ror * (ur*ur + vr*vr) );
   hr  = (qr[3] + pr)/ror;
   cr  = sqrt(gamma*pr/ror);
   //
   // primitive state differences
   //
   dro = ror - rol;
   du  =  ur - ul;
   dv  =  vr - vl;
   dp  =  pr - pl;
   dc  =  cr - cl;
   //
   // face normal velocities (contravariant)
   //
   uconr = ur*nx + vr*ny;
   uconl = ul*nx + vl*ny;
   //
   // Roe average state
   //
   fact = sqrt(ror/rol);
   A    = 1.0 /(1.0 + fact);
   B    = fact/(1.0 + fact);
   
   robar   = rol*fact;
   ubar    = ul*A + ur*B;
   vbar    = vl*A + vr*B;
   hbar    = hl*A + hr*B;
   cbar    = gm1*(hbar - 0.5*(ubar*ubar + vbar*vbar));
   cbar    = sqrt(cbar);
   uconbar = ubar*nx + vbar*ny;
   //
   // Eigenvalues
   //
   eig1 = fabs(uconbar);
   eig2 = fabs(uconbar + cbar);
   eig3 = fabs(uconbar - cbar);

   *specRadius = eig1 + cbar;
   
   if ( (*specRadius != *specRadius))
   {  
      tracef(eig1);
      tracef(cbar);
      tracef(ql[0]);
      tracef(qr[0]);
      exit(1);
   }


   //
   //
   // ROE FLUX CALCULATIONS
   //
   //

   term1 = -eig1 + 0.5*(eig2 + eig3);
   term2 = 0.5*(eig2 - eig3);
   del1  = term1*dp/cbar/cbar + term2*robar*(uconr - uconl)/cbar;
   del2  = term1*(uconr - uconl)*robar + term2*dp/cbar;

   for (i = 0; i < 5; i++)
      fm[i] = eig1*(qr[i]-ql[i]);

   fm[0] =  fm[0] + del1;
   fm[1] =  fm[1] + del1*ubar + del2*nx;
   fm[2] =  fm[2] + del1*vbar + del2*ny;
   fm[3] =  fm[3] + del1*hbar + del2*uconbar;
   //
   // native flux for left and right states 
   //
   fl[0] = rol*uconl;
   fr[0] = ror*uconr;

   fl[1] = rol*uconl*ul + nx*pl;
   fr[1] = ror*uconr*ur + nx*pr;

   fl[2] = rol*uconl*vl + ny*pl;
   fr[2] = ror*uconr*vr + ny*pr;

   fl[3] = ( ql[3] + pl ) * uconl;
   fr[3] = ( qr[3] + pr ) * uconr;

   //
   // total flux
   //
   for (i = 0; i < 4; i++)
   {
      flux[i] = 0.5*( fl[i] + fr[i] - fm[i] );
      flux[i] *= area;
   }
   *specRadius *= area;
}

// ##############################################################################
//
// THREE-DIMENSIONAL ROX FLUX
//
// ##############################################################################
void SOLVER::fluxRoe3D(double *nxyz, // normal vector of face
                       double *ql,   // left states
                       double *qr,   // right states
                       double *flux, // 
                       double *specRadius,
                       double gamma)
{
   //
   // Local variable declaration
   //
   int    i;
   double gm1;
   double nxd,nyd,nzd,area,invArea;
   double rol,ul,vl,wl,pl,hl,cl;
   double ror,ur,vr,wr,pr,hr,cr;
   double uconl,uconr;
   double nx,ny,nz;
   double fm[5];
   double robar,srl,srr,denom;
   double ubar,vbar,wbar,hbar,uconbar,uvwbar,cbar;
   double dp,dro,du,dv,dw,ducon,duvw,dc;
   double eig1,eig2,eig3;
   double term;
   double lamL,lamR,lamB;
   
   double fact,A,B,term1,term2,del1,del2;
   double eps1,eps2,eps3;
   double eig1L,eig2L,eig3L;
   double eig1R,eig2R,eig3R;
   double fl[5],fr[5];

   gm1 = gamma - 1.0;

   nxd = nxyz[0];
   nyd = nxyz[1];
   nzd = nxyz[2];
    
   area    = sqrt(nxd*nxd + nyd*nyd + nzd*nzd);
   invArea = 1./area;

   nx = nxd*invArea;
   ny = nyd*invArea;
   nz = nzd*invArea;
   //
   // back calculate left and right primitive states
   //
   rol = ql[0];
   ul  = ql[1]/ql[0];
   vl  = ql[2]/ql[0];
   wl  = ql[3]/ql[0];
   pl  = gm1*( ql[4] - 0.5 * rol * (ul*ul + vl*vl + wl*wl) );
   hl  = (ql[4] + pl)/rol;
   cl  = sqrt(gamma*pl/rol);
    
   ror = qr[0];
   ur  = qr[1]/qr[0];
   vr  = qr[2]/qr[0];
   wr  = qr[3]/qr[0];
   pr  = gm1*( qr[4] - 0.5 * ror * (ur*ur + vr*vr + wr*wr) );
   hr  = (qr[4] + pr)/ror;
   cr  = sqrt(gamma*pr/ror);
   //
   // primitive state differences
   //
   dro = ror - rol;
   du  =  ur - ul;
   dv  =  vr - vl;
   dw  =  wr - wl;
   dp  =  pr - pl;
   dc  =  cr - cl;
   //
   // face normal velocities (contravariant)
   //
   uconr = ur*nx + vr*ny + wr*nz;
   uconl = ul*nx + vl*ny + wl*nz;
   //
   // Roe average state
   //
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
   //
   // Eigenvalues
   //
   eig1 = fabs(uconbar);
   eig2 = fabs(uconbar + cbar);
   eig3 = fabs(uconbar - cbar);

   *specRadius = eig1+cbar;

   //
   //
   // ROE FLUX CALCULATIONS
   //
   //

   term1 = -eig1 + 0.5*(eig2 + eig3);
   term2 = 0.5*(eig2 - eig3);
   del1  = term1*dp/cbar/cbar + term2*robar*(uconr - uconl)/cbar;
   del2  = term1*(uconr - uconl)*robar + term2*dp/cbar;

   for (i = 0; i < 5; i++)
      fm[i] = eig1*(qr[i]-ql[i]);

   fm[0] =  fm[0] + del1;
   fm[1] =  fm[1] + del1*ubar + del2*nx;
   fm[2] =  fm[2] + del1*vbar + del2*ny;
   fm[3] =  fm[3] + del1*wbar + del2*nz;
   fm[4] =  fm[4] + del1*hbar + del2*uconbar;
   //
   // native flux for left and right states 
   //
   fl[0] = rol*uconl;
   fr[0] = ror*uconr;

   fl[1] = rol*uconl*ul + nx*pl;
   fr[1] = ror*uconr*ur + nx*pr;

   fl[2] = rol*uconl*vl + ny*pl;
   fr[2] = ror*uconr*vr + ny*pr;

   fl[3] = rol*uconl*wl + nz*pl;
   fr[3] = ror*uconr*wr + nz*pr;

   fl[4] = ( ql[4] + pl ) * uconl;
   fr[4] = ( qr[4] + pr ) * uconr;

   //
   // total flux
   //
   for (i = 0; i < 5; i++)
      flux[i] = 0.5*( fl[i] + fr[i] - fm[i] );


}
// ##############################################################################
// END OF FILE
// ##############################################################################
