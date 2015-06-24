// ##################################################################
//
// reconstruction.cpp
//
// ##################################################################
#include "reconstruction.h"

namespace RECONSTRUCTION
{

// ##################################################################
//
// MUSCL RECONSTRUCTION
//
// ##################################################################
   void MUSCL(double **flux,
               double **ql,
               double **qr,
               double **f2,
               int is,
               int ie,
               double th,
               double qt,
               double eps,
               int imax,
               int nq)
   {  
      int    i,n,im;
      double thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt;

      im = imax-1;
      //
      // First order in space
      //
      if (qt == 0.0)
      {
         for(n=0;n<nq;n++)
         {
            for(i=is;i<=ie;i++)
            {
               ql[i][n]=flux[i][n];
               qr[i][n]=flux[i][n];
            }
         }
         return;
      }
      //
      // Third order in space
      //
      thm = 1.0 - th;
      thp = 1.0 + th;
      
      for(n = 0; n < nq; n++)
      {
         for(i=0;i<im;i++)
            f2[i][n] = flux[i+1][n]-flux[i][n];

         for(i=is;i<=ie;i++)
         {
            f2i      = f2[i][n];
            f2i1     = f2[i-1][n];
            a1       = 3.0*(f2i*f2i1+eps);
            a2       = 2.0*(f2i-f2i1)*(f2i-f2i1) + a1;
            f3i      =  a1/a2 ;
            f3qt     = qt*f3i;
            ql[i][n] = flux[i][n]+f3qt*( thm*f2i1 + thp*f2i );
            qr[i][n] = flux[i][n]-f3qt*( thp*f2i1 + thm*f2i );
         } // i loop
      } // n loop
   }

// ##################################################################
//
// WENO 5 RECONSTRUCTION
//
// ##################################################################

   double WENO5(double a,  //fluxes at each edge
                double b, // left states
                double c, // right states
                double d,   // start index of chain
                double e,   // end index of chain
                double epsw)
   {
      int    p;
      double b1,b2,djm1,ejm1,dj,ej,djp1,ejp1;
      double dis0,dis1,dis2,q30,q31,q32,d01,d02,a1ba0,a2ba0;
      double w0,w1,w2;
      double sol;
     

      p     = 1;
      
      b1    = 13./12.;
      b2    = 1./6.;
      djm1  = a-2.*b+c;
      ejm1  = a-4.*b+3.*c;
      dj    = b-2.*c+d;
      ej    = b-d;
      djp1  = c-2.*d+e;
      ejp1  = 3.*c-4.*d+e;
      dis0  = b1*djm1*djm1+0.25*ejm1*ejm1+epsw;
      dis1  = b1*dj*dj+0.25*ej*ej+epsw;
      dis2  = b1*djp1*djp1+0.25*ejp1*ejp1+epsw;
      
      dis0  = pow(dis0,p);
      dis1  = pow(dis1,p);
      dis2  = pow(dis2,p);

      q30   = 2.*a-7.*b+11.*c;
      q31   = -b+5.*c+2.*d;
      q32   = 2.*c+5.*d-e;
      d01   = dis0/dis1;
      d02   = dis0/dis2;
      a1ba0 = 6.*d01;
      a2ba0 = 3.*d02;
      w0    = 1.0/(1.0+a1ba0+a2ba0);
      w1    = a1ba0*w0;
      w2    = 1.-w0-w1;
      sol   = b2*(w0*q30+w1*q31+w2*q32);  //output
      
      return sol;
           
   } // end function


// ##################################################################
//
// WENO RECONSTRUCTION (WRAPPER FUNCTION)
//   
// ##################################################################

   void WENO(double **f,  //fluxes at each edge
             double **ql, // left states
             double **qr, // right states
             int    is,   // start index of chain
             int    ie,   // end index of chain
             double epsw,  // epsilon in limiter
             int    imax, // chainSize
             int    nq)   // number of variables
   {
      int    i,n,im,k;
      double f0[imax-1],f1[imax-1],f2[imax-1];

      epsw = 1e-15;
      im   = imax-2;

      for (n = 0; n < nq; n++) // number of flow variable
      {
         for (i = 0; i <= ie; i++) // loop through the chain
            f0[i] = f[i][n];

         for (i = is; i < im; i++)
         {
            ql[i][n] = WENO5(f0[i-2],f0[i-1],f0[i],f0[i+1],f0[i+2],epsw);
            qr[i][n] = WENO5(f0[i+2],f0[i+1],f0[i],f0[i-1],f0[i-2],epsw);
         }
      } // n loop

   } //end function

// ##################################################################
//
// MUSCL RECONSTRUCTION DERIVATIVE
//
// ##################################################################

   void MUSCL_DERIV(double **f,
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
                    int      nq)
   {  
      //
      // local variable declaration
      //
      int    i,n,im;
      double thm,thp,f2i,f2i1,a1,a2,epsf,f3i,f3qt;
      double df2i,df2i1,da1,da2,df3i,df3qt;
      //
      // 1st order 
      //
      im = imax-1;
      if (qt == 0.0)
      {
         for(n=0;n<nq;n++)
         {
           for(i=is;i<=ie;i++)
           {
             ql[i][n]=dq[i][n];
             qr[i][n]=dq[i][n];
           }
      }
      return;
      }
      //
      // 3rd order 
      //
      thm    = 1.0 - th;
      thp    = 1.0 + th;
      //
      // find the tangent
      //
      for(n = 0; n < nq; n++)
      {
         for(i = 0; i < im; i++) 
            f2[i][n] = f[i+1][n]-f[i][n];

         for(i = is; i <= ie; i++)
         {
            f2i    = f2[i][n];
            f2i1   = f2[i-1][n];
            df2i   = dq[i+1][n] - dq[i  ][n];
            df2i1  = dq[i  ][n] - dq[i-1][n];
 
            a1     = 3.0*(f2i*f2i1+eps);
            da1    = 3.0*(df2i*f2i1+f2i*df2i1);
            a2     = 2.0*(f2i-f2i1)*(f2i-f2i1) + a1;
            da2    = 2.0*2.0*(f2i-f2i1)*(df2i-df2i1)+da1;
            f3i    =  a1/a2 ;
            df3i   =  da1/a2 -a1*da2/(a2*a2);
            f3qt   = qt*f3i;
            df3qt  = qt*df3i;
 
            ql[i][n]= dq[i][n]+df3qt*( thm*f2i1 + thp*f2i ) + f3qt*(thm*df2i1+thp*df2i);
            qr[i][n]= dq[i][n]-df3qt*( thp*f2i1 + thm*f2i ) - f3qt*(thp*df2i1+thm*df2i);
         } // i loop
      }// n loop
  } // end function

//##########################################################
//
// WENO RECONSTRUCTION DERIVATIVE 
//
//##########################################################
  double WENO_5(double a,  //fluxes at each edge
                double b, // left states
                double c, // right states
                double d,   // start index of chain
                double e,   // end index of chain
                double a_d,
                double b_d, // left states
                double c_d, // right states
                double d_d,   // start index of chain
                double e_d,   // end index of chain
                double epsw)
  {
    double b1,b2,djm1,ejm1,dj,ej,djp1,ejp1;
    double dis0,dis1,dis2,q30,q31,q32,d01,d02,a1ba0,a2ba0;
    double w0,w1,w2;
    double sol;
    
    double djm1_d,ejm1_d,dj_d,ej_d,djp1_d,ejp1_d;
    double dis0_d,dis1_d,dis2_d,q30_d,q31_d,q32_d;
    double term1,term2,term3,dbeta0,dbeta1,dbeta2,beta0,beta1,beta2;
    double dw0_up_1,dw0_up_2,dw1_up_1,dw1_up_2,dw2_up_1,dw2_up_2;
    double dw0,dw1,dw2,dw0_down,dw1_down,dw2_down;

    //original term
    b1   = 13./12.;
    b2   = 1./6.;
    djm1 = a-2.*b+c; //beta0_0
    ejm1 = a-4.*b+3.*c; //beta0_1
    dj   = b-2.*c+d; //beta1_0
    ej   = b-d; //beta1_1
    djp1 = c-2.*d+e; //beta2_0
    ejp1 = 3.*c-4.*d+e; //beta2_1
    dis0 = b1*djm1*djm1+0.25*ejm1*ejm1+epsw;
    dis1 = b1*dj*dj+0.25*ej*ej+epsw;
    dis2 = b1*djp1*djp1+0.25*ejp1*ejp1+epsw;


    q30   = 2.*a-7.*b+11.*c;
    q31   = -b+5.*c+2.*d;
    q32   = 2.*c+5.*d-e;
    d01   = dis0/dis1;
    d02   = dis0/dis2;
    a1ba0 = 6.*d01;
    a2ba0 = 3.*d02;
    w0    = 1.0/(1.0+a1ba0+a2ba0);
    w1    = a1ba0*w0;
    w2    = 1.-w0-w1;


    //deriveative term
    djm1_d = a_d-2.*b_d+c_d; //dbeta 0_0
    ejm1_d = a_d-4.*b_d+3.*c_d; //dbeta 0_1
    dj_d   = b_d-2.*c_d+d_d; //dbeta 1_0
    ej_d   = b_d - d_d; //dbeta1_1
    djp1_d = c_d -2.*d_d+e_d; //dbeta2_0
    ejp1_d = 3.*c_d-4.*d_d+e_d; //dbeta2_1
    dis0_d = b1*2.*djm1*djm1_d+0.25*2.*ejm1*ejm1_d;
    dis1_d = b1*2.*dj*dj_d + 0.25*2.*ej*ej_d;
    dis2_d = b1*2.*djp1*djp1_d + 0.25*2.*ejp1*ejp1_d;
    q30_d  = 2.*a_d - 7.*b_d + 11.*c_d;
    q31_d  = -b_d+5.*c_d+2.*d_d;
    q32_d  = 2.*c_d + 5.*d_d -e_d;


    //let's find!!
    term1 = b2*w0*q30_d;
    term2 = b2*w1*q31_d;
    term3 = b2*w2*q32_d;

    //dbeta
    dbeta0 = b1*2*djm1*djm1_d + 0.25*2*ejm1*ejm1_d;
    dbeta1 = b1*2*dj*dj_d     + 0.25*2*ej*ej_d;
    dbeta2 = b1*2*djp1*djp1_d + 0.25*2*ejp1*ejp1_d;

    beta0 = dis0;
    beta1 = dis1;
    beta2 = dis2;

    //dw0
    dw0_down = beta1*beta2+6*beta0*beta2+3*beta0*beta1;
    dw0_up_1 = (dbeta1*beta2 + beta1*dbeta2)*dw0_down;
    dw0_up_2 = beta1*beta2*(dbeta1*beta2+beta1*dbeta2+6*(dbeta0*beta2+beta0*dbeta2)+3*(dbeta0*beta1+beta0*dbeta1));

    dw0 = (dw0_up_1 - dw0_up_2)/(dw0_down*dw0_down);

    //dw1
    dw1_down = dw0_down;
    dw1_up_1 = 6*(dbeta0*beta2 + beta0*dbeta2)*dw1_down;
    dw1_up_2 = 6*dw0_up_2/beta1*beta0;

    dw1 = (dw1_up_1 - dw1_up_2)/(dw1_down*dw1_down);

    //dw2
    dw2_down = dw0_down;
    dw2_up_1 = 3*(dbeta0*beta1 + beta0*dbeta1)*dw2_down;
    dw2_up_2 = 3*dw0_up_2/beta2*beta0;

    dw2 = (dw2_up_1 - dw2_up_2)/(dw2_down*dw2_down);

    sol = term1 + term2 +term3 + b2*dw0*q30 + b2*dw1*q31 + b2*dw2*q32; 

    return sol;
          
          
  } // end function
// ##################################################################
//
// WENO RECONSTRUCTION DERIVATIVE (WRAPPER)
//
// ##################################################################
  void WENO_DERIV(double **f,
                  double **ql,
                  double **qr,
                  double **dq,
                  int      is,
                  int      ie,
                  double   eps,
                  int      imax,
                  int      nq)
  {  
    int    i,n,im;
    double epsf;
    double f0[imax],df0[imax];

    im = imax-2;
    for (n=0;n<nq;n++)
    {
      for(i=0;i<=ie;i++) // loop through the chain
      {
        f0[i] = f[i][n];
        df0[i] = dq[i][n];
      }

      for (i=is;i<im;i++)
      {
        ql[i][n]=WENO_5( f0[i-2], f0[i-1], f0[i], f0[i+1], f0[i+2],
                        df0[i-2],df0[i-1],df0[i],df0[i+1],df0[i+2],eps);

        qr[i][n]=WENO_5( f0[i+2], f0[i+1], f0[i], f0[i-1], f0[i-2],
                        df0[i+2],df0[i+1],df0[i],df0[i-1],df0[i-2],eps);
       
      }
    }
  } // end function



};
// ##################################################################
// END OF FILE
// ##################################################################
