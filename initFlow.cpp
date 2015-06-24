// ##################################################################
//
// initFlow.cpp
//
// contains routines to initialize the flow based on the
// test condition and BC conditions
//
// ##################################################################
#include "solver.h"

using namespace CODEVARS;
// ##################################################################

void SOLVER::initFlow(void)
{

   //
   // local variable initialization
   //
   int i,j,m,workSize;
   //
   // basic property initialization
   //
   Gamma    = 1.4;
   gm1      = Gamma - 1.;
   invGamma = 1./Gamma;
   pr       = 0.72;
   prtr     = THIRD;
   rinf     = ONE;
   pinf     = ONE/Gamma;

   rey = rey/Mach;

   // set number of ghost cells
   nGhost = 2;                    // MUSCL/WENO3
   if(order == 5) nGhost = 3; // WENO5

   // error checking
   if (test != 0)
   {
      cout << " HAMSTR: Particular test condition("<<test<<
         ") not yet implemented.\n";
      exit(1);
   }   
   // Mach scaled velocities
   if (test == 0)
   {
#ifdef Dim3
      uinf = Mach*cos(alpha*DEG2RAD)*cos(beta*DEG2RAD);
      vinf = Mach*cos(alpha*DEG2RAD)*sin(beta*DEG2RAD);
      winf = Mach*sin(alpha*DEG2RAD);
#else
      uinf = Mach*cos(alpha*DEG2RAD);
      vinf = Mach*sin(alpha*DEG2RAD);
#endif
   }      
   // energy
   einf = pinf/(gm1) + HALF*rinf*
              (uinf*uinf + vinf*vinf);
#ifdef Dim3   
   einf = pinf/(gm1) + HALF*rinf*(uinf*uinf + 
              vinf*vinf + winf*winf);
#endif
   
  
   //
   // Far-field BC
   //
   m = 0;
   for (i = 0; i < nCell; i++)
   {
      q[m]  = rinf;            
      qt[m] = q[m];
      m++;
      q[m]  = rinf*uinf;
      qt[m] = q[m];
      m++;
      q[m]  = rinf*vinf;
      qt[m] = q[m];
      m++;

#ifdef Dim3
      q[m]  = rinf*winf;
      qt[m] = q[m];
      m++;      
#endif

      // E*rho
      q[m]  = einf;
      qt[m] = q[m];
      m++;

   }

}

// ##################################################################
// END OF FILE
// ##################################################################

