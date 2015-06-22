// ##################################################################
//
// initFlow.cpp
//
// contains routines to initialize the flow based on the
// test condition and BC conditions
//
// ##################################################################
#include "solver.h"
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
   sb->gamma    = 1.4;
   sb->gm1      = sb->gamma - 1.;
   sb->invGamma = 1./sb->gamma;
   sb->pr       = 0.72;
   sb->prtr     = THIRD;
   sb->rinf     = ONE;
   sb->pinf     = ONE/sb->gamma;

   sb->rey = sb->rey/sb->Mach;

   // set number of ghost cells
   sb->nGhost = 2;                    // MUSCL/WENO3
   if(sb->order == 5) sb->nGhost = 3; // WENO5

   // error checking
   if (sb->test != 0)
   {
      cout << " HAMSTR: Particular test condition("<<sb->test<<
         ") not yet implemented.\n";
      exit(1);
   }   
   // Mach scaled velocities
   if (sb->test == 0)
   {
#ifdef Dim3
      sb->uinf = sb->Mach*cos(sb->alpha*DEG2RAD)*cos(sb->beta*DEG2RAD);
      sb->vinf = sb->Mach*cos(sb->alpha*DEG2RAD)*sin(sb->beta*DEG2RAD);
      sb->winf = sb->Mach*sin(sb->alpha*DEG2RAD);
#else
      sb->uinf = sb->Mach*cos(sb->alpha*DEG2RAD);
      sb->vinf = sb->Mach*sin(sb->alpha*DEG2RAD);
#endif
   }      
   // energy
   sb->einf = sb->pinf/(sb->gm1) + HALF*sb->rinf*
              (sb->uinf*sb->uinf + sb->vinf*sb->vinf);
#ifdef Dim3   
   sb->einf = sb->pinf/(sb->gm1) + HALF*sb->rinf*(sb->uinf*sb->uinf + 
              sb->vinf*sb->vinf + sb->winf*sb->winf);
#endif
   
  
   //
   // Far-field BC
   //
   m = 0;
   for (i = 0; i < mb->nCell; i++)
   {
      sb->q[m]  = sb->rinf;            
      sb->qt[m] = sb->q[m];
      m++;
      sb->q[m]  = sb->rinf*sb->uinf;
      sb->qt[m] = sb->q[m];
      m++;
      sb->q[m]  = sb->rinf*sb->vinf;
      sb->qt[m] = sb->q[m];
      m++;

#ifdef Dim3
      sb->q[m]  = sb->rinf*sb->winf;
      sb->qt[m] = sb->q[m];
      m++;      
#endif

      // E*rho
      sb->q[m]  = sb->einf;
      sb->qt[m] = sb->q[m];
      m++;

   }

}

// ##################################################################
// END OF FILE
// ##################################################################

