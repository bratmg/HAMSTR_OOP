// ##################################################################
//
// march.cpp
//
// ##################################################################
#include "solver.h"
#include <time.h>
// ##################################################################

// ##################################################################
//
// march in time (wrapper function)
//
// ##################################################################
void SOLVER::marchInTime(void)
{
   int n;
   clock_t start,stop;
   double  cpuTime = ZERO;

   // march one step in time
   for (n = 0; n < sb->nsteps; n++)
   {
      // start the clock
      start = clock();

      // step through solution for one time step
      stepSolution();

      // stop the clock
      stop  = clock();

      // compute CPU time used
      cpuTime += (((double) (stop - start)) / CLOCKS_PER_SEC);

      //
      // print to screen
      //
      if(n==0)
      {
         printf("\n --------------------------------------------------------------\n");
         printf(" |   Iter  | L_2 (density) | L_inf (density) | CPU Time (sec) |\n");
         printf(" --------------------------------------------------------------\n");
      }

      printf(" |  %*d  |  %e |   %e  |  %e  |\n",
         5,n,sb->L2Norm,sb->LInfNorm,cpuTime );

      if (n == sb->nsteps-1)
      {
         printf(" --------------------------------------------------------------\n");
         printf(" |   Iter  | L_2 (density) | L_inf (density) | CPU Time (sec) |\n");
         printf(" --------------------------------------------------------------\n");
      }

   }


}

// ##################################################################
//
// step solution
//
// ##################################################################
void SOLVER::stepSolution(void)
{
   double coef;

   //
   // Euler explicit
   // 
   if (strcmp(stepType,"euler") == 0 )
   {
      if (sb->visc) 
      {
         cout << "Not yet implemented. exit.\n";
         exit(1);
         // computeRHSv;
      }
      else
         computeRHS("explicit");

      coef = 1.;
      updateSolution(sb->q,sb->q,coef);
   }
   //
   // Runge - Kutta 3
   //
   else if (strcmp(stepType,"rk3") == 0 )
   {
      // RK step 1 
      computeRHS("explicit");
      coef=0.25;
      updateSolution(sb->q,sb->qt,coef);
      coef=8.0/15;
      updateSolution(sb->q,sb->q,coef);

      // RK step 2 
      computeRHS("explicit");
      coef=5.0/12;
      updateSolution(sb->qt,sb->q,coef);

      // RK step 3 
      computeRHS("explicit");
      coef=3.0/4.0;
      updateSolution(sb->qt,sb->q,coef);
   }
   //
   // ADI ( Alternating Direction Implicit)
   //
   else if (strcmp(stepType,"ADI")==0)
   {
      // if(sb->idual==0)
      // {
         computeRHS("implicit");
         ADI();
      // }
      // else //dual time stepping
      // {
      //    printf("Error. Dual Time stepping not yet implemented.\n");
      //    exit(1);
      //    // for (i = 0; i < NQ*mb->nCell; i++) sb->pq[i] = sb->q[i];

      //    // for(k = 0; k < sb->ndual; k++)
      //    // {
      //    //    DualcomputeRHSk(g,s,l2norm,s->cflnum);
      //    //    DualADI(g,s,s->cflnum,dt);
      //    // }

      //    // for (i=0;i<NVAR*g->ncells;i++) s->q[i] = s->pq[i];
      // }
   }
}  

// ##################################################################
//
// update solution
//
// using 1st order Euler explicit for now
// ##################################################################
void SOLVER::updateSolution(double *qold, double *qnew, double coef)
{
   int    i,j,m;
   double dtfac;
   //
   // 1st order euler explicit for now
   //
   m = 0;
   //
   for(i = 0 ; i < mb->nCell; i++)
   {
      for(j = 0; j < NQ; j++)
      {
         dtfac    = coef*sb->CFL/sb->sigma[i];

         if(sb->r[m] != sb->r[m])
            printf("oops: %d %d %d %lf\n",i,j,mb->nCell,sb->sigma[i] );

         qnew[m] = qold[m] + dtfac*sb->r[m];
         m++;
      } // j loop
   } // i loop
}
// ##################################################################
// END OF FILE
// ##################################################################
