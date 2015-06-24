// ##################################################################
//
// solver.h
//
//
// ##################################################################
#ifndef SOLVER_H
#define SOLVER_H

#include "meshBlock.h"
#include "solnBlock.h"
#include "codetypes.h"
// ##################################################################
class SOLVER
{

   private:
      double     dt;
      char      *stepType;
      // MESHBLOCK *mb;
      // SOLNBLOCK *sb;
   public:
      // basic constructor
      SOLVER(); 
      // {
      //    mb = mbin;sb=sbin;
      //    stepType = sb->scheme;
      // };

      // basic destructor
      ~SOLVER() {};

      // initialize flow
      void initFlow(void);
      void marchInTime(void);
      void stepSolution(void);
      void computeRHS(char *);
      void computeLinearRHS(void);
      void updateSolution(double *, double *, double);
      void basicScreenOutput(void);
      // computeRHSv;

      // Roe's flux (Make function pointer??)
      void fluxRoe2D(double *,double *,double *,double *,double *,double);
      void fluxRoe3D(double *,double *,double *,double *,double *,double);

      void jacobianRoe2D(double *,double *,double *,double **,double **,double);
      void jacobianRoe3D(double *,double *,double *,double **,double **,double);
      
      // matrix inversion schemes
      void ADI(void);

};
#endif
// ##################################################################
// END OF FILE
// ##################################################################
