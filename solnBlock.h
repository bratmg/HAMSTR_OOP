// ##################################################################
//
// solnBlock.h
//
// ##################################################################
#ifndef SOLNBLOCK_H
#define SOLNBLOCK_H

#include "utils.h"
#include "codetypes.h"


using namespace CODEVARS;
// ##################################################################
class SOLNBLOCK
{
   // friend class SOLVER;
   // // only the getVariable function of MESHBLOCK class has access
   // // friend int MESHBLOCK::getViscosity(void) {return visc};

   // private:
   //    int    order;
   //    int    timeacc;
   //    int    nsteps;
   //    int    nwrite;
   //    int    msweep;
   //    int    test;
   //    int    outform;
   //    double Mach;
   //    double alpha;
   //    double beta;
   //    double rey;
   //    double CFL;
   //    double dt;
   //    char   timeInteg[20];
   //    char   scheme[20];

   //    double *q;     // q -variables [rho rho*u rho*v e]
   //    double *qt;
   //    double *qtt;
   //    double *dq;    // \delta q -variables [rho rho*u rho*v e]
   //    double *ddq;   // \delta\delta q - variables [rho rho*u rho*v e]
   //    double *ddqb;  // \delta\delta q - variables [rho rho*u rho*v e]
   //    double *ddqf;  // \delta\delta q - variables [rho rho*u rho*v e]
   //    double *r;     // solution residual at cell centroids
   //    double *r0;    // solution residual at cell centroids
   //    double *sigma; // line integral of spectral radius per cell
   //    double *dtac;  // time scaling - nominally equaly to dt/dv but can be 
   //                   // combination of pseudo time step as well
   //    // int *itag;
   //    int    nGhost;
   //    double uinf,vinf,winf,einf,pinf,rinf; // inf primitive variables
   //    double gamma,invGamma,gm1,c2b,rgas,pr,prtr; // fluid thermodynamic properties
   //    double LInfNorm,L2Norm;
   //    // 
   //    // work arrays for processing fluxes
   //    // 
   //    double  **f;   // flux vector
   //    double  **fv;  // viscous flux vector
   //    double  **ql;  // left state (cons. variables)
   //    double  **qr;  // right state (cons. variables)
   //    double  **dqr; // change in left state
   //    double  **dql; // change in right state
   //    double  **df;  // change in flux
   //    double  **f2;
   //    double  **F;
   //    double  **Q;
   //    double ***A;
   //    double ***B;
   //    double ***C;
   //    double ***D;   // diagonal matrix
   //    FACEMAT  *ff;  // NQ x NQ matrix at a particular face (Jacobian)
   //    //
   //    // variables from MESHBLOCK
   //    //
   //    int    visc;
      

   public:
      // basic contructor
      SOLNBLOCK() {};

      // basic destructor
      ~SOLNBLOCK() {};

      // read input from files
      void readSolutionInputs(void);

      // initialize the flow
      void initFlow(void);
      

};
#endif
// ##################################################################
// END OF FILE
// ##################################################################
