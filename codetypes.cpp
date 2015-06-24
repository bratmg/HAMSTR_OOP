// ##################################################################
//
// codetypes.cpp
//
// ##################################################################
#include "codetypes.h"

namespace CODEVARS
{

// ==================================================================  
// MESH RELATED VARIABLES
// ==================================================================
      int      myid;
      //
      // quantity variables   
      //
      int      nNode; // number of nodes
      int      nCell; // number of cells
      int      nFace; // number of faces  
      int      nBface; // number of boundary faces  
      int      nColor; // loop colours
      int      nChain; // number of chains/loops
      int      nMaxChain; // maximum length of chain
      int      nChainFace; // faces on all chains
      int      nStrand; // number of strands

      // variables from SOLNBLOCK
      int      sol_visc;
      //
      // array variables
      //
      int      *conn; // connectivity array
      int      *faces; // face 
      int      *chainsPerColor; // loops organized per colour
      int      *faceStartPerChain;
      int      *chainConn; // connectivity along each chain
      int      *neig;
      int      *c2f;
      int      *c2chain;
      int      *cindx;
      int      *bFaces;
      double   *vol;
      double    volmax,volmin;
      double   *x; // position array      
      double  **normVec; // normal vector
      double ***refMtx; // reflection matrix

// ==================================================================  
// SOLVER RELATED VARIABLES
// ==================================================================
      int    order;
      int    timeacc;
      int    nsteps;
      int    nwrite;
      int    msweep;
      int    test;
      int    outform;
      double Mach;
      double alpha;
      double beta;
      double rey;
      double CFL;
      double dt;
      char   timeInteg[20];
      char   scheme[20];

      double *q;     // q -variables [rho rho*u rho*v e]
      double *qt;
      double *qtt;
      double *dq;    // \delta q -variables [rho rho*u rho*v e]
      double *ddq;   // \delta\delta q - variables [rho rho*u rho*v e]
      double *ddqb;  // \delta\delta q - variables [rho rho*u rho*v e]
      double *ddqf;  // \delta\delta q - variables [rho rho*u rho*v e]
      double *r;     // solution residual at cell centroids
      double *r0;    // solution residual at cell centroids
      double *sigma; // line integral of spectral radius per cell
      double *dtac;  // time scaling - nominally equaly to dt/dv but can be 
                     // combination of pseudo time step as well
      // int *itag;
      int    nGhost;
      double uinf,vinf,winf,einf,pinf,rinf; // inf primitive variables
      double Gamma,invGamma,gm1,c2b,rgas,pr,prtr; // fluid thermodynamic properties
      double LInfNorm,L2Norm;
      // 
      // work arrays for processing fluxes
      // 
      double  **flux;   // flux vector
      double  **fv;  // viscous flux vector
      double  **ql;  // left state (cons. variables)
      double  **qr;  // right state (cons. variables)
      double  **dqr; // change in left state
      double  **dql; // change in right state
      double  **df;  // change in flux
      double  **flux2;
      double  **F;
      double  **Q;
      double ***A;
      double ***B;
      double ***C;
      double ***D;   // diagonal matrix
      FACEMAT  *ff;  // NQ x NQ matrix at a particular face (Jacobian)
      //
      // variables from MESHBLOCK
      //
      int    visc;

};

// ##################################################################
// END OF FILE
// ##################################################################
