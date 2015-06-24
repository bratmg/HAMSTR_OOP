// ##################################################################
//
// codetypes.h
//
// ##################################################################
#ifndef CODETYPES_H
#define CODETYPES_H
// ##################################################################
#include "utils.h"

typedef struct FACEMAT; //forward declaration

namespace CODEVARS
{

// ==================================================================  
// MESH RELATED VARIABLES
// ==================================================================
      extern int      myid;
      //
      // quantity variables   
      //
      extern int      nNode; // number of nodes
      extern int      nCell; // number of cells
      extern int      nFace; // number of faces  
      extern int      nBface; // number of boundary faces  
      extern int      nColor; // loop colours
      extern int      nChain; // number of chains/loops
      extern int      nMaxChain; // maximum length of chain
      extern int      nChainFace; // faces on all chains
      extern int      nStrand; // number of strands
      //
      // array variables
      //
      extern int      *conn; // connectivity array
      extern int      *faces; // face 
      extern int      *chainsPerColor; // loops organized per colour
      extern int      *faceStartPerChain;
      extern int      *chainConn; // connectivity along each chain
      extern int      *neig;
      extern int      *c2f;
      extern int      *c2chain;
      extern int      *cindx;
      extern int      *bFaces;
      extern double   *vol;
      extern double    volmax,volmin;
      extern double   *x; // position array      
      extern double  **normVec; // normal vector
      extern double ***refMtx; // reflection matrix

// ==================================================================  
// SOLVER RELATED VARIABLES
// ==================================================================
      extern int    order;
      extern int    timeacc;
      extern int    nsteps;
      extern int    nwrite;
      extern int    msweep;
      extern int    test;
      extern int    outform;
      extern double Mach;
      extern double alpha;
      extern double beta;
      extern double rey;
      extern double CFL;
      extern double dt;
      extern char   timeInteg[20];
      extern char   scheme[20];

      extern double *q;     // q -variables [rho rho*u rho*v e]
      extern double *qt;
      extern double *qtt;
      extern double *dq;    // \delta q -variables [rho rho*u rho*v e]
      extern double *ddq;   // \delta\delta q - variables [rho rho*u rho*v e]
      extern double *ddqb;  // \delta\delta q - variables [rho rho*u rho*v e]
      extern double *ddqf;  // \delta\delta q - variables [rho rho*u rho*v e]
      extern double *r;     // solution residual at cell centroids
      extern double *r0;    // solution residual at cell centroids
      extern double *sigma; // line integral of spectral radius per cell
      extern double *dtac;  // time scaling - nominally equaly to dt/dv but can be 
                     // combination of pseudo time step as well
      // int *itag;
      extern int    nGhost;
      extern double uinf,vinf,winf,einf,pinf,rinf; // inf primitive variables
      extern double Gamma,invGamma,gm1,c2b,rgas,pr,prtr; // fluid thermodynamic properties
      extern double LInfNorm,L2Norm;
      // 
      // work arrays for processing fluxes
      // 
      extern double  **flux;   // flux vector
      extern double  **fv;  // viscous flux vector
      extern double  **ql;  // left state (cons. variables)
      extern double  **qr;  // right state (cons. variables)
      extern double  **dqr; // change in left state
      extern double  **dql; // change in right state
      extern double  **df;  // change in flux
      extern double  **flux2;
      extern double  **F;
      extern double  **Q;
      extern double ***A;
      extern double ***B;
      extern double ***C;
      extern double ***D;   // diagonal matrix
      extern FACEMAT  *ff;  // NQ x NQ matrix at a particular face (Jacobian)
      //
      // variables from MESHBLOCK
      //
      extern int    visc;

};


// ##################################################################
typedef struct HOLEMAP
{
  int existWall;
  int nx[3];
  int *samLocal;
  int *sam;
  double extents[6];
} HOLEMAP;

// ##################################################################
typedef struct OBB
{
  double xc[3];
  double dxc[3];
  double vec[3][3];
}OBB;

// ##################################################################
typedef struct DONORLIST
{
  int donorData[3];
  double donorRes;
  struct DONORLIST *next;
} DONORLIST;

// ##################################################################
typedef struct PACKET
{
  int nints;
  int nreals;
  int *intData;
  REAL *realData;
} PACKET;

// ##################################################################
typedef struct INTERPLIST
{
  int cancel;
  int nweights;
  int receptorInfo[2];
  int *inode;
  double *weights;
} INTERPLIST;

// ##################################################################
typedef struct INTEGERLIST
{
  int inode;
  struct INTEGERLIST *next;
} INTEGERLIST;

// ##################################################################
typedef struct FACEMAT
{
  double lmat[NQ][NQ];
  double rmat[NQ][NQ];
} FACEMAT;

#endif
// ##################################################################
// END OF FILE
// ##################################################################
