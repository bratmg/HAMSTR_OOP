// ##################################################################
//
// meshBlock.h
//
// contains routines to set up the initial mesh and re-organizing
// the strcutures to make is cache friendly
// ##################################################################
#ifndef MESHBLOCK_H
#define MESHBLOCK_H

#include "utils.h"
#include "codetypes.h"

// ##################################################################
class MESHBLOCK
{
   // solver needs access to private variables
   friend class SOLVER;

   private:
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

   public:
      // basic contructor and destructor
      MESHBLOCK() {};

      ~MESHBLOCK();

      // read data from files
      void readMeshData(int myid);

      // preprocess data 
      void preprocess(void);      

      //
      // utility functions
      //
      void setViscosity(int visc) {sol_visc = visc;};

};
#endif
// ##################################################################
// END OF FILE
// ##################################################################
