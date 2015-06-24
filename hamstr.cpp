// ##################################################################
//
// hamstr.cpp
//
// Treated as a signal-caller. Tells other routines what to do.
// ##################################################################
#include "hamstr.h"
// ##################################################################

// ==================================================================
// constructor
// ==================================================================

HAMSTR::HAMSTR(const int nprocs, const int myid)
{

   // set ngrid
   nGrid  = nprocs;
   procID = myid;

   // initialize
   mb   = new MESHBLOCK[1];
   sb   = new SOLNBLOCK[1];
   
}


// ==================================================================
// read inputs hamstr
// read the input files and data files for mesh
// ==================================================================

void HAMSTR::readInputs(int myid)
{  
   
   cout << endl;
   cout << " HAMSTR: Read input files ... ";
   mb->readMeshData(myid);
   sb->readSolutionInputs();

   // create instace for SOLVER class here
   // not done in HAMSTR constructor as variables
   // necessary for memory allocation within SOLVERs
   // constructor are not available prior to this point
   solv = new SOLVER[1];

   cout << " [done]\n";
}

// // ==================================================================
// // preprocess
// // ==================================================================

void HAMSTR::preprocessor(void)
{
   cout << " HAMSTR: Preprocess data  ... ";
   mb->preprocess();
   cout << " [done]\n";

   cout << " HAMSTR: Initialize flow  ... ";
   solv->initFlow();
   cout << " [done]\n";

   solv->basicScreenOutput();

}

// // ==================================================================
// // step solution for marching in time
// // ==================================================================

void HAMSTR::marchInTime(void)
{
   cout << " HAMSTR: Marching in time ... \n";
   solv->marchInTime();
   // cout << " [done]\n";  
}

// ##################################################################
// END OF FILE
// ##################################################################
