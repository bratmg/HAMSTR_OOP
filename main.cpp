// ##################################################################
//
// Driver file for HAMSTR
//
// ##################################################################
#include "hamstr.h"


//using namespace Reconstruction;

// ##################################################################
int main(int argc, char **argv)
{
   //
   // Initialize MPI
   //
   int nprocs,myid,ierr;
#ifdef useMPI   
   ierr = MPI_Init(&argc, &argv);
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#else
   nprocs = 1;
   myid   = 0;
#endif   
   //
   // create member for HAMSTR class
   //
   HAMSTR myHam(nprocs,myid);
   myHam.welcome();
   //
   // read inputs and data
   // 
   myHam.readInputs(myid);
   //
   // preprocess and reorganize data
   //
   myHam.preprocessor();
   //
   // march in time
   //
   myHam.marchInTime();
   //
   // thanks and finalize
   //
   myHam.thanks();

#ifdef useMPI   
   MPI_Finalize();
#endif
   
   return 0;
}
// ##################################################################
// END OF FILE
// ##################################################################
