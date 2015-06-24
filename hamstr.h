// ##################################################################
//
// hamstr.h
//
// Header files for HAMSTR
// ##################################################################
#ifndef HAMSTR_H
#define HAMSTR_H

#include "meshBlock.h"
#include "solnBlock.h"
#include "solver.h"
#include "codetypes.h"
//#include "commBlock.h"

// ##################################################################
class HAMSTR
{
   private:
      int        nGrid;
      int        procID;
      MESHBLOCK *mb;
      SOLNBLOCK *sb;
      SOLVER    *solv;
      //COMMBLOCK *cb;

   public:
      // basic constructor
      HAMSTR(const int , const int);

      // basic destructor
      ~HAMSTR(){};// {delete mb; delete sb;};

      // read inputs
      void readInputs(int);
    
      // reorganize data
      void preprocessor(void);

      // communication routines

      // march one step
      void marchInTime(void);

      // welcome and thanks
      void welcome(void);
      void thanks(void);


};
#endif
// ##################################################################
// END OF FILE
// ##################################################################
