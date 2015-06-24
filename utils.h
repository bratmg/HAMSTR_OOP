// ##################################################################
//
// utils.h
//
// ##################################################################
#ifndef UTILS_H
#define UTILS_H



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <string>


// if 2D
#define NDIM  2
#define NFACE 4
#define NNODE 4
#define FNODE 2
#define NQ    4

// if 3D
#ifdef Dim3

#define NDIM  3
#define NFACE 6
#define NNODE 8
#define FNODE 4
#define NQ    5

#endif
using namespace std;
/*====================================================================*/
/*  Floating point definition                                         */
/*====================================================================*/
# define REAL double
/*====================================================================*/
/*  Base for indexing (0 or 1)                                        */
/*====================================================================*/
# define BASE 0
/*====================================================================*/
/*  Define arithmetic constants                                       */
/*====================================================================*/
#define ZERO               0.0e+00
#define ONE                1.0e+00
#define TWO                2.0e+00
#define THREE              3.0e+00
#define FOUR               4.0e+00
#define HALF               0.5e+00
#define THIRD              0.3333333333333333e+00
#define QUARTER            0.25e+00
#define FIFTH              0.2e+00
#define PI                 3.1415926535897932e+00
#define RAD2DEG            (180.0/PI)
#define DEG2RAD            (PI/180.0)
#define BIGVALUE           1.0e+15
#define BIGINT             2147483647
#define TOL                1.0e-10
#define NFRINGE            3
#define NVAR               6
/*==================================================================*/
/* inline debugging tools                                             */
/*==================================================================*/
# define tracei(x)  printf("\n#hamstr:\t"#x" = %d\n",x);
# define tracef(x)  printf("\n#hamstr:\t"#x" = %.16e\n",x);
# define traces(x)  printf("\n#hamstr:\t"#x"\n");
# define min(x,y)  (x) < (y) ? (x) : (y)
# define max(x,y)  (x) > (y) ? (x) : (y)
# define debug(x,y)  printf("#\nhamstr:\t"#x"=%d,"#y"=%d\n",x,y);
# define stdwrite(x) if (myid==0) printf("\n#hamstr:\t"#x"\n");
# define dstr(x) printf("\n#hamstr:\t"#x"\n");
# define ditch(x,y) {dstr(x); tracei(y); MPI_Abort(MPI_COMM_WORLD,ierr);}
/*====================================================================*/
/*  Numerical Tools                                                   */
/*====================================================================*/
#define Sign(a1,a2)\
        (((a2) < ZERO)? - fabs(a1): fabs(a1))
#define MAX(a1,a2)\
        (((a1) >= (a2))? (a1): (a2))
#define MIN(a1,a2)\
        (((a1) <= (a2))? (a1): (a2))
#define Abs(aa)\
        (((aa) >= 0)? (aa): -(aa))
#define Round(aa)\
        (int) ((fabs((aa) - floor(aa)) >= HALF)? ceil(aa): floor(aa))
#define swap(a,b) { a=a+b;b=a-b;a=a-b;}
/*====================================================================*/
/*  Boundary condition IDs                                            */
/*====================================================================*/
#define BC_FREESTREAM  -1
#define BC_WALL_INVISC -2
#define BC_WALL_INSC   -3
#define BC_EXTRAPOL    -4
#define BC_INTERDOMAIN -5

#endif
// ##################################################################
// END OF FILE
// ##################################################################
