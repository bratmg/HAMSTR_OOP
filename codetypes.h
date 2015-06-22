// ##################################################################
//
// codetypes.h
//
// ##################################################################
#ifndef CODETYPES_H
#define CODETYPES_H
// ##################################################################
#include "utils.h"

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
