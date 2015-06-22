// ##################################################################
//
// preprocess.cpp
//
// 
// ##################################################################
#include "hamstr.h"

#define BOUNDARYLIMIT 3.0
// ##################################################################

//
// find mesh properties
//
void MESHBLOCK::preprocess(void)
{   
   //
   // local variable declaration
   //
   int i,j,m,i1,i2,i3,i4;
   int leftCell,rightCell;
   int leftFaceIdx,rightFaceIdx;
   int f,f1,f2;
   int icell,iface;

   double x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,dvedge;
   double r_mid[3],s_vec[NDIM],r1,r2;
   double xa,xb,ya,yb,za,zb,dsnorm;
   // 
   // variable allocation
   //
   neig    = new int [NFACE*nCell];
   c2f     = new int [NFACE*nCell];
   c2chain = new int [NDIM*nCell];
   bFaces  = new int [nFace];
   cindx   = new int [nMaxChain+5];
   vol     = new double [nCell];
   normVec = new double* [nFace];
   refMtx  = new double** [nFace];   
   for (i = 0; i < nFace; i++) 
   {
      normVec[i] = new double [NDIM];
      refMtx[i]  = new double* [NDIM];
      for (j = 0; j < NDIM; j++)
         refMtx[i][j] = new double [NDIM];
   }

   //
   for (i = 0; i < nCell; i++) vol[i] = ZERO; 
   
// ==================================================================
// loop through all the faces
// ==================================================================

   m = 0;
   for (i = 0; i < nFace; i++)
   {
      leftCell     = faces[(FNODE+4)*i + FNODE     ];
      rightCell    = faces[(FNODE+4)*i + FNODE + 2 ];
      leftFaceIdx  = faces[(FNODE+4)*i + FNODE + 1 ];
      rightFaceIdx = faces[(FNODE+4)*i + FNODE + 3 ];

      // make a list of faces for a given cell
      c2f[(FNODE+2)*leftCell + leftFaceIdx] = i;
      if (rightCell > -1 ) 
         c2f[(FNODE+2)*rightCell + rightFaceIdx] = i;
      // 
      // collect neighbours
      //
      c2f[(FNODE+2)*leftCell + leftFaceIdx] = rightCell;
      if (rightCell > -1 ) 
         c2f[(FNODE+2)*rightCell + rightFaceIdx] = leftCell;
      //
      // indices of edge nodes
      //
      i1 = faces[(FNODE+4)*i  ];
      i2 = faces[(FNODE+4)*i+1];
      //
      // positions
      x1 = x[NDIM*i1  ];
      y1 = x[NDIM*i1+1];

      x2 = x[NDIM*i2  ];
      y2 = x[NDIM*i2+1];

      normVec[i][0] =  (y2-y1);
      normVec[i][1] = -(x2-x1);

      dvedge = HALF*(x1*y2-x2*y1);

      dsnorm = normVec[i][0]*normVec[i][0]
             + normVec[i][1]*normVec[i][1];
   
      dsnorm = 1./dsnorm;             

      // create reflection matrix
      refMtx[i][0][0] = 1. - 2.*normVec[i][0]*normVec[i][0]*dsnorm;
      refMtx[i][0][1] =    - 2.*normVec[i][0]*normVec[i][1]*dsnorm;

      refMtx[i][1][0] =    - 2.*normVec[i][1]*normVec[i][0]*dsnorm;
      refMtx[i][1][1] = 1. - 2.*normVec[i][1]*normVec[i][1]*dsnorm;

#ifdef Dim3 /* three-dimensional space */
      //
      // If in 3D (only the excess variables are included)
      //
      i3 = faces[(FNODE+4)*i+2];
      i4 = faces[(FNODE+4)*i+3];

      z1 = x[NDIM*i1+2];
      z2 = x[NDIM*i2+2];

      x3 = x[NDIM*i3  ];
      y3 = x[NDIM*i3+1];
      z3 = x[NDIM*i3+2];

      x4 = x[NDIM*i4  ];
      y4 = x[NDIM*i4+1];
      z4 = x[NDIM*i4+2];

      r_mid[0] = QUARTER*(x1+x2+x3+x4);
      r_mid[1] = QUARTER*(y1+y2+y3+y4);
      r_mid[2] = QUARTER*(z1+z2+z3+z4);

      xa = x3-x1;
      xb = x2-x4;
      ya = y3-y1;
      yb = y2-y4;
      za = z3-z1;
      zb = z2-z4;

      normVec[i][0] = HALF*(za*yb-ya*zb);
      normVec[i][1] = HALF*(xa*zb-za*xb);
      normVec[i][2] = HALF*(ya*xb-xa*yb);
     
      dvedge = THIRD*(r_mid[0]*normVec[i][0]
                     +r_mid[1]*normVec[i][1]
                     +r_mid[2]*normVec[i][2]);

      dsnorm = normVec[i][0]*normVec[i][0]
             + normVec[i][1]*normVec[i][1]
             + normVec[i][2]*normVec[i][2];
   
      dsnorm = 1./dsnorm;             

      // create reflection matrix
      refMtx[i][0][0] = 1. - 2.*normVec[i][0]*normVec[i][0]*dsnorm;
      refMtx[i][0][1] =    - 2.*normVec[i][0]*normVec[i][1]*dsnorm;
      refMtx[i][0][2] =    - 2.*normVec[i][0]*normVec[i][2]*dsnorm;

      refMtx[i][1][0] =    - 2.*normVec[i][1]*normVec[i][0]*dsnorm;
      refMtx[i][1][1] = 1. - 2.*normVec[i][1]*normVec[i][1]*dsnorm;
      refMtx[i][1][2] =    - 2.*normVec[i][1]*normVec[i][2]*dsnorm;

      refMtx[i][2][0] =    - 2.*normVec[i][2]*normVec[i][0]*dsnorm;
      refMtx[i][2][1] =    - 2.*normVec[i][2]*normVec[i][1]*dsnorm;
      refMtx[i][2][2] = 1. - 2.*normVec[i][2]*normVec[i][2]*dsnorm;

#endif      
      //
      // compute cell volumes
      //
      vol[leftCell] += dvedge;

      if(rightCell > -1) vol[rightCell] -= dvedge;
      //
      // Identify boundary face and implement loop boundary conditions
      //
      if (faces[(2+NFACE)*i + FNODE + 2] == -1)
      {

         i1 = faces[(2+NFACE)*i    ];
         i2 = faces[(2+NFACE)*i + 1];

         r1 = sqrt(x[2*i1]*x[2*i1] + x[2*i1+1]*x[2*i1+1]);
         r2 = sqrt(x[2*i2]*x[2*i2] + x[2*i2+1]*x[2*i2+1]);

         if (0.5*(r1+r2) < BOUNDARYLIMIT)
         {
            bFaces[m] = i;
            m++;
            faces[(2+NFACE)*i + FNODE + 2] = -(sol_visc + 2);
            icell                          = faces[(2+NFACE)*i+FNODE  ];
            iface                          = faces[(2+NFACE)*i+FNODE+1];
#ifdef Dim3       
            neig[6*icell+iface]            = -(sol_visc + 2);
#else
            neig[3*icell+iface]            = -(sol_visc + 2);
#endif
         }
      }

   } // i loop

// ==================================================================
// compute min and max volumes
// ==================================================================

   volmin = BIGVALUE;
   volmax = ZERO;

   // evaluate the minimum and maximum cell volume
   for (i = 0; i < nCell; i++)
   {
      volmin = (volmin < vol[i]) ? volmin : vol[i]; 
      volmax = (volmax > vol[i]) ? volmax : vol[i];
   }
   //  
   // Find cell to chain connectivity
   //
   for (i = 0; i < nCell; i++)
   {
      c2chain[NDIM*i  ] = -1;
      c2chain[NDIM*i+1] = -1;
#ifdef Dim3      
      c2chain[NDIM*i+2] = -1;
#endif
   }
   //
   // loop through all the loops
   //
   for (i = 0; i < nChain; i++)
   {
      f1 = faceStartPerChain[i  ];
      f2 = faceStartPerChain[i+1];

      // loop through the chain indices
      for(f = f1; f < f2-1; f++)
      {
        iface    = chainConn[f]; // index of the face/edge
        leftCell = faces[(FNODE+4)*iface + FNODE]; // corresponding left cell 

        if (c2chain[NDIM*leftCell + (NDIM-2)] > -1)
           c2chain[NDIM*leftCell + (NDIM-1)] = i;
        else
           c2chain[NDIM*leftCell + (NDIM-2)] = i;
      } // f loop
   } // i loop


}
// ##################################################################
// END OF FILE
// ##################################################################

