// ##################################################################
//
// io.cpp
//
// contains input output routines 
// ##################################################################
#include "meshBlock.h"
#include "solnBlock.h"
#include "hamstr.h"
#include "solver.h"
// ##################################################################

void SOLNBLOCK::readSolutionInputs(void)
{

   char c;
   FILE * fptr;

   if ( (fptr = fopen("input.hamstr","r")) == NULL)
   {
      printf("Error: Input file 'input.hamstr' missing.\n");
      exit(1);
   }
   else
   {
      fscanf(fptr,"Mach=%lf\n",&(Mach));
      fscanf(fptr,"alpha=%lf\n",&(alpha));
      fscanf(fptr,"beta=%lf\n",&(beta));
      fscanf(fptr,"rey=%lf\n",&(rey));
      fscanf(fptr,"scheme=%s\n",scheme);
      fscanf(fptr,"time integration=%s\n",timeInteg);
      fscanf(fptr,"order=%d\n",&(order));
      fscanf(fptr,"timeacc=%d\n",&(timeacc));
      fscanf(fptr,"nsteps=%d\n",&(nsteps));
      fscanf(fptr,"nwrite=%d\n",&(nwrite));
      fscanf(fptr,"dt=%lf\n",&(dt));
      fscanf(fptr,"CFL=%lf\n",&(CFL));
      fscanf(fptr,"msweep=%d\n",&(msweep));
      fscanf(fptr,"visc=%d\n",&(visc));
      fscanf(fptr,"testcase=%d\n",&(test));
      fscanf(fptr,"output=%d\n",&(outform));
      fclose(fptr);
   }

   if(test!=0)
   {
      cout << "\nERROR: Code not implemented for all test cases\n";
      exit(1);
   }

}

// ##################################################################
//
// MESHBLOCK::READDATA
//
// ##################################################################
void MESHBLOCK::readMeshData(int myid)
{
   
   
   //
   // local variables
   //
   int i,j;
   char fCoord[40], fConn[40],fOface[40],fQloop[40],fIqloop[40],fNcolor[40];
   char fCommu[80],fDomain[80],fStrand[50];
   FILE * fptr;
   //
   // open file based on domain ID (proc ID)
   //
   if(myid < 10)
   {
      sprintf(fCoord, "./input/domain00%i/coord.dat",myid);
      sprintf(fConn,  "./input/domain00%i/conn.dat",myid);
      sprintf(fNcolor,"./input/domain00%i/ncolors.dat",myid);
      sprintf(fIqloop,"./input/domain00%i/iqloops.dat",myid);
      sprintf(fQloop, "./input/domain00%i/qloops.dat",myid);

      sprintf(fCommu, "./input/domain00%i/quadCommu00%i.dat",myid,myid);
      sprintf(fDomain, "./input/domain00%i/subdomain00%i.dat",myid,myid);

#ifdef Dim3
      sprintf(fOface, "./input/domain00%i/ofaces.dat",myid);
      sprintf(fStrand, "./input/domain00%i/nstrands.dat",myid);
#else
      sprintf(fOface, "./input/domain00%i/qedges.dat",myid);
#endif      

   }
   else if(myid < 100)
   {
      sprintf(fCoord, "./input/domain0%i/coord.dat",myid);
      sprintf(fConn,  "./input/domain0%i/conn.dat",myid);
      sprintf(fNcolor,"./input/domain0%i/ncolors.dat",myid);
      sprintf(fIqloop,"./input/domain0%i/iqloops.dat",myid);
      sprintf(fQloop, "./input/domain0%i/qloops.dat",myid);

      sprintf(fCommu, "./input/domain0%i/quadCommu00%i.dat",myid,myid);
      sprintf(fDomain, "./input/domain0%i/subdomain00%i.dat",myid,myid);

#ifdef Dim3
      sprintf(fOface, "./input/domain0%i/ofaces.dat",myid);
      sprintf(fStrand, "./input/domain0%i/nstrands.dat",myid);
#else
      sprintf(fOface, "./input/domain0%i/qedges.dat",myid);
#endif 

   }
   else if (myid < 1000)
   {
      sprintf(fCoord, "./input/domain%i/coord.dat",myid);
      sprintf(fConn,  "./input/domain%i/conn.dat",myid);
      sprintf(fNcolor,"./input/domain%i/ncolors.dat",myid);
      sprintf(fIqloop,"./input/domain%i/iqloops.dat",myid);
      sprintf(fQloop, "./input/domain%i/qloops.dat",myid);

      sprintf(fCommu, "./input/domain%i/quadCommu00%i.dat",myid,myid);
      sprintf(fDomain, "./input/domain%i/subdomain00%i.dat",myid,myid);

#ifdef Dim3
      sprintf(fOface, "./input/domain%i/ofaces.dat",myid);
      sprintf(fStrand, "./input/domain%i/nstrands.dat",myid);
#else
      sprintf(fOface, "./input/domain%i/qedges.dat",myid);
#endif 

   } 
   else
   {
      printf("\nError. Code needs to be modified for beyond 1,000 procs\n");
   }
      
// ==================================================================
// read coordinate file
// ==================================================================

   if( (fptr = fopen(fCoord,"r")) == NULL)
   {
      printf("\n File missing - %s. Stopping\n", fCoord);exit(1);
   }
   else
   {
      fscanf(fptr, "%d", &(nNode));
      x = new double [NDIM*nNode];

#ifdef Dim3
      for (i = 0; i < nNode; i++)
         fscanf(fptr,"%lf %lf %lf\n",
         &(x[NDIM*i]),&(x[NDIM*i+1]),&(x[NDIM*i+2]));
#else
      double temp;
      for (i = 0; i < nNode; i++)
         fscanf(fptr,"%lf %lf %lf\n",
            &(x[NDIM*i]),&(x[NDIM*i+1]),&temp);

#endif      
   }
   fclose(fptr);

// ==================================================================
// read connectivity file
// ==================================================================

   if((fptr=fopen(fConn,"r")) == NULL)
   {
      printf("\n File missing - %s. Stopping\n", fConn);exit(1);  
   }
   else
   {
      fscanf(fptr,"%d",&(nCell));

      conn = new int [NNODE*nCell];

#ifdef Dim3
      for(i = 0; i < nCell; i++)    
         fscanf(fptr,"%d %d %d %d %d %d %d %d\n",
         &(conn[NNODE*i  ]),&(conn[NNODE*i+1]),
         &(conn[NNODE*i+2]),&(conn[NNODE*i+3]),
         &(conn[NNODE*i+4]),&(conn[NNODE*i+5]),
         &(conn[NNODE*i+6]),&(conn[NNODE*i+7]));
#else
      for(i = 0; i < nCell; i++)   
      {
         fscanf(fptr,"%d %d %d %d\n",
         &(conn[NNODE*i  ]),&(conn[NNODE*i+1]),
         &(conn[NNODE*i+2]),&(conn[NNODE*i+3]));
      }
#endif      

      fclose(fptr);    
   } // if file exists

      
// ==================================================================
// Read colours
// ==================================================================

   if ( (fptr = fopen(fNcolor,"r")) == NULL)
   {
      printf("\n File missing - %s. Stopping\n", fNcolor);exit(1);    
   }
   else
   {
      fscanf(fptr,"%d",&(nColor));
      chainsPerColor = new int [nColor];
      for(i = 0; i < nColor; i++)
         fscanf(fptr,"%d\n",&(chainsPerColor[i]));

      fclose(fptr);
   }

// ==================================================================
// read chain information (chain size)
// ==================================================================

   if ((fptr = fopen(fIqloop,"r")) == NULL)
   {
      printf("\n File missing - %s. Stopping\n", fIqloop);exit(1);    
   }
   else
   {
      fscanf(fptr,"%d\n",&(nChain));

      faceStartPerChain = new int [nChain];
      nMaxChain         = 0;

      j = 0;
      for(i = 0; i < nChain; i++)
      {      
         fscanf(fptr,"%d\n",&(faceStartPerChain[j]));
         //
         // this routine is for delete the duplicated one
         // at every layers (maybe rectify this issue in meshgen)
         //
         if (i > 0) 
         {
            if(faceStartPerChain[j]==faceStartPerChain[j-1]) j--;
            
            nMaxChain = max(nMaxChain,
               faceStartPerChain[j]-faceStartPerChain[j-1]);
         }
         j++;
      }
      nChain = j-1;
      fclose(fptr);
   }

// ==================================================================
// read chain information (chain connectivity)
// ==================================================================

   if ((fptr = fopen(fQloop,"r")) == NULL)
   {
      printf("\n File missing - %s. Stopping\n", fQloop);exit(1);    
   }
   else
   {
      fscanf(fptr,"%d",&(nChainFace));

      chainConn = new int [nChainFace];

      for(i = 0; i < nChainFace; i++)
         fscanf(fptr,"%d",&(chainConn[i]));

      fclose(fptr);

   }  

#ifdef Dim3   

// ==================================================================
// read faces file (THREE DIMENSIONS)
// ==================================================================

   if ( (fptr = fopen(fOface,"r")) == NULL)
   {
      printf("\n File missing - %s. Stopping\n", fOface);exit(1);    
   }
   else
   {
      fscanf(fptr,"%d",&(nFace));
      faces = new int [8*nFace];

      for(i = 0; i < nFace; i++)
      {
         fscanf(fptr, "%d %d %d %d %d %d %d %d\n",
            &(faces[8*i  ]),  //n1
            &(faces[8*i+1]),  //n2
            &(faces[8*i+2]),  //n3
            &(faces[8*i+3]),  //n4
            &(faces[8*i+4]),  //c1
            &(faces[8*i+6]),  //c2
            &(faces[8*i+5]),  //e1
            &(faces[8*i+7])); //e2

         // boundary cell should be -1 at upperlayers
         // this is due to just add the total number of faces at each layer
         // to the previous layer face index
         if(faces[8*i+7] == BC_FREESTREAM) 
            faces[8*i+6] = BC_FREESTREAM;
         //
         // swap if left cell = -1  
         //
         if (faces[8*i+4] == BC_FREESTREAM || 
             faces[8*i+4] == BC_WALL_INVIS || 
             faces[8*i+4] == BC_INTERDOMAIN) 
         {
            swap(faces[8*i  ] , faces[8*i+1]);//n2,n4 
            swap(faces[8*i+2] , faces[8*i+3]);//n2,n4 
            swap(faces[8*i+4] , faces[8*i+6]);//cell index
            swap(faces[8*i+5] , faces[8*i+7]);//element index  
         }
      }
      fclose(fptr);
   }

// ==================================================================
// read chain information (chain connectivity)
// ==================================================================

   if ((fptr = fopen(fStrand,"r")) == NULL)
   {
      printf("\n File missing - %s. Stopping\n", fStrand);exit(1);    
   }
   else
   {
      fscanf(fptr,"%d",&(nStrand));
      fclose(fptr);
   }  

#else

// ==================================================================
// read faces file (TWO DIMENSIONS)
// ==================================================================

   if ( (fptr = fopen(fOface,"r")) == NULL)
   {
      printf("\n File missing - %s. Stopping\n", fOface);exit(1);    
   }
   else
   {
      fscanf(fptr,"%d",&(nFace));
      faces = new int [6*nFace];

      for(i = 0; i < nFace; i++)
      {
         fscanf(fptr, "%d %d %d %d %d %d\n",
            &(faces[6*i  ]),  //n1
            &(faces[6*i+1]),  //n2
            &(faces[6*i+2]),  //c1
            &(faces[6*i+4]),  //c2
            &(faces[6*i+3]),  //e1
            &(faces[6*i+5])); //e2

         //
         // swap if left cell = -1  
         //
         if (faces[6*i+2] == BC_FREESTREAM) 
         {
            swap(faces[6*i  ] , faces[6*i+1]);//n2,n4 
            swap(faces[6*i+2] , faces[6*i+4]);//n2,n4 
            swap(faces[6*i+3] , faces[6*i+5]);//cell index
         }
      }
      fclose(fptr);
   }
#endif
      
   
}

// ##################################################################
//
// HAMSTR::WELCOME
//
// ##################################################################
void HAMSTR::welcome(void)
{
   cout << endl;
   cout << "######################################################################\n";
   cout << "#                                                                    #\n";
   cout << "#                            HAMSTR                                  #\n";
   cout << "#                                                                    #\n";
   cout << "#              HAMILTONIAN-STRAND FINITE VOLUME SOLVER               #\n";
   cout << "#                                                                    #\n";
   cout << "######################################################################\n";
}


// ##################################################################
//
// HAMSTR::THANKS
//
// ##################################################################
void HAMSTR::thanks(void)
{
   cout << endl;
   cout << "######################################################################\n";
   cout << "#                                                                    #\n";
   cout << "#                       END OF SIMULATION                            #\n";
   cout << "#                                                                    #\n";
   cout << "######################################################################\n";
}

// ##################################################################
//
// HAMSTR::BASICSCREENOUTPUT
//
// ##################################################################
void SOLVER::basicScreenOutput(void)
{
   printf("======================================================================\n");
   printf("                       GRID RELATED OUTPUTS                           \n");
   printf("======================================================================\n");
   printf("\n");
   printf(" Number of nodes            %d\n",mb->nNode);
   printf(" Number of cells            %d\n",mb->nCell);
   printf(" Number of faces            %d\n",mb->nFace);
   printf(" Number of colours          %d\n",mb->nColor);
   printf(" Number of chains           %d\n",mb->nChain);
   printf(" Number of chain faces      %d\n",mb->nChainFace);
   printf(" Maximum length of chain    %d\n",mb->nMaxChain);
   printf(" Minimum cell volume        %0.8e\n",mb->volmin);
   printf(" Maximum cell volume        %0.8e\n",mb->volmax);
   printf("\n");
   printf("======================================================================\n");
   printf("                  SOLVER AND FLOW RELATED OUTPUTS                     \n");
   printf("======================================================================\n");
   printf("\n");
   printf(" Mach number                %0.8e\n",sb->Mach);
   printf(" Alpha                      %0.8e\n",sb->alpha);
   printf(" Beta                       %0.8e\n",sb->beta);
   printf(" U infinity                 %0.8e\n",sb->uinf);
   printf(" V infinity                 %0.8e\n",sb->vinf);
   printf(" Beta                       %0.8e\n",sb->beta);
   printf(" CFL number                 %0.8e\n",sb->CFL);
   printf(" Number of steps            %d\n",sb->nsteps);
   printf(" Time marching scheme       %s\n",sb->scheme);
   printf(" Time integration scheme    %s\n",sb->timeInteg);
   printf(" Reconstruction accuracy    %d\n",sb->order);

   printf("\n");
   printf("======================================================================\n");
}
// ##################################################################
// END OF FILE
// ##################################################################

