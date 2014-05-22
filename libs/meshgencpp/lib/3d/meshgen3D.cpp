#include "meshdefs.h"
#include <cstdlib>
#include <cstdio>
#include <time.h>

int meshgen3D(int argc, char**argv)
{  
  // Output title information
  printf("\n");
  printf("   -------------------------------------------------   \n");
  printf("   |  MeshGenC++: An Unstructured Grid Generator   |   \n");
  printf("   |    C++ Version of Per-Olof Persson and        |   \n");
  printf("   |    Gilbert Strang's DistMesh algorithm for    |   \n");
  printf("   |    unstructured grid generation               |   \n");
  printf("   |                                               |   \n");
  printf("   |                   Written by                  |   \n");
  printf("   |    James A. Rossmanith and the DoGPack Team   |   \n");
  printf("   -------------------------------------------------   \n");
  printf("\n"); 

  // Needed functions
  double SignedDistance(point);
  double GridSpacing(point);
  void ParseArguments(int argc,
		      char**argv,
		      char* outputdir);
  void MeshInputData(double&,
		     point&,
		     point&,
		     int&,
		     point*&,
		     int&,
		     char*&); 
  void MeshCreateInitial(point xyzmin,
			 point xyzmax,
			 double h0,
			 double geps,
			 int numfixpts, 
			 point* fixpt,
			 int& numpts,
			 point*& p,
			 double (*SignedDistance)(point),
			 double (*GridSpacing)(point));
  void MeshIterate(int maxiters,
		   double geps,
		   double ttol,
		   double dptol,
  		   double Fscale,
		   double deltat,
		   double deps,
		   double h0,
  		   int numfixpts,
		   int& numpts,
		   int& numtet,
  		   point*& p,
		   tetra*& t,
		   double (*SD)(point),
  		   double (*GS)(point)); 
  void MeshCleanup(const char* GridType, 
		   double h0, 
		   double geps, 
		   int& numpts, 
		   int& numtet, 
		   point*& p, 
		   tetra*& t, 
		   double*& volume, 
		   double*& cdual,
		   double (*SignedDistance)(point));
  void MeshCreateStructured(double h0,
			    point xyzmin,
			    point xyzmax,
			    int& numpts, 
			    int& numtet,
			    point*& p,
			    tetra*& t);
  void MeshBoundaryNodes(int numpts,
			 int numtet,
			 point p[],
			 tetra t[],
			 double volume[],
			 int& numbnd,
			 int*& bnd_node,
			 double (*SD)(point));
  void MeshOutput(int numtet,
		  int numtet_phys,
		  int numghost,
		  int numpts,
  		  int numpts_phys,
		  int numbnd,
		  int numedges,
		  point p[], 
  		  tetra t[],
		  int bnd_node[],
		  double area[], 
		  double cdual[],
		  char* outputdir);

  // Get data from input file
  double h0;
  point xyzmin;
  point xyzmax;
  int numfixpts;
  int maxiters;
  point *fixpt;
  fixpt = new point;
  char* GridType = new char[12];
  int numpts,numtet;
  point* p;  
  tetra* t;
  
  // Parse arguments -- sets directory to which output will be sent,
  //                    the default is "output"  
  char outputdir[] = "mesh_output";
  ParseArguments(argc,argv,outputdir);

  // Get data from input file
  MeshInputData(h0,xyzmin,xyzmax,numfixpts,fixpt,maxiters,GridType);
  
  const double geps   = 0.1 * h0;
  const double ttol   = 0.1;
  const double dptol  = 0.001;
  const double deltat = 0.1;
  const double Fscale = 1.1;
  const double deps   = 1.4901161119384766e-8 * h0;
  
  if (GridType[0]=='U')
    {
      // Create other needed constants
      srand((unsigned int)time(0));	// seed for random number generation      
      
      // Create initial mesh      
      MeshCreateInitial(xyzmin,xyzmax,h0,geps,numfixpts,fixpt,numpts,p,
			&SignedDistance,&GridSpacing);
  
      // Move points and remesh as necessary to create equilibrium mesh      
      MeshIterate(maxiters,geps,ttol,dptol,Fscale,deltat,deps,h0,
      		  numfixpts,numpts,numtet,p,t,&SignedDistance,&GridSpacing);
    }
  else if (GridType[0]=='S')
    {
      // Create structured mesh
      MeshCreateStructured(h0,xyzmin,xyzmax,numpts,numtet,p,t);
    }
  else
    {
      printf("\n");
      printf("  ERROR in meshgen3D.cpp:  GridType = %s", GridType);
      printf(" is not supported \n\n");
      exit(1);
    }

  // Cleanup mesh (remove duplicate points, correctly orient triangles, and
  //               compute triangle areas)
  double* volume;
  double* cdual;
  volume = new double[numtet];
  cdual  = new double[numpts];
  MeshCleanup(GridType,h0,geps,numpts,numtet,p,t,volume,cdual,&SignedDistance);

  // Find all nodes that live on the boundary
  int* bnd_node;
  int numbnd;
  MeshBoundaryNodes(numpts,numtet,p,t,volume,numbnd,
		    bnd_node,&SignedDistance);

  // Output grid
  int numtet_phys = numtet;
  int numghost = 0;
  int numpts_phys = numpts;
  int numfaces = 0;
  MeshOutput(numtet,numtet_phys,numghost,numpts,numpts_phys,numbnd,
	     numfaces,p,t,bnd_node,volume,cdual,outputdir);
  
  delete fixpt;
  delete GridType;
  delete p;  
  delete t;
  delete volume;
  delete cdual;
  delete bnd_node;
      
  return 0;
}
