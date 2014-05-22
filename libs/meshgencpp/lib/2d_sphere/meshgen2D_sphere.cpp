#include "meshdefs.h"
#include "mesh.h"

int meshgen2D_sphere(int argc, char**argv)
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
  void ParseArguments(int argc,
		      char**argv,
		      char* outputdir);
  void MeshInputData(double&, int&);
  void MeshCreateIcosa(const double radius, int& numpts, int& numtri,
		       point*& p, triangle*& t);
  void MeshSubdivide(const int level, const double radius, 
		     int& numpts, int& numtri, point*& p, triangle*& t);
  void MeshSphereCoords(const int, const point* p, point*& psphere);
  void MeshEdgeData(const int& numpts, const int& numtri, const point p[],
		    const point psphere[], const triangle t[], const double area[],
		    int& numedges, dTensor2*& edge, iTensor2*& tedge,
		    iTensor2*& eelem);
  void MeshOrientEdge(mesh& Mesh);
  void MeshComputeAreas(const int numtri, const int numpts, 
			const point p[], triangle t[], 
			double*& area, double*& cdual);
  void ScreenOutput(const mesh& Mesh);
  void MeshStore(point p[], triangle t[], int bnd_node[], 
		 int ghost_link[], double area[], double cdual[], 
		 dTensor2*& edge, iTensor2*& tedge, 
		 iTensor2*& eelem, mesh& Mesh);
  
  // Parse arguments -- sets directory to which output will be sent,
  //                    the default is "output"  
  char outputdir[] = "mesh_output";
  ParseArguments(argc,argv,outputdir);

  // Get data from input file
  double radius;
  int level;
  MeshInputData(radius,level);

  // Create icosahedral mesh
  point* p;
  triangle* t;
  int numpts,numtri;
  MeshCreateIcosa(radius,numpts,numtri,p,t);

  // Subdivide mesh "level" number of times
  if (level>0)
    {
      MeshSubdivide(level,radius,numpts,numtri,p,t);
    }

  // Compute element areas and dual-element areas
  double* area = new double[numtri];
  double* cdual = new double[numpts];
  MeshComputeAreas(numtri,numpts,p,t,area,cdual);

  // Store spherical coordinates of p in psphere
  point* psphere = new point[numpts];
  MeshSphereCoords(numpts,p,psphere);

  // Create information about edges
  int numedges = 1;
  dTensor2* edge;
  iTensor2* tedge;
  iTensor2* eelem;
  int* bnd_node = new int[1];
  int* ghost_link = new int[1];
  edge = new dTensor2(numedges,4);
  tedge = new iTensor2(numtri,3);
  eelem = new iTensor2(numedges,2);
  MeshEdgeData(numpts,numtri,p,psphere,t,area,numedges,
               edge,tedge,eelem);

  // Store all mesh info 
  mesh Mesh(numtri,numtri,numpts,numpts,0,numedges,0);
  MeshStore(psphere,t,bnd_node,ghost_link,area,cdual,edge,tedge,eelem,Mesh);

  // Adjust edge information so that the element to the "left"
  // of the edge is always "upwind" of the unit normal to the
  // edge and the element to the "right" of the edge is always
  // "downwind" of the unit normal to the edge
  MeshOrientEdge(Mesh);

  // Output grid
  ScreenOutput(Mesh);
  Mesh.OutputMesh(outputdir);
  
  return 0;
}
