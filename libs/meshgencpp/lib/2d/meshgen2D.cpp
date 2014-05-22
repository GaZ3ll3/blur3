#include "meshdefs.h"
#include "mesh.h"
#include <time.h>
#include "meshgen2D.h"

int meshgen2D(int argc, char**argv)
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
    
  // Parse arguments -- sets directory to which output will be sent,
  //                    the default is "mesh_output"  
  char* outputdir = new char[12];
  outputdir[0] ='m';
  outputdir[1] ='e';
  outputdir[2] ='s';
  outputdir[3] ='h';
  outputdir[4] ='_';
  outputdir[5] ='o';
  outputdir[6] ='u';
  outputdir[7] ='t';
  outputdir[8] ='p';
  outputdir[9] ='u';
  outputdir[10]='t';
  outputdir[11]='\0';
  
  ParseArguments(argc,argv,outputdir);
  RunStartScript(outputdir);  
  
  // Get data from input file
  double h0;
  point xymin;
  point xymax;
  int numfixpts;
  int maxiters;
  int sub_factor;
  point *fixpt;
  char* GridType = new char[12];
  fixpt = new point;
  int numpts,numtri;
  point* p;
  triangle* t;
  cart2d cart2dinfo;
  MeshInputData(h0,xymin,xymax,numfixpts,fixpt,maxiters,sub_factor,GridType);
  
  // Pre-process input data if necessary
  MeshPreProcess(GridType,h0,xymin,xymax,numfixpts,fixpt);    

  // Create needed constants
  srand((unsigned)time(0));	// seed for random number generation
  const double geps   = 0.001 * h0;
  const double ttol   = 0.1;
  const double dptol  = 0.001;
  const double Fscale = 1.2;
  const double deltat = 0.2;
  const double deps   = 1.49012 * pow(10.0, -8) * h0;

  if (GridType[0]=='U')
    {
      // Create initial mesh      
      MeshCreateInitial(xymin,xymax,h0,geps,numfixpts,fixpt,numpts,p,
			&SignedDistance,&GridSpacing);
      
      // Move points and remesh as necessary to create equilibrium mesh      
      MeshIterate(maxiters,geps,ttol,dptol,Fscale,deltat,deps,h0,
		  numfixpts,numpts,numtri,p,t,&SignedDistance,&GridSpacing);
    }
  else if (GridType[0]=='S')
    {
      // Create structured mesh      
      MeshCreateStructured(h0,xymin,xymax,cart2dinfo,numpts,numtri,p,t);
    }
  else
    {
      printf("\n");
      printf("  ERROR in meshgen2D.cpp:  GridType = %s \n",GridType);
      printf(" is not supported \n");
      printf("\n");
      exit(1);
    }

  // Cleanup mesh (remove duplicate points, correctly orient triangles, and
  //               compute triangle areas)
  double* area;
  double* cdual;
  area = new double[numtri];
  cdual = new double[numpts];
  MeshCleanup(GridType,h0,geps,numpts,numtri,
	      p,t,area,cdual,xymin,xymax,
	      &SignedDistance);

  // Post-process mesh if necessary
  MeshPostProcess1(GridType,h0,numpts,numtri,p,t,area,cdual,
		   &SignedDistance);
  
  // Find all nodes that live on the boundary
  int* bnd_node;
  int numbnd;
  MeshBoundaryNodes(numpts,numtri,p,t,area,numbnd,
		    bnd_node,&SignedDistance);
  
  // Add ghost cells
  int numpts_phys = numpts;
  int numtri_phys = numtri;
  int numghost;  
  int num_ext_node;
  int* ghost_link;
  int* proper_ghostcell;
  int* ext_node_link;
  MeshAddGhostCells(GridType,cart2dinfo,numpts,numtri,numghost,num_ext_node,p,t,
		    area,cdual,ghost_link,proper_ghostcell,ext_node_link,
		    &SignedDistance);

  // Create information about edges
  int numedges;
  int numbndedges;
  dTensor2* edge;
  iTensor2* tedge;
  iTensor2* eelem;
  iTensor2* enode;
  iTensor2* tedge_orientation;
  iTensor1* bnd_edge;  
  tedge_orientation = new iTensor2(numtri,3);
  for (int i=1; i<=numtri; i++)
    for (int j=1; j<=3; j++)
      {  tedge_orientation->set(i,j, 0 );  }
  numedges = 1;
  numbndedges = 1;
  
  MeshEdgeData(h0,numpts,numtri,numghost,p,t,area,proper_ghostcell,
	       numedges,numbndedges,edge,tedge,eelem,enode,bnd_edge,
	       &SignedDistance);
  
  // Store all mesh info
  mesh Mesh(numtri,numtri_phys,numpts,numpts_phys,numbnd,numedges,numbndedges);
  MeshStore(p,t,bnd_node,ghost_link,ext_node_link,area,cdual,edge,tedge,
            tedge_orientation,eelem,enode,bnd_edge,Mesh);  
  delete fixpt;
  delete p;
  delete t;
  delete area;
  delete cdual;
  delete bnd_node;
  delete ghost_link;
  delete proper_ghostcell;
  delete ext_node_link;
  delete edge;
  delete tedge;
  delete eelem;
  delete enode;
  delete tedge_orientation;
  delete bnd_edge;
  
  // Adjust edge information so that the element to the "left"
  // of the edge is always "upwind" of the unit normal to the
  // edge and the element to the "right" of the edge is always
  // "downwind" of the unit normal to the edge
  MeshOrientEdge(Mesh);
  
  // Post-process mesh if necessary
  MeshPostProcess2(GridType,Mesh,&SignedDistance);
  
  // Compute number elements attached to each node
  Mesh.Compute_NumElemsPerNode();
  
  // Compute Jacobian matrix for each element
  Mesh.ComputeJacobian();
  
  // Find all elements that are adjacent to current element
  Mesh.SetAdjacency();

  // Find all elements that are adjacent to current element
  Mesh.ComputeBoundingBox();
  
  // Create a sub-mesh if necessary
  if (sub_factor>1)
    {  
      Mesh.CreateSubMesh(sub_factor,
			 xymin.x,xymax.x,
			 xymin.y,xymax.y,
			 deps,&SignedDistance); 
      
    }
  
  // Post-process mesh if necessary
  MeshPostProcess3(sub_factor,GridType,outputdir,Mesh);
  
  // Output grid
  ScreenOutput(Mesh);
  Mesh.OutputMesh(outputdir);
  mesh AnotherMesh(numtri,numtri_phys,numpts,numpts_phys,numbnd,numedges,numbndedges);
  AnotherMesh.InputMesh(outputdir);

  bool is_same = Mesh.Compare(AnotherMesh,true);
  
  return is_same;
}
