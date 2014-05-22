#ifndef _MESHGEN2D_H_
#define _MESHGEN2D_H_

// Needed functions
double SignedDistance(point);
double GridSpacing(point);
void ParseArguments(int argc,char**argv,char*& outputdir);
void RunStartScript(char*& outputdir);
void MeshInputData(double&,point&,point&,int&,point*&,int&,int&,char*&);
void MeshCreateInitial(point,point,double,double,int,point*,int&,point*&,
		       double (*SD)(point),double (*GS)(point));
void MeshIterate(int maxiters,double geps,double ttol,double dptol,
		 double Fscale,double deltat,double deps,double h0,
		 int numfixpts,int& numpts,int& numtri,
		 point*& p,triangle*& t,double (*SD)(point),
		 double (*GS)(point));
void MeshCreateStructured(double h0, point xymin, point xymax, cart2d& meshinfo,
			  int& numpts, int& numtri, point*& p, triangle*& t);
void MeshCleanup(char*& GridType, double h0, double geps, 
		 int& numpts, int& numtri, point*& p, triangle*& t,		   
		 double*& area, double*& cdual, point xymin, point xymax,
		 double (*SD)(point));
void MeshBoundaryNodes(int numpts, int numtri, point p[], triangle t[],
		       double area[], int& numbnd, int*& bnd_node,
		       double (*SD)(point));
void MeshAddGhostCells(char*& GridType, 
		       cart2d& cart2dinfo,
		       int& numpts, 
		       int& numtri, 
		       int& numghost, 
		       int& num_ext_node,
		       point*& p, 
		       triangle*& t, 
		       double*& area, 
		       double*& cdual, 
		       int*& ghost_link, 
		       int*& proper_ghostcell,
		       int*& ext_node_link,
		       double (*SignedDistance)(point));
void MeshEdgeData(const double& h0,
		  const int& numpts, 
		  const int& numtri, 
		  const int& numghost, 
		  const point p[], 
		  const triangle t[], 
		  const double area[],
		  const int proper_ghostcell[],
		  int& numedges, 
		  int& numbndedges,
		  dTensor2*& edge,
		  iTensor2*& tedge,
		  iTensor2*& eelem,
		  iTensor2*& enode,
		  iTensor1*& bnd_edge,
		  double (*SignedDistance)(point));
void MeshOrientEdge(mesh& Mesh);
void ScreenOutput(const mesh& Mesh);
void MeshStore(point p[], 
	       triangle t[], 
	       int bnd_node[], 
	       int ghost_link[],
	       int ext_node_link[],
	       double area[], 
	       double cdual[], 
	       dTensor2*& edge, 
	       iTensor2*& tedge,
	       iTensor2*& tedge_orientation,
	       iTensor2*& eelem, 
	       iTensor2*& enode,
	       iTensor1*& bnd_edge,
	       mesh& Mesh);
void MeshPreProcess(char*& GridType, 
		    double& h0, 
		    point& xymin, 
		    point& xymax, 
		    int& numfixpts, 
		    point*& fixpt);
void MeshPostProcess1(char*&, double& h0, int& numpts, int& numtri, 
		      point*& p, triangle*& t, double*& area, 
		      double*& cdual, double (*SD)(point));
void MeshPostProcess2(char*&, mesh& Mesh, double (*SD)(point));
void MeshPostProcess3(const int,char*&,char*&, mesh& Mesh);

#endif
