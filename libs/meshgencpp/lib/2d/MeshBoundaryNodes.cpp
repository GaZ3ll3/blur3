#include "meshdefs.h"

//
// Find all nodes that lie directly on the boundary
//
void MeshBoundaryNodes(int numpts, 
		       int numtri, 
		       point p[], 
		       triangle t[],
		       double area[], 
		       int& numbnd, 
		       int*& bnd_node,
		       double (*SignedDistance)(point))
{
  // Sorting functions
  void QuickSort_int(int*& a, int*& index, int lo, int hi);
  void Unique_int(int*& vec, int*& index, int& size);

  const double hmin = 1.0e-14;
  int k = 0;
  int* bnd_node_tmp = new int[3*numtri];
  
  for (int i=0; i<numtri; i++)
    {
      point ptmp;
      ptmp.x = p[t[i].n1].x;
      ptmp.y = p[t[i].n1].y;	
      ptmp.z = 0.0;
      double dtmp = SignedDistance(ptmp);
      if (fabs(dtmp) <= hmin)
	{  
	  bnd_node_tmp[k] = t[i].n1; 
	  k = k+1; 
	}
      
      ptmp.x = p[t[i].n2].x;
      ptmp.y = p[t[i].n2].y;	   
      dtmp = SignedDistance(ptmp);
      if (fabs(dtmp) <= hmin)
	{
	  bnd_node_tmp[k] = t[i].n2;
	  k = k+1; 
	}
      
      ptmp.x = p[t[i].n3].x;
      ptmp.y = p[t[i].n3].y;	   
      dtmp = SignedDistance(ptmp);
      if (fabs(dtmp) <= hmin)
	{
	  bnd_node_tmp[k] = t[i].n3;
	  k = k+1;
	}
    }
  
  numbnd = k;
  int* index = new int[numbnd];
  QuickSort_int(bnd_node_tmp,index,0,numbnd-1);
  Unique_int(bnd_node_tmp,index,numbnd);
  
  bnd_node = new int[numbnd];
  for (int j=0; j<numbnd; j++)
    {
      bnd_node[j] = bnd_node_tmp[j];
    }    
  delete[] bnd_node_tmp;    
  delete[] index;
}
