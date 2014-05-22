#include "meshdefs.h"
#include <cstdio>
#include <cstdlib>

//
// Create information about edges
//
void MeshEdgeData(const int& numpts, 
		  const int& numtri,
		  const point p[],
		  const point psphere[],
		  const triangle t[], 
		  const double area[],
		  int& numedges, 
		  dTensor2*& edge, 
		  iTensor2*& tedge,
		  iTensor2*& eelem)
{
  // Needed sorting functions
  void QuickSort(double*& a, int*& index, int lo, int hi);
  void Unique(double*& vec, int*& index, int*& index_reverse, int& size);

  // Create helper variables and arrays
  dTensor2 edge_tmp(3*numtri,6);
  double* edge_label;
  edge_label = new double[3*numtri];

  // Loop over each element
  int k = 0;
  for (int i=0; i<numtri; i++)
    {
      double xmid,ymid,zmid;
      int tp1 = t[i].n1;
      int tp2 = t[i].n2;
      int tp3 = t[i].n3;
    
      // look at first edge
      k = k+1;
      edge_tmp.set(k,1, p[tp1].x );
      edge_tmp.set(k,2, p[tp1].y );
      edge_tmp.set(k,3, p[tp1].z );
      edge_tmp.set(k,4, p[tp2].x );
      edge_tmp.set(k,5, p[tp2].y );
      edge_tmp.set(k,6, p[tp2].z );

      xmid = 0.5*( p[tp1].x + p[tp2].x );
      ymid = 0.5*( p[tp1].y + p[tp2].y );
      zmid = 0.5*( p[tp1].z + p[tp2].z );
      edge_label[k-1] = xmid + 1e4*ymid + 1e8*zmid;

      // look at second edge
      k = k+1;
      edge_tmp.set(k,1, p[tp1].x );
      edge_tmp.set(k,2, p[tp1].y );
      edge_tmp.set(k,3, p[tp1].z );
      edge_tmp.set(k,4, p[tp3].x );
      edge_tmp.set(k,5, p[tp3].y );
      edge_tmp.set(k,6, p[tp3].z );

      xmid = 0.5*( p[tp1].x + p[tp3].x );
      ymid = 0.5*( p[tp1].y + p[tp3].y );
      zmid = 0.5*( p[tp1].z + p[tp3].z );
      edge_label[k-1] = xmid + 1e4*ymid + 1e8*zmid;
    
      // look at third edge
      k = k+1;
      edge_tmp.set(k,1, p[tp2].x );
      edge_tmp.set(k,2, p[tp2].y );
      edge_tmp.set(k,3, p[tp2].z );
      edge_tmp.set(k,4, p[tp3].x );
      edge_tmp.set(k,5, p[tp3].y );
      edge_tmp.set(k,6, p[tp3].z );

      xmid = 0.5*( p[tp2].x + p[tp3].x );
      ymid = 0.5*( p[tp2].y + p[tp3].y );
      zmid = 0.5*( p[tp2].z + p[tp3].z );
      edge_label[k-1] = xmid + 1e4*ymid + 1e8*zmid;
    }

  // Keep only unique edges
  int numedges_tmp = 3*numtri;
  int* index;
  index = new int[3*numtri];
  int* index_reverse;
  for (int i=0; i<(3*numtri); i++)
    {  index[i] = i;  }

  QuickSort(edge_label,index,0,3*numtri-1);  
  Unique(edge_label,index,index_reverse,numedges_tmp);

  // Store information up to this point
  numedges = numedges_tmp;
  edge  = new dTensor2(numedges,4);
  tedge = new iTensor2(numtri,3);
  eelem = new iTensor2(numedges,2);

  // edge coordinates
  for (int i=1; i<=numedges; i++)
    {
      int mu = index[i-1]+1;
      
      double   x1 = edge_tmp.get(mu,1);
      double   y1 = edge_tmp.get(mu,2);
      double   z1 = edge_tmp.get(mu,3);
      double lon1 = atan2(y1,x1) - pi;
      double lat1 = 0.5*pi - acos(z1/sqrt(x1*x1+y1*y1+z1*z1));

      double   x2 = edge_tmp.get(mu,4);
      double   y2 = edge_tmp.get(mu,5);
      double   z2 = edge_tmp.get(mu,6);
      double lon2 = atan2(y2,x2) - pi;
      double lat2 = 0.5*pi - acos(z2/sqrt(x2*x2+y2*y2+z2*z2));

      edge->set(i,1, lon1 );
      edge->set(i,2, lat1 );
      edge->set(i,3, lon2 );
      edge->set(i,4, lat2 );	
    }

  // edges attached to element i (physical elements)
  for (int i=1; i<=numtri; i++)
    {
      tedge->set(i,1, index_reverse[3*i-3]+1 );
      tedge->set(i,2, index_reverse[3*i-2]+1 );
      tedge->set(i,3, index_reverse[3*i-1]+1 );
    }
  
  // find the elements on either side of each edge
  for (int i=1; i<=numedges; i++)
    {
      int mfound = 0;
      int k = 1;
      while (mfound<2 && k<=numtri)
	{
	  for (int n=1; n<=3; n++)
	    {
	      int j = tedge->get(k,n);

	      if (i==j)
		{
		  mfound = mfound+1;
		  eelem->set(i,mfound, k);
		}
	    }
	  k=k+1;
	}

      if (mfound!=2)
	{
	  printf(" ERROR in MeshEdgeData.cpp:  mfound!=2. Should never get here ... \n");
	  printf("      i = %8i\n",i);
	  printf("      k = %8i\n",k);
	  printf(" mfound = %8i\n",mfound);
	  printf("\n");
	  exit(1);
	}

    }

  delete edge_label;
  delete index_reverse;
  delete index;

}
