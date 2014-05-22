#include "meshdefs.h"
#include <cstdlib>
#include <cstdio>

//
// Create information about edges
//
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
		  double (*SignedDistance)(point))
{
  // Needed sorting functions
  void QuickSort(double*& a, int*& index, int lo, int hi);
  void Unique(double*& vec, int*& index, int*& index_reverse, int& size);

  // Create helper variables and arrays
  const int numtri_phys = numtri-numghost;
  dTensor2 edge_tmp(3*numtri_phys + numghost,4);
  dTensor2 enode_tmp(3*numtri_phys + numghost,2);
  double* edge_label;
  edge_label = new double[3*numtri_phys + numghost];
  const double hmin = 0.01*h0;
  const double bignumber = 1.0/hmin;

  // Find xmin and ymax
  void XminYmin(const int& numpts, 
		const point p[],
		double& Xmin,
		double& Ymin);
  double Xmin,Ymin;
  XminYmin(numpts,p,Xmin,Ymin);

  // Loop over each physical element
  int k = 0;
  for (int i=0; i<numtri_phys; i++)
    {
      const int tp1 = t[i].n1;
      const int tp2 = t[i].n2;
      const int tp3 = t[i].n3;
    
      // look at first edge
      k = k+1;
      edge_tmp.set(k,1, p[tp1].x );
      edge_tmp.set(k,2, p[tp1].y );
      edge_tmp.set(k,3, p[tp2].x );
      edge_tmp.set(k,4, p[tp2].y );
      
      enode_tmp.set(k,1, tp1+1 );
      enode_tmp.set(k,2, tp2+1 );

      {
	const double xmid = 0.5*( p[tp1].x + p[tp2].x );
	const double ymid = 0.5*( p[tp1].y + p[tp2].y );
	edge_label[k-1] = (bignumber)*(xmid-Xmin+10.0e0) + (ymid-Ymin+10.0e0);
      }

      // look at second edge
      k = k+1;
      edge_tmp.set(k,1, p[tp1].x );
      edge_tmp.set(k,2, p[tp1].y );
      edge_tmp.set(k,3, p[tp3].x );
      edge_tmp.set(k,4, p[tp3].y );

      enode_tmp.set(k,1, tp1+1 );
      enode_tmp.set(k,2, tp3+1 );
      
      {
	const double xmid = 0.5*( p[tp1].x + p[tp3].x );
	const double ymid = 0.5*( p[tp1].y + p[tp3].y );
	edge_label[k-1] = (bignumber)*(xmid-Xmin+10.0e0) + (ymid-Ymin+10.0e0);
      }
    
      // look at third edge
      k = k+1;
      edge_tmp.set(k,1, p[tp2].x );
      edge_tmp.set(k,2, p[tp2].y );
      edge_tmp.set(k,3, p[tp3].x );
      edge_tmp.set(k,4, p[tp3].y );
      
      enode_tmp.set(k,1, tp2+1 );
      enode_tmp.set(k,2, tp3+1 );

      {
	const double xmid = 0.5*( p[tp2].x + p[tp3].x );
	const double ymid = 0.5*( p[tp2].y + p[tp3].y );
	edge_label[k-1] = (bignumber)*(xmid-Xmin+10.0e0) + (ymid-Ymin+10.0e0);
      }
    }

  // loop over each ghost element
  for (int i=numtri_phys; i<numtri; i++)
    {
      const int tp1 = t[i].n1;
      const int tp2 = t[i].n2;
      const int tp3 = t[i].n3;      
    
      // look at first edge
      point ptmp;
      ptmp.x = p[tp3].x;
      ptmp.y = p[tp3].y;
      ptmp.z = 0.0;
      if (SignedDistance(ptmp)>hmin)
	{
	  k = k+1;
	  edge_tmp.set(k,1, p[tp1].x );
	  edge_tmp.set(k,2, p[tp1].y );
	  edge_tmp.set(k,3, p[tp2].x );
	  edge_tmp.set(k,4, p[tp2].y );

	  enode_tmp.set(k,1, tp1+1 );
	  enode_tmp.set(k,2, tp2+1 );
	  
	  const double xmid = 0.5*( p[tp1].x + p[tp2].x );
	  const double ymid = 0.5*( p[tp1].y + p[tp2].y );
	  edge_label[k-1] = (bignumber)*(xmid-Xmin+10.0e0) + (ymid-Ymin+10.0e0);
	}
      
      // look at second edge
      ptmp.x = p[tp2].x;
      ptmp.y = p[tp2].y;
      if (SignedDistance(ptmp)>hmin)
	{	  
	  k = k+1;
	  edge_tmp.set(k,1, p[tp1].x );
	  edge_tmp.set(k,2, p[tp1].y );
	  edge_tmp.set(k,3, p[tp3].x );
	  edge_tmp.set(k,4, p[tp3].y );
	  
	  enode_tmp.set(k,1, tp1+1 );
	  enode_tmp.set(k,2, tp3+1 );
	  
	  const double xmid = 0.5*( p[tp1].x + p[tp3].x );
	  const double ymid = 0.5*( p[tp1].y + p[tp3].y );
	  edge_label[k-1] = (bignumber)*(xmid-Xmin+10.0e0) + (ymid-Ymin+10.0e0);
	}
	 
      // look at third edge
      ptmp.x = p[tp1].x;
      ptmp.y = p[tp1].y;
      if (SignedDistance(ptmp)>hmin)
	{	 
	  k = k+1;
	  edge_tmp.set(k,1, p[tp2].x );
	  edge_tmp.set(k,2, p[tp2].y );
	  edge_tmp.set(k,3, p[tp3].x );
	  edge_tmp.set(k,4, p[tp3].y );

	  enode_tmp.set(k,1, tp2+1 );
	  enode_tmp.set(k,2, tp3+1 );

	  const double xmid = 0.5*( p[tp2].x + p[tp3].x );
	  const double ymid = 0.5*( p[tp2].y + p[tp3].y );
	  edge_label[k-1] = (bignumber)*(xmid-Xmin+10.0e0) + (ymid-Ymin+10.0e0);
	}
    }

  // Keep only unique edges
  int numedges_tmp = (3*numtri_phys + numghost);
  int* index;
  index = new int[3*numtri_phys + numghost];
  int* index_reverse;
  for (int i=0; i<(3*numtri_phys + numghost); i++)
    {  index[i] = i;  }

  QuickSort(edge_label,index,0,3*numtri_phys+numghost-1);  
  Unique(edge_label,index,index_reverse,numedges_tmp);

  // Store information up to this point
  numedges = numedges_tmp;
  edge  = new dTensor2(numedges,4);
  tedge = new iTensor2(numtri,3);
  eelem = new iTensor2(numedges,2);
  enode = new iTensor2(numedges,2);

  // edge coordinates
  for (int i=1; i<=numedges; i++)
    {
      for (int k=1; k<=2; k++)
	{
	  enode->set(i,k,enode_tmp.get(index[i-1]+1,k));
	}
      for (int k=1; k<=4; k++)
	{
	  edge->set(i,k, edge_tmp.get(index[i-1]+1,k) );
	}
    }

  // initialize tedge
  for (int i=1; i<=numtri; i++)
    for (int k=1; k<=3; k++)
      {
	tedge->set(i,k, -1 );
      }

  // edges attached to element i (physical elements)
  for (int i=1; i<=numtri_phys; i++)
    {
      tedge->set(i,1, index_reverse[3*i-3]+1 );
      tedge->set(i,2, index_reverse[3*i-2]+1 );
      tedge->set(i,3, index_reverse[3*i-1]+1 );
    }

  // edges attached to element i (ghost elements)
  k=3*numtri_phys;
  for (int i=(numtri_phys+1); i<=numtri; i++)
    {
      k = k+1;
      tedge->set(i,1, index_reverse[k-1]+1 );
      tedge->set(i,2, -1 );
      tedge->set(i,3, -1 );
    }
  
  // find the elements on either side of each edge
  for (int i=1; i<=numedges; i++)
    {
      int mfound = 0;
      k = 1;
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
	  printf("      i = %i\n",i);
	  printf("      k = %i\n",k);
	  printf(" mfound = %i\n",mfound);
	  printf("\n");
	  exit(1);
	}

    }
  delete index;
  delete index_reverse;
  delete edge_label;
  
  // Find all edges that live on the boundary of the mesh
  numbndedges = 0;
  iTensor1 bnd_edge_tmp(numedges);
  for (int i=1; i<=numedges; i++)
    {
      int t1 = eelem->get(i,1);
      int t2 = eelem->get(i,2);

      if (t1>(numtri-numghost) && t2<=(numtri-numghost))
	{
	  numbndedges = numbndedges + 1;
	  bnd_edge_tmp.set( numbndedges, i );
	}
      else if (t2>(numtri-numghost) && t1<=(numtri-numghost))
	{
	  numbndedges = numbndedges + 1;
	  bnd_edge_tmp.set( numbndedges, i );
	}
    }
  bnd_edge = new iTensor1(numbndedges);
  for (int i=1; i<=numbndedges; i++)
    {
      bnd_edge->set(i, bnd_edge_tmp.get(i) );
    }

}


void XminYmin(const int& numpts, 
	      const point p[],
	      double& Xmin,
	      double& Ymin)
{
  Xmin = p[0].x;
  Ymin = p[0].y;

  for (int i=1; i<numpts; i++)
    {
      double xx = p[i].x;
      if (xx < Xmin)
	{  Xmin = xx;  }
      double yy = p[i].y;
      if (yy < Ymin)
	{  Ymin = yy;  }
    }
}
