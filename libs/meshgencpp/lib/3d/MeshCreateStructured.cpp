#include "meshdefs.h"

void MeshCreateStructured(double h0, point xyzmin, point xyzmax, int& numpts, 
			  int& numtet, point*& p, tetra*& t)
{
  // Find structured grid parameters mx, my, mz, dx, dy, and dz
  int mx = round((xyzmax.x - xyzmin.x)/h0)+1;
  int my = round((xyzmax.y - xyzmin.y)/h0)+1;
  int mz = round((xyzmax.z - xyzmin.z)/h0)+1;
  
  double dx = (xyzmax.x-xyzmin.x)/double(mx-1);
  double dy = (xyzmax.y-xyzmin.y)/double(my-1);
  double dz = (xyzmax.z-xyzmin.z)/double(mz-1);
  
  // Dimension p and t arrays
  numtet = 6*(mx-1)*(my-1)*(mz-1);
  numpts = mx*my*mz;
  p = new point[numpts];
  t = new tetra[numtet];

  // Set node information
  int m = 0; 
  for (int k=0; k<=(mz-1); k++)
    for (int j=0; j<=(my-1); j++)
      for (int i=0; i<=(mx-1); i++)
	{
	  p[m].x = xyzmin.x + double(i)*dx;
	  p[m].y = xyzmin.y + double(j)*dy;
	  p[m].z = xyzmin.z + double(k)*dz;
	  m = m+1;
	}
  
  // Set element information
  m = 0;
  for (int k=1; k<=(mz-1); k++)
    for (int j=1; j<=(my-1); j++)
      for (int i=1; i<=(mx-1); i++)	
	{
	  int p1 = i + (j-1)*mx + (k-1)*mx*my - 1;
	  int p2 = p1 + 1;
	  int p3 = p1 + mx;
	  int p4 = p3 + 1;
	  int p5 = p1 + mx*my;
	  int p6 = p5 + 1;
	  int p7 = p3 + mx*my;
	  int p8 = p7 + 1;

	  t[m].n1 = p2;
	  t[m].n2 = p6;
	  t[m].n3 = p7;
	  t[m].n4 = p8;
	  m = m+1;
	  
	  t[m].n1 = p2;
	  t[m].n2 = p4;
	  t[m].n3 = p7;
	  t[m].n4 = p8;
	  m = m+1;
	  
	  t[m].n1 = p2;
	  t[m].n2 = p3;
	  t[m].n3 = p4;
	  t[m].n4 = p7;
	  m = m+1;
	  
	  t[m].n1 = p2;
	  t[m].n2 = p5;
	  t[m].n3 = p6;
	  t[m].n4 = p7;
	  m = m+1;
	      
	  t[m].n1 = p1;
	  t[m].n2 = p2;
	  t[m].n3 = p5;
	  t[m].n4 = p7;
	  m = m+1;

	  t[m].n1 = p1;
	  t[m].n2 = p2;
	  t[m].n3 = p3;
	  t[m].n4 = p7;
	  m = m+1;
	}

}
