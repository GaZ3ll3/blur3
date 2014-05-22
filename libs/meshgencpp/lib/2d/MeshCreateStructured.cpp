#include "meshdefs.h"

void MeshCreateStructured(double h0, 
			  point xymin, 
			  point xymax, 
			  cart2d& cart2dinfo,
			  int& numpts, 
			  int& numtri, 
			  point*& p, 
			  triangle*& t)
{
  // Find structured grid parameters mx, my, dx, and dy
  int mx = round((xymax.x - xymin.x)/h0)+1;
  int my = round((xymax.y - xymin.y)/h0)+1;
  
  double dx = (xymax.x-xymin.x)/double(mx-1);
  double dy = (xymax.y-xymin.y)/double(my-1);

  cart2dinfo.mx = mx;
  cart2dinfo.my = my;
  cart2dinfo.dx = dx;
  cart2dinfo.dy = dy;
  cart2dinfo.xlow  = xymin.x;
  cart2dinfo.ylow  = xymin.y;
  cart2dinfo.xhigh = xymax.x;
  cart2dinfo.yhigh = xymax.y;
  
  // Dimension p and t arrays
  numtri = 2*(mx-1)*(my-1);
  numpts = mx*my;
  p = new point[numpts];
  t = new triangle[numtri];
  
  // Set node information
  int k = 0;      
  for (int j=0; j<=(my-1); j++)
    for (int i=0; i<=(mx-1); i++)
      {
	p[k].x = xymin.x + double(i)*dx;
	p[k].y = xymin.y + double(j)*dy;
	k = k+1;
      }

  // Set element information
  k = 0;
  for (int j=1; j<=(my-1); j++)
    for (int i=1; i<=(mx-1); i++)	
      {	
	t[k].n1 = i       + (j-1)*mx - 1;
	t[k].n2 = i+1+mx  + (j-1)*mx - 1;
	t[k].n3 = i+mx    + (j-1)*mx - 1;
	k = k+1;
	
	t[k].n1 = i      + (j-1)*mx - 1;
	t[k].n2 = i+1    + (j-1)*mx - 1;
	t[k].n3 = i+1+mx + (j-1)*mx - 1;
	k = k+1;      
      }
}
