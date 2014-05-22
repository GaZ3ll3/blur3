#include <cstdlib>
#include "meshdefs.h"

void MeshCreateInitial(point xymin, 
		       point xymax, 
		       double h0, 
		       double geps, 
		       int numfixpts, 
		       point* fixpt, 
		       int& numpts, point*& p,
		       double (*SignedDistance)(point),
		       double (*GridSpacing)(point))
{
  int i,j;
  int m = static_cast<int>(ceil((xymax.x-xymin.x)/h0));
  int n = static_cast<int>(ceil((xymax.y-xymin.y)/h0/sqrt(3.0)*2));
  int maxpoints = m*n;

  // create needed dynamic arrays
  point* temp1;
  point* temp2;
  double* r0;
  double* random;
  
  temp1 = new point[maxpoints];
  temp2 = new point[maxpoints];
  p = new point[maxpoints + numfixpts];
  r0 = new double[maxpoints];
  random = new double[maxpoints];
  
  // Generate uniform grid
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      {
	temp1[i*n+j].x = xymin.x + i*h0;
	temp1[i*n+j].y = xymin.y + j*h0*sqrt(3.0)/2;
      }

  // Shift even rows
  for (i=0; i<m; i++)
    for (j=1; j<n; j+=2)
      {  temp1[i*n+j].x += h0/2;  }

  // Delete points outside region
  int temp2elts = 0;
  for (i=0; i<m; i++)
    for (j=0; j<n; j++)
      if (SignedDistance(temp1[i*n+j]) < geps)
	{
	  temp2[temp2elts] = temp1[i*n+j];
	  temp2elts++;
	}
  
  // Probability to keep point
  for (i=0; i<temp2elts; i++)
    {  r0[i]=(1.0 / (GridSpacing(temp2[i])) / (GridSpacing(temp2[i])));  }

  // find max(r0)
  double maxr0 = r0[0];
  for (i=1; i<temp2elts; i++)
    if (r0[i] > maxr0)
      {  maxr0 = r0[i];  }
	
  // Rejection method
  for (i=0; i<numfixpts; i++)			// add fixed points
    {  p[i]=fixpt[i];  }
  
  numpts = numfixpts;					// number of elements in p
  
  for (i=0; i<temp2elts; i++)			
    {
      random[i]=((RAND_MAX - rand()) / static_cast<double>(RAND_MAX));
      if (random[i] < (r0[i]/maxr0))
	{
	  p[numpts] = temp2[i];
	  numpts++;
	}
    }

  delete [] temp1;
  delete [] r0;
  delete [] random;
  delete [] temp2;
}
