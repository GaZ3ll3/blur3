#include <cstdlib>
#include "meshdefs.h"

void MeshCreateInitial(point xyzmin, point xyzmax, double h0, double geps, int numfixpts, 
                       point* fixpt, int& numpts, point*& p,
                       double (*SignedDistance)(point),
                       double (*GridSpacing)(point))
{
  // create integers to iterate for loops
  int i,j,k;
  
  // determine number of points in initial grid
  int m = static_cast<int>(ceil((xyzmax.x-xyzmin.x)/h0) + 1);
  int n = static_cast<int>(ceil((xyzmax.y-xyzmin.y)/h0) + 1);
  int o = static_cast<int>(ceil((xyzmax.z-xyzmin.z)/h0) + 1);
  int maxpoints = m*n*o;

  // create needed dynamic arrays
  point* temp1;
  point* temp2;
  double* r0;
  double* random;
  
  temp1  = new  point[maxpoints];
  temp2  = new  point[maxpoints];
  p      = new  point[maxpoints + numfixpts];
  r0     = new double[maxpoints];
  random = new double[maxpoints];
  
  // Generate uniform grid
  for (k=0; k<o; k++)
    for (j=0; j<n; j++)
      for (i=0; i<m; i++)
	{
	  temp1[(k*n*m)+(j*m)+i].x = xyzmin.x + i*h0;
	  temp1[(k*n*m)+(j*m)+i].y = xyzmin.y + j*h0;
	  temp1[(k*n*m)+(j*m)+i].z = xyzmin.z + k*h0;
	}
  
  // Delete points outside region
  int temp2elts = 0;
  for (k=0; k<o; k++)
    for (j=0; j<n; j++)
      for (i=0; i<m; i++)
	if (SignedDistance(temp1[(k*n*m)+(j*m)+i]) < geps)
	  {
	    temp2[temp2elts].x = temp1[(k*n*m)+(j*m)+i].x;
	    temp2[temp2elts].y = temp1[(k*n*m)+(j*m)+i].y;
	    temp2[temp2elts].z = temp1[(k*n*m)+(j*m)+i].z;
	    temp2elts++;	    
	  }

  // Probability to keep point
  for (i=0; i<temp2elts; i++)
    {  r0[i]= GridSpacing(temp2[i]);  }

  // find min(r0)
  double minr0 = r0[0];
  for (i=1; i<temp2elts; i++)
    {
      if (r0[i] < minr0)
	{  minr0 = r0[i];  }
    }

  // Rejection method
  if (numfixpts>0)
    {
      for (i=0; i<numfixpts; i++) // add fixed points
	{  p[i]=fixpt[i];  }
    }

  numpts = numfixpts;  // number of elements in p
  for (i=0; i<temp2elts; i++)			
    {
      random[i]=((RAND_MAX - rand()) / static_cast<double>(RAND_MAX));
      if (random[i] <= (minr0*minr0*minr0)/(r0[i]*r0[i]*r0[i]))
	{
	  p[numpts].x = temp2[i].x;
	  p[numpts].y = temp2[i].y;
	  p[numpts].z = temp2[i].z;
	  numpts++;
	}
    }

  delete [] temp1;
  delete [] r0;
  delete [] random;
  delete [] temp2;
}
