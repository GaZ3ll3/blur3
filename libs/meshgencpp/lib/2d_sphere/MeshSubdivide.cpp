#include "meshdefs.h"
#include <cstdio>

void MeshSubdivide(const int level, 
		   const double radius, 
		   int& numpts, 
		   int& numtri, 
		   point*& p, 
		   triangle*& t)
{
  void QuickSort_bar(bar* b, int* key, int i, int j);
  
  for (int ell=1; ell<=level; ell++)
    {
      bar* bars;
      bar* bartemp;
      int* priority;
      bars = new bar[3*numtri];
      bartemp = new bar[3*numtri];
      priority = new int[3*numtri];
      
      int numbar = 0;
      int temp;
  
      // find unique listing of bars
      for (int i=0; i<numtri; i++)
	{
	  bartemp[3*i].n1   = t[i].n1;
	  bartemp[3*i].n2   = t[i].n2;
	  bartemp[3*i+1].n1 = t[i].n1;
	  bartemp[3*i+1].n2 = t[i].n3;
	  bartemp[3*i+2].n1 = t[i].n2;
	  bartemp[3*i+2].n2 = t[i].n3;
	}

      for (int i=0; i<3*numtri; i++)
	{
	  // sort each bar so that n1 < n2
	  if (bartemp[i].n1 > bartemp[i].n2)
	    {
	      temp = bartemp[i].n1;
	      bartemp[i].n1 = bartemp[i].n2;
	      bartemp[i].n2 = temp;
	    }

	  // assign priority for sorting list of bars
	  priority[i] = numpts*bartemp[i].n1 + bartemp[i].n2;
	}

      // sort the bars according to priority      
      QuickSort_bar(bartemp,priority,0,3*numtri-1);

      // keep only unique bars
      int value = 0;
      for (int i=0; i<3*numtri; i++)
	if (priority[i] != value)
	  {
	    bars[numbar].n1 = bartemp[i].n1;
	    bars[numbar].n2 = bartemp[i].n2;
	    value = priority[i];
	    numbar++;
	  }

      // store all old points and newly created points
      // (located at bar midpoints) in pnew
      point* pnew;
      pnew = new point[numpts+numbar];
      for (int i=0; i<numpts; i++)
	{
	  pnew[i].x = p[i].x;
	  pnew[i].y = p[i].y;
	  pnew[i].z = p[i].z;
	}

      for (int i=0; i<numbar; i++)
	{
	  pnew[i+numpts].x = .5 * (p[bars[i].n1].x + p[bars[i].n2].x);
	  pnew[i+numpts].y = .5 * (p[bars[i].n1].y + p[bars[i].n2].y);
	  pnew[i+numpts].z = .5 * (p[bars[i].n1].z + p[bars[i].n2].z);
	  bars[i].bt = i+numpts;
	}
      numpts = numpts+numbar;
  
      // create new triangle list
      triangle* tnew;
      tnew = new triangle[4*numtri];      

      for (int i=0; i<numtri; i++)
	{
	  int mstop = 0;
	  int j = 1;
	  int ind1,ind2,ind3;
	  int bin1,bin2,tin1,tin2,tin3;
	  
	  tin1 = t[i].n1;
	  tin2 = t[i].n2;
	  tin3 = t[i].n3;
      
	  while(j<=numbar && mstop<3)
	    {
	      bin1 = bars[numbar-j].n1;
	      bin2 = bars[numbar-j].n2;
	      
	      if ((tin1 == bin1) && (tin2 == bin2))
		{ 
		  ind1 = numbar-j; 
		  mstop = mstop + 1;
		}
	      else if ((tin1 == bin1) && (tin3 == bin2))
		{
		  ind2 = numbar-j;
		  mstop = mstop + 1;
		}
	      else if ((tin2 == bin1) && (tin3 == bin2))
		{
		  ind3 = numbar-j;
		  mstop = mstop + 1;
		}
	      j = j+1;
	    }
	  tnew[4*i+0].n1 = t[i].n1;
	  tnew[4*i+0].n2 = bars[ind1].bt;
	  tnew[4*i+0].n3 = bars[ind2].bt;
	  tnew[4*i+1].n1 = t[i].n2;
	  tnew[4*i+1].n2 = bars[ind1].bt;
	  tnew[4*i+1].n3 = bars[ind3].bt;
	  tnew[4*i+2].n1 = t[i].n3;
	  tnew[4*i+2].n2 = bars[ind2].bt;
	  tnew[4*i+2].n3 = bars[ind3].bt;
	  tnew[4*i+3].n1 = bars[ind1].bt;
	  tnew[4*i+3].n2 = bars[ind2].bt;
	  tnew[4*i+3].n3 = bars[ind3].bt;
	}

      delete []t;
      delete []p;

      numtri = 4*numtri;
      
      t = tnew;
      p = pnew;
      
      // project new points to sphere surface
      for (int i=0; i<numpts; i++)
	{
	  // Get (x,y,z) Cartesian coordinates
	  double x = p[i].x;
	  double y = p[i].y;
	  double z = p[i].z;
	  
	  // Compute spherical coordinates
	  double   rad = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	  double theta = acos(z/rad);
	  double   phi = atan2(y,x);
	  
	  // Project point onto a sphere of radius "radius"
	  p[i].x = radius * sin(theta) * cos(phi);
	  p[i].y = radius * sin(theta) * sin(phi);
	  p[i].z = radius * cos(theta);
	}

      // Report that current level is finished
      printf("  Finished level:  %5i\n",ell);      

      delete bars;
      delete bartemp;
      delete priority;
    }
  
}


