#include "meshdefs.h"
#include <cstdio>
#include <cstdlib>

//   Cleanup mesh (remove duplicate points, correctly orient triangles, and
//                 compute triangle areas)
//
void MeshCleanup(char*& GridType, 
		 double h0, 
		 double geps, 
		 int& numpts, 
		 int& numtri, 
		 point*& p, 
		 triangle*& t, 		 
		 double*& area, 
		 double*& cdual,
		 point xymin, 
		 point xymax,
		 double (*SignedDistance)(point))
{
  // Sorting functions
  void QuickSort(double*& a, int*& index, int lo, int hi);

  // ************************************************************************
  // ********** Remove any points that may appear more than once ************
  // ************************************************************************
  double* plist = new double[numpts];
  int* index = new int[numpts];
  
  if (GridType[0]=='U')
    {
      for (int i=0; i<numpts; i++)
        {
	  plist[i] = ((p[i].x-xymin.x)/(xymax.x-xymin.x)) 
	    + (600.0/h0)*(1.0+(p[i].y-xymin.y)/(xymax.y-xymin.y));
	  index[i] = i;
        }

      QuickSort(plist,index,0,numpts-1);
      
      point* psorted = new point[numpts];
      for (int i=0; i<numpts; i++)
        {
	  int ktmp = index[i];
	  psorted[i].x = p[ktmp].x;
	  psorted[i].y = p[ktmp].y;
        }   
      
      // Remove points that are too close (i.e., the same)
      point* punique = new point[numpts];
      punique[0].x = psorted[0].x;
      punique[0].y = psorted[0].y;
      int k = 0;
      for (int i=1; i<numpts; i++)
        {
	  double dij = sqrt(pow(psorted[i].x-punique[k].x,2)+
			    pow(psorted[i].y-punique[k].y,2));
	  if (dij>=1.0e-8)
            {
	      k = k+1;
	      punique[k].x = psorted[i].x;
	      punique[k].y = psorted[i].y;
            }
        }

      numpts = k+1;
      delete[] p;
      p = new point[numpts];

      for (int i=0; i<numpts; i++)
        {            
	  p[i].x = punique[i].x;
	  p[i].y = punique[i].y;
        }   
      delete[] plist;
      delete[] index;
      delete[] psorted;
      delete[] punique;

      FILE* pointfile = fopen("pointfile.txt","w");
      fprintf(pointfile,"%8i\n%8i\n",2,numpts);
      for (int i=0; i<numpts; i++)
        { fprintf(pointfile,"%24.16e %24.16e\n",p[i].x,p[i].y); }
      fclose(pointfile);

      // System call for Delaunay triangulation
      system("$QHULL/bin/qdelaunay Qt i  < ./pointfile.txt > ./trifile.txt");

      // Read-in triangulation from qhull
      int numtritemp;
      FILE* trifile = fopen("trifile.txt","r");
      fscanf(trifile,"%i",&numtritemp);
      triangle* ttemp = new triangle[numtritemp];
      point* pmid = new point[numtritemp];

      for (int i=0; i<numtritemp; i++)
        {  fscanf(trifile,"%i %i %i",&ttemp[i].n1,&ttemp[i].n2,&ttemp[i].n3);  }
      fclose(trifile);
      
      // Remove files
      system("rm -f pointfile.txt trifile.txt");

      numtri = 0;
      for (int i=0; i<numtritemp; i++)
        {
	  pmid[i].x = (p[ttemp[i].n1].x + p[ttemp[i].n2].x + p[ttemp[i].n3].x) / 3.0;
	  pmid[i].y = (p[ttemp[i].n1].y + p[ttemp[i].n2].y + p[ttemp[i].n3].y) / 3.0;
	  if (SignedDistance(pmid[i]) < -geps)
            {  numtri = numtri + 1;  }
        }
      
      k=0;
      delete[] t;
      t = new triangle[numtri];
      for (int i=0; i<numtritemp; i++)
        {
	  if (SignedDistance(pmid[i]) < -geps)
            {
	      t[k].n1 = ttemp[i].n1;
	      t[k].n2 = ttemp[i].n2;
	      t[k].n3 = ttemp[i].n3;
	      k = k + 1;
            }
        }
      
      delete[] pmid;
      delete[] ttemp;
    }

  // ************************************************************************
  // ************* find triangle areas and dual mesh areas ******************
  // ************************************************************************
  delete[] cdual;
  delete[] area;
  cdual = new double[numpts];
  area  = new double[numtri];
  for (int i=0; i<numpts; i++)
    {  cdual[i] = 0.0;  }

  // compute triangle areas, force CCW orientation, compute dual mesh areas
  for (int i=0; i<numtri; i++)
    {
      area[i] = .5 * (p[t[i].n1].x * (p[t[i].n2].y-p[t[i].n3].y) + 
		      p[t[i].n2].x * (p[t[i].n3].y-p[t[i].n1].y) +
		      p[t[i].n3].x * (p[t[i].n1].y-p[t[i].n2].y));
      if (area[i] < 0.0)
        {
	  double temp = t[i].n2;
	  t[i].n2 = t[i].n3;
	  t[i].n3 = temp;
	  area[i] = fabs(area[i]);
        }
      cdual[t[i].n1] += (area[i]/3.0);
      cdual[t[i].n2] += (area[i]/3.0);
      cdual[t[i].n3] += (area[i]/3.0);
    }

}
