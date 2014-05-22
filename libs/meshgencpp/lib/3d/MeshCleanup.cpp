#include "meshdefs.h"
#include <cstdlib>
#include <cstdio>

//
//   Cleanup mesh (remove duplicate points, correctly orient triangles, and
//                 compute triangle areas)
//
void MeshCleanup(const char* GridType, 
		 double h0, 
		 double geps, 
		 int& numpts, 
		 int& numtet, 
		 point*& p, 
		 tetra*& t, 
		 double*& volume, 
		 double*& cdual,
		 double (*SignedDistance)(point))
{
  // Sorting functions
  void QuickSort(double*& a, int*& index, int lo, int hi);
  
  // ************************************************************************
  // ********** Remove any points that may appear more than once ************
  // ************************************************************************
  if (GridType[0]=='U')
    {  
      double maxp = -1.e8;
      double minp =  1.e8;
      double snap;      
      double tmp1,tmp2;
      double* plist = new double[numpts];
      int* index = new int[numpts];
      int ktmp;
      point* pnew = new point[numpts];
      int numtettemp = numtet;
      tetra* ttemp = new tetra[numtettemp];
      point*  pmid = new point[numtettemp];
      for (int i=0; i<=numtettemp; i++)
	{
	  ttemp[i].n1 = t[i].n1;
	  ttemp[i].n2 = t[i].n2;
	  ttemp[i].n3 = t[i].n3;
	  ttemp[i].n4 = t[i].n4;
	}      
      
      for (int i=0; i<numpts; i++)
	{
	  if (p[i].x > maxp)
	    { maxp = p[i].x; }
	  
	  if (p[i].y > maxp)
	    { maxp = p[i].y; }
	  
	  if (p[i].z > maxp)
	    { maxp = p[i].z; }
	  
	  if (p[i].x < minp)
	    { minp = p[i].x; }
	  
	  if (p[i].y < minp)
	    { minp = p[i].y; }
	  
	  if (p[i].z < minp)
	    { minp = p[i].z; }
	}        
      
      snap = (maxp-minp)*1024.e0*pow(2.0,-52);
  
      for (int i=0; i<numpts; i++)
	{
	  tmp1 = p[i].x/snap;
	  tmp2 = floor(tmp1);
	  if ( (tmp1-tmp2) > 0.5e0 )
	    { tmp2 = tmp2 + 1.0e0; } 
	  p[i].x = tmp2*snap;
	  
	  tmp1 = p[i].y/snap;
	  tmp2 = floor(tmp1);
	  if ( (tmp1-tmp2) > 0.5e0 )
	    { tmp2 = tmp2 + 1.0e0; } 
	  p[i].y = tmp2*snap;
	  
	  tmp1 = p[i].z/snap;
	  tmp2 = floor(tmp1);
	  if ( (tmp1-tmp2) > 0.5e0 )
	    { tmp2 = tmp2 + 1.0e0; } 
	  p[i].z = tmp2*snap;
	}
      
      for (int i=0; i<numpts; i++)
	{
	  plist[i] = (36000.e0/h0)*p[i].x + (600.e0/h0)*p[i].y + p[i].z;
	  index[i] = i;
	}
      
      QuickSort(plist,index,0,numpts-1);
      
      point* psorted = new point[numpts];
      for (int i=0; i<numpts; i++)
	{
	  ktmp = index[i];
	  psorted[i].x = p[ktmp].x;
	  psorted[i].y = p[ktmp].y;
	  psorted[i].z = p[ktmp].z;
	}        
      
      // Remove points that are too close (i.e., the same)
      point* punique = new point[numpts];
      punique[0].x = psorted[0].x;
      punique[0].y = psorted[0].y;
      punique[0].z = psorted[0].z;
      int k = 0;
      double dij;
      for (int i=1; i<numpts; i++)
	{
	  dij = sqrt(pow(psorted[i].x-punique[k].x,2)+
		     pow(psorted[i].y-punique[k].y,2)+
		     pow(psorted[i].z-punique[k].z,2));
	  if (dij>=1.0e-8)
	    {
	      k = k+1;
	      punique[k].x = psorted[i].x;
	      punique[k].y = psorted[i].y;
	      punique[k].z = psorted[i].z;
	    }
	}
      
      numpts = k+1;
      delete[] p;
      p = new point[numpts];
      
      for (int i=0; i<numpts; i++)
	{            
	  p[i].x = punique[i].x;
	  p[i].y = punique[i].y;
	  p[i].z = punique[i].z;
	}         

      delete psorted;
      delete punique;
      
      FILE* pointfile = fopen("pointfile.txt","w");
      fprintf(pointfile,"%8i\n%8i\n",3,numpts);
      for (int i=0; i<numpts; i++)
	{ fprintf(pointfile,"%24.16e %24.16e %24.16e\n",p[i].x,p[i].y,p[i].z); }
      fclose(pointfile);
	  
      // System call for Delaunay triangulation
      system ("$QHULL/bin/qdelaunay Qt i  < ./pointfile.txt > ./tetfile.txt");

      // Read-in triangulation from qhull
      FILE* tetfile = fopen("tetfile.txt","r");
      fscanf(tetfile,"%i",&numtettemp);
      delete ttemp;
      delete pmid;
      ttemp = new tetra[numtettemp];
      pmid = new point[numtettemp];

      for (int i=0; i<numtettemp; i++)
	{  fscanf(tetfile,"%i %i %i %i",&ttemp[i].n1,&ttemp[i].n2,&ttemp[i].n3,&ttemp[i].n4);  }
      fclose(tetfile);
      
      // Remove files
      system("rm -f pointfile.txt tetfile.txt");

      numtet = 0;    
      for (int i=0; i<numtettemp; i++)
	{
	  pmid[i].x = (p[ttemp[i].n1].x + p[ttemp[i].n2].x + p[ttemp[i].n3].x + p[ttemp[i].n4].x) / 4.0e0;
	  pmid[i].y = (p[ttemp[i].n1].y + p[ttemp[i].n2].y + p[ttemp[i].n3].y + p[ttemp[i].n4].y) / 4.0e0;
	  pmid[i].z = (p[ttemp[i].n1].z + p[ttemp[i].n2].z + p[ttemp[i].n3].z + p[ttemp[i].n4].z) / 4.0e0;
	  if (SignedDistance(pmid[i]) < -geps)
	    {	  
	      numtet++;
	    }
	}  
      
      delete[] t;
      t = new tetra[numtet];
      ktmp = 0;
      for (int i=0; i<numtettemp; i++)
	{
	  pmid[i].x = (p[ttemp[i].n1].x + p[ttemp[i].n2].x + p[ttemp[i].n3].x + p[ttemp[i].n4].x) / 4.0e0;
	  pmid[i].y = (p[ttemp[i].n1].y + p[ttemp[i].n2].y + p[ttemp[i].n3].y + p[ttemp[i].n4].y) / 4.0e0;
	  pmid[i].z = (p[ttemp[i].n1].z + p[ttemp[i].n2].z + p[ttemp[i].n3].z + p[ttemp[i].n4].z) / 4.0e0;
	  if (SignedDistance(pmid[i]) < -geps)
	    {
	      t[ktmp].n1 = ttemp[i].n1;
	      t[ktmp].n2 = ttemp[i].n2;
	      t[ktmp].n3 = ttemp[i].n3;
	      t[ktmp].n4 = ttemp[i].n4;
	      ktmp++;	  
	    }
	}
       
      delete plist;
      delete index;
      delete pnew;
      delete ttemp;
      delete pmid;

    }

  // ************************************************************************
  // *************** get triangle areas  ( == .5 * |A x B| ) ****************
  // ************************* and dual mesh areas **************************
  // ************************************************************************
  volume = new double[numtet];
  cdual  = new double[numpts];
  
  double x1,x2,x3,x4;
  double y1,y2,y3,y4;
  double z1,z2,z3,z4;
  
  for (int i=0; i<numtet; i++)
    {      
      x1 = p[t[i].n1].x;
      y1 = p[t[i].n1].y;
      z1 = p[t[i].n1].z;
      
      x2 = p[t[i].n2].x;
      y2 = p[t[i].n2].y;
      z2 = p[t[i].n2].z;
      
      x3 = p[t[i].n3].x;
      y3 = p[t[i].n3].y;
      z3 = p[t[i].n3].z;
      
      x4 = p[t[i].n4].x;
      y4 = p[t[i].n4].y;
      z4 = p[t[i].n4].z;
      
      volume[i] = fabs(  x1*y2*z3 - x1*y2*z4 - x1*y3*z2 + x1*y3*z4 + x1*y4*z2 - x1*y4*z3 
			 - x2*y1*z3 + x2*y1*z4 + x2*y3*z1 - x2*y3*z4 - x2*y4*z1 + x2*y4*z3
			 + x3*y1*z2 - x3*y1*z4 - x3*y2*z1 + x3*y2*z4 + x3*y4*z1 - x3*y4*z2
			 - x4*y1*z2 + x4*y1*z3 + x4*y2*z1 - x4*y2*z3 - x4*y3*z1 + x4*y3*z2 )/6.0e0;
    }
  
  for (int i=0; i<numtet; i++)
    {
      cdual[t[i].n1] += volume[i]/4.0e0;
      cdual[t[i].n2] += volume[i]/4.0e0;
      cdual[t[i].n3] += volume[i]/4.0e0;
      cdual[t[i].n4] += volume[i]/4.0e0;
    }
}
