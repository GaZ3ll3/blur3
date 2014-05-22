#include "meshdefs.h"
#include <cstdlib>
#include <cstdio>

//
//  Takes initial point distribution "p", triangulates with "t", and then
//  iterates on the point distribution and triangulation until an equilibrium
//  mesh is created (or until the number of iterations exceeds "maxiters").
//
void MeshIterate(int maxiters, 
		 double geps, 
		 double ttol, 
		 double dptol, 
		 double Fscale, 
                 double deltat, 
		 double deps, 
		 double h0, 
		 int numfixpts,
		 int& numpts, 
		 int& numtet, 
                 point*& p, 
		 tetra*& t,
		 double (*SignedDistance)(point),
                 double (*GridSpacing)(point))
{
  // Sorting function
  void QuickSort_bar(bar* b, int* key, int i, int j);
  
  // Set up all needed variables and dynamic arrays for while-loop iteration
  bool retriang = true;		// force initial retriangulation
  bool terminate = false;		// end when termination criterion satisfied

  // create and allocate dynamic arrays
  point* pold;
  tetra* ttemp;
  point* pmid;
  bar* bartemp;
  bar* bars;
  int* priority;
  point* barvec;
  double* L;
  point* midpt;
  double* hbars;
  double* L0;
  double* F;
  point* Fvec;
  point* Ftot;
  double* dist;
  
  int maxtet = 9*numpts;
  int maxbar = 6*maxtet;
  pold = new point[numpts];
  ttemp = new tetra[maxtet];
  pmid = new point[maxtet];
  t = new tetra[maxtet]; 
  bartemp = new bar[maxbar];
  bars = new bar[maxbar];
  priority = new int[maxbar];
  barvec = new point[maxbar];
  L = new double[maxbar];
  midpt = new point[maxbar];
  hbars = new double[maxbar];
  L0 = new double[maxbar];
  F = new double[maxbar];
  Fvec = new point[maxbar];
  Ftot = new point[numpts];
  dist = new double[numpts];
  
  //*************************************************************************
  //*** Iterate until maxiterations reached or termination criterion met ****
  //*************************************************************************
  int counter = 0;  // ***************************
  int retris = 0;
  int numbar = 0;
  while (terminate == false  && counter < maxiters) //***********
    {							
      counter++;    // ***************************
      
      printf("  Iteration %i\n",counter);
      
      // retriangulate if necessary
      if (retriang == true)
	{
	  retris++; // ***************************
	  // if retriangulating, save old positions
	  for (int i=0; i<numpts; i++)
	    {
	      pold[i].x = p[i].x;
	      pold[i].y = p[i].y;
	      pold[i].z = p[i].z;
	    }
	  
	  FILE* pointfile = fopen("pointfile.txt","w");
	  fprintf(pointfile,"%8i\n%8i\n",3,numpts);
	  for (int i=0; i<numpts; i++)
	    { fprintf(pointfile,"%24.16e %24.16e %24.16e\n",p[i].x,p[i].y,p[i].z); }
	  fclose(pointfile);
	  
	  // System call for Delaunay triangulation
	  printf("    Retriangulating ... \n\n");
	  system ("$QHULL/bin/qdelaunay Qt i  < ./pointfile.txt > ./tetfile.txt");

	  // Read-in triangulation from qhull
	  FILE* tetfile = fopen("tetfile.txt","r");
	  int numtettemp;
	  fscanf(tetfile,"%i",&numtettemp);
	  tetra* ttemp = new tetra[numtettemp];
	  point* pmid = new point[numtettemp];

	  for (int i=0; i<numtettemp; i++)
	    {  fscanf(tetfile,"%i %i %i %i",&ttemp[i].n1,&ttemp[i].n2,&ttemp[i].n3,&ttemp[i].n4);  }
	  fclose(tetfile);

	  // Remove files
	  system("rm -f pointfile.txt tetfile.txt");
	  
	  retriang = false;
	  numtet = 0;
	  for (int i=0; i<numtettemp; i++)
	    {
	      pmid[i].x = (p[ttemp[i].n1].x + p[ttemp[i].n2].x + p[ttemp[i].n3].x + p[ttemp[i].n4].x) / 4.0e0;
	      pmid[i].y = (p[ttemp[i].n1].y + p[ttemp[i].n2].y + p[ttemp[i].n3].y + p[ttemp[i].n4].y) / 4.0e0;
	      pmid[i].z = (p[ttemp[i].n1].z + p[ttemp[i].n2].z + p[ttemp[i].n3].z + p[ttemp[i].n4].z) / 4.0e0;
	      if (SignedDistance(pmid[i]) < -geps)
		{
		  t[numtet].n1 = ttemp[i].n1;
		  t[numtet].n2 = ttemp[i].n2;
		  t[numtet].n3 = ttemp[i].n3;
		  t[numtet].n4 = ttemp[i].n4;
		  numtet++;
		}
	    }
	  
	  // find unique listing of bars
	  numbar = 0;
	  for (int i=0; i<numtet; i++)
	    {
	      bartemp[6*i].n1 = t[i].n1;
	      bartemp[6*i].n2 = t[i].n2;
	      bartemp[6*i+1].n1 = t[i].n1;
	      bartemp[6*i+1].n2 = t[i].n3;
	      bartemp[6*i+2].n1 = t[i].n1;
	      bartemp[6*i+2].n2 = t[i].n4;
	      bartemp[6*i+3].n1 = t[i].n2;
	      bartemp[6*i+3].n2 = t[i].n3;
	      bartemp[6*i+4].n1 = t[i].n2;
	      bartemp[6*i+4].n2 = t[i].n4;
	      bartemp[6*i+5].n1 = t[i].n3;
	      bartemp[6*i+5].n2 = t[i].n4;
	    }
	  for (int i=0; i<6*numtet; i++)
	    {
	      // sort each bar so that n1 < n2
	      if (bartemp[i].n1 > bartemp[i].n2)
		{
		  int temp = bartemp[i].n1;
		  bartemp[i].n1 = bartemp[i].n2;
		  bartemp[i].n2 = temp;
		}
	      // assign priority for sorting list of bars
	      priority[i] = numpts*bartemp[i].n1 + bartemp[i].n2;
	    }
	  QuickSort_bar(bartemp, priority, 0, 6*numtet - 1);
	  // keep only unique bars
	  int value = 0;
	  for (int i=0; i<6*numtet; i++)
	    if (priority[i] != value)
	      {
		bars[numbar].n1 = bartemp[i].n1;
		bars[numbar].n2 = bartemp[i].n2;
		value = priority[i];
		numbar++;
	      }
	}
      
      // ********************************************************************
      // **************** Move mesh points by computed forces ***************
      // ********************************************************************
      // create list of bar vectors "barvec" and bar lengths "L"
      double sumL3 = 0.0;
      double sumhbars3 = 0.0;
      for (int i=0; i<numbar; i++)
	{
	  barvec[i].x = p[bars[i].n1].x - p[bars[i].n2].x;
	  barvec[i].y = p[bars[i].n1].y - p[bars[i].n2].y;
	  barvec[i].z = p[bars[i].n1].z - p[bars[i].n2].z;
	  L[i] = sqrt((barvec[i].x*barvec[i].x + barvec[i].y*barvec[i].y
		       +barvec[i].z*barvec[i].z));
	  midpt[i].x = (p[bars[i].n1].x + p[bars[i].n2].x)/2.0;
	  midpt[i].y = (p[bars[i].n1].y + p[bars[i].n2].y)/2.0;
	  midpt[i].z = (p[bars[i].n1].z + p[bars[i].n2].z)/2.0;
	  hbars[i] = GridSpacing(midpt[i]);
	  sumL3 += (L[i]*L[i]*L[i]);
	  sumhbars3 += (hbars[i]*hbars[i]*hbars[i]);
	}
      for (int i=0; i<numbar; i++)
	{
	  L0[i] = hbars[i] * Fscale * pow((sumL3 / sumhbars3),(1.0/3.0));
	  F[i] = L0[i] - L[i];
	  if (F[i] < 0.0e0)
	    { F[i] = 0.0e0; }
	  Fvec[i].x = F[i]/L[i] * barvec[i].x;
	  Fvec[i].y = F[i]/L[i] * barvec[i].y;
	  Fvec[i].z = F[i]/L[i] * barvec[i].z;
	}	
      // find forces on each point (reset Ftot to 0, then add Fvec components)
      for (int i=0; i<numpts; i++)
        { 
	  Ftot[i].x = 0.0e0;
	  Ftot[i].y = 0.0e0;
	  Ftot[i].z = 0.0e0;
        }
      for (int i=0; i<numbar; i++)
	{
	  Ftot[bars[i].n1].x += Fvec[i].x;
	  Ftot[bars[i].n1].y += Fvec[i].y;
	  Ftot[bars[i].n1].z += Fvec[i].z;
	  Ftot[bars[i].n2].x -= Fvec[i].x;
	  Ftot[bars[i].n2].y -= Fvec[i].y;
	  Ftot[bars[i].n2].z -= Fvec[i].z;
	}
      // set forces on fixed points to 0
      if (numfixpts>0)
        {
	  for (int i=0; i<numfixpts; i++)
            {  
	      Ftot[i].x = 0.0e0;
	      Ftot[i].y = 0.0e0;
	      Ftot[i].z = 0.0e0;
            }
        }
      // update node positions
      for (int i=0; i<numpts; i++)
	{
	  p[i].x += deltat * Ftot[i].x;
	  p[i].y += deltat * Ftot[i].y;
	  p[i].z += deltat * Ftot[i].z;
	}

      // move outside points back to the boundary
      
      for (int i=0; i<numpts; i++)
	{
	  dist[i] = SignedDistance(p[i]);
	  if (dist[i] > 0.0e0)
	    {
	      point ptemp = p[i];
	      ptemp.x += deps;
	      double dgradx = (SignedDistance(ptemp)-SignedDistance(p[i])) / deps;
	      ptemp = p[i];
	      ptemp.y += deps;
	      double dgrady = (SignedDistance(ptemp)-SignedDistance(p[i])) / deps;
	      ptemp = p[i];
	      ptemp.z += deps;
	      double dgradz = (SignedDistance(ptemp)-SignedDistance(p[i])) / deps;
	      p[i].x -= dist[i]*dgradx;
	      p[i].y -= dist[i]*dgrady;
	      p[i].z -= dist[i]*dgradz;
	    }
	}
      
      // check to retriangulate
      double bigmove = 0.0;
      double bigmove_interior = 0.0;
      for (int i=0; i<numpts; i++)
	{
	  double pmove = sqrt(pow(p[i].x-pold[i].x,2) +  
			      pow(p[i].y-pold[i].y,2) + 
			      pow(p[i].z-pold[i].z,2));

	  if (pmove > bigmove)
	    { 
	      bigmove = pmove; 
	    }

	  if (pmove > bigmove_interior && dist[i]<-geps)
	    {
	      bigmove_interior = pmove;
	    }
	}
      retriang = false;
      if (bigmove > (ttol*h0))
	{
	  retriang = true;
	}
      
      // check to terminate      
      if (bigmove_interior < (dptol*h0))
	{
	  terminate = true;
	}
      
    }

  delete pold;
  delete ttemp;
  delete pmid;
  delete bartemp;
  delete bars;
  delete priority;
  delete barvec;
  delete L;
  delete midpt;
  delete hbars;
  delete L0;
  delete F;
  delete Fvec;
  delete Ftot;
  delete dist;
}
