#include "meshdefs.h"
#include <cstdio>
#include <cstdlib>

//
//  Takes initial point distribution "p", triangulates with "t", and then
//  iterates on the point distribution and triangulation until an equilibrium
//  mesh is created (or until the number of iterations exceeds "maxiters").
//
void MeshIterate(int maxiters, double geps, double ttol, double dptol, double Fscale, 
		 double deltat, double deps, double h0, int numfixpts, int& numpts, int& numtri, 
		 point*& p, triangle*& t,double (*SignedDistance)(point),
		 double (*GridSpacing)(point))
{
  // Sorting function
  void QuickSort_bar(bar* b, int* key, int i, int j);
  
  // Set up all needed variables and dynamic arrays for while-loop iteration
  bool retriang = true;	      // force initial retriangulation
  bool terminate = false;     // end when termination criterion satisfied
  
  // create and allocate dynamic arrays
  point* pold;
  triangle* ttemp;
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

  pold = new point[numpts];
  ttemp = new triangle[2*numpts];
  pmid = new point[2*numpts];
  t = new triangle[2*numpts];
  bartemp = new bar[6*numpts];
  bars = new bar[6*numpts];
  priority = new int[6*numpts];
  barvec = new point[6*numpts];
  L = new double[6*numpts];
  midpt = new point[6*numpts];
  hbars = new double[6*numpts];
  L0 = new double[6*numpts];
  F = new double[6*numpts];
  Fvec = new point[6*numpts];
  Ftot = new point[numpts];
  dist = new double[numpts];

  //*************************************************************************
  //*** Iterate until maxiterations reached or termination criterion met ****
  //*************************************************************************
  int counter = 0;
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
	    }	  
	  	  
	  FILE* pointfile = fopen("pointfile.txt","w");
	  fprintf(pointfile,"%8i\n%8i\n",2,numpts);
	  for (int i=0; i<numpts; i++)
	    { fprintf(pointfile,"%24.16e %24.16e\n",p[i].x,p[i].y); }
	  fclose(pointfile);
	  
	  // System call for Delaunay triangulation
	  printf("    Retriangulating ... \n\n");
	  system ("$QHULL/bin/qdelaunay Qt i  < ./pointfile.txt > ./trifile.txt");
	  
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

	  retriang = false;
	  numtri = 0;
	  for (int i=0; i<numtritemp; i++)
	    {
	      pmid[i].x = (p[ttemp[i].n1].x + p[ttemp[i].n2].x + p[ttemp[i].n3].x) / 3;
	      pmid[i].y = (p[ttemp[i].n1].y + p[ttemp[i].n2].y + p[ttemp[i].n3].y) / 3;
	      if (SignedDistance(pmid[i]) < -geps)
		{
		  t[numtri].n1 = ttemp[i].n1;
		  t[numtri].n2 = ttemp[i].n2;
		  t[numtri].n3 = ttemp[i].n3;
		  numtri++;
		}
	    }
	  
	  // find unique listing of bars
	  numbar = 0;
	  for (int i=0; i<numtri; i++)
	    {
	      bartemp[3*i].n1 = t[i].n1;
	      bartemp[3*i].n2 = t[i].n2;
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
		  int temp = bartemp[i].n1;
		  bartemp[i].n1 = bartemp[i].n2;
		  bartemp[i].n2 = temp;
		}
	      // assign priority for sorting list of bars
	      priority[i] = numpts*bartemp[i].n1 + bartemp[i].n2;
	    }
	  
	  QuickSort_bar(bartemp, priority, 0, 3*numtri - 1);

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
	}
      
      // ********************************************************************
      // **************** Move mesh points by computed forces ***************
      // ********************************************************************
      // create list of bar vectors "barvec" and bar lengths "L"
      double sumL2 = 0.0;
      double sumhbars2 = 0.0;
      for (int i=0; i<numbar; i++)
	{
	  barvec[i].x = p[bars[i].n1].x - p[bars[i].n2].x;
	  barvec[i].y = p[bars[i].n1].y - p[bars[i].n2].y;
	  L[i] = sqrt(static_cast<double>(barvec[i].x*barvec[i].x + barvec[i].y*barvec[i].y));
	  midpt[i].x = (p[bars[i].n1].x + p[bars[i].n2].x)/2.0;
	  midpt[i].y = (p[bars[i].n1].y + p[bars[i].n2].y)/2.0;
	  hbars[i] = GridSpacing(midpt[i]);
	  sumL2 += (L[i]*L[i]);
	  sumhbars2 += (hbars[i]*hbars[i]);
	}
      for (int i=0; i<numbar; i++)
	{
	  L0[i] = hbars[i] * Fscale * sqrt(static_cast<double>(sumL2 / sumhbars2));
	  F[i] = L0[i] - L[i];
	  if (F[i] < 0)
	    F[i] = 0;
	  Fvec[i].x = F[i]/L[i] * barvec[i].x;
	  Fvec[i].y = F[i]/L[i] * barvec[i].y;
	}	

      // find forces on each point (reset Ftot to 0, then add Fvec components)
      for (int i=0; i<numpts; i++)
	Ftot[i].x = Ftot[i].y = 0;
      for (int i=0; i<numbar; i++)
	{
	  Ftot[bars[i].n1].x += Fvec[i].x;
	  Ftot[bars[i].n1].y += Fvec[i].y;
	  Ftot[bars[i].n2].x -= Fvec[i].x;
	  Ftot[bars[i].n2].y -= Fvec[i].y;
	}
      // set forces on fixed points to 0
      for (int i=0; i<numfixpts; i++)
	Ftot[i].x = Ftot[i].y = 0;
      // update node positions
      for (int i=0; i<numpts; i++)
	{
	  p[i].x += deltat * Ftot[i].x;
	  p[i].y += deltat * Ftot[i].y;
	}
      
      // move outside points back to the boundary
      for (int i=0; i<numpts; i++)
	{
	  dist[i] = SignedDistance(p[i]);
	  if (dist[i] > 0)
	    {
	      int mm=0;
	      double tmp = fabs(dist[i]);
	      while (tmp>1.0e-14)
		{
		  point ptemp = p[i];
		  ptemp.x += deps;
		  double dgradx = (SignedDistance(ptemp)-SignedDistance(p[i])) / deps;
		  ptemp = p[i];
		  ptemp.y += deps;
		  double dgrady = (SignedDistance(ptemp)-SignedDistance(p[i])) / deps;
		  double norm2 = pow(dgradx,2)+pow(dgrady,2);
		  p[i].x -= dist[i]*dgradx/norm2;
		  p[i].y -= dist[i]*dgrady/norm2;
		  dist[i] = SignedDistance(p[i]);
		  tmp = fabs(dist[i]);
		  mm=mm+1;
		}
	    }	  
	}

      // check to retriangulate
      double bigmove = 0.0;
      for (int i=0; i<numpts; i++)
	{
	  double pmove = sqrt(static_cast<double>((p[i].x-pold[i].x)*(p[i].x-pold[i].x) 
						  + (p[i].y-pold[i].y)*(p[i].y-pold[i].y)));
	  if (pmove > bigmove)
	    bigmove = pmove;
	}
      if ((bigmove / h0) > ttol)
	retriang = true;
      
      
      // check to terminate
      bigmove = 0.0;
      for (int i=0; i<numpts; i++)
	if (dist[i] < -geps)
	  {
	    double pmove = sqrt(static_cast<double>((deltat*Ftot[i].x *Ftot[i].x)
						    + (deltat*Ftot[i].y *Ftot[i].y)));
	    if (pmove > bigmove)
	      bigmove = pmove;
	  }
      if ((bigmove / h0) < dptol)
	terminate = true;
      
    }

   delete [] pmid;
   delete [] priority;
   delete [] bartemp;
   delete [] midpt;
   delete [] barvec;
   delete [] L;
   delete [] hbars;
   delete [] L0;
   delete [] F;
   delete [] bars;
   delete [] Fvec;
   delete [] Ftot;
   delete [] dist;

}
