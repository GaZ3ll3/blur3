#include "meshdefs.h"
#include <cstdio>

// Read-in input data
void MeshInputData(double& h0, 
		   point& xymin, 
		   point& xymax, 
		   int& numfixpts, 
		   point*& fixpt, 
		   int& maxiters, 
		   int& sub_factor, 
		   char*& GridType)
{
  FILE* infile = fopen("input2D.data","r");  
  char buffer[256];

  fscanf(infile,"%s",GridType);  
  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%lf",&h0);
  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%i",&maxiters);
  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%i",&sub_factor);
  fgets(buffer, sizeof buffer, infile);
  fgets(buffer, sizeof buffer, infile);
  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%lf",&xymin.x);
  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%lf",&xymax.x);
  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%lf",&xymin.y);
  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%lf",&xymax.y);
  fgets(buffer, sizeof buffer, infile);
  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%i",&numfixpts);
  fgets(buffer, sizeof buffer, infile);
  fgets(buffer, sizeof buffer, infile);
  fgets(buffer, sizeof buffer, infile);

  fixpt = new point[numfixpts];
  for (int i = 0; i<numfixpts; i++)
    {  
      fscanf(infile,"%lf %lf",&fixpt[i].x,&fixpt[i].y);
    }

  fclose(infile);
}
