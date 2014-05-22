#include "meshdefs.h"
#include <cstdio>

// Read-in input data
void MeshInputData(double& h0, 
		   point& xyzmin, 
		   point& xyzmax, 
		   int& numfixpts, 
		   point*& fixpt,
		   int& maxiters,
		   char*& GridType)
{
  FILE* infile = fopen("input3D.data","r");
  char buffer[256];

  fscanf(infile,"%s",GridType);  
  fgets(buffer, sizeof buffer, infile);
  fscanf(infile,"%lf",&h0); 
  fgets(buffer, sizeof buffer, infile);
  fscanf(infile,"%i",&maxiters);
  fgets(buffer, sizeof buffer, infile);

  fgets(buffer, sizeof buffer, infile);
  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%lf",&xyzmin.x);
  fgets(buffer, sizeof buffer, infile);
  fscanf(infile,"%lf",&xyzmax.x);
  fgets(buffer, sizeof buffer, infile);
  fscanf(infile,"%lf",&xyzmin.y);
  fgets(buffer, sizeof buffer, infile);
  fscanf(infile,"%lf",&xyzmax.y);
  fgets(buffer, sizeof buffer, infile);
  fscanf(infile,"%lf",&xyzmin.z);
  fgets(buffer, sizeof buffer, infile);
  fscanf(infile,"%lf",&xyzmax.z);
  fgets(buffer, sizeof buffer, infile);

  fgets(buffer, sizeof buffer, infile);

  fscanf(infile,"%i",&numfixpts);
  fgets(buffer, sizeof buffer, infile);

  fgets(buffer, sizeof buffer, infile);
  fgets(buffer, sizeof buffer, infile);

  fixpt = new point[numfixpts];
  for (int i = 0; i<numfixpts; i++)
    {
      fscanf(infile,"%lf %lf %lf",&fixpt[i].x,&fixpt[i].y,&fixpt[i].z);
    }
  fclose(infile);

  return;
}
