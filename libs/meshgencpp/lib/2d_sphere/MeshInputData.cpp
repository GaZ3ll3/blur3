#include "meshdefs.h"
#include <cstdio>

// Read-in input data
void MeshInputData(double& radius, 
		   int& level)
{
  FILE* read = fopen("input.data","r");
  char buffer[256];

  fscanf(read,"%lf",&radius);  
  fgets(buffer, sizeof buffer, read);
  fscanf(read,"%i",&level);  
  fgets(buffer, sizeof buffer, read);
  fclose(read);

  return;
}
