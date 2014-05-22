#include "meshdefs.h"

// Store spherical coordinates of p in psphere
void MeshSphereCoords(const int numpts, const point* p, point*& psphere)
{

  for (int i=0; i<numpts; i++)
    {
      double x = p[i].x;
      double y = p[i].y;
      double z = p[i].z;

      double   rad = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
      double theta = acos(z/rad);
      double   phi = atan2(y,x);

      // Longitude
      psphere[i].x = phi - pi;

      // Latitude
      psphere[i].y = 0.5*pi - theta;
    }

}
