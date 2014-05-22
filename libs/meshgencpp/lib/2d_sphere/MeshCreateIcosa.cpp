#include "meshdefs.h"

void MeshCreateIcosa(const double radius, int& numpts, int& numtri,
		     point*& p, triangle*& t)
{
  p = new point[12];
  t = new triangle[20];
  numpts = 12;
  numtri = 20;

  // Some constants
  const double GR = (1.0+sqrt(5.0))/2.0;
  double scale = 2.0*radius/sqrt(10.0+2.0*sqrt(5.0));

  // Points and elements for a regular icosahedron
  p[0].x = p[1].x = p[2].x = p[3].x = p[4].z = p[5].z = 0;
  p[6].z = p[7].z = p[8].y = p[9].y = p[10].y= p[11].y= 0;
  p[0].z = p[2].z = p[4].y = p[6].y = p[8].x = p[9].x = GR;
  p[1].z = p[3].z = p[5].y = p[7].y = p[10].x= p[11].x= -1*GR;
  p[0].y = p[1].y = p[4].x = p[5].x = p[8].z = p[10].z= 1;
  p[2].y = p[3].y = p[6].x = p[7].x = p[9].z = p[11].z= -1;
  
  t[0].n1  = t[1].n1  = t[2].n1  = t[3].n1  = t[4].n1  = 0;
  t[7].n1  = t[8].n1  = t[17].n1 = t[18].n1 = t[19].n1 = 1;
  t[0].n2  = t[4].n2  = t[12].n1 = t[13].n1 = t[14].n1 = 2;
  t[5].n1  = t[6].n1  = t[7].n2  = t[8].n2  = t[9].n1  = 3;
  t[1].n2  = t[2].n2  = t[16].n1 = t[17].n2 = t[18].n2 = 4;
  t[6].n2  = t[13].n2 = t[14].n2 = t[15].n1 = t[5].n2  = 5;
  t[2].n3  = t[3].n2  = t[10].n1 = t[18].n3 = t[19].n2 = 6;
  t[5].n3  = t[9].n2  = t[11].n1 = t[12].n2 = t[13].n3 = 7;
  t[0].n3  = t[1].n3  = t[14].n3 = t[15].n2 = t[16].n2 = 8;
  t[6].n3  = t[7].n3  = t[15].n3 = t[16].n3 = t[17].n3 = 9;
  t[3].n3  = t[4].n3  = t[10].n2 = t[11].n2 = t[12].n3 = 10;
  t[8].n3  = t[9].n3  = t[10].n3 = t[11].n3 = t[19].n3 = 11;

  // Scale points so that we have a sphere of radius "radius"
  for (int i=0; i<12; i++)
    {
      p[i].x = scale * p[i].x;
      p[i].y = scale * p[i].y;
      p[i].z = scale * p[i].z;
    }

  // Rotate icosahedron so that there is one point at each of the poles
  double ynew;
  double znew;
  double angle = asin(-radius*p[0].y/(pow(p[0].y,2) + pow(p[0].z,2)));

  for (int i=0; i<12; i++)
    {
      ynew =  p[i].y * cos(angle) + p[i].z * sin(angle);
      znew = -p[i].y * sin(angle) + p[i].z * cos(angle);
      p[i].y = ynew;
      p[i].z = znew;
    }

}
