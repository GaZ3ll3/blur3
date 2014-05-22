#include "meshdefs.h"

void MeshComputeAreas(const int numtri, const int numpts, 
		      const point p[], triangle t[], 
		      double*& area, double*& cdual)
{
  // **************************************************************************
  // ***************** get triangle areas and dual mesh areas *****************
  // **************************************************************************
  for (int i=0; i<numpts; i++)
    {  cdual[i] = 0;  }

  double A,B,C;  
  double xt,yt,zt,dot_prod;
  int tmp_switch;
  
  for (int i=0; i<numtri; i++)
    {
      A = p[t[i].n1].z*(p[t[i].n2].y-p[t[i].n3].y)+p[t[i].n2].z*(p[t[i].n3].y-p[t[i].n1].y)
	+p[t[i].n3].z*(p[t[i].n1].y-p[t[i].n2].y);

      B = p[t[i].n1].x*(p[t[i].n2].z-p[t[i].n3].z)+p[t[i].n2].x*(p[t[i].n3].z-p[t[i].n1].z)
	+p[t[i].n3].x*(p[t[i].n1].z-p[t[i].n2].z);

      C = p[t[i].n1].y*(p[t[i].n2].x-p[t[i].n3].x)+p[t[i].n2].y*(p[t[i].n3].x-p[t[i].n1].x)
	+p[t[i].n3].y*(p[t[i].n1].x-p[t[i].n2].x);
      
      area[i] = 0.5 * sqrt(A*A + B*B + C*C);
      cdual[t[i].n1] += (area[i]/3.0);
      cdual[t[i].n2] += (area[i]/3.0);
      cdual[t[i].n3] += (area[i]/3.0);
      
      xt = (p[t[i].n1].x + p[t[i].n2].x + p[t[i].n3].x)/3.0;
      yt = (p[t[i].n1].y + p[t[i].n2].y + p[t[i].n3].y)/3.0;
      zt = (p[t[i].n1].z + p[t[i].n2].z + p[t[i].n3].z)/3.0;
     
      dot_prod = xt*A + yt*B + zt*C;

      if (dot_prod > 0.0) 
	{ 
	  tmp_switch = t[i].n2;
	  t[i].n2 = t[i].n3;
	  t[i].n3 = tmp_switch;
	}      
    }

}
