#include "meshdefs.h"
#include "dog_math.h"
#include <cstdlib>
#include <cstdio>
#include <assert.h>

//
// Output mesh information
//
void MeshOutput(int numtet, 
		int numtet_phys,
		int numghost,
		int numpts,
		int numpts_phys,
		int numbnd,
		int numfaces,
		point p[],
		tetra t[], 
		int bnd_node[],
		double volume[],
		double cdual[],
		char* outputdir)
{
  // run startscript
  // to create output directory if it does not exist
  // and copy data files to output directory
  char command_str[1024];
  int numchars = snprintf(command_str,1024,
			  "if test -f startscript && test -x startscript;\n"
			  "then ./startscript %s %d\n"
			  "else ${MESHGENCPP}/lib/3d/startscript %s %d\n"
			  "fi", outputdir, 3, outputdir, 3);
  assert(numchars<1023);
  assert(numchars>0);
  int exit_status = 0;
  exit_status = system(command_str);
  
  double minAngle = 180.0e0;
  double minvolume = volume[0];
  double maxvolume = volume[0];
  double totalvolume = volume[0];
  double angle1,angle2,angle3;
  int xv[4][3];
  double x1,y1,z1,x2,y2,z2,x3,y3,z3;
  double vax,vay,vaz,vam;
  double vbx,vby,vbz,vbm;
  double vcx,vcy,vcz,vcm;
  double vab,vac;
    
  // Find minimum angle
  for (int i=0; i<numtet; i++)
    {      
      xv[0][0] = t[i].n1;
      xv[0][1] = t[i].n2;
      xv[0][2] = t[i].n3;
      
      xv[1][0] = t[i].n1;
      xv[1][1] = t[i].n2;
      xv[1][2] = t[i].n4;
      
      xv[2][0] = t[i].n2;
      xv[2][1] = t[i].n3;
      xv[2][2] = t[i].n4;
      
      xv[3][0] = t[i].n1;
      xv[3][1] = t[i].n3;
      xv[3][2] = t[i].n4;
      
      for (int k=0; k<4; k++)
        {
	  x1 = p[xv[k][0]].x;
	  y1 = p[xv[k][0]].y;
	  z1 = p[xv[k][0]].z;
          
	  x2 = p[xv[k][1]].x;
	  y2 = p[xv[k][1]].y;
	  z2 = p[xv[k][1]].z;
          
	  x3 = p[xv[k][2]].x;
	  y3 = p[xv[k][2]].y;
	  z3 = p[xv[k][2]].z;
            
	  vax = x2-x1;
	  vay = y2-y1;
	  vaz = z2-z1;
          
	  vbx = x3-x1;
	  vby = y3-y1;
	  vbz = z3-z1;
            
	  vcx = x2-x3;
	  vcy = y2-y3;
	  vcz = z2-z3;
                    
	  vam = sqrt(vax*vax + vay*vay + vaz*vaz);
	  vbm = sqrt(vbx*vbx + vby*vby + vbz*vbz);
	  vcm = sqrt(vcx*vcx + vcy*vcy + vcz*vcz);
	  vab = fabs(vax*vbx + vay*vby + vaz*vbz);
	  vac = fabs(vax*vcx + vay*vcy + vaz*vcz);
          
	  angle1 = (180.0/pi)*acos(vab/(vam*vbm));
	  angle2 = (180.0/pi)*acos(vac/(vam*vcm));
	  angle3 = (180.0 - angle1 - angle2);         
	  minAngle = Min( minAngle, Min(Min(angle1,angle2),angle3));
        }
    }

  // Find min and max volumes
  for (int i=1; i<numtet; i++)
    {
      double tmp = volume[i];

      totalvolume = totalvolume + tmp;

      if (tmp>maxvolume)
	{ maxvolume = tmp; }
      else if (tmp<minvolume)
	{ minvolume = tmp; }
    }

  // Output data to files
  char fname[1024];
  snprintf(fname,1024,"%s/mesh_node.dat",outputdir);
  FILE* write1 = fopen(fname,"w");
  for (int i=0; i<(numpts-1); i++)
    {
      fprintf(write1,"%24.16e %24.16e %24.16e\n",
	      p[i].x,
	      p[i].y,
	      p[i].z);
    }
  fprintf(write1,"%24.16e %24.16e %24.16e",
	  p[numpts-1].x,
	  p[numpts-1].y,
	  p[numpts-1].z);
  fclose(write1);

  snprintf(fname,1024,"%s/mesh_tnode.dat",outputdir);
  FILE* write2 = fopen(fname,"w");
  for (int i=0; i<(numtet-1); i++)
    {
      fprintf(write2,"%8i %8i %8i %8i\n",
	      t[i].n1+1,
	      t[i].n2+1,
	      t[i].n3+1,
	      t[i].n4+1);
    }
  fprintf(write2,"%8i %8i %8i %8i",
	  t[numtet-1].n1+1,
	  t[numtet-1].n2+1,
	  t[numtet-1].n3+1,
	  t[numtet-1].n4+1);
  fclose(write2);
  
  snprintf(fname,1024,"%s/mesh_bnd_node.dat",outputdir);
  FILE* write3 = fopen(fname,"w");
  for (int i=0; i<(numbnd-1); i++)
    {
      fprintf(write3,"%8i\n",bnd_node[i]+1);
    }
  fprintf(write3,"%8i",bnd_node[numbnd-1]+1);
  fclose(write3);

  // Output summary of results to screen
  printf("\n");
  printf("  SUMMARY OF RESULTS:\n");
  printf("  -------------------\n");
  printf("          Number of Elements:  %8i\n",numtet);
  printf(" Number of Physical Elements:  %8i\n",numtet_phys);
  printf("    Number of Ghost Elements:  %8i\n",numghost);
  printf("             Number of Nodes:  %8i\n",numpts);
  printf("    Number of Physical Nodes:  %8i\n",numpts_phys);
  printf("    Number of Boundary Nodes:  %8i\n",numbnd);
  printf("             Number of Faces:  %8i\n",numfaces);
  printf("\n");
  printf("        Total Volume Covered:  %24.16e\n",totalvolume);
  printf("   Volume Ratio: small/large:  %24.16e\n",minvolume/maxvolume);
  printf("    Angle Ratio: minAngle/60:  %24.16e\n",minAngle/60.0);
  printf("\n");
}
