#include "meshdefs.h"
#include "mesh.h"

//
// Add ghost cells
//
void MeshAddGhostCells(char*& GridType, 
		       cart2d& cart2dinfo,
		       int& numpts, 
		       int& numtri, 
		       int& numghost, 
		       int& num_ext_node,
		       point*& p, 
		       triangle*& t, 
		       double*& area, 
		       double*& cdual, 
		       int*& ghost_link, 
		       int*& proper_ghostcell,
		       int*& ext_node_link,
		       double (*SignedDistance)(point))
{
  // Temporary storage
  point pnew[2*numpts];
  triangle tnew[2*numtri];
  double area_new[2*numtri];
  double cdual_new[2*numpts];
  int ghost_link_new[2*numtri];
  int numghost_new = 0;
  int numpts_new = numpts;
  int numtri_new = numtri;
  int ext_new[2*numtri];
  int num_ext = 0;
  
  // Compute the minimal distance parameter hmin
  double min_parea = 1.0e10;
  for (int i=0; i<numtri; i++)
    {
      if (min_parea > area[i])
        { min_parea = area[i]; }
    }
  double hmin = sqrt(2.0*min_parea)/5.0e0;
  
  // Store old t,area,p,cdual in new versions of each
  for (int i=0; i<numtri; i++)
    {
      tnew[i].n1 = t[i].n1;
      tnew[i].n2 = t[i].n2;
      tnew[i].n3 = t[i].n3;
      
      area_new[i] = area[i];
    }
  
  for (int i=0; i<numpts; i++)
    {
      pnew[i].x = p[i].x;
      pnew[i].y = p[i].y;
      
      cdual_new[i] = cdual[i];
    }
  
  // Loop over every element in order to determine
  // if a ghost cell should be created
  for (int i=0; i<numtri; i++)
    {
      point A,B,C,Anew;
      int mfound = 0;
      int jstore[4];
      int jopp;
      jstore[1] = 0;
      jstore[2] = 0;
      jstore[3] = 0;
      
      // Look at ith element
      int tmp[4];
      tmp[1] = t[i].n1;
      tmp[2] = t[i].n2;
      tmp[3] = t[i].n3;
      
      // Check how many nodes of ith element lie on the boundary (either 0, 1, 2, or 3)
      // If there are 0 or 1 nodes on the bnd, then do nothing and proceed to next element
      for (int j=1; j<=3; j++)
        {
	  point ptmp;
	  ptmp.x = p[tmp[j]].x;
	  ptmp.y = p[tmp[j]].y;
	  
	  if (fabs(SignedDistance(ptmp))<=hmin)
            {
	      mfound = mfound + 1;
	      jstore[mfound] = j;
            }
	  else
            {
	      A.x = ptmp.x;
	      A.y = ptmp.y;
	      jopp = tmp[j];
            }
        }    
      
      // Exactly three nodes of a triangle are on the boundary
      // NOTE: 3 nodes on boundary <==> exactly two edges are boundary edges
      //   (assuming the mesh has more than 1 element and these elements are connected)
      if (mfound==3)
        {
	  int kcheck = 0;
	  
	  // check if edge 12 is on the boundary
	  point ptmp;
	  ptmp.x = 0.5*( p[tmp[1]].x + p[tmp[2]].x );
	  ptmp.y = 0.5*( p[tmp[1]].y + p[tmp[2]].y );

	  if (fabs(SignedDistance(ptmp))<=hmin)
            {
	      kcheck=kcheck+1;
	      
	      jstore[1] = 1;
	      jstore[2] = 2;
	      jstore[3] = 3;
	      
	      A.x = p[tmp[jstore[3]]].x;
	      A.y = p[tmp[jstore[3]]].y;
	      
	      B.x = p[tmp[jstore[1]]].x;
	      B.y = p[tmp[jstore[1]]].y;

	      C.x = p[tmp[jstore[2]]].x;
	      C.y = p[tmp[jstore[2]]].y;
	      	      
	      Anew.x = B.x + C.x - A.x;
	      Anew.y = B.y + C.y - A.y;
	      
	      // Add new node and corresponding new element
	      numtri_new = numtri_new+1;
	      numpts_new = numpts_new+1;
	      
	      pnew[numpts_new-1].x = Anew.x;
	      pnew[numpts_new-1].y = Anew.y;
	      
	      tnew[numtri_new-1].n1 = tmp[jstore[1]];
	      tnew[numtri_new-1].n2 = tmp[jstore[2]];
	      tnew[numtri_new-1].n3 = numpts_new-1;
	      
	      // New node link to old node
	      num_ext = num_ext+1;
	      ext_new[num_ext-1] = tmp[jstore[3]];

	      // Add ghost index;
	      numghost_new = numghost_new+1;
	      ghost_link_new[numghost_new-1] = i;
	      
	      // Compute signed element area, change orientation if necessary
	      area_new[numtri_new-1] = .5*(pnew[tnew[numtri_new-1].n1].x * 
					   (pnew[tnew[numtri_new-1].n2].y
					    -pnew[tnew[numtri_new-1].n3].y) + 
					   pnew[tnew[numtri_new-1].n2].x * 
					   (pnew[tnew[numtri_new-1].n3].y
					    -pnew[tnew[numtri_new-1].n1].y) +
					   pnew[tnew[numtri_new-1].n3].x * 
					   (pnew[tnew[numtri_new-1].n1].y
					    -pnew[tnew[numtri_new-1].n2].y));
	      
	      if (area_new[numtri_new-1] < 0.0)
                {
		  int temp = tnew[numtri_new-1].n2;
		  tnew[numtri_new-1].n2 = tnew[numtri_new-1].n3;
		  tnew[numtri_new-1].n3 = temp;
		  area_new[numtri_new-1] = fabs(area_new[numtri_new-1]);
                }
	      
	      // Add element area to the appropriate dual cell area
	      cdual_new[tnew[numtri_new-1].n1] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n2] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n3] += (area_new[numtri_new-1]/3.0);
            }
	  
	  // check if edge 13 is on the boundary
	  ptmp.x = 0.5*( p[tmp[1]].x + p[tmp[3]].x );
	  ptmp.y = 0.5*( p[tmp[1]].y + p[tmp[3]].y );
	  
	  if (fabs(SignedDistance(ptmp))<=hmin)
            {
	      kcheck=kcheck+1;
	      
	      jstore[1] = 3;
	      jstore[2] = 1;
	      jstore[3] = 2;
	      
	      A.x = p[tmp[jstore[3]]].x;
	      A.y = p[tmp[jstore[3]]].y;
	      
	      B.x = p[tmp[jstore[1]]].x;
	      B.y = p[tmp[jstore[1]]].y;
	      
	      C.x = p[tmp[jstore[2]]].x;
	      C.y = p[tmp[jstore[2]]].y;

	      Anew.x = B.x + C.x - A.x;
	      Anew.y = B.y + C.y - A.y;
	      
	      // Add new node and corresponding new element
	      numtri_new = numtri_new+1;
	      numpts_new = numpts_new+1;
	      
	      pnew[numpts_new-1].x = Anew.x;
	      pnew[numpts_new-1].y = Anew.y;
	      
	      tnew[numtri_new-1].n1 = tmp[jstore[1]];
	      tnew[numtri_new-1].n2 = tmp[jstore[2]];
	      tnew[numtri_new-1].n3 = numpts_new-1;
	      
	      // New node link to old node
	      num_ext = num_ext+1;
	      ext_new[num_ext-1] = tmp[jstore[3]];
	      
	      // Add ghost index;
	      numghost_new = numghost_new+1;
	      ghost_link_new[numghost_new-1] = i;
	      
	      // Compute signed element area, change orientation if necessary
	      area_new[numtri_new-1] = .5*(pnew[tnew[numtri_new-1].n1].x * 
					   (pnew[tnew[numtri_new-1].n2].y
					    -pnew[tnew[numtri_new-1].n3].y) + 
					   pnew[tnew[numtri_new-1].n2].x * 
					   (pnew[tnew[numtri_new-1].n3].y
					    -pnew[tnew[numtri_new-1].n1].y) +
					   pnew[tnew[numtri_new-1].n3].x * 
					   (pnew[tnew[numtri_new-1].n1].y
					    -pnew[tnew[numtri_new-1].n2].y));
	      
	      if (area_new[numtri_new-1] < 0.0)
                {
		  int temp = tnew[numtri_new-1].n2;
		  tnew[numtri_new-1].n2 = tnew[numtri_new-1].n3;
		  tnew[numtri_new-1].n3 = temp;
		  area_new[numtri_new-1] = fabs(area_new[numtri_new-1]);
                }
	      
	      // Add element area to the appropriate dual cell area
	      cdual_new[tnew[numtri_new-1].n1] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n2] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n3] += (area_new[numtri_new-1]/3.0);
            }
	  
	  // check if edge 23 is on the boundary
	  ptmp.x = 0.5*( p[tmp[2]].x + p[tmp[3]].x );
	  ptmp.y = 0.5*( p[tmp[2]].y + p[tmp[3]].y );
	  
	  if (fabs(SignedDistance(ptmp))<=hmin)
            {
	      kcheck=kcheck+1;
	      
	      jstore[1] = 2;
	      jstore[2] = 3;
	      jstore[3] = 1;

	      A.x = p[tmp[jstore[3]]].x;
	      A.y = p[tmp[jstore[3]]].y;
	      
	      B.x = p[tmp[jstore[1]]].x;
	      B.y = p[tmp[jstore[1]]].y;
	      
	      C.x = p[tmp[jstore[2]]].x;
	      C.y = p[tmp[jstore[2]]].y;

	      Anew.x = B.x + C.x - A.x;
	      Anew.y = B.y + C.y - A.y;
	      
	      // Add new node and corresponding new element
	      numtri_new = numtri_new+1;
	      numpts_new = numpts_new+1;
	      
	      pnew[numpts_new-1].x = Anew.x;
	      pnew[numpts_new-1].y = Anew.y;
	      
	      tnew[numtri_new-1].n1 = tmp[jstore[1]];
	      tnew[numtri_new-1].n2 = tmp[jstore[2]];
	      tnew[numtri_new-1].n3 = numpts_new-1;
	      
	      // New node link to old node
	      num_ext = num_ext+1;
	      ext_new[num_ext-1] = tmp[jstore[3]];
	      
	      // Add ghost index;
	      numghost_new = numghost_new+1;
	      ghost_link_new[numghost_new-1] = i;
	      
	      // Compute signed element area, change orientation if necessary
	      area_new[numtri_new-1] = .5*(pnew[tnew[numtri_new-1].n1].x * 
					   (pnew[tnew[numtri_new-1].n2].y
					    -pnew[tnew[numtri_new-1].n3].y) + 
					   pnew[tnew[numtri_new-1].n2].x * 
					   (pnew[tnew[numtri_new-1].n3].y
					    -pnew[tnew[numtri_new-1].n1].y) +
					   pnew[tnew[numtri_new-1].n3].x * 
					   (pnew[tnew[numtri_new-1].n1].y
					    -pnew[tnew[numtri_new-1].n2].y));
	      
	      if (area_new[numtri_new-1] < 0.0)
                {
		  int temp = tnew[numtri_new-1].n2;
		  tnew[numtri_new-1].n2 = tnew[numtri_new-1].n3;
		  tnew[numtri_new-1].n3 = temp;
		  area_new[numtri_new-1] = fabs(area_new[numtri_new-1]);
                }

	      // Add element area to the appropriate dual cell area
	      cdual_new[tnew[numtri_new-1].n1] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n2] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n3] += (area_new[numtri_new-1]/3.0);
            }
	  
	  if (kcheck!=2)
            {
	      printf("\n");
	      printf(" ERROR in MeshAddGhostCells.cpp: \n");
	      printf("       Found a triangle with 3 nodes on the boundary.\n");
	      printf("       This triangle should have exactly 2 edges on the boundary,\n");
	      printf("         but in fact it has %i\n",kcheck);
	      printf("\n");
	      exit(1);
            }
        }
      // Exactly two nodes of a triangle are on the boundary
      // NOTE: 2 nodes on boundary  <==> either exactly one edge is a boundary edge  -or-
      //                                 near a corner and exactly zero edges are a bnd edge
      else if (mfound==2)
        {
	  
	  // Compute edge midpoint
	  point ptmp;
	  ptmp.x = 0.5*( p[tmp[jstore[1]]].x + p[tmp[jstore[2]]].x );
	  ptmp.y = 0.5*( p[tmp[jstore[1]]].y + p[tmp[jstore[2]]].y );
	  
	  // Proceed if edge midpoint is also on the bnd <==> exactly one edge is a bnd edge
	  if (fabs(SignedDistance(ptmp))<=hmin)
            {
	      if (jstore[1]>jstore[2])
                {
		  int tmps = jstore[2];
		  jstore[2] = jstore[1];
		  jstore[1] = tmps;
                }
	      	      
	      B.x = p[tmp[jstore[1]]].x;
	      B.y = p[tmp[jstore[1]]].y;
	      
	      C.x = p[tmp[jstore[2]]].x;
	      C.y = p[tmp[jstore[2]]].y;

	      Anew.x = B.x + C.x - A.x;
	      Anew.y = B.y + C.y - A.y;
	      
	      // Add new node and corresponding new element
	      numtri_new = numtri_new+1;
	      numpts_new = numpts_new+1;
	      
	      pnew[numpts_new-1].x = Anew.x;
	      pnew[numpts_new-1].y = Anew.y;
	      
	      tnew[numtri_new-1].n1 = tmp[jstore[1]];
	      tnew[numtri_new-1].n2 = tmp[jstore[2]];
	      tnew[numtri_new-1].n3 = numpts_new-1;
	      
	      // New node link to old node
	      num_ext = num_ext+1;
	      ext_new[num_ext-1] = jopp;// tmp[jstore[3]];
	      
	      // Add ghost index;
	      numghost_new = numghost_new+1;
	      ghost_link_new[numghost_new-1] = i;

	      // Compute signed element area, change orientation if necessary
	      area_new[numtri_new-1] = .5*(pnew[tnew[numtri_new-1].n1].x * 
					   (pnew[tnew[numtri_new-1].n2].y
					    -pnew[tnew[numtri_new-1].n3].y) + 
					   pnew[tnew[numtri_new-1].n2].x * 
					   (pnew[tnew[numtri_new-1].n3].y
					    -pnew[tnew[numtri_new-1].n1].y) +
					   pnew[tnew[numtri_new-1].n3].x * 
					   (pnew[tnew[numtri_new-1].n1].y
					    -pnew[tnew[numtri_new-1].n2].y));
	      
	      if (area_new[numtri_new-1] < 0.0)
                {
		  int temp = tnew[numtri_new-1].n2;
		  tnew[numtri_new-1].n2 = tnew[numtri_new-1].n3;
		  tnew[numtri_new-1].n3 = temp;
		  area_new[numtri_new-1] = fabs(area_new[numtri_new-1]);
                }
	      
	      // Add element area to the appropriate dual cell area
	      cdual_new[tnew[numtri_new-1].n1] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n2] += (area_new[numtri_new-1]/3.0);
	      cdual_new[tnew[numtri_new-1].n3] += (area_new[numtri_new-1]/3.0);
            }
        }
    }
  
  // Resize all appropriate vectors
  delete[] p;
  delete[] t;
  delete[] area;
  delete[] cdual;
  
  numpts = numpts_new;
  numtri = numtri_new;
  numghost = numghost_new;
  num_ext_node = num_ext;
  
  p = new point[numpts];
  t = new triangle[numtri];
  area = new double[numtri];
  cdual = new double[numpts];
  ghost_link = new int[numghost];
  proper_ghostcell = new int[numghost];
  ext_node_link = new int[num_ext];
  
  for (int i=0; i<numtri; i++)
    {
      t[i].n1 = tnew[i].n1;
      t[i].n2 = tnew[i].n2;
      t[i].n3 = tnew[i].n3;
      
      area[i] = area_new[i];
    }
  
  for (int i=0; i<numpts; i++)
    {
      p[i].x = pnew[i].x;
      p[i].y = pnew[i].y;
      
      cdual[i] = cdual_new[i];
    }
  
  for (int i=0; i<numghost; i++)
    {
      ghost_link[i] = ghost_link_new[i];
    }
  
  for (int i=0; i<numghost; i++)
    {
      proper_ghostcell[i] = 1;
    }
  
  for (int i=0; i<num_ext; i++)
    {
      ext_node_link[i] = ext_new[i];
    }
}
