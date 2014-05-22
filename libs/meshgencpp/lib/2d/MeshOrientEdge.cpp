#include "meshdefs.h"
#include "mesh.h"

//
// Adjust edge information so that the element to the "left"
// of the edge is always "upwind" of the unit normal to the
// edge and the element to the "right" of the edge is always
// "downwind" of the unit normal to the edge
//
void MeshOrientEdge(mesh& Mesh)
{
  // Reshuffle edges as they are labeled for each element
  // I want the following setup for each element:
  //
  //                3
  //                /\
  //               /  \
  //              /    \
  //          2  /      \ 1
  //            /        \
  //           /    i     \
  //          /            \
  //       1 /--------------\ 2
  //                3
  //
  for (int i=1; i<=Mesh.get_NumPhysElems(); i++)
    {
      int e1_new=-1;
      int e2_new=-1;
      int e3_new=-1;

      int e1 = Mesh.get_tedge(i,1);
      int e2 = Mesh.get_tedge(i,2);
      int e3 = Mesh.get_tedge(i,3);

      int n1_1 = Mesh.get_enode(e1,1);
      int n1_2 = Mesh.get_enode(e1,2);

      int n2_1 = Mesh.get_enode(e2,1);
      int n2_2 = Mesh.get_enode(e2,2);

      int n3_1 = Mesh.get_enode(e3,1);
      int n3_2 = Mesh.get_enode(e3,2);

      int node1 = Mesh.get_tnode(i,1);
      int node2 = Mesh.get_tnode(i,2);
      int node3 = Mesh.get_tnode(i,3);

      int numcase = 10000;
      if (node1!=n1_1 && node1!=n1_2)
	{ numcase = 1; }
      else
	{
	  if (node2!=n1_1 && node2!=n1_2)
	    { numcase = 2; }
	  else
	    {
	      if (node3!=n1_1 && node3!=n1_2)
		{ numcase = 3; } 
	    }
	}

      switch(numcase)
	{
	case 1:

	  e1_new = e1;
	  if (node2!=n2_1 && node2!=n2_2)
	    {
	      e2_new = e2;
	      e3_new = e3;
	    }
	  else
	    {
	      e2_new = e3;
	      e3_new = e2;
	    }
	  
	  break;

	case 2:

	  e1_new = e2;
	  if (node1!=n2_1 && node1!=n2_2)
	    {
	      e2_new = e1;
	      e3_new = e3;
	    }
	  else
	    {
	      e2_new = e3;
	      e3_new = e1;
	    }

	  break;

	case 3:

	  e1_new = e3;
	  if (node1!=n2_1 && node1!=n2_2)
	    {
	      e2_new = e1;
	      e3_new = e2;
	    }
	  else
	    {
	      e2_new = e2;
	      e3_new = e1;
	    }

	  break;

	default:
	  printf("\n");
	  printf(" Error in MeshOrientEdge.cpp \n");
	  printf("    e1 = %i \n",e1);
	  printf("    e2 = %i \n",e2);
	  printf("    e3 = %i \n",e3);
	  printf("  n1_1 = %i \n",n1_1);
	  printf("  n1_2 = %i \n",n1_2);
	  printf("  n2_1 = %i \n",n2_1);
	  printf("  n2_2 = %i \n",n2_2);
	  printf("  n3_1 = %i \n",n3_1);
	  printf("  n3_2 = %i \n",n3_2);
	  printf(" node1 = %i \n",node1);
	  printf(" node2 = %i \n",node2);
	  printf(" node3 = %i \n",node3);
	  printf(" \n");
	  exit(1);
	}
      
      Mesh.set_tedge(i,1, e1_new );
      Mesh.set_tedge(i,2, e2_new );
      Mesh.set_tedge(i,3, e3_new );
  
    }

  int count1 = 0;
  int count2 = 0;

  // Loop over all Edges
  for (int i=1; i<=Mesh.get_NumEdges(); i++)
    {
      // Edge coordinates
      double x1 = Mesh.get_edge(i,1);
      double y1 = Mesh.get_edge(i,2);
      double x2 = Mesh.get_edge(i,3);
      double y2 = Mesh.get_edge(i,4);
      
      // Edge length
      double L = sqrt(pow(x2-x1,2) + pow(y2-y1,2));

      // Unit normal to edge
      dTensor1 nhat(2);      
      nhat.set(1, (y2-y1)/L );
      nhat.set(2, (x1-x2)/L );

      // Nodes on either side of edge
      int ileft  = Mesh.get_eelem(i,1);
      int iright = Mesh.get_eelem(i,2);
      
      double xc_l = (Mesh.get_node(Mesh.get_tnode(ileft,1),1)
		   + Mesh.get_node(Mesh.get_tnode(ileft,2),1)
		   + Mesh.get_node(Mesh.get_tnode(ileft,3),1))/3.0;
      
      double yc_l = (Mesh.get_node(Mesh.get_tnode(ileft,1),2)
		   + Mesh.get_node(Mesh.get_tnode(ileft,2),2)
		   + Mesh.get_node(Mesh.get_tnode(ileft,3),2))/3.0;

      double xc_r = (Mesh.get_node(Mesh.get_tnode(iright,1),1)
		   + Mesh.get_node(Mesh.get_tnode(iright,2),1)
		   + Mesh.get_node(Mesh.get_tnode(iright,3),1))/3.0;
      
      double yc_r = (Mesh.get_node(Mesh.get_tnode(iright,1),2)
		   + Mesh.get_node(Mesh.get_tnode(iright,2),2)
		   + Mesh.get_node(Mesh.get_tnode(iright,3),2))/3.0;

      dTensor1 rhat(2);      
      rhat.set(1, xc_r - xc_l );
      rhat.set(2, yc_r - yc_l );

      double ddot = nhat.get(1)*rhat.get(1) + nhat.get(2)*rhat.get(2);

      if (ddot<0.0) // Then need to interchange nodes
	{
	  count2 = count2 + 1;

	  // Switch element order
	  Mesh.set_eelem(i,1,iright);
	  Mesh.set_eelem(i,2,ileft);

	  // Then repeat test to make sure logic is sound

	  // Nodes on either side of edge
	  ileft  = Mesh.get_eelem(i,1);
	  iright = Mesh.get_eelem(i,2);
	  
	  xc_l = (Mesh.get_node(Mesh.get_tnode(ileft,1),1)
		  + Mesh.get_node(Mesh.get_tnode(ileft,2),1)
		  + Mesh.get_node(Mesh.get_tnode(ileft,3),1))/3.0;
	  
	  yc_l = (Mesh.get_node(Mesh.get_tnode(ileft,1),2)
		  + Mesh.get_node(Mesh.get_tnode(ileft,2),2)
		  + Mesh.get_node(Mesh.get_tnode(ileft,3),2))/3.0;
	  
	  xc_r = (Mesh.get_node(Mesh.get_tnode(iright,1),1)
		  + Mesh.get_node(Mesh.get_tnode(iright,2),1)
		  + Mesh.get_node(Mesh.get_tnode(iright,3),1))/3.0;
	  
	  yc_r = (Mesh.get_node(Mesh.get_tnode(iright,1),2)
		  + Mesh.get_node(Mesh.get_tnode(iright,2),2)
		  + Mesh.get_node(Mesh.get_tnode(iright,3),2))/3.0;
	  
	  rhat.set(1, xc_r - xc_l );
	  rhat.set(2, yc_r - yc_l );

	  double ddot_new = nhat.get(1)*rhat.get(1) + nhat.get(2)*rhat.get(2);
	  
	  if (ddot_new<0.0) // this should never happen
	    {
	      printf("\n");
	      printf(" ERROR in MeshOrientEdge.cpp: ddot_new < 0 \n");
	      printf("\n");
	      exit(1);
	    }		
	}
      else
	{  count1 = count1 + 1;  }
    }

  // Store mesh orientation:
  //    if element "j" is on the upwind side of the
  //    unit normal vector to edge "e", then
  //        tedge_orientation(j) = 1
  //    else
  //        tedge_orientation(j) = -1
  //
  for (int j=1; j<=Mesh.get_NumPhysElems(); j++)
    for (int k=1; k<=3; k++)
      {
	int e = Mesh.get_tedge(j,k);
	
	if (j==Mesh.get_eelem(e,1))
	  {
	    Mesh.set_tedge_orientation(j,k, 1 );
	  }
	else
	  {
	    if (j==Mesh.get_eelem(e,2))
	      {
		Mesh.set_tedge_orientation(j,k, -1 );
	      }
	    else
	      {
		printf("\n");
		printf(" Error in MeshOrientEdge.cpp\n");
		printf("                   j = %i\n",j);
		printf("                   k = %i\n",k);
		printf("                   e = %i\n",e);
		printf(" Mesh.get_eelem(e,1) =  %i\n",Mesh.get_eelem(e,1));
		printf(" Mesh.get_eelem(e,2) =  %i\n",Mesh.get_eelem(e,2));
		printf("\n");
	      }
	  }
      }

}
