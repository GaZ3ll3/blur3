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

      // Check periodicity condition
      double LX = x1-x2;
      double LY = y2-y1;
      if ( LX>=pi )
	{  LX = LX - 2*pi;  }
      else if (LX<=-pi)
	{  LX = LX + 2*pi;  }

      if ( LY>=pi )
	{  LY = LY - 2*pi;  }
      else if (LY<=-pi)
	{  LY = LY + 2*pi;  }

      if (fabs(y1-0.5*pi)<=1.0e-12 || fabs(y1+0.5*pi)<=1.0e-12 ||
	  fabs(y2-0.5*pi)<=1.0e-12 || fabs(y2+0.5*pi)<=1.0e-12)
	{
	  LX = 0.0e0;
	}
	
      // Edge length
      double L = sqrt(pow(LX,2) + pow(LY,2));

      // Unit normal to edge
      dTensor1 nhat(2);      
      nhat.set(1, LY/L );
      nhat.set(2, LX/L );

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
      double RX = xc_r - xc_l;
      double RY = yc_r - yc_l;
      
      if ( RX>=pi )
	{  RX = RX - 2*pi;  }
      else if (RX<=-pi)
	{  RX = RX + 2*pi;  }

      if ( RY>=pi )
	{  RY = RY - 2*pi;  }
      else if (RY<=-pi)
	{  RY = RY + 2*pi;  }

      if (fabs(yc_r-0.5*pi)<=1.0e-12 || fabs(yc_r+0.5*pi)<=1.0e-12 ||
	  fabs(yc_l-0.5*pi)<=1.0e-12 || fabs(yc_l+0.5*pi)<=1.0e-12)
	{
	  RX = 0.0e0;
	}
      
      rhat.set(1, RX );
      rhat.set(2, RY );

      double ddot = nhat.get(1)*rhat.get(1) + nhat.get(2)*rhat.get(2);

      if (ddot<0.0) // Then need to interchange nodes
	{
	  count2 = count2 + 1;

	  // Switch element order
	  Mesh.set_eelem(i,1,iright);
	  Mesh.set_eelem(i,2,ileft);		
	}
      else
	{  count1 = count1 + 1;  }
    }

  //cout << endl;
  //cout << " NumEdges = " << Mesh.get_NumEdges() << endl;
  //cout << "   count1 = " << count1 << endl;
  //cout << "   count2 = " << count2 << endl;
  //cout << endl;

}
