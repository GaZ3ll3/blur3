#include "meshdefs.h"
#include "mesh.h"

void MeshStore(point p[], 
	       triangle t[], 
	       int bnd_node[], 
	       int ghost_link[],
	       double area[], 
	       double cdual[], 
	       dTensor2*& edge, 
	       iTensor2*& tedge, 
	       iTensor2*& eelem, 
	       mesh& Mesh)
{
  int NumElems      = Mesh.get_NumElems();
  int NumGhostElems = Mesh.get_NumGhostElems();
  int NumNodes      = Mesh.get_NumNodes();
  int NumBndNodes   = Mesh.get_NumBndNodes();
  int NumEdges      = Mesh.get_NumEdges();

  // node
  for (int i=1; i<=NumNodes; i++)
    {
      Mesh.set_node(i,1, p[i-1].x);
      Mesh.set_node(i,2, p[i-1].y);
    }

  // tnode
  for (int i=1; i<=NumElems; i++)
    {
      Mesh.set_tnode(i,1, t[i-1].n1 + 1 );
      Mesh.set_tnode(i,2, t[i-1].n2 + 1 );
      Mesh.set_tnode(i,3, t[i-1].n3 + 1 );
    }
  
  // bnd_node
  for (int i=1; i<=NumBndNodes; i++)
    {
      Mesh.set_bnd_node(i, bnd_node[i-1] + 1 );
    }

  // ghost_link
  for (int i=1; i<=NumGhostElems; i++)
    {
      Mesh.set_ghost_link(i, ghost_link[i-1] + 1 );
    }


  // area_prim
  for (int i=1; i<=NumElems; i++)
    {
      Mesh.set_area_prim(i, area[i-1] );
    }

  // area_dual
  for (int i=1; i<=NumNodes; i++)
    {
      Mesh.set_area_dual(i, cdual[i-1] );
    }

  // edge
  for (int i=1; i<=NumEdges; i++)
    for (int k=1; k<=4; k++)
      {
	Mesh.set_edge(i,k, edge->get(i,k) );
      }
  
  // tedge
  for (int i=1; i<=NumElems; i++)
    for (int k=1; k<=3; k++)
      {
	Mesh.set_tedge(i,k, tedge->get(i,k) );
      }
  
  // eelem
  for (int i=1; i<=NumEdges; i++)
    for (int k=1; k<=2; k++)
      {
	Mesh.set_eelem(i,k, eelem->get(i,k) );
      }
}
