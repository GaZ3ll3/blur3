#include "meshdefs.h"
#include "mesh.h"

void MeshStore(point p[], 
	       triangle t[], 
	       int bnd_node[], 
	       int ghost_link[], 
	       int ext_node_link[],
	       double area[], 
	       double cdual[], 
	       dTensor2*& edge, 
	       iTensor2*& tedge,
	       iTensor2*& tedge_orientation,
	       iTensor2*& eelem, 
	       iTensor2*& enode,
	       iTensor1*& bnd_edge,
	       mesh& Mesh)
{
  const int NumElems      = Mesh.get_NumElems();
  const int NumGhostElems = Mesh.get_NumGhostElems();
  const int NumNodes      = Mesh.get_NumNodes();
  const int NumPhysNodes  = Mesh.get_NumPhysNodes();
  const int NumBndNodes   = Mesh.get_NumBndNodes();
  const int NumEdges      = Mesh.get_NumEdges();
  const int NumBndEdges   = Mesh.get_NumBndEdges();

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

  // ext_node_link
  for (int i=1; i<=(NumNodes-NumPhysNodes); i++)
    {
      Mesh.set_ext_node_link(i, ext_node_link[i-1] + 1 );
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
  
  // tedge_orientation
  for (int i=1; i<=NumElems; i++)
    for (int k=1; k<=3; k++)
      {
	Mesh.set_tedge_orientation(i,k, tedge_orientation->get(i,k) );
      }

  // eelem
  for (int i=1; i<=NumEdges; i++)
    for (int k=1; k<=2; k++)
      {
	Mesh.set_eelem(i,k, eelem->get(i,k) );
      }

  // enode
  for (int i=1; i<=NumEdges; i++)
    for (int k=1; k<=2; k++)
      {
	Mesh.set_enode(i,k, enode->get(i,k) );
      }

  // bnd_edge
  for (int i=1; i<=NumBndEdges; i++)
    {
      Mesh.set_bnd_edge(i, bnd_edge->get(i) );
    }

  // Compute min & max edge lengths
  Mesh.compute_min_max_edge_length();
}
