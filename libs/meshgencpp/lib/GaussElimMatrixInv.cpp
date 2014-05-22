#include "meshdefs.h"
#include <stdlib.h>
#include <stdio.h>

void GaussElimMatrixInv(dTensor2 inMat,
			dTensor2& outMat)
//
//    # Gaussian Elimination with Partial Pivoting Routine.
//    # Produces the inverse matrix outMat from the matrix inMat.
//    # Both of these matrices must be square and the same size.
{
  const int n1 = inMat.getsize(1);
  const int m1 = inMat.getsize(2);
  const int n2 = outMat.getsize(1);
  const int m2 = outMat.getsize(2);
  
  if (n1!=m1 || n2!=m2) 
    { 
      printf(" ERROR in GaussElim.cpp:  matrices are not square \n \n");
      exit(1);
    }
  
  if (n1!=n2 || m1!=m2) 
    { 
      printf(" ERROR in GaussElim.cpp:  input and output matrices are of different size  \n \n");
      exit(1);
    }
  
  const int n = n1;
  
  // Initialize outMat to the identity matrix
  for (int i=1; i<=n; i++)
    {
      for (int j=1; j<=n; j++)
	{
	  if (i==j)
	    { outMat.set(i,j,1.0e0); }
	  else
	    { outMat.set(i,j,0.0e0); }
	}
    }
  
  for (int i=1; i<=n; i++)
    {
      // Select pivot
      double maxA = -100.0e0;
      int p;
      for (int j=i; j<=n; j++)
	{
	  if ( fabs(inMat.get(j,i)) > maxA )
	    {
	      p = j;
	      maxA = fabs(inMat.get(j,i));
	    }
	}
      
      // See if matrix is singular
      if (maxA <= 1.0e-14)
	{
	  printf("  Error in GaussElim.cpp: cannot invert system \n \n");
	  exit(1);
	}
	
      // Pivot matrices
      if (p!=i)
	{
	  for (int j=1; j<=n; j++)
	    {
	      double tmp = inMat.get(i,j);
	      inMat.set(i,j, inMat.get(p,j));
	      inMat.set(p,j, tmp);
	      
	      tmp = outMat.get(i,j);
	      outMat.set(i,j, outMat.get(p,j));
	      outMat.set(p,j, tmp);
	    }
	}
	
      // Eliminate below diagonal
      for (int j=(i+1); j<=n; j++)
	{
	  double dm = inMat.get(j,i)/inMat.get(i,i);
	  for (int k=1; k<=n; k++)
	    {
	      inMat.set(j,k, (inMat.get(j,k) - dm*inMat.get(i,k)) );
	      outMat.set(j,k, (outMat.get(j,k) - dm*outMat.get(i,k)) );
	    }
	}
	
      // Eliminate above diagonal
      for (int j=1; j<=(i-1); j++)
	{
	  double dm = inMat.get(j,i)/inMat.get(i,i);
	  for (int k=1; k<=n; k++)
	    {
	      inMat.set(j,k, (inMat.get(j,k) - dm*inMat.get(i,k)) );
	      outMat.set(j,k, (outMat.get(j,k) - dm*outMat.get(i,k)) );
	    }
	}
      
      // Rescale current row
      double tmp = inMat.get(i,i);
      for (int k=1; k<=n; k++)
	{
 	  inMat.set(i,k, (inMat.get(i,k)/tmp) );
	  outMat.set(i,k, (outMat.get(i,k)/tmp) );
	}
      
    }
  
}
