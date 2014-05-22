#include "meshdefs.h"

void Unique_int(int*& vec, int*& index, int& size)
{
  int i;
  int tmp = vec[0];
  int k = 0;
  int* other_vec = new int[size];
  other_vec[0] = vec[0];
  
  for (i=1; i<size; i++)
    {
      if (vec[i]!=tmp)
        {
	  k = k+1;
	  other_vec[k] = vec[i];
	  tmp = vec[i];
	  index[k] = index[i];
        }
    }
  
  size = k+1;
  for (i=0; i<size; i++)
    {
      vec[i] = other_vec[i];
    }
  
  delete[] other_vec;
}      
