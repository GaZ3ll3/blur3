#include "meshdefs.h"

void QuickSort(double*& a, int*& index, int lo, int hi)
{
  //  lo is the lower index, hi is the upper index
  //  hi the region of array a that is to be sorted
  
  int i = lo;
  int j = hi;
  double x=a[(i+j)/2];
  double h;
  int itmp;
  
  //  partition
  while(i<=j) 
    {           
      while (a[i]<x) 
        {i++;} 
      
      while (a[j]>x) 
        {j--;}
      
      if (i<=j)
        {
	  h=a[i]; 
	  a[i]=a[j];
	  a[j]=h;
          
	  itmp = index[i];
	  index[i] = index[j];
	  index[j] = itmp;
          
	  i++; j--;
        }
    }
  
  //  recursion
  if (lo<j) QuickSort(a, index, lo, j);
  if (i<hi) QuickSort(a, index, i, hi);
}