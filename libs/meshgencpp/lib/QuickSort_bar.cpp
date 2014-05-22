#include "meshdefs.h"

// Function to sort a list of "bars"
void QuickSort_bar(bar* b, int* key, int i, int j)
{
  // Swap function
  void Swap_bar(bar* b, int* key, int i, int j);

  // Partition function
  int Partition_bar(bar* b, int* key, int l, int r, int pivot);
  
  // Quicksort
  int pivotindex = (i+j)/2;
  Swap_bar(b, key, j, pivotindex);

  int k = Partition_bar(b, key, i-1, j, key[j]);
  Swap_bar(b, key, j, k);

  if ((k-i) > 1) QuickSort_bar(b, key, i, k-1);
  if ((j-k) > 1) QuickSort_bar(b, key, k+1, j);

  return;
}

// Function to swap bars i and j
void Swap_bar(bar* b, int* key, int i, int j)
{
  int temp;
  temp = key[i];
  key[i] = key[j];
  key[j] = temp;
  temp = b[i].n1;
  b[i].n1 = b[j].n1;
  b[j].n1 = temp;
  temp = b[i].n2;
  b[i].n2 = b[j].n2;
  b[j].n2 = temp;
  
  return;
}

// Function to partition bar according to key
int Partition_bar(bar* b, int* key, int l, int r, int pivot)
{
  void Swap_bar(bar*,int*,int,int);

  do
    {
      while (key[++l] < pivot);
      while (r && (key[--r]) > pivot);
      Swap_bar(b, key, l, r);
    } while (l < r);
  Swap_bar(b, key, l, r);
  return(l);
}
