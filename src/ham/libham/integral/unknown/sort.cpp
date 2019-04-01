//
// Sort the data in ascend order: The simplest sorting
//
// 
// This implements the bubble sort algorithm. 
// Works for smaller N, say N<1000
//

#include <cstdlib>
#include <iostream>

using namespace std;

#define NUM 10

void BubbleSort(int N, double * coeff);
void randVec(int N, double * c0);

void main(){
  int     i, n;
  double  coeff[NUM];

  n=NUM;
  randVec(n,coeff);
  cout << " Coeff before the sorting:\n";
  for (i=0;i<n;i++)
    cout << coeff[i]<<' ';
  cout<<endl;

  BubbleSort(n,coeff);
  cout << " Coeff After the sorting:\n";
  for (i=0;i<n;i++)
    cout << coeff[i]<<' ';
  cout<<endl;
}



void BubbleSort(int N, double * coeff){

  double tmp;
  int i, j;

  for ( i=N-1;i>=0; i--)
    for (j=1;j<=i;j++){
      if (coeff[j-1]<coeff[j]){
	tmp = coeff[j-1];
	coeff[j-1]=coeff[j];
	coeff[j]=tmp;
      }
    }
}


void randVec(int N, double * c0){

  for (int i=0;i<N;i++)
    c0[i]= (double) rand();
}

//
// coeff has been sorted 
//
int getSize(int n, int * nSize, double ** coeff, double cutoff){

  int totalSize, i;
  double tmp;
  double ind[n];

  for (i=0;i<n;i++)
    ind[i]=0;

  totalSize=0;
  while (true){
    tmp=1.0;
    for (i=0;i<n;i++){
      tmp *= coeff[ind[i]];
      if (tmp<cutoff) break;      
    }
    
    if (tmp>cutoff)
      totalSize++;

    for (j=i+1;j<n;j++)
      ind[j]=0;

    ind[i]++;

    for (j=i;j>0;j--){
      if (ind[i]>nSize[i]){
	ind[i]=0; ind[i-1]++;
      } 
    }

    if (ind[0]>nSize[0]) break;
    
  }

}
