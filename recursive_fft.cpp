#include <iostream>
#include <complex>
#include <cstdlib>
#include "dados.h"
#define MAX 65536

using namespace std;

void FFT(complex<double>* f, int N)
{
  if(N > 1){
    complex<double> evenVec[N/2];
    complex<double> oddVec[N/2];
    int index = 0;

    for(int i = 0;i < N - 1;i += 2){
      evenVec[index] = f[i];
      oddVec[index] = f[i + 1];
      index++; 
    }


    FFT(evenVec,N/2);
    FFT(oddVec,N/2);

    complex<double> taxaW = polar(1., -2. * M_PI / N);
    complex<double> W = 1;

    for(int i = 0; i < N/2;i++){
      f[i] = evenVec[i] + (W*oddVec[i]);
      f[i + (N/2)] = evenVec[i] - (W*oddVec[i]);
      W *= taxaW;
    }
  }
}

int main()
{
    
  complex<double> vec[MAX];
  for(int i = 0; i < n; i++) {
    vec[i] = entrada65536[i];
  }
  FFT(vec, n);
  cout << "...printing the FFT of the array specified" << endl;
  for(int j = 0; j < n; j++)
    cout << vec[j] << endl;
  return 0;
}