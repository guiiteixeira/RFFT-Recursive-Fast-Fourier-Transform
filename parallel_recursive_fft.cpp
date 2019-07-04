#include <iostream>
#include <complex>
#include <cstdlib>
#include "mpi.h"
#include <omp.h>
#include "dados.h"
#define MAX 65536

using namespace std;

int log2(int N)    /*function to calculate the log2(.) of int numbers*/
{
  int k = N, i = 0;
  while(k) {
    k >>= 1;
    i++;
  }
  return i - 1;
}


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

    complex<double> W[N/2];
    W[0] = 1;
    W[1] = polar(1., -2. * M_PI / N);

    #pragma omp parallel for
    for(int i = 2;i < (N/2);i++){
      W[i] = pow(W[1],i);
    }

    #pragma omp parallel for shared(W)
    for(int i = 0; i < N/2;i++){
      f[i] = evenVec[i] + (W[i]*oddVec[i]);
      f[i + (N/2)] = evenVec[i] - (W[i]*oddVec[i]);
    }
  }
}

int main(int argc, char** argv)
{
	int size,rank;

  MPI_Init(&argc,&argv);
  MPI_Comm_size (MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status status;


  double d = 1.0;
  //cout << "specify sampling step" << endl; //just write 1 in order to have the same results of matlab fft(.)
  //cin >> d;
  complex<double> vec[MAX];
  complex<double> oddVec[MAX/2];
  complex<double> evenVec[MAX/2];
  int index = 0;
  //cout << "specify the array" << endl;
  if(rank == 0){
    for(int i = 0; i < (n - 1); i += 2) {
      evenVec[index] = entrada65536[i];
      oddVec[index]= entrada65536[i+1];
      index++;
    }
     
    MPI_Send(&oddVec, ((n/2) * sizeof(complex<double>)), MPI_BYTE, 1,0, MPI_COMM_WORLD);

    FFT(evenVec,(n/2));

    MPI_Send(&evenVec, ((n/2) * sizeof(complex<double>)), MPI_BYTE, 1,0, MPI_COMM_WORLD);
    
  }
  else if(rank == 1){
     MPI_Recv(&oddVec, ((n/2) * sizeof(complex<double>)), MPI_BYTE, 0,0, MPI_COMM_WORLD, &status);
     
     FFT(oddVec,(n/2));

     MPI_Recv(&evenVec, ((n/2) * sizeof(complex<double>)), MPI_BYTE, 0,0, MPI_COMM_WORLD, &status);

    complex<double> W[n/2];
    W[0] = 1;
    W[1] = polar(1., -2. * M_PI / n);

    omp_set_num_threads(4);

    #pragma omp parallel for
    for(int i = 2;i < (n/2);i++){
     W[i] = pow(W[1],i);
    }
    #pragma omp parallel for shared(W)
    for(int i = 0;i < (n/2);i++){
      vec[i] = evenVec[i] + (W[i]*oddVec[i]);
      vec[i + (n/2)] = evenVec[i] - (W[i]*oddVec[i]);
    }
  }

  if(rank == 1){

    cout << "done.." << endl;
    for(int j = 0; j < n; j++)
      cout << vec[j] << endl;
  }
  MPI_Finalize();
  return 0;
}