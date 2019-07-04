# RFFT - Recursive Fast Fourier Transform
Sequential and Parallel implementation of the Recursive Fast Fourier Transform using openMP and MPI with C++.
Concurrent and Parallel Programming - GCC177. Federal University of Lavras - Brazil.

Compiling sequential algorithm:

    -> g++ recursive_fft.cpp -o fftexe
    
Compiling parallel algorithm:

    -> mpic++ -fopenmp parallel_recursive_fft.cpp -o fftexe.o
    
Running parallel algorithm:

    -> mpirun -np 2 fftexe.o 
