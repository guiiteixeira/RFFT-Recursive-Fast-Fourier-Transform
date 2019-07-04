# RFFT - Recursive Fast Fourier Transform
Sequential and Parallel implementation of the Recursive Fast Fourier Transform using openMP and MPI with C++.
Concurrent and Parallel Programming - GCC177. Federal University of Lavras - Brazil.

Compiling sequential algorithm:

    -> g++ recursive_fft.cpp -o fftexe
    
Compiling parallel algorithm:

    -> mpic++ -fopenmp parallel_recursive_fft.cpp -o fftexe.o
    
Running sequential algorithm:

    -> ./fftexe

Running parallel algorithm:

    -> mpirun -np 2 fftexe.o 

There is two datasets, one has a size of 8 and the other has a size of 65536. To switch between datasets you must replace 'entrada65536' for 'entrada8' in both files, 'recursive_fft.cpp' and 'parallel_recursive_fft.cpp', and set n = 8 in 'dados.h'. The datasets are in the file 'dados.h'.
