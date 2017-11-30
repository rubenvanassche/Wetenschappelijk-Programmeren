#! /bin/bash
export MKLROOT=/opt/intel/mkl
g++-7  -DMKL_ILP64 -m64 -fopenmp -I${MKLROOT}/include -c main.cpp
g++-7  ${MKLROOT}/lib/libmkl_intel_ilp64.a ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl main.o -o main
rm main.o
time OMP_NUM_THREADS=4 ./main
