#include <iostream>
#include <math.h>
#include <list>
#include "mkl_vsl.h"
#include "omp.h"
#include <stdio.h>


long seed = -67867890;

bool inFigure(double x, double y, double z){
  if(exp(x) <= y){
    if(y*sin(z) >= 0){
      return true;
    }
  }

  return false;
}

int calculate(int N = 1000, int LOOPS = 1000){
    int threadID;
    int numberOfThreads = omp_get_max_threads();

    std::cout << "Amount of threads using:" << numberOfThreads << std::endl;
    VSLStreamStatePtr stream[numberOfThreads];

    for(int i = 0; i < numberOfThreads; i++){
        int errorCode = 0;

        errorCode = vslNewStream(&stream[i], VSL_BRNG_MCG59, seed);
        if(errorCode){
            printf("vslNewStream failed\n");
            return 1;
        }

        errorCode = vslLeapfrogStream(stream[i], i, numberOfThreads);
        if(errorCode) {
            printf("vslLeapfrogStream failed\n");
            return 1;
        }
    }

    long randomNumbersPerThread = 3 * N;
    int figure_tk = 0;
    int figure_tn = 0;

    #pragma omp parallel private(threadID) num_threads(numberOfThreads) reduction( + : figure_tk, figure_tn )
    {
        threadID = omp_get_thread_num();
        double randomNumbers[randomNumbersPerThread];

        for (int j = 0; j < LOOPS;++j){
            vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream[threadID], randomNumbersPerThread, randomNumbers, 0, 1);

            int n = 0;
            int figure_k = 0;
            int figure_n = 0;

            for(int i = 0; i < N; ++i ){
                double x = randomNumbers[n++];
                double y = randomNumbers[n++];
                double z = randomNumbers[n++];

                // Manipulate the numbers
                x = 2*x;
                y = y + 1;
                z = 4*z - 1;

                /// check for points in cone or sphere respectively
                if(inFigure(x, y, z) == true){
                  figure_k++;
                }

                figure_n++;
            }

            figure_tk += figure_k;
            figure_tn += figure_n;
        }

    }   // end of thread


    double figure_tv = 8*((double)figure_tk/(double)figure_tn);
    std::cout << "Total Points: " << figure_tn << std::endl;
    std::cout << "Volume: " << figure_tv << std::endl;


    for(int i = 0; i < numberOfThreads; i++){
        vslDeleteStream(&stream[i]);
    }
}

int main(){
  for(int i = 10;i < pow(10, 100);i = i*2){
    calculate(1000, i);
  }

}
