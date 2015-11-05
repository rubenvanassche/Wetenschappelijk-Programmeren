#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <math.h>
#include "util.h"

gsl_matrix* calculateTriangularMatrix(int size){
  // Create the triangular matrix
  int matrix_size = size*size;
  double data[matrix_size];

  // initialize data zero
  for (int i = 0; i < size; i++){  /* OUT OF RANGE ERROR */
    for (int j = 0; j <  size; j++){
      data[size*i + j] = 0;
    }
  }

  // Upper diagonal
  int upper_counter = 1;
  while(true){
    if(upper_counter >= matrix_size){
      break;
    }

    data[upper_counter] = -1.0;
    upper_counter += size + 1;
  }
  // Middle diagonal
  int middle_counter = 0;
  while(true){
    if(middle_counter >= matrix_size){
      break;
    }

    int index_on_diagonal = middle_counter % size;
    if((index_on_diagonal + 2) % 3 == 0){
      data[middle_counter] = (2.0*(double(index_on_diagonal) + 2.0))/3.0;
    }else{
      data[middle_counter] = 1;
    }

    middle_counter += size + 1;
  }
  // Lower diagonal
  int lower_counter = size;
  while(true){
    if(lower_counter >= matrix_size){
      break;
    }

    data[lower_counter] = 1.0;
    lower_counter += size + 1;
  }

  gsl_matrix * m =  gsl_matrix_alloc (size, size);
  for (int i = 0; i < size; i++){  /* OUT OF RANGE ERROR */
    for (int j = 0; j <  size; j++){
      gsl_matrix_set (m, i, j, data[i*size + j]);
    }
  }
  return m;
}

gsl_vector* calculateYVector(int size){
  double data[size];

  for(int i = 0;i < size;i++){
      data[i] = 0;
  }

  data[0] = 1; // First one should be 1
  gsl_vector * v = gsl_vector_alloc(size);
  for(int i = 0; i < size;i++){
    gsl_vector_set(v, i, data[i]);
  }

  return v;
}

double calculatePrecision(double result){
  return fabs((M_E-2.0-result)/(M_E));
}

double process(int size){
  gsl_matrix* m = calculateTriangularMatrix(size);
  gsl_vector* v = calculateYVector(size);

  // Allocate the vector to store the x values
  gsl_vector *x = gsl_vector_alloc (size);

  // GSL requirements
  int s;
  gsl_permutation * p = gsl_permutation_alloc(size);

  // Do an LU decomposition(GSL uses Gaussian elimination with partial pivoting)
  gsl_linalg_LU_decomp (m, p, &s);

  // Now solve this matrix using the y values and store the x values in x
  gsl_linalg_LU_solve (m, p, v, x);

  // get x1
  double result = gsl_vector_get(x, 0);

  // cleaning
  gsl_permutation_free (p);
  gsl_vector_free (x);
  gsl_vector_free (v);
  gsl_matrix_free (m);
  return result;
}

int main (void){
  double precisionToReach = pow(10, -10);

  std::cout << "Find e-2" << std::endl;
  std::cout << "------------------" << std::endl;
  std::cout << "Precision to reach: " << precisionToReach << std::endl << std::endl;

  int i = 1;
  while(true){
    double result = process(i);
    double precision = calculatePrecision(result);

    std::cout << "n: " << i + 1 << "  result: " << result << "    precision:" << precision << std::endl;

    if(precision <= precisionToReach){
      std::cout << std::endl << std::endl;
      std::cout << "------------------" << std::endl;
      std::cout << "Needed " << i + 1 << "x" << i + 1 << " matrix" << std::endl;
      std::cout << "Result: " << result << std::endl;
      std::cout << "Precision:" << precision << std::endl;
      std::cout << "Absolute fault:" << fabs((M_E - 2.0) - result) << std::endl;

      break;
    }

    i++;
  }
}
