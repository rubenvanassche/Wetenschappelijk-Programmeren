#include "util.h"

void print_matrix(gsl_matrix * m){
  int rows = m->size1;
  int columns = m->size2;

  for (int i = 0; i < rows; i++){  /* OUT OF RANGE ERROR */
    for (int j = 0; j <  columns; j++){
      std::cout << gsl_matrix_get (m, i, j) << "   ";
    }
    std::cout << std::endl;
  }
}

void print_vector(gsl_vector * v){
  int rows = v->size;

  for (int i = 0; i < rows; i++){  /* OUT OF RANGE ERROR */
      std::cout << gsl_vector_get (v, i) << std::endl;
  }
}
