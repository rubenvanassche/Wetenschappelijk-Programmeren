//
// Created by Ruben Van Assche on 31/10/17.
//
#include <gsl/gsl_linalg.h>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector_double.h>

#ifndef SYSTEM_OF_LINEAR_EQUATIONS_UTILS_H
#define SYSTEM_OF_LINEAR_EQUATIONS_UTILS_H

void printMatrix(gsl_matrix * m, bool latex = false){
    int rows = m->size1;
    int columns = m->size2;

    for (int i = 0; i < rows; i++){
        for (int j = 0; j <  columns; j++){
            if(latex == true){
                std::cout << gsl_matrix_get (m, i, j) << " ";
                if(j != columns - 1){
                    std::cout << "& ";
                }
            }else{
                std::cout << gsl_matrix_get (m, i, j) << "  ";
            }

        }
        if(latex == true){
            std::cout << "\\\\" << std::endl;
        }else{
            std::cout << std::endl;
        }

    }
}

void printVector(gsl_vector * v, bool latex = false){
    int rows = v->size;

    for (int i = 0; i < rows; i++){
        if(latex == true){
            std::cout << gsl_vector_get (v, i) << " \\\\" << std::endl;
        }else{
            std::cout << gsl_vector_get (v, i) << std::endl;
        }

    }
}

double conditionNumber(gsl_matrix* m){
    if(m->size1 != m->size2){
        throw new std::runtime_error("Can only calculate the condition number of an square matrix");
    }

    int size = m->size1;

    gsl_vector* v = gsl_vector_alloc(size);
    gsl_matrix* temp = gsl_matrix_alloc(size,size);
    gsl_vector* w = gsl_vector_alloc(size);

    gsl_linalg_SV_decomp(m, temp, v, w);

    double min = std::numeric_limits<double>::max();
    double max = -std::numeric_limits<double>::min();

    for(int i = 0; i < size;i++){
        double result = gsl_vector_get(v, i);
        if(result > max){
            max = result;
        }

        if(result < min){
            min = result;
        }
    }

    gsl_vector_free(v);
    gsl_matrix_free(temp);
    gsl_vector_free(w);

    return max/min;
}


#endif //SYSTEM_OF_LINEAR_EQUATIONS_UTILS_H

