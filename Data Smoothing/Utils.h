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

void residualVector(gsl_matrix* A, gsl_vector* x, gsl_vector* b, gsl_vector* residual){
    // Ax
    for(int i = 0;i < A->size1;i++){
        double result = 0;
        for(int j = 0;j < A->size2;j++){
            result += gsl_vector_get(x, j) * gsl_matrix_get(A, i, j);
        }
        gsl_vector_set(residual, i, result);
    }

    // y - Ax
    for(int i = 0;i < A->size1;i++){
        double result = gsl_vector_get(b, i) - gsl_vector_get(residual, i);
        gsl_vector_set(residual, i, result);
    }
}

double conditionNumber(gsl_matrix* m){
    if(m->size1 != m->size2){
        //throw new std::runtime_error("Can only calculate the condition number of an square matrix");
    }

    int size1 = m->size1;
    int size2 = m->size2;

    gsl_vector* v = gsl_vector_alloc(size1);
    gsl_matrix* temp = gsl_matrix_alloc(size2,size1);
    gsl_vector* w = gsl_vector_alloc(size1);

    gsl_linalg_SV_decomp(m, temp, v, w);

    double min = std::numeric_limits<double>::max();
    double max = -std::numeric_limits<double>::min();

    for(int i = 0; i < std::min(size1, size2);i++){
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

