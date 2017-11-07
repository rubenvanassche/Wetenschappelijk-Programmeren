#include <iostream>
#include <math.h>
#include <list>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include "PointsWriter.h"
#include "Utils.h"

double combinate(const unsigned int n, const unsigned int k){
    return ( gsl_sf_fact(n) )/(gsl_sf_fact(k) * gsl_sf_fact(n - k));
}

void solve(gsl_matrix* A, gsl_vector* x, gsl_vector* b){
    int s;
    gsl_permutation * p = gsl_permutation_alloc (A->size1);

    gsl_linalg_LU_decomp (A, p, &s);
    gsl_linalg_LU_solve (A, p, b, x);

    gsl_permutation_free (p);
}

/**
 * Calculates an error vector
 *
 * @param x1 The exact computed value
 * @param x2 The linear computed value
 * @param error Vector with error values
 */
void errorVector(gsl_vector* x1, gsl_vector* x2, gsl_vector* error){
    int size = x1->size;

    for(int i = 0;i < size;i++){
        double result = gsl_vector_get(x1, i) - gsl_vector_get(x2, i);
        gsl_vector_set(error, i, result);
    }
}

/**
 * Calculates the residual vector
 *
 * @param A
 * @param x
 * @param b
 * @param residual
 */
void residualVector(gsl_matrix* A, gsl_vector* x, gsl_vector* b, gsl_vector* residual){
    int size = A->size1;

    // Ax
    for(int i = 0;i < size;i++){
        for(int j = 0;j < size;j++){
            double result = gsl_vector_get(x, j) * gsl_matrix_get(A, i, j);
            gsl_vector_set(residual, j, result);
        }
    }

    // y - Ax
    for(int i = 0;i < size;i++){
        double result = gsl_vector_get(b, i) - gsl_vector_get(residual, i);
        gsl_vector_set(residual, i, result);
    }
}


void build(unsigned  int n){
    // Allocate vectors and matrices
    gsl_matrix* A = gsl_matrix_alloc(n, n);
    gsl_vector* x1 = gsl_vector_alloc(n); // xj = 1
    gsl_vector* x2 = gsl_vector_alloc(n); // xj = unkown
    gsl_vector* residual = gsl_vector_alloc(n);
    gsl_vector* error = gsl_vector_alloc(n);
    gsl_vector* b = gsl_vector_alloc(n);

    // Set the values
    for(unsigned  int i = 1;i <= n;i++){
        for(unsigned  int j = 1;j <= n;j++){
            gsl_matrix_set(A, i - 1, j - 1, combinate(i + j - 2, j - 1));
        }

        gsl_vector_set(x1, i - 1, 1);
        gsl_vector_set(b, i - 1, combinate(n + i - 1, i));
    }
    std::cout << "A: " << std::endl <<  std::endl;
    printMatrix(A, false);
    std::cout << std::endl;


    // Calculate condition number
    double kA = conditionNumber(A);
    std::cout << "Condition Number : " <<  kA << std::endl;
    std::cout << "Normal rounding error : " << kA*pow(2,-52) << std::endl;

    // Solve the system
    solve(A, x2, b);

    std::cout << "x Vector : " << std::endl << std::endl;
    printVector(x2, true);

    // Error vector
    errorVector(x1, x2, error);
    std::cout << std::endl << "Error Vector : " << std::endl << std::endl;
    printVector(error);
    std::cout << "Norm : " << gsl_blas_dnrm2(error) << std::endl;

    // Residual vector
    residualVector(A, x2, b, residual);
    std::cout << std::endl << "Residual Vector : " << std::endl << std::endl;
    printVector(residual);
    std::cout << "Norm : " << gsl_blas_dnrm2(residual) << std::endl;

    // Deallocation
    gsl_vector_free(b);
    gsl_vector_free(x1);
    gsl_vector_free(x2);
    gsl_vector_free(error);
    gsl_vector_free(residual);
    gsl_matrix_free(A);
}



int main(){
    std::cout.precision(4);

    for(int i = 1;i <= 4;i++){
        std::cout << "n = " << i*3 << std::endl;
        std::cout << "-----------------------" << std::endl;

        build(i*3);

        std::cout << std::endl;
    }

    std::cout << "ULP double : " << std::numeric_limits<double>::epsilon() << std::endl;
    std::cout << "ULP required : " << pow(2, -52) << std::endl;
}