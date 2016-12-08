#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_blas.h>
#include <sstream>
#include <iomanip>
#include "util.h"
#include <gsl/gsl_chebyshev.h>

double pi = 3.141592653589793;

double chebychevZero(int j){
    int i = 11;
    return cos(((2.0*j-1.0)*pi)/(2.0*i));
}

std::vector< std::pair<double, double> > calculateChebychevPoints(){
    std::vector< std::pair<double, double> > points;

    int i = 11;
    for(int j = 1;j <= i;j++){
        double x = chebychevZero(j);

        points.push_back(std::make_pair(x, f(x)));
    }

    return points;
}



gsl_vector* buildYVector(int m){
    gsl_vector* v = gsl_vector_alloc(m);
    for(int i = 1;i <= m;i++){
        gsl_vector_set(v, i - 1, f(chebychevZero(i)));
    }

    return v;
}

gsl_matrix* buildMatrixA(int n, int m){
    gsl_matrix* A = gsl_matrix_alloc(m, n);

    for(int i = 0;i < m; i++){
        for(int j = 0;j < n;j++){
            double x = chebychevZero(i+1);

            gsl_matrix_set(A, i, j, pow(x, j));
        }
    }

    return A;
}

gsl_vector* generateChebychevLeastSquares(){
    int n = 5;
    int m = 11;

    gsl_vector* c = gsl_vector_alloc(n);
    gsl_vector* y = buildYVector(m);
    gsl_vector* r = gsl_vector_alloc(m);
    gsl_matrix* A = buildMatrixA(n, m);
    gsl_matrix* cov = gsl_matrix_alloc(n, n);
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (m, n);
    double chisq;

    gsl_multifit_linear(A, y, c, cov, &chisq, work);
    gsl_multifit_linear_residuals(A, y, c, r);


    std::cout << "Residu: " << gsl_blas_dnrm2(r) << std::endl;

    gsl_multifit_linear_free(work);
    gsl_matrix_free(cov);
    gsl_matrix_free(A);
    gsl_vector_free(y);
    gsl_vector_free(r);


    return c;
}

gsl_vector* generateChebychevLeastSquares2(){
    int n = 5;
    int m = 11;

    gsl_vector* c = gsl_vector_alloc(n);
    gsl_vector* y = buildYVector(m);

    for(int j = 0;j <= n-1;j++){
        double top = 0.0;
        double bottom = 0.0;

        for(int i = 1;i <= m;i++){
            double fj = pow(chebychevZero(i), j);
            top += fj*gsl_vector_get(y, i-1);
            bottom += pow(fj, 2);

            std::cout << fj << std::endl;
            std::cout << "(i:" << i << ",j:" << j << ")  top: " << top << "   bottom: " << bottom << std::endl;
        }
        std::cout << "\n";

        gsl_vector_set(c, j, top/bottom);
    }

    gsl_vector_free(y);

    return c;
}

double calculateChebychevLeastSquares(gsl_vector* c, double x){
    double output = 0.0;

    for(int i = 0;i < c->size;i++){
        double value = gsl_vector_get(c, i);
        double xValue = pow(x, i);

        output += value*xValue;
    }

    return output;
}

void fitChebychevLeastSquares(gsl_vector* c){
    std::vector< std::pair<double, double> >* points = new std::vector< std::pair<double, double> >;

    for(double x = -1.0;x <= 1.0;x += 0.001){
        double y = calculateChebychevLeastSquares(c, x);

        std::pair<double, double> point = std::make_pair(x,y);
        points->push_back(point);
    }

    writePointsFile(*points, "ChebyChevLeastSquares.dat");
}

int main (void){
    // Set COUT Precision
    std::cout.precision(20);

    calculateBasicPoints(-1.0, 1, "BasicPoints1.dat");

    calculateChebychevPoints();
    gsl_vector* c = generateChebychevLeastSquares();
    gsl_vector* c2 = generateChebychevLeastSquares2();
    print_vector(c);
    std::cout << "" << std::endl;
    print_vector(c2);
    fitChebychevLeastSquares(c);
}