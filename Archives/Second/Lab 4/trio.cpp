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
#include <gsl/gsl_vector_double.h>

double pi = 3.141592653589793;

double calculateA(gsl_vector* ti, gsl_vector* xi, int j){
    int n = ti->size;

    double total = 0.0;
    for(int k = 0;k < n;k++){
        double xk = gsl_vector_get(xi, k);
        double tk = gsl_vector_get(ti, k);
        total += xk*cos(j*tk);
    }

    total *= 2.0/n;
    return total;
}

double calculateB(gsl_vector* ti, gsl_vector* xi, int j){
    int n = ti->size;

    double total = 0.0;
    for(int k = 0;k < n;k++){
        double xk = gsl_vector_get(xi, k);
        double tk = gsl_vector_get(ti, k);
        total += xk*sin(j*tk);
    }

    total *= 2.0/n;
    return total;
}

double calculatePoint(std::vector<double>* a, std::vector<double>* b, double x){
    double y = 0.0;

    // A0
    y += a->at(0)/2.0;

    for(int i = 1;i < a->size();i++){
        y += a->at(i)*cos(i*x);
        y += b->at(i)*sin(i*x);
    }

    return y;
}

void fit(int m){
    gsl_vector* ti = gsl_vector_alloc(11);
    for(int i = 0;i <= 10;i++){
        gsl_vector_set(ti, i, -pi+((2*pi*i)/(11.0)));
    }

    gsl_vector* xi = gsl_vector_alloc(11);
    for(int i = 0;i <= 10;i++){
        double t = gsl_vector_get(ti, i);
        gsl_vector_set(xi, i, f(t));
    }


    std::vector<double>* a = new std::vector<double>;
    std::vector<double>* b = new std::vector<double>;

    for(int i = 0;i <= m;i++){
        double ai = calculateA(ti, xi, i);
        double bi = calculateB(ti, xi, i);
        a->push_back(ai);
        b->push_back(bi);

        std::cout << "a" << i << ": " << ai << ",    b" << i << ": " << bi << std::endl;
    }

    std::vector< std::pair<double, double> > points;
    for(double x = -pi;x < pi;x += 0.01){
        points.push_back(std::make_pair(x, calculatePoint(a, b, x)));
    }

    writePointsFile(points, "TrioLeastSquares.dat");

}

int main (void){
    // Set COUT Precision
    std::cout.precision(20);
    calculateBasicPoints(-pi, pi, "BasicPoints2.dat");
    fit(2);
}