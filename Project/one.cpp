#include <iostream>
#include <list>
#include <vector>
#include "utils.h"

// GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit.h>

// n >= 2
std::list<std::pair<double, double>> getDataPoints(int n){
    n = n - 1;
    std::list<std::pair<double, double>> points;

    double a = 0;
    double b = Utils::PI;
    double h = (a+b)/n;

    for(int i = 0;i <= n;i++){
        double x = a + i*h;
        points.push_back(std::make_pair(x, sin(x)));
    }

    return points;
}

void interpolate(int n, std::string type, const gsl_interp_type * T){
    auto points = getDataPoints(n);
    n = points.size();

    std::cout << type << " interpolation, points: " <<  n  << std::endl;

    double x[points.size()];
    double y[points.size()];

    // Load Points
    int counter = 0;
    for(auto point : points){
        x[counter] = point.first;
        y[counter] = point.second;
        counter++;
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *interpolation = gsl_spline_alloc(T, n);
    gsl_spline_init(interpolation, x, y, n);

    double totalError = 0;
    counter = 0;

    std::list<std::pair<double, double>> plotFunction;
    std::list<std::pair<double, double>> interpolatedFunction;
    std::list<std::pair<double, double>> errorFunction;

    for(double xi = 0;xi <= Utils::PI;xi += 0.01){
        // Calculate y by using the interpolation
        double yi = gsl_spline_eval(interpolation, xi, acc);
        double fi = sin(xi);

        // Calculate the error
        double error = fabs(fi - yi);
        totalError += error;
        if(error > pow(10, -10)){
            // Try with a higher amount of samples
            interpolate(n + 1, type,  T);
            return;
        }

        plotFunction.push_back(std::make_pair(xi, fi));
        interpolatedFunction.push_back(std::make_pair(xi, yi));
        errorFunction.push_back(std::make_pair(xi, error));

        counter++;
    }

    // Write points
    Utils::writePointsFile(plotFunction, "Onegraphs/Sin-Plot.dat");
    Utils::writePointsFile(interpolatedFunction, "Onegraphs/" + type + "-interpolated.dat");
    Utils::writePointsFile(errorFunction, "Onegraphs/" + type + "-error.dat");


    gsl_spline_free(interpolation);
    gsl_interp_accel_free(acc);

    totalError = totalError/counter;
    std::cout << "Total Error : " << totalError << std::endl;
}

void leastSquares(int n){
    double chisq;

    gsl_matrix* X = gsl_matrix_alloc(n, 3);
    gsl_vector* y = gsl_vector_alloc(n);

    gsl_vector* q = gsl_vector_alloc(3);
    gsl_matrix* cov = gsl_matrix_alloc(3, 3);

    auto points = getDataPoints(n);

    int counter = 0;
    for(auto point : points){
        // Set matrix
        gsl_matrix_set(X, counter, 0, 1.0);
        gsl_matrix_set(X, counter, 1, pow(point.first, 1));
        gsl_matrix_set(X, counter, 2, pow(point.first, 2));

        // Set y Vector
        gsl_vector_set (y, counter, point.second);

        counter++;
    }

    Utils::print_matrix(X);

    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(n, 3);
    gsl_multifit_linear(X, y, q, cov, &chisq, work);
    gsl_multifit_linear_free (work);

    double a = gsl_vector_get(q,0);
    double b = gsl_vector_get(q,1);
    double c = gsl_vector_get(q,2);

    std::cout << "Least Squares: a+bx+cx^2" << std::endl;
    std::cout << "a : " << a << std::endl;
    std::cout << "b : " << b << std::endl;
    std::cout << "c : " << c << std::endl;


    std::list<std::pair<double, double>> interpolatedFunction;
    std::list<std::pair<double, double>> errorFunction;

    for(double xi = 0;xi <= Utils::PI;xi += 0.01){
        // Calculate y by using the interpolation
        double yi = a + pow(xi, 1)*b + pow(xi, 2)*c;
        double fi = sin(xi);

        // Calculate the error
        double error = fabs(fi - yi);

        interpolatedFunction.push_back(std::make_pair(xi, yi));
        errorFunction.push_back(std::make_pair(xi, error));

        counter++;
    }

    // Write points

    Utils::writePointsFile(interpolatedFunction, "Onegraphs/LeastSquares-interpolated.dat");
    Utils::writePointsFile(errorFunction, "Onegraphs/LeastSquares-error.dat");

    gsl_matrix_free(cov);
    gsl_vector_free(q);

    gsl_vector_free(y);
    gsl_matrix_free(X);

}

int main() {
    std::cout.precision(20);
    interpolate(3, "Polynomial",  gsl_interp_polynomial);
    interpolate(3, "Spline",  gsl_interp_cspline);
    leastSquares(13);
}