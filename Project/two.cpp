#include <iostream>
#include <list>
#include <vector>
#include "utils.h"

// GSL
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit.h>


// n >= 120
std::vector<std::pair<double, double>> getDataPoints(int n){
    if(n < 120){
        throw std::runtime_error("n should be more then 120");
    }

    std::vector<std::pair<double, double>> points;

    double a = -Utils::PI;
    double b = Utils::PI;
    double h = (b-a)/n;

    for(int i = 0;i <= n;i++){
        double x = a + i*h;
        if(-Utils::PI - 1 <= x and x < 0.0){
            points.push_back(std::make_pair(x, -1));
        }else if(0.0 <= x and x < Utils::PI + 1){
            points.push_back(std::make_pair(x, 1));
        }else{
            throw std::runtime_error("Unkown number for range");
        }
    }

    return points;
}

void polynomialInterpolant(std::vector<std::pair<double, double>> generatedPoints, int degree,  bool equidistant){
    int n = degree + 1; // We need degree + 1 points for polynomial interpolation of degree x

    std::vector<std::pair<double, double>> points;
    if(equidistant == true){
        double a = 0;
        double b = generatedPoints.size();
        int h = floor((b-a)/n);

        for(int i = 0;i < n;i++){
            int x = a + i*h;
            std::pair<double, double> pair = generatedPoints.at(x);
            points.push_back(pair);
        }
    }else{
        for(int i = 0;i < n;i++){
            auto pair = Utils::select_randomly(generatedPoints.begin(), generatedPoints.end());
            points.push_back(*pair);
            generatedPoints.erase(pair);
        }

        std::sort(points.begin(), points.end());
    }

    double x[n];
    double y[n];

    // Load Points
    int counter = 0;
    for(auto point : points){
        x[counter] = point.first;
        y[counter] = point.second;
        counter++;

        std::cout << point.first << ", " << point.second << std::endl;
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *interpolation = gsl_spline_alloc(gsl_interp_polynomial, n);
    gsl_spline_init(interpolation, x, y, n);
}

int main() {
    std::cout.precision(20);

    polynomialInterpolant(getDataPoints(120), 10, false);
}