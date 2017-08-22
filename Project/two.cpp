#include <iostream>
#include <list>
#include <vector>
#include "utils.h"

// GSL
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_chebyshev.h>

enum TYPE{RANDOM, EQUIDISTANT};

double function(double x){
    if(-Utils::PI - 1 <= x and x < 0.0){
        return -1.0;
    }else if(0.0 <= x and x < Utils::PI + 1){
        return 1.0;
    }else{
        throw std::runtime_error("Unkown number for range");
    }
}


// n >= 120 gets one more
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
        points.push_back(std::make_pair(x, function(x)));
    }

    return points;
}

// Returns the error
std::vector<std::pair<double, double> > polynomialInterpolant(std::vector<std::pair<double, double>> generatedPoints, int degree,  TYPE equidistant, double *error, bool *stopped){
    int n = degree + 1; // We need degree + 1 points for polynomial interpolation of degree x

    std::vector<std::pair<double, double>> points;
    if(equidistant == EQUIDISTANT){
        double a = 0;
        double b = generatedPoints.size();
        int h = ceil((b-a)/n);

        for(int i = 0;i < n;i++){
            int x = a + i*h;

            if(x >= generatedPoints.size()){
                // Couln't make this polynom
                *stopped = true;
                return points;
            }

            std::pair<double, double> pair = generatedPoints.at(x);
            points.push_back(pair);
        }
    }else if(equidistant == RANDOM){
        for(int i = 0;i < n;i++){
            auto pair = Utils::select_randomly(generatedPoints.begin(), generatedPoints.end());
            points.push_back(*pair);
            generatedPoints.erase(pair);
        }

        std::sort(points.begin(), points.end());
    }else{

    }

    double x[n];
    double y[n];

    // Load Points
    int counter = 0;
    for(auto point : points){
        x[counter] = point.first;
        y[counter] = point.second;
        counter++;
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *interpolation = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(interpolation, x, y, n);

    *error = 0.0;

    for(double xi = points.front().first;xi <= points.back().first;xi += 0.01){
        double yi = gsl_spline_eval(interpolation, xi, acc);


        *error += fabs(yi - function(xi));
    }


    if(equidistant == EQUIDISTANT){
        //std::cout << "Degree(" << degree <<") Polynomial Interpolation Equidsitant error : " << *error << std::endl;
    }else if(equidistant == RANDOM){
        //std::cout << "Degree(" << degree <<") Polynomial Interpolation Random error : " << *error << std::endl;
    }else{

    }

    return points;
}

void generatePointsFileFromPolynomialInterpolated(std::vector<std::pair<double, double>> points, std::string name){
    int n = points.size();

    double x[n];
    double y[n];

    // Load Points
    int counter = 0;
    for(auto point : points){
        x[counter] = point.first;
        y[counter] = point.second;
        counter++;
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *interpolation = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline_init(interpolation, x, y, n);

    std::list<std::pair<double, double>> interpolatedFunction;


    for(double xi = points.front().first;xi <= points.back().first;xi += 0.01){
        double yi = gsl_spline_eval(interpolation, xi, acc);

        interpolatedFunction.push_back(std::make_pair(xi, yi));
    }

    Utils::writePointsFile(interpolatedFunction, "Twographs/" + name + "-interpolated.dat");
}

double cheb(int grade, double x){
    if(grade == 0){
        return 1.0;
    }else if(grade == 1){
        return x;
    }else if(grade > 1){
        return 2.0 * x * cheb(grade - 1, x) - cheb(grade - 2, x);
    }else{

    }
}

double polynomialLeastSquares(int n, double given){
    /*
    auto points = getDataPoints(120);
    int m = points.size();
    double a = -Utils::PI;
    double b = Utils::PI;

    // Scale to chebychev, not nesairy
    std::vector<double> scalers;
    for(int k = 1;k <= m;k++){
        double scaler = -cos(((2.0*k - 1.0)/(2.0*m))*Utils::PI);
        scalers.push_back(scaler);
    }

    // Scale to interval [–1, 1]
    for(int k = 1;k <= m;k++){
        double x = points.at(k-1).first;
        //double scaled =  -1 + 2*((x + Utils::PI)/(2*Utils::PI));
        double scaled = (scalers.at(k-1) + 1.0)*((b - a)/2.0) + a;
        points.at(k-1).first = scaled;
        points.at(k-1).second = function(scaled);
    }


    std::vector<double> coeff;
    for(int i = 0;i <= n;i++){
        double upper = 0.0;
        double lower = 0.0;

        for(int k = 1;k <= m;k++){
            upper += (points.at(k-1).second * cheb(i, scalers.at(k-1)));
            lower += (pow(cheb(i, scalers.at(k-1)), 2));
        }

        std::cout << upper << " , " << lower << std::endl;

        coeff.push_back(upper/lower);
    }

    double result;
    // Let's calculate
    for(int i = 0;i <= n;i++){
        //result += coeff.at(i) * cheb(i, 2*((given-a)/(b-a) - 1));
        result += coeff.at(i)*given;
    }
    */

    std::vector<double> cj;
    for(int j = 0;j <= n - 1;j++){
        double result = 0.0;

        for(int k = 1;k <= n;k++){
            double xk = cos( (Utils::PI * (k - 0.5))/(n) );
            double xkj = cos( (j*Utils::PI * (k - 0.5))/(n) );
            result += function(xk)*xkj;
        }

        result = result/2.0;

        cj.push_back(result);
    }

    double result;
    for(int k = 0;k <= n-1;k++){
        result += cj.at(k)*cheb(k, given);
    }

    result -= 0.5*cj.at(0);

    return result;
}

void calculatePolynomialLeastSquares(){
    double error = 0.0;
    int n = 3;

    std::list<std::pair<double, double>> interpolatedFunction;

    for(double xi = -Utils::PI;xi <= Utils::PI;xi += 0.1){
        double yi = polynomialLeastSquares(n, xi);

        interpolatedFunction.push_back(std::make_pair(xi, yi));

        error += fabs(yi - function(xi));
    }

    Utils::writePointsFile(interpolatedFunction, "Twographs/PolynomialLS-interpolated.dat");

    std::cout << "Polynomial Least Squares n: " << n << ", error: " << error << std::endl;
}

void calculatePolynomialInterpolant(){
    double minError = pow(10,100000);
    int degree = 0;
    std::vector<std::pair<double, double> > equidistantPoints;
    for(int i = 3;i <= 240;i++){
        double error = 0.0;
        bool stopped = false;
        auto temp = polynomialInterpolant(getDataPoints(240), i, EQUIDISTANT, &error, &stopped);
        if(error < minError and stopped == false){
            minError = error;
            degree = i;
            equidistantPoints = temp;
        }
    }

    std::cout << "Polynomial Equidistant function with degree: " << degree << " and error: " << minError << std::endl;


    minError = pow(10,100000);
    degree = 0;
    std::vector<std::pair<double, double> > randomPoints;
    for(int i = 3;i <= 240;i++){
        double error = 0.0;
        bool stopped = false;
        auto temp = polynomialInterpolant(getDataPoints(240), i, RANDOM, &error, &stopped);
        if(error < minError and stopped == false){
            minError = error;
            degree = i;
            randomPoints = temp;
        }
    }

    // Generate file
    generatePointsFileFromPolynomialInterpolated(equidistantPoints, "PolynomialEquidistant");
    generatePointsFileFromPolynomialInterpolated(randomPoints, "PolynomialRandom");

    std::cout << "Polynomial Random function with degree: " << degree << " and error: " << minError << std::endl;
}



double aj(std::vector<std::pair<double, double>> points, int j){
    double result = 0.0;
    int n = points.size();

    for(int k = 0;k < n;k++){
        result += points.at(k).second*cos(j * points.at(k).first);
    }

    result = result * (2.0/n);

    return result;

}



double bj(std::vector<std::pair<double, double>> points, int j){
    double result = 0.0;
    int n = points.size();

    for(int k = 0;k < n;k++){ // till n-1 otherwise periodic
        result += points.at(k).second*sin(j * points.at(k).first);
    }

    result = result * (2.0/n);

    return result;

}

double a0(std::vector<std::pair<double, double>> points){
    return aj(points, 0);
}

// m -> amount of coef
// x -> x value for generated function
double triGoniometric(std::vector<std::pair<double, double>> points, int m, double x){
    int n = points.size();


    double approx = a0(points)/2.0;

    if(n == 2*m){
        for(int i = 1;i <= m - 1;i++){
            approx += aj(points, i)*cos(i * x);
            approx += bj(points, i)*sin(i * x);
        }
        approx += aj(points, m)*0.5*cos(m * x);
    }else{
        for(int i = 1;i <= m;i++){
            approx += aj(points, i)*cos(i * x);
            approx += bj(points, i)*sin(i * x);
        }
    }

    return approx;
}

void calculateTriogoniometricLeastSquares(){
    // Least Squares -> n > 2m + 1

    auto points = getDataPoints(120);
    int m = 25; // 121 > 41
    double error = 0.0;

    std::list<std::pair<double, double>> interpolatedFunction;

    for(double xi = -Utils::PI;xi <= Utils::PI;xi += 0.01){
        double yi = triGoniometric(points, m, xi);

        error += fabs(yi - function(xi));

        interpolatedFunction.push_back(std::make_pair(xi, yi));
    }

    Utils::writePointsFile(interpolatedFunction, "Twographs/trioLS-interpolated.dat");

    std::cout << "Triogoniometric Least Squares with m : " << m << ", n : " << points.size() << ", error : " << error << std::endl;
}

void calculateTriogoniometricInterpolant(){
    // Least Squares -> n > 2m + 1

    auto points = getDataPoints(121);
    int m = 61; // 121 > 41
    double error = 0.0;

    std::list<std::pair<double, double>> interpolatedFunction;

    for(double xi = -Utils::PI;xi <= Utils::PI;xi += 0.01){
        double yi = triGoniometric(points, m, xi);

        error += fabs(yi - function(xi));

        interpolatedFunction.push_back(std::make_pair(xi, yi));
    }

    Utils::writePointsFile(interpolatedFunction, "Twographs/trioI-interpolated.dat");

    std::cout << "Triogoniometric Interpolant with m : " << m << ", n : " << points.size() << ", error : " << error << std::endl;
}

int main() {
    std::cout.precision(20);

    calculatePolynomialInterpolant();
    calculateTriogoniometricLeastSquares();
    calculateTriogoniometricInterpolant();
    calculatePolynomialLeastSquares();
}
