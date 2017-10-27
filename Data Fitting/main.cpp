#include <iostream>
#include <math.h>
#include <list>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "PointsWriter.h"

const double pi = 3.14159265359;

double f(double x){
    return atan(x);
}

// With i in [0,16]
double x(int i){
    return -1.0 + (2.0/16.0)*i;
}

// With i in [0,16]
double t(int i){
    return cos( ( (2*i + 1)*pi )/(34) );
}

int getAmountOfEquidistantGridPoints(double spacing = 0.1){
    return (int)(2/spacing) + 1;
}

double* equidistantGrid(double spacing = 0.1){
    int points = getAmountOfEquidistantGridPoints(spacing);
    double* x = new double[points];

    for(int i = 0;i < points;i++){
        double value = -1 + (i*spacing);
        x[i] = value;
    }

    return x;
}


/**
 * Polynomial
 */
void P(){
    int points = 17;
    double xx[points];
    double yy[points];

    // Calculate Points
    for(int i = 0;i <= points - 1;i++){
        xx[i] = x(i);
        yy[i] = f(xx[i]);
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_polynomial, points);
    gsl_spline_init(spline, xx, yy, points);

    PointsWriter function("data/polynomial-function.data");
    PointsWriter interpolation("data/polynomial-interpolation.data");
    PointsWriter error("data/polynomial-error.data");

    function.header("X", "Y");
    interpolation.header("X", "Y");
    error.header("X", "Y");

    for (double xi = -1;xi <= 1; xi += 0.01) {
        double fi = f(xi);
        double yi = gsl_spline_eval(spline, xi, acc);
        double err = (yi - fi)/(fi);

        function.line(xi, fi);
        interpolation.line(xi, yi);
        error.line(xi, err);
    }

    PointsWriter approx("data/polynomial-approx.data");
    approx.header("X", "Y");

    for(int i = 0;i < points;i++){
        double yi = yy[i];
        double xi = xx[i];
        double calculatedYi = gsl_spline_eval(spline, xi, acc);

        approx.line(xi, fabs(yi - calculatedYi));
    }


    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}


/**
 * Mock ChebyShev
 */
double Q(double spacing = 0.01){
    int points = 17;
    double xx[points];
    double yy[points];

    // Generate chebyshev points
    double chebPoints[points];
    for(int i = 0;i < points;i++){
        chebPoints[i] = t(i);
    }

    // Generate an equisidstant grid
    double* grid = equidistantGrid(spacing);
    int gridPoints = getAmountOfEquidistantGridPoints(spacing);

    //PointsWriter gridWriter("data/mockcheb-grid.data");
    //PointsWriter chebWriter("data/mockcheb-cheb.data");

    //gridWriter.header("X", "Y");
    //chebWriter.header("X", "Y");

    // Write grid points to file
    for(int i = 0;i < gridPoints;i++){
        //gridWriter.line(grid[i], 0);
    }

    // Write Cbehyshev points to file
    for(int i = 0;i < points;i++){
        //chebWriter.line(chebPoints[i], 0);
    }

    // Now calculate the mock chebyshev points
    for(int i = 0;i < points;i++){
        int pointIndex = 0;
        double minDistance = 10;

        for(int j = 0;j < gridPoints;j++){
            double gridX = grid[j];
            double chebX = chebPoints[i];
            double distance = fabs(gridX - chebX);

            if(distance < minDistance){
                pointIndex = j;
                minDistance = distance;
            }
        }

        xx[i] = grid[pointIndex];
        grid[pointIndex] = 10; // Set it to a distance so it would not be used again
    }

    // Sort the points
    std::sort(xx, xx + points);

    for(int i = 0;i < points;i++){
        yy[i] = f(xx[i]);
    }


    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_polynomial, points);
    gsl_spline_init(spline, xx, yy, points);

    /*
    PointsWriter function("data/mockcheb-function.data");
    PointsWriter interpolation("data/mockcheb-interpolation.data");
    PointsWriter error("data/mockcheb-error.data");

    function.header("X", "Y");
    interpolation.header("X", "Y");
    error.header("X", "Y");
     */

    for (double xi = xx[0];xi <= xx[16]; xi += 0.01) {
        double fi = f(xi);
        double yi = gsl_spline_eval(spline, xi, acc);
        double err = (yi - fi)/(fi);

        /*
        function.line(xi, fi);
        interpolation.line(xi, yi);
        error.line(xi, err);
         */
    }

    /*
    PointsWriter approx("data/mockcheb-approx.data");
    approx.header("X", "Y");
     */

    double goodness = 0;
    for(int i = 0;i < points;i++){
        double yi = yy[i];
        double xi = xx[i];
        double calculatedYi = gsl_spline_eval(spline, xi, acc);

        goodness += fabs(calculatedYi);

        //approx.line(xi, fabs(yi - calculatedYi));
    }


    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);

    return goodness;
}

/**
 * Natural Cubic Spline
 */
void S(){
    int points = 17;
    double xx[points];
    double yy[points];

    // Calculate Points
    for(int i = 0;i <= points - 1;i++){
        xx[i] = x(i);
        yy[i] = f(xx[i]);
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, points);
    gsl_spline_init(spline, xx, yy, points);

    PointsWriter function("data/spline-function.data");
    PointsWriter interpolation("data/spline-interpolation.data");
    PointsWriter error("data/spline-error.data");

    function.header("X", "Y");
    interpolation.header("X", "Y");
    error.header("X", "Y");

    for (double xi = -1;xi <= 1; xi += 0.01) {
        double fi = f(xi);
        double yi = gsl_spline_eval(spline, xi, acc);
        double err = (yi - fi)/(fi);

        function.line(xi, fi);
        interpolation.line(xi, yi);
        error.line(xi, err);
    }


    PointsWriter approx("data/spline-approx.data");
    approx.header("X", "Y");

    for(int i = 0;i < points;i++){
        double yi = yy[i];
        double xi = xx[i];
        double calculatedYi = gsl_spline_eval(spline, xi, acc);

        approx.line(xi, fabs(yi - calculatedYi));
    }

    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
}

void calculateGoodnessMockChebychev(){
    double spacing = 0;
    double goodness = 3000;
    for(double i = 0.01;i < 0.1;i += 0.01){
        std::cout << spacing << " & " << goodness << "\\\\" << std::endl;

        double calculation = Q(i);
        if(calculation < goodness){
            spacing = i;
            goodness = calculation;
        }
    }
    std::cout << "Spacing : " << spacing << "  Goodness : "  << goodness << std::endl;

}


int main() {
    std::cout.precision(20);

    P();
    Q(0.095300000000001702793);
    S();
}