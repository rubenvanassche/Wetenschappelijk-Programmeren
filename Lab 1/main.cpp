#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <string>

/*
int main (void){
  int i;
  double xi, yi, x[10], y[10];

  printf ("#m=0,S=2\n");

  for (i = 0; i < 10; i++)
    {
      x[i] = i + 0.5 * sin (i);
      y[i] = i + cos (i * i);
      printf ("%g %g\n", x[i], y[i]);
    }

  printf ("#m=1,S=0\n");

  {
    gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
    gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_cspline, 10);

    gsl_spline_init (spline, x, y, 10);

    for (xi = x[0]; xi < x[9]; xi += 0.01)
      {
        yi = gsl_spline_eval (spline, xi, acc);
        printf ("%g %g\n", xi, yi);
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
  }
  return 0;
}
*/

void writePointsFile(std::vector< std::pair<double, double> >& points, std::string filename){
    std::ofstream file;
    file.open(filename);

    file << "#  X   Y" << std::endl;
    for(auto i: points){
        file << i.first << "   " << i.second << std::endl;
    }

    file.close();
}

void cubicSplineInterpolation(double x[], double y[], int size){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline_periodic, size);
    gsl_spline_init(spline, x, y, size);

    std::vector< std::pair<double, double> > points;

    for(double xi = x[0];xi <= x[size-1];xi += 0.01){
        double yi = gsl_spline_eval(spline, xi, acc);

        auto pair = std::make_pair(xi,yi);
        points.push_back(pair);
    }

    writePointsFile(points, "CubicSpline.dat");

}

void naturalSplineInterpolation(double x[], double y[], int size){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, size);
    gsl_spline_init(spline, x, y, size);

    std::vector< std::pair<double, double> > points;

    for(double xi = x[0];xi <= x[size-1];xi += 0.01){
        double yi = gsl_spline_eval(spline, xi, acc);

        auto pair = std::make_pair(xi,yi);
        points.push_back(pair);
    }

    writePointsFile(points, "naturalSpline.dat");

}

void polynomialInterpolation(double x[], double y[], int size){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_interp* interpolant = gsl_interp_alloc (gsl_interp_polynomial, size);
    gsl_interp_init(interpolant, x, y ,size);

    std::vector< std::pair<double, double> > points;

    for(double xi = x[0];xi <= x[size-1];xi += 0.01){
        double yi = gsl_interp_eval(interpolant, x, y, xi, acc);

        auto pair = std::make_pair(xi,yi);
        points.push_back(pair);
    }

    writePointsFile(points, "Polynomial.dat");
}

int main (void){
    // Set COUT Precision
    std::cout.precision(20);

    std::cout << "Building Points..." << std::endl;

    double x[7] = {0,23,37,54,74,88,130};
    double y[7] = {2.8, 3.6, 4.4, 5.5, 6.4, 7.2, 8.3};
    int size = 7;

    std::vector< std::pair<double, double> > points;

    for(int i = 0;i < 7;i++){
        auto pair = std::make_pair(x[i],y[i]);
        points.push_back(pair);

        std::cout << "x: " << x[i] << ", y: " << y[i] << std::endl;
    }

    writePointsFile(points, "Points.dat");


    std::cout << "Building Interpolants ..." << std::endl;
    polynomialInterpolation(x, y, size);
    naturalSplineInterpolation(x, y, size);
    cubicSplineInterpolation(x, y, size);



}
