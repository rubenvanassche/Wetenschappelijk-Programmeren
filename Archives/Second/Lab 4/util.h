#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <gsl/gsl_multifit.h>
#include <sstream>
#include <iomanip>

template <typename T> std::string to_pstring(const T a_value, const int n);

void print_matrix(gsl_matrix * m);

void print_vector(gsl_vector * v);

void writePointsFile(std::vector< std::pair<double, double> >& points, std::string filename);

double f(double x);

void calculateBasicPoints(double a, double b, std::string filename);