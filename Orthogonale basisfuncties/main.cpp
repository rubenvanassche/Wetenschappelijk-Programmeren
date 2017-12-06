#include <iostream>
#include <math.h>
#include <list>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "PointsWriter.h"

double f(double x){
    return cos(11*x);
}

int main() {
    std::cout.precision(20);

    std::cout << "Hello world" << std::endl;
}