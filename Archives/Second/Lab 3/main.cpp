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
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_blas.h>
#include <sstream>
#include <iomanip>

template <typename T>
std::string to_pstring(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}

void print_matrix(gsl_matrix * m){
    int rows = m->size1;
    int columns = m->size2;

    for (int i = 0; i < rows; i++){  /* OUT OF RANGE ERROR */
        for (int j = 0; j <  columns; j++){
            std::cout << gsl_matrix_get (m, i, j) << "   ";
        }
        std::cout << std::endl;
    }
}

void print_vector(gsl_vector * v){
    int rows = v->size;

    for (int i = 0; i < rows; i++){  /* OUT OF RANGE ERROR */
        std::cout << gsl_vector_get (v, i) << std::endl;
    }
}


void writePointsFile(std::vector< std::pair<double, double> >& points, std::string filename){
    std::ofstream file;
    file.open(filename);

    file << "#  X   Y" << std::endl;
    for(auto i: points){
        file << i.first << "   " << i.second << std::endl;
    }

    file.close();
}

std::vector< std::pair<double, double> >* calculatePoints(){
    std::vector< std::pair<double, double> >* points = new std::vector< std::pair<double, double> >;


    for(double x = -1.0;x <= 1.0;x += 0.125){
        double y = 1.0/(1.0+25*pow(x, 2));


        std::pair<double, double> point = std::make_pair(x,y);
        points->push_back(point);
    }

    return points;
}

gsl_matrix* buildMatrixA(int n, int m){
    gsl_matrix* A = gsl_matrix_alloc(m, n);
    double step = 2.0/(m-1);

    for(int i = 0;i < m; i++){
        for(int j = 0;j < n;j++){
            double x = -1 + i*step;
            double y = 1.0/(1.0+25*pow(x, 2));

            gsl_matrix_set(A, i, j, pow(x, j));
        }
    }

    return A;
}

gsl_vector* buildYVector(int m){
    gsl_vector* Y = gsl_vector_alloc(m);
    double step = 2.0/(m-1);

    for(int j = 0;j < m;j++){
        double x = -1 + j*step;
        double y = 1.0/(1.0+25*pow(x, 2));

        gsl_vector_set(Y, j, y);
    }

    return Y;

}

// Generate the coefficients with m points and of grade n
gsl_vector* fit(int n, int m){
    n = n + 1; // Grade 6 needs 7 terms

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

std::string buildEquation(gsl_vector* c){
    std::string output = "Y = ";

    for(int i = 0;i < c->size;i++){
        double value = gsl_vector_get(c, i);

        // Add the sign
        if(i != 0){
            if(value >= 0.0){
                output += " + ";
            }
        }

        // Add value
        output += to_pstring(value, 10);
        output += " ";

        // Add x^i
        if(i > 1){
            output += "X^{";
            output += std::to_string(i);
            output += "} ";
        }else if(i == 1){
            output += "X ";
        }else{}
    }

    return output;
}

double fittedFunction(gsl_vector* c, double x){
    double output = 0.0;

    for(int i = 0;i < c->size;i++){
        double value = gsl_vector_get(c, i);
        double xValue = pow(x, i);

        output += value*xValue;
    }

    return output;
}

std::vector< std::pair<double, double> >* calculateFittedPoints(gsl_vector* c){
    std::vector< std::pair<double, double> >* points = new std::vector< std::pair<double, double> >;

    for(double x = -1.0;x <= 1.0;x += 0.001){
        double y = fittedFunction(c, x);

        std::pair<double, double> point = std::make_pair(x,y);
        points->push_back(point);
    }

    return points;
}



int main (void){
    // Set COUT Precision
    std::cout.precision(20);

    // Write original points to file
    auto originalPoints = calculatePoints();
    writePointsFile(*originalPoints, "points.dat");

    std::vector<int> indices = {6,7,8,10,11,12,13,14,15};
    for(auto i : indices){
        std::cout << "------------------------" << std::endl;
        std::cout << "n = " << i << std::endl;
        std::cout << "------------------------" << std::endl;
        gsl_vector* c = fit(i, 17);

        std::string equation = buildEquation(c);
        std::cout << equation << std::endl;

        std::string filename = "fittedpointsN";
        filename += std::to_string(i);
        filename += ".dat";

        auto points = calculateFittedPoints(c);
        writePointsFile(*points, filename);

        gsl_vector_free(c);
    }



}
