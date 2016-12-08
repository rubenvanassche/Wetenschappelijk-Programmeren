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



int main (void){
    // Set COUT Precision
    std::cout.precision(20);

    std::cout << "xxx" << std::endl;

}
