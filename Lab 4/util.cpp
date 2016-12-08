#include "util.h"


template <typename T>
std::string to_pstring(const T a_value, const int n)
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

double f(double x){
    double e = 2.718281828459045;
    double pi = 3.141592653589793;

    return pow(e, -x)*sin(pi*x);
}

void calculateBasicPoints(double a, double b, std::string filename){
    std::vector< std::pair<double, double> > points;

    for(double i = a;i <= b;i += 0.01){
        points.push_back(std::make_pair(i, f(i)));
    }

    writePointsFile(points, filename);
}