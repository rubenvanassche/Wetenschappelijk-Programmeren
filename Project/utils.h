//
// Created by Ruben Van Assche on 20/08/17.
//

#ifndef PROJECT_UTILS_H_H
#define PROJECT_UTILS_H_H

#include <math.h>
#include <list>
#include <fstream>
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include  <random>
#include  <iterator>
#include <string>
#include <sstream>

namespace Utils{
    const double PI = 3.141592653589793238463;

    void writePointsFile(std::list< std::pair<double, double> >& points, std::string filename){
        std::ofstream file;
        file.open(filename);

        file << "#  X   Y" << std::endl;
        for(auto i: points){
            file << i.first << "   " << i.second << std::endl;
        }

        file.close();
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

    // https://stackoverflow.com/questions/6942273/get-random-element-from-container
    template<typename Iter, typename RandomGenerator>
    Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
        std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
        std::advance(start, dis(g));
        return start;
    }

    template<typename Iter>
    Iter select_randomly(Iter start, Iter end) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        return select_randomly(start, end, gen);
    }

    template <typename T>
    std::string to_pstring(const T a_value, const int n = 6)
    {

        std::ostringstream out;
        out.precision(6);
        out  << a_value;
        return out.str();
    }
};


#endif //PROJECT_UTILS_H_H
