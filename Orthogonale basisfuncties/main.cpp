#include <iostream>
#include <math.h>
#include <list>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include "PointsWriter.h"
#include <vector>
#include <string>

struct Coefficients{
    int index = 0;
    double a = 0;
    double b = 0;
};

double f(double x){
    return cos(11*x);
}

void calculateFunction(){
    PointsWriter original("data/original.dat");
    original.header("X", "Y");

    for(double xi = 0;xi <= M_PI;xi += 0.001){
        original.line(xi, f(xi));
    }
}

Coefficients coefficients(std::vector<double> x, std::vector<double> y, unsigned int j){
    unsigned long int n = x.size();

    Coefficients c;
    c.index = j;

    for(unsigned long int k = 0;k < n;k++){
        c.a += y.at(k)*cos(j*x.at(k));
        c.b += y.at(k)*sin(j*x.at(k));
    }

    c.a = (2.0/n)*c.a;
    c.b = (2.0/n)*c.b;

    return c;
}

// Evaluates the interpolant in x
double eval(std::vector<Coefficients>& coeff, int n, double x){
    double output = 0;
    unsigned long int m = coeff.size() - 1;

    output += coeff.at(0).a*(1.0/2.0);
    if(m == 0){
        return output; // Case only a0 should be calculate
    }

    if(2*m + 1 == n or 2*m + 1 < n or 2*m == n){
        // Aproximation(2m + 1 < n) and interpolation(2m + 1 == n) (2*m == n)
        for(int i = 1;i <= m;i++){
            if(2*m == n and i == m){
                output += coeff.at(i).a*cos(i*x)*(1.0/2.0);
                output += 0; // no b following cursus
                continue;
            }

            output += coeff.at(i).a*cos(i*x);
            output += coeff.at(i).b*sin(i*x);
        }
    }else{
        std::cout << "Something went wrong! no interpolation or approximation"  << std::endl;
        return 0;
    }


    return output;
}

std::string interpolant(std::vector<Coefficients>& coeff, int n){
    std::string output = "";
    unsigned long int m = coeff.size() - 1;
    double border = 0.00000001;

    output += std::to_string(coeff.at(0).a*(1.0/2.0));
    output += " ";
    if(m == 0){
        return output; // Case only a0 should be calculate
    }

    // Aproximation(2m + 1 < n) and interpolation(2m + 1 == n) (2*m == n)
    for(int i = 1;i <= m;i++){
        if(2*m == n and i == m){
            double a =coeff.at(i).a*(1.0/2.0);

            if(a > border){
                output += " + ";
                output += std::to_string(a);
                output += "cos(" ;
                output += std::to_string(i);
                output +=  "x)";
            }

            continue;
        }

        double a = coeff.at(i).a;
        double b = coeff.at(i).b;

        if(a > border){
            output += " + ";
            output += std::to_string(a);
            output += "cos(" ;
            output += std::to_string(i);
            output +=  "x)";
        }

        if(b > border){
            output += " + ";
            output += std::to_string(b);
            output += "sin(" ;
            output += std::to_string(i);
            output +=  "x)";
        }
    }


    return output;
}

std::string coefficients(std::vector<Coefficients> coeff){
    std::string output;

    for(auto it : coeff){
        output += "---------------\n";
        output += "a" + std::to_string(it.index) + " : " + std::to_string(it.a) + "\n";
        output += "b" + std::to_string(it.index) + " : " + std::to_string(it.b) + "\n";
    }

    return output;
}

void rescalePoints(std::vector<double>& points){
    double scale = double(points.size() - 1)/(double)(points.size());
    scale *= 2;

    for(int i = 0;i < points.size();i++){
        points[i] = points.at(i)*scale;
    }
}

void fit(){
    // Calculate points
    // Calculate examples from cursus
    //std::vector<double> x = {0, (2.0*M_PI)/5.0,(4.0*M_PI)/5.0,(6.0*M_PI)/5.0, (8.0*M_PI)/5.0 };
    //std::vector<double> y = {1,3,2,0,-1};
    //std::vector<double> x = {0, (2.0*M_PI*1.0)/8.0, (2.0*M_PI*2.0)/8.0, (2.0*M_PI*3.0)/8.0, (2.0*M_PI*4.0)/8.0, (2.0*M_PI*5.0)/8.0, (2.0*M_PI*6.0)/8.0, (2.0*M_PI*7.0)/8.0};
    //std::vector<double> y = {1, 1, 1, 1, 0, 0, 0, 0};
    // Original points
    std::vector<double> x = {0, M_PI/2, M_PI};
    std::vector<double> y = {f(0), f(M_PI/2), f(M_PI)};
    // Calculate more then 3 points
    /*
    std::vector<double> x;
    std::vector<double> y;
    int size = 23;
    for(int i = 0;i < size;i++){
        double xi = 0 + ((M_PI)/(double)(size-1))*i;

        x.push_back(xi);
        y.push_back(f(xi));
    }
     */
    // Calculate 4 points in the original function
    //std::vector<double> x = {0, M_PI/2, M_PI, (3.0/2.0)*M_PI};
    //std::vector<double> y = {f(0), f(M_PI/2), f(M_PI), f((3.0/2.0)*M_PI)};

    // Transform the points to [0, 2pi]
    //rescalePoints(x);

    int n = x.size();

    // Scaling factor to scale interpolant back to original interval
    double scale = double(n)/(double)(n-1);
    scale *= 1.0/2.0;

    PointsWriter points("data/points.dat");
    points.header("X", "Y");

    for(int i = 0;i < x.size();i++){
        points.line(x.at(i), y.at(i));

        // Function when rescaling points back to original interval
        //points.line(scale*x.at(i), y.at(i));
    }

    // Calculate coefficients
    std::vector<Coefficients> coeff;
    int m = 1;


    for(int i = 0;i <= m;i++){
        coeff.push_back(coefficients(x, y, i));
    }

    PointsWriter interpolated("data/interpolated.dat");
    interpolated.header("X", "Y");

    for(double xi = 0;xi <= 2*M_PI;xi += 0.001){
        interpolated.line(xi, eval(coeff, n, xi));

        // Function when rescaling points back to original interval
        //interpolated.line(scale*xi, eval(coeff, n, xi));
    }

    std::cout << interpolant(coeff, n) << std::endl;

    std::cout << coefficients(coeff) << std::endl;

}

int main() {
    std::cout.precision(20);

    calculateFunction();
    fit();

}