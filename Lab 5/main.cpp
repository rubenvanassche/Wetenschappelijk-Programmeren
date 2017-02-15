#include <functional>
#include <iostream>
#include <math.h>
#include <map>
#include <vector>
#include <limits>
#include <list>



double trapezium(std::function<double(double)> f, int k, double a, double b) {
    int n = pow(2, k);
    double h = (b - a) / n;

    double sum = 0.0;
    double x = a;

    for (int i = 1; i < n; i++) {
        x += h;
        sum += h * f(x);
    }

    sum += (h / 2) * f(a + h);
    sum += (h / 2) * f(b - h);

    return sum;
}


double calculateSum(std::function<double(double)> f, int k, double a, double b) {
    if (k == 1) {
        double n = pow(2, k);
        double h = (b - a) / n;

        return h * f(a + h);
    }

    double nprev = pow(2, k - 1);
    double n = nprev * 2;
    double h = (b - a) / n;

    double incremental = h * 2;
    double toCalc = (nprev - 2) / 2;

    double sum;
    double posa = a + h;
    double posb = b - h;
    sum += h * f(posa);
    sum += h * f(posb);
    for (int i = 0; i < toCalc; i++) {
        posa += incremental;
        posb -= incremental;
        sum += h * f(posa);
        sum += h * f(posb);
    }

    return sum;
}

double trapeziume(std::function<double(double)> f, int k, double a, double b) {
    double sum = 0.0;

    for (int i = 1; i <= k; i++) {
        double n = pow(2, i);

        if (i == 1) {
            double h = (b - a) / n;

            sum += (h/2) * (f(a) + f(b));
            sum += calculateSum(f, i, a, b);
        } else {
            sum /= 2;
            sum += calculateSum(f, i, a, b);
        }
    }


    return sum;
}

// Calculates the trapeoid rule till the relative error is smaller then maxError
int trapeziumewitherror(std::function<double(double)> f, double a, double b, double maxError){
    int k = 1;
    double prevsum = 0.0;
    double sum = 0.0;

    while(true){
        // Set prevsum equal to the previous sum
        prevsum = sum;

        double n = pow(2, k);

        if (k == 1) {
            double h = (b - a) / n;

            sum += (h/2) * (f(a) + f(b));
            sum += calculateSum(f, k, a, b);
        } else {
            sum /= 2;
            sum += calculateSum(f, k, a, b);
        }

        // Calculate error
        double error = (prevsum-sum)/(sum);

        if(fabs(error) <= maxError){
            std::cout << "Found integration using " << n << " intervals(n = " << k;
            std::cout << ") with an error of " << error << " and solution: " << sum << std::endl;
            break;
        }

        // Raise intervals
        std::cout << "Using " << n << " intervals: " << sum << " error: " << error << std::endl;
        k++;
    }

    return k;
}

int main() {
    std::cout.precision(20);

    int a = 1;
    int b = 2;
    std::function<double(double)> f = [](double x) {
        return 1.0 / x;
    };

    int k = trapeziumewitherror(f, 1, 2, pow(2, -40));

    std::cout << "Calculate with n = " << k << " using basic trapezoid rule" << std::endl;
    std::cout << trapezium(f, k, 1, 2);
}
