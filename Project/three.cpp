#include <iostream>
#include <list>
#include <vector>
#include "utils.h"

// GSL
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>

double func (double x[], size_t dim, void * p){
    if (dim != 2)
    {
        fprintf (stderr, "error: dim != 2");
        abort ();
    }

    if(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2) <= 0.25){
        return fabs(log(x[0] + x[1]));
    }else{
        return 0.0;
    }
}

// x 0 -> 1
// y 0 -> 1
// MAPLE
// with(VectorCalculus);
// int(abs(ln(x+y)), [x, y] = Circle(⟨0.5,0.5⟩, .5));
// .2550197403

double calculate(int samples){
    double res, err;

    double xl[2] = {0.0, 0.0};
    double xu[2] = {1.0, 1.0};

    const gsl_rng_type *T;
    gsl_rng *r;

    gsl_monte_function G = { &func, 2, 0 };

    gsl_rng_env_setup ();

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_plain_state *s = gsl_monte_plain_alloc (2);
    gsl_monte_plain_integrate (&G, xl, xu, 2, samples, r, s, &res, &err);
    gsl_monte_plain_free (s);
    gsl_rng_free (r);

    return res;
}

int main() {
    std::cout.precision(20);

    int samples = 10;
    while(true){
        double result = calculate(samples);
        std::cout << "Samples: " << samples << ", result: " << result << std::endl;

        int test = result*100;
        if(test == 25){
            break;
        }else{
            samples += 10;
        }
    }

}