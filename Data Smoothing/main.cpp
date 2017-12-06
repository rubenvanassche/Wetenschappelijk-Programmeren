#include <iostream>
#include <math.h>
#include <list>
#include "PointsWriter.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector_double.h>
#include "Utils.h"

const double pi = 3.14159265359;

enum Basis{MONOMIAL, CHEBYSHEV, LEGENDRE};

// Should be sampled in -1, 1
double yFunction(double x){
    return 1 + 2*x + pow(x, 2);
}

double eval(gsl_vector* v, double x, Basis basis){
    double out = 0;

    if(basis == MONOMIAL) {
        out += gsl_vector_get(v, 0)*1;
        out += gsl_vector_get(v, 1)*x;
        out += gsl_vector_get(v, 2)*pow(x, 2);
    }else if(basis == LEGENDRE){
        out += gsl_vector_get(v, 0)*1;
        out += gsl_vector_get(v, 1)*gsl_sf_legendre_P1(x);
        out += gsl_vector_get(v, 2)*gsl_sf_legendre_P2(x);
    }else if(basis == CHEBYSHEV){
        out += gsl_vector_get(v, 0)*1;
        out += gsl_vector_get(v, 1)*x;
        out += gsl_vector_get(v, 2)*((2*pow(x, 2)) - 1);
    }

    return out;
}

double* getXPoints(){
    double* points = new double[201];

    for(int i = 0;i < 201;i++){
        points[i] = -1.0 +((double)i/100.0);
    }

    return points;
}

double* getYPoints(const double* xPoints){
    double* yPoints = new double[201];

    for(int i = 0;i < 201;i++){
        yPoints[i] = yFunction(xPoints[i]);
    }

    return yPoints;
}

double* getNoisyPoints(const double* points){
    gsl_rng_env_setup();

    const gsl_rng_type * T = gsl_rng_mt19937;
    gsl_rng * r = gsl_rng_alloc(T);

    double* noisyPoints = new double[201];

    for(int i = 0;i < 201;i++){
        double noise = gsl_rng_uniform(r)*2 - 1;
        noisyPoints[i] = points[i] + noise;
    }

    gsl_rng_free(r);

    return noisyPoints;
}

gsl_vector* smooth(const double* originalPoints, const double* xPoints, const double* yPoints, Basis basis){
    gsl_matrix* A = gsl_matrix_alloc(201, 3);
    gsl_matrix* A2 = gsl_matrix_alloc(201, 3);
    gsl_vector* b = gsl_vector_alloc(3);
    gsl_vector* tau = gsl_vector_alloc(3);
    gsl_vector* y = gsl_vector_alloc(201);
    gsl_vector* original = gsl_vector_alloc(201);
    gsl_vector* residual = gsl_vector_alloc(201);
    gsl_vector* residual2 = gsl_vector_alloc(201);
    gsl_vector* error = gsl_vector_alloc(201);

    // Fill in the values
    for(int i = 0;i < 201;i++){
        // Set matrix
        if(basis == MONOMIAL) {
            gsl_matrix_set(A, i, 0, 1);
            gsl_matrix_set(A, i, 1, xPoints[i]);
            gsl_matrix_set(A, i, 2, pow(xPoints[i], 2));
        }else if(basis == LEGENDRE){
            gsl_matrix_set(A, i, 0, 1);
            gsl_matrix_set(A, i, 1, gsl_sf_legendre_P1(xPoints[i]));
            gsl_matrix_set(A, i, 2, gsl_sf_legendre_P2(xPoints[i]));
        }else if(basis == CHEBYSHEV){
            gsl_matrix_set(A, i, 0, 1);
            gsl_matrix_set(A, i, 1, xPoints[i]);
            gsl_matrix_set(A, i, 2, 2*pow(xPoints[i], 2) - 1);
        }

        // Set y vector
        gsl_vector_set(original, i, originalPoints[i]);
        gsl_vector_set(y, i, yPoints[i]);
    }

    gsl_matrix_memcpy(A2,A);

    // Check the condition number
    std::cout << "Condition Number : " << conditionNumber(A) << std::endl;

    gsl_linalg_QR_decomp(A, tau);
    gsl_linalg_QR_lssolve(A, tau, y, b, residual);

    if(basis == MONOMIAL) {
        std::cout << gsl_vector_get(b, 0) << " + " << gsl_vector_get(b, 1) << "x + " << gsl_vector_get(b, 2) << "x^2" << std::endl;
    }else if(basis == LEGENDRE){
        std::cout << gsl_vector_get(b, 0) << " + " << gsl_vector_get(b, 1) <<  "x + " << gsl_vector_get(b, 2) << "*0.5*(3x^2-1)" << std::endl;
    }else if(basis == CHEBYSHEV){
        std::cout << gsl_vector_get(b, 0) << " + " << gsl_vector_get(b, 1) << "x + " << gsl_vector_get(b, 2) << "*(2x^2-1)" << std::endl;
    }

    std::cout << "Euclidean norm residual(with noise) vector: " << gsl_blas_dnrm2(residual) << std::endl;



    // Check the error
    for(int i = 0; i < 201;i++){
        double xi = xPoints[i];
        double yi = yFunction(xi);
        double ei = eval(b, xi, basis);

        gsl_vector_set(error, i, fabs(yi-ei));
    }

    std::cout << "Euclidean norm error vector: " << gsl_blas_dnrm2(error) << std::endl;

    // Calculate residual(without noise)
    residualVector(A2, b, original, residual2);
    std::cout << "Euclidean norm residual(without noise) vector: " << gsl_blas_dnrm2(residual2) << std::endl;


    gsl_vector_free(error);
    gsl_vector_free(residual);
    gsl_vector_free(residual2);
    gsl_vector_free(y);
    gsl_vector_free(original);
    gsl_vector_free(tau);
    gsl_matrix_free(A2);
    gsl_matrix_free(A);

    return b;
}





int main() {
    std::cout.precision(10);

    const double* xPoints = getXPoints();
    const double* yPoints = getYPoints(xPoints);
    const double* noisyPoints = getNoisyPoints(yPoints);

    // Write original points and noisy points
    PointsWriter originalPW("data/original.dat");
    originalPW.header("x", "y");

    PointsWriter noisyPW("data/noisy.dat");
    noisyPW.header("x", "y");
    for(int i = 0;i < 201;i++){
        originalPW.line(xPoints[i], yPoints[i]);
        noisyPW.line(xPoints[i], noisyPoints[i]);
    }

    // Now let's smooth
    std::cout << "---------------------" << std::endl << "Original Points" << std::endl << "---------------------" << std::endl;
    gsl_vector* oSmoothed = smooth(yPoints, xPoints, yPoints, MONOMIAL);
    std::cout << "---------------------" << std::endl << "Noised Points" << std::endl << "---------------------" << std::endl;
    gsl_vector* mSmoothed = smooth(yPoints, xPoints, noisyPoints, MONOMIAL);
    std::cout << "---------------------" << std::endl << "Noised Points Legendre" << std::endl << "---------------------" << std::endl;
    gsl_vector* lSmoothed = smooth(yPoints, xPoints, noisyPoints, LEGENDRE);
    std::cout << "---------------------" << std::endl << "Noised Points Chebyshev" << std::endl << "---------------------" << std::endl;
    gsl_vector* cSmoothed = smooth(yPoints, xPoints, noisyPoints, CHEBYSHEV);

    // Calculate points
    PointsWriter OLSPW("data/OLS.dat");
    OLSPW.header("x", "y");

    PointsWriter MnoisyLSPW("data/MnoisyLS.dat");
    MnoisyLSPW.header("x", "y");

    PointsWriter LnoisyLSPW("data/LnoisyLS.dat");
    LnoisyLSPW.header("x", "y");

    PointsWriter CnoisyLSPW("data/CnoisyLS.dat");
    CnoisyLSPW.header("x", "y");


    for(double xi = -1;xi < 1;xi += 0.001){
        OLSPW.line(xi, eval(oSmoothed, xi, MONOMIAL));
        MnoisyLSPW.line(xi, eval(mSmoothed, xi, MONOMIAL));
        LnoisyLSPW.line(xi, eval(lSmoothed, xi, LEGENDRE));
        CnoisyLSPW.line(xi, eval(cSmoothed, xi, CHEBYSHEV));
    }


    gsl_vector_free(oSmoothed);
    gsl_vector_free(mSmoothed);
    gsl_vector_free(lSmoothed);
    gsl_vector_free(cSmoothed);
    delete[] noisyPoints;
    delete[] yPoints;
    delete[] xPoints;
}
