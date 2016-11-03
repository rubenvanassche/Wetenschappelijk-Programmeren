#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include <vector>
#include <utility>
#include <fstream>
#include <limits>
#include <string>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_double.h>

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

gsl_matrix* vectorToMatrix(gsl_vector* v){
    gsl_matrix* m = gsl_matrix_alloc(v->size, 1);

    for(int i = 0;i < v->size;i++){
        gsl_matrix_set(m, i, 0, gsl_vector_get(v, i));
    }

    return m;
}

gsl_vector* matrixToVector(gsl_matrix* m){
    gsl_vector* v = gsl_vector_alloc(m->size1);


    for(int i = 0;i < m->size1;i++){
        gsl_vector_set(v, i, gsl_matrix_get(m, i, 0));
    }

    return v;
}

class HilbertMatrix{
public:
    int n = 0;
    int s = 0;
    double c = 0.0;
    double determinant = 0.0;
    double conditionNumber = 0.0;
    double mNorm = 0.0;
    gsl_matrix* m = nullptr;
    gsl_matrix* mLU = nullptr;
    gsl_matrix* mLUS = nullptr;
    gsl_permutation* p = nullptr;
    gsl_vector* x = nullptr;
    gsl_vector* y = nullptr;


    HilbertMatrix(int n){
        this->n = n;

        // Initiate Hilbert Matrix
        this->m = gsl_matrix_alloc (n, n);
        for(int i = 1;i <= n;i++){
            for(int j = 1;j <= n;j++){
                double hilbert = 1.0/(i+j-1);
                gsl_matrix_set(this->m, i - 1, j - 1, hilbert);
            }
        }

        // Do an LU Decomposition
        this->mLU = gsl_matrix_alloc (n, n);

        gsl_matrix_memcpy(this->mLU, this->m);

        // Do an LU decompostion
        this->p = gsl_permutation_alloc(n);
        gsl_linalg_LU_decomp(this->mLU, this->p, &(this->s));
    }

    double getDeterminant(){
        this->determinant = gsl_linalg_LU_det(this->mLU, this->s);
        return this->determinant;
    }

    double getConditionNumber(){
        gsl_vector* vector = gsl_vector_alloc(this->n);
        gsl_matrix* tempmatrix = gsl_matrix_alloc(this->n,this->n);
        gsl_matrix* tempmatrix2 = gsl_matrix_alloc(this->n,this->n);
        gsl_vector* work = gsl_vector_alloc(this->n);

        gsl_matrix_memcpy(tempmatrix, this->m);

        gsl_linalg_SV_decomp(tempmatrix, tempmatrix2, vector, work);

        double min = std::numeric_limits<double>::max();
        double max = -std::numeric_limits<double>::min();
        for(int i = 0; i < this->n;i++){
            double result = gsl_vector_get(vector, i);
            if(result > max){
                max = result;
            }

            if(result < min){
                min = result;
            }
        }

        this->mNorm = max;
        this->conditionNumber = max/min;
        return this->conditionNumber;
    }

    double getC(){
        double ULP = std::numeric_limits<double>::epsilon();

        // x* - x
        gsl_vector* xx = gsl_vector_alloc(this->n);
        gsl_vector_memcpy(xx, this->solveViaMaple());
        gsl_vector_sub(xx, this->solveViaGEPP());
        double xxNorm = gsl_blas_dnrm2(xx);

        // K(A)*x**ULP
        double xNorm = gsl_blas_dnrm2(this->solveViaGEPP());

        this->c = xxNorm/(this->conditionNumber*xNorm*ULP);


        std::cout << xxNorm << " <= " << this->c*this->conditionNumber*xNorm*ULP << std::endl;

        return this->c;
    }

    double getResidu(){
        double ULP = std::numeric_limits<double>::epsilon();

        // Ax*
        gsl_matrix* axm = gsl_matrix_alloc(this->n, 1);
        gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1.0, this->m, vectorToMatrix(this->x), 0.0, axm);
        gsl_vector* ax  = matrixToVector(axm);

        // y - AX*
        gsl_vector* yy = gsl_vector_alloc(this->n);
        gsl_vector_memcpy(yy, this->y);
        gsl_vector_sub(yy, ax);


        double yax = gsl_blas_dnrm2(yy);

        // c*A*x**ULP
        double caxu = this->c*this->mNorm*gsl_blas_dnrm2(this->x)*ULP;

        std::cout << yax << " <= " << caxu << std::endl;

        return yax;
    }


    gsl_vector* solveViaGEPP(){
        // Initiate vectors
        this->y = gsl_vector_alloc(this->n);
        this->x = gsl_vector_alloc(this->n);

        // Inititiate MLUS matrix
        this->mLUS = gsl_matrix_alloc (this->n, this->n);

        for(int i = 0;i < this->n;i++){
            if(i == 0){
                gsl_vector_set(this->y, i, 1.0);
            }else{
                gsl_vector_set(this->y, i, 0.0);
            }
        }

        // Copy the original LU matrix to an temporary matrix
        gsl_matrix_memcpy(this->mLUS, this->mLU);


        // Solve the temporary LU matrix
        gsl_linalg_LU_solve (this->mLUS, this->p, this->y, x);

        return x;
    }

    gsl_vector* solveViaMaple(){
        /*
         * This vectors are pre calculated with Maple
         * using the following commands:
         * m := HilbertMatrix(3, 3)
         * b := Vector[column]([1, 0, 0])
         * x := LinearSolve(m, b)
         */

        // Initiate vectors
        this->x = gsl_vector_alloc(this->n);
        std::vector<double> values;

        if(this->n == 3){
            values.push_back(9);
            values.push_back(-36);
            values.push_back(30);
        }else if(this->n == 6){
            values.push_back(36);
            values.push_back(-360);
            values.push_back(3360);
            values.push_back(-7560);
            values.push_back(7560);
            values.push_back(-2772);
        }else if(this->n == 9){
            values.push_back(81);
            values.push_back(-3240);
            values.push_back(41580);
            values.push_back(-249480);
            values.push_back(810810);
            values.push_back(-1513512);
            values.push_back(1621620);
            values.push_back(-926640);
            values.push_back(218790);
        }else if(this->n == 12){
            values.push_back(144);
            values.push_back(-10296);
            values.push_back(240240);
            values.push_back(-2702700);
            values.push_back(17297280);
            values.push_back(-68612544);
            values.push_back(176432256);
            values.push_back(-299304720);
            values.push_back(332560800);
            values.push_back(-232792560);
            values.push_back(93117024);
            values.push_back(-16224936);
        }else{
            std::cout << "Couldn't find pre calculated vector, using linear solution" << std::endl;
            return this->solveViaGEPP();
        }

        for(int i = 0;i < this->n;i++){
            gsl_vector_set(x, i, values.at(i));
        }

        return x;
    }

    void printMatrix(){
        int rows = this->m->size1;
        int columns = this->m->size2;

        for (int i = 0; i < rows; i++){  /* OUT OF RANGE ERROR */
            for (int j = 0; j <  columns; j++){
                std::cout << gsl_matrix_get (this->m, i, j) << "   ";
            }
            std::cout << std::endl;
        }
    }

    void calc(){
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Hilbert(" << this->n << ")" << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
        this->printMatrix();
        this->solveViaGEPP();
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << "Determinant: " << this->getDeterminant() << std::endl;
        std::cout << "Condition Number: " << this->getConditionNumber() << std::endl;
        std::cout << "--------------- Y Vector ----------------" << std::endl;
        print_vector(this->y);
        std::cout << "--------------- X Vector ----------------" << std::endl;
        print_vector(this->x);
        std::cout << "-----------------------------------------" << std::endl;
        this->getC();
        std::cout << "c: "  << this->c << std::endl;
        std::cout << "Residu:" << std::endl;
        this->getResidu();
        std::cout << std::endl << std::endl << std::endl;
    }
};


int main (void){
    // Set COUT Precision
    std::cout.precision(4);


    for(int i = 3; i < 15; i += 3){
        auto h = HilbertMatrix(i);
        h.calc();
    }
}
