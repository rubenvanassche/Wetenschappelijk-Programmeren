#include <stdio.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <iostream>
#include "util.h"
#include <math.h>
#include <vector>
#include <utility>
#include <exception>
#include <limits>

class Point{
public:
  double x;
  double y;

  Point(double fx, double fy){
    this->x = fx;
    this->y = fy;
  };
};

typedef std::vector<Point> ValueList;

double generateTcoefficientbase2(int i, double x){
  if(i == 0){
    return 1.0;
  }else if(i == 1){
    return x;
  }else if(i == 2){
    return pow(x, 2);
  }else if(i == 3){
    return pow(x, 3);
  }else{
    throw std::runtime_error("Invalid i");
  }
}

double generateTcoefficientbase1(int i, double x){
  if(i == 0){
    return 1.0;
  }else if(i == 1){
    return x;
  }else if(i == 2){
    return 2*pow(x, 2) - 1;
  }else if(i == 3){
    return 4*pow(x, 3) - 2*pow(x, 2) - 2*x;
  }else{
    throw std::runtime_error("Invalid i");
  }
}

gsl_matrix* generatebase1matrix(ValueList* values){
  gsl_matrix* matrix = gsl_matrix_alloc (20, 4);

  for(int i = 0; i < 20;i++){
    for(int j = 0; j < 4;j++){
      gsl_matrix_set(matrix, i, j, generateTcoefficientbase1(j, values->at(i).x));
    }
  }

  return matrix;
}

gsl_matrix* generatebase2matrix(ValueList* values){
  gsl_matrix* matrix = gsl_matrix_alloc (20, 4);

  for(int i = 0; i < 20;i++){
    for(int j = 0; j < 4;j++){
      gsl_matrix_set(matrix, i, j, generateTcoefficientbase2(j, values->at(i).x));
    }
  }

  return matrix;
}

gsl_vector* generateyvector(ValueList* values){
  gsl_vector* vector = gsl_vector_alloc(20);

  for(int i = 0;i < 20; i++){
    gsl_vector_set(vector, i, values->at(i).y);
  }

  return vector;
}

double calculateConditionNumber(gsl_matrix* matrix){
  gsl_vector* vector = gsl_vector_alloc(4);
  gsl_matrix* tempmatrix = gsl_matrix_alloc(4,4);
  gsl_vector* work = gsl_vector_alloc(4);

  gsl_linalg_SV_decomp(matrix, tempmatrix, vector, work);

  double min = std::numeric_limits<double>::max();
  double max = -std::numeric_limits<double>::min();
  for(int i = 0; i < 4;i++){
    double result = gsl_vector_get(vector, i);
    if(result > max){
      max = result;
    }

    if(result < min){
      min = result;
    }
  }

  gsl_vector_free(vector);
  gsl_matrix_free(tempmatrix);
  gsl_vector_free(work);

  return max/min;
}

ValueList*  initializeValues(){
  ValueList* values = new ValueList();
  values->push_back(Point(0.0, -0.8 ));
  values->push_back(Point(0.6, -0.34));
  values->push_back(Point(1.5, 0.59));
  values->push_back(Point(1.7, 0.59));
  values->push_back(Point(1.9, 0.23));
  values->push_back(Point(2.1, 0.1));
  values->push_back(Point(2.3, 0.28));
  values->push_back(Point(2.6, 1.03));
  values->push_back(Point(2.8, 1.5));
  values->push_back(Point(3.0, 1.44));
  values->push_back(Point(3.6, 0.74));
  values->push_back(Point(4.7, -0.82));
  values->push_back(Point(5.2, -1.27));
  values->push_back(Point(5.7, -0.92));
  values->push_back(Point(5.8, -0.92));
  values->push_back(Point(6.0, -1.04));
  values->push_back(Point(6.4, -0.79));
  values->push_back(Point(6.9, -0.06));
  values->push_back(Point(7.6, 1.00));
  values->push_back(Point(8.0, 0.00));

  return values;
}

void scaleValues(ValueList* values){
  for(int i = 0;i < 20; i++){
    values->at(i).x = (values->at(i).x/4.0) - 1;
  }
}

void base1(){
  std::cout << "--------------" << std::endl;
  std::cout << "Base 1" << std::endl;
  std::cout << "--------------" << std::endl;

  ValueList* values = initializeValues();

  gsl_matrix* matrixcondition = generatebase1matrix(values);
  std::cout << "Condition number: " << calculateConditionNumber(matrixcondition) << std::endl;

  gsl_matrix* matrix = generatebase1matrix(values);
  gsl_vector* tau = gsl_vector_alloc(4);
  gsl_vector* y = generateyvector(values);
  gsl_vector* t = gsl_vector_alloc(4);
  gsl_vector* residual = gsl_vector_alloc(20);
  gsl_linalg_QR_decomp(matrix, tau);
  print_vector(tau);
  gsl_linalg_QR_lssolve (matrix, tau, y, t, residual);

  std::cout << "t vector:" << std::endl;
  print_vector(t);

  std::cout << "Minimal Residual vector:" << std::endl;
  print_vector(residual);

  // Cleanup
  delete values;
  gsl_matrix_free(matrixcondition);
  gsl_matrix_free(matrix);
  gsl_vector_free(tau);
  gsl_vector_free(y);
  gsl_vector_free(t);
  gsl_vector_free(residual);
}

void base2(){
  std::cout << "--------------" << std::endl;
  std::cout << "Base 2" << std::endl;
  std::cout << "--------------" << std::endl;

  ValueList* values = initializeValues();

  gsl_matrix* matrixcondition = generatebase2matrix(values);
  std::cout << "Condition number: " << calculateConditionNumber(matrixcondition) << std::endl;

  gsl_matrix* matrix = generatebase2matrix(values);
  gsl_vector* tau = gsl_vector_alloc(4);
  gsl_vector* y = generateyvector(values);
  gsl_vector* t = gsl_vector_alloc(4);
  gsl_vector* residual = gsl_vector_alloc(20);
  gsl_linalg_QR_decomp(matrix, tau);
  print_vector(tau);
  gsl_linalg_QR_lssolve (matrix, tau, y, t, residual);

  std::cout << "t vector:" << std::endl;
  print_vector(t);

  std::cout << "Minimal Residual vector:" << std::endl;
  print_vector(residual);

  // Cleanup
  delete values;
  gsl_matrix_free(matrixcondition);
  gsl_matrix_free(matrix);
  gsl_vector_free(tau);
  gsl_vector_free(y);
  gsl_vector_free(t);
  gsl_vector_free(residual);
}

void base1scaled(){
  std::cout << "--------------" << std::endl;
  std::cout << "Base 1 Scaled" << std::endl;
  std::cout << "--------------" << std::endl;

  ValueList* values = initializeValues();
  scaleValues(values);

  gsl_matrix* matrixcondition = generatebase1matrix(values);
  std::cout << "Condition number: " << calculateConditionNumber(matrixcondition) << std::endl;

  gsl_matrix* matrix = generatebase1matrix(values);
  gsl_vector* tau = gsl_vector_alloc(4);
  gsl_vector* y = generateyvector(values);
  gsl_vector* t = gsl_vector_alloc(4);
  gsl_vector* residual = gsl_vector_alloc(20);
  gsl_linalg_QR_decomp(matrix, tau);
  gsl_linalg_QR_lssolve (matrix, tau, y, t, residual);

  std::cout << "t vector:" << std::endl;
  print_vector(t);

  std::cout << "Minimal Residual vector:" << std::endl;
  print_vector(residual);

  // Cleanup
  delete values;
  gsl_matrix_free(matrixcondition);
  gsl_matrix_free(matrix);
  gsl_vector_free(tau);
  gsl_vector_free(y);
  gsl_vector_free(t);
  gsl_vector_free(residual);
}

void base2scaled(){
  std::cout << "--------------" << std::endl;
  std::cout << "Base 2 Scaled" << std::endl;
  std::cout << "--------------" << std::endl;

  ValueList* values = initializeValues();
  scaleValues(values);

  gsl_matrix* matrixcondition = generatebase2matrix(values);
  std::cout << "Condition number: " << calculateConditionNumber(matrixcondition) << std::endl;

  gsl_matrix* matrix = generatebase2matrix(values);
  gsl_vector* tau = gsl_vector_alloc(4);
  gsl_vector* y = generateyvector(values);
  gsl_vector* t = gsl_vector_alloc(4);
  gsl_vector* residual = gsl_vector_alloc(20);
  gsl_linalg_QR_decomp(matrix, tau);
  gsl_linalg_QR_lssolve (matrix, tau, y, t, residual);

  std::cout << "t vector:" << std::endl;
  print_vector(t);

  std::cout << "Minimal Residual vector:" << std::endl;
  print_vector(residual);

  // Cleanup
  delete values;
  gsl_matrix_free(matrixcondition);
  gsl_matrix_free(matrix);
  gsl_vector_free(tau);
  gsl_vector_free(y);
  gsl_vector_free(t);
  gsl_vector_free(residual);
}

int main(){
  // set preicison
  std::cout.precision(30);
  base1();
  base2();
  base1scaled();
  base2scaled();
}
