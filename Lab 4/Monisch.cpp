#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <algorithm>


// Caclulates the y value of the chebychev polynom
double calculatePoint(double point){
  /*
  double first = point;
  double third = first*point*point;
  double fifth = third*point*point;
  double seventh = fifth*point*point;

  return 64*seventh - 112*fifth + 56*third - 7*first;
  */

 return fabs(cos(7.0*acos(point))/64);
}

// Calculate the extrema of T7, index should be in range 0...7
double calculateExtrema(int index){
  return cos( (index*1.0*M_PI)/7 );
}


int main (void){
  // calculate extremas and write them to a file, also store them in a vector
  std::vector<double> extremas;

  std::cout << std::endl;
  std::cout << "EXTREMA" << std::endl;
  std::cout << "-------" << std::endl;
  // Calculate zeros
  for(int i = 0; i <= 7;i++){
    double extrema = calculateExtrema(i);
    extremas.push_back(extrema);
    std::cout << "T7(" << extrema << ") = " << calculatePoint(extrema) << "  (should be : " <<  1.0/64.0  << ")"<< std::endl;;
  }


  return 0;
}
