#include <stdio.h>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <algorithm>

class FileWriter{
public:
  std::string filename;
  std::vector< std::pair<double, double>> buffer;

  FileWriter(std::string name){
    this->filename = name;
    // Clean file
    std::ofstream file;
    file.open(this->filename, std::ofstream::out | std::ofstream::trunc);
    file.close();
  }

  void write(double x, double y){
    std::pair<double, double> coordinates = std::make_pair(x,y);
    buffer.push_back(coordinates);
  }

  ~FileWriter(){
    std::ofstream file;
    file.open(this->filename);
    file << "#   X   Y\n";

    for(int i = 0; i < buffer.size();i++){
      file << buffer[i].first << "   " << buffer[i].second << "\n";
    }

    file.close();
  }
};

// Caclulates the y value of the chebychev polynom
double calculatePoint(double point){
  /*
  double first = point;
  double third = first*point*point;
  double fifth = third*point*point;
  double seventh = fifth*point*point;

  return 64*seventh - 112*fifth + 56*third - 7*first;
  */

 return cos(7.0*acos(point));
}

// Caclulate the zeros of T7, index should be in range 1...7
double calculateZero(int index){
  return cos( ((2.0*index - 1.0)*M_PI)/14.0 );
}

// Calculate the extrema of T7, index should be in range 0...7
double calculateExtrema(int index){
  return cos( (index*1.0*M_PI)/7 );
}


int main (void){
  // calculate zeros and write them to a file, also store them in a vector

  FileWriter writerZeros("T7_zeros.dat");
  std::vector<double> zeros;

  std::cout << "ZEROS" << std::endl;
  std::cout << "-----" << std::endl;
  // Calculate zeros
  for(int i = 1; i <= 7;i++){
    double zero = calculateZero(i);
    zeros.push_back(zero);
    writerZeros.write(zero, 0);
    std::cout << "T7(" << zero << ") = " << calculatePoint(zero) << std::endl;
  }

  // calculate extremas and write them to a file, also store them in a vector
  FileWriter writerExtremas("T7_extremas.dat");
  std::vector<double> extremas;

  std::cout << std::endl;
  std::cout << "EXTREMA" << std::endl;
  std::cout << "-------" << std::endl;
  // Calculate zeros
  for(int i = 0; i <= 7;i++){
    double extrema = calculateExtrema(i);
    extremas.push_back(extrema);
    writerExtremas.write(extrema, pow(-1, i));
    std::cout << "T7(" << extrema << ") = " << calculatePoint(extrema) << "  (should be : " <<  pow(-1, i)  << ")"<< std::endl;;
  }

  std::cout << std::endl;
  std::cout << "CONCLUSION" << std::endl;
  std::cout << "---------" << std::endl;

  // Concatenate zeros and extremas
  std::vector<double> points;
  points.insert( points.end(), zeros.begin(), zeros.end() );
  points.insert( points.end(), extremas.begin(), extremas.end() );

  // Find minimum and maximum of zeros and extremas combined
  auto minElement = std::min_element(std::begin(points), std::end(points));
  auto maxElement = std::max_element(std::begin(points), std::end(points));
  double min = *minElement;
  double max = *maxElement;

  FileWriter writerPoints("T7.dat");
  double range = max - min;

  int n = 1000;
  double step = range/1000;

  double x = min;
  for(int i = 0; i <= n; i++){
    double y = calculatePoint(x);
    writerPoints.write(x, y);

    x += step;
  }


  std::cout << "Plot from " << min << " till " << max << std::endl;
  std::cout << n << " points written to T7.dat" << std::endl;
  std::cout << "To generate graph, run: gnuplot T7.gnuplot" << std::endl;


  return 0;
}
