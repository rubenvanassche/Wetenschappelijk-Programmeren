#include <functional>
#include <iostream>
#include <math.h>
#include <map>
#include <vector>
#include <limits>
#include <list>


double trapezium(std::function<double(double)> f, int k, double a,  double b){
  int n = pow(2, k);
  double h = (b-a)/n;

  double sum = 0.0;
  double x = a;

  for(int i = 1;i < n;i++){
    x += h;
    sum += h*f(x);
  }

  sum += (h/2)*f(a+h);
  sum += (h/2)*f(b-h);

  return sum;
}


double calculateSum(std::function<double(double)> f, int k, double a, double b){
  if(k == 1){
    double n = pow(2, k);
    double h = (b-a)/n;

    return h*f(0);
  }

   double nprev = pow(2, k-1);
   double n = nprev*2;
   double h = (b-a)/n;

   double incremental = h*2;
   double toCalc = (nprev-2)/2;

   double sum;
   double pos = a+h;
   sum += h*f(pos);
   sum += h*f(-pos);
   for(int i = 0;i < toCalc;i++){
     pos += incremental;
     sum += h*f(pos);
     sum += h*f(-pos);
   }

   return sum;
}

double trapeziume(std::function<double(double)> f, int k, double a,  double b){
  double sum = 0.0;
  for(int i = 1;i <= k;i++){
    double n = pow(2, i);

    if(i == 1){
      double h = (b-a)/n;

      sum += (h/2)*(f(a+h)+f(b-h));
      sum += calculateSum(f, i, a, b);
    }else{
      sum /= 2;
      sum += calculateSum(f, i, a, b);
    }
  }

  std::cout << sum << std::endl;

  return sum;
}

int main(){
  std::cout.precision(30);

  int a = -1;
  int b = 1;
  std::function<double(double)> f = [](double x){
    return -log(1 + x)*log(1 - x);
  };

  double check = pow(2, -23);

  double prev = trapeziume(f, 1, a, b);
  int k = 1;
  /*
  std::cout << "CALCULATE integral(-log(1 + x)*log(1 - x))" << std::endl;
  while(true){
    std::cout << "k: " << k << "  n: " << pow(2, k) << "  -> " << prev << std::endl;
    double now = trapeziume(f, k + 1, a, b);

    double fault = fabs((prev-now)/now);
    std::cout << "fault: " << fault << std::endl;
    if(fault <= check){
      std::cout << "SOLUTION:" << now << std::endl;
      return now;
    }else{
      prev = now;
      k++;
    }
  }*/
  trapeziume(f, 27, a, b);
}
