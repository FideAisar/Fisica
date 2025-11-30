#include "mylibrary.h"

double to(double base, int pot) {

  double num = 1.;
  for(int i=0; i<pot; i++) {
    num = num * base;
  }

  return num;
}


double ln(double x, int iterations = 100) {
  double y = (x - 1) / (x + 1);
  double y_squared = y * y;
  double result = 0.0;
  double term = y;
    
  for (int n = 1; n <= iterations; n += 2) {
    result += term / n;
    term *= y_squared;
  }
    
  return 2 * result;
}
