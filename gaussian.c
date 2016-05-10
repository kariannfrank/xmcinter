#include <math.h>

//double gaussian(double mean, double stddev, x)
// args = [x, mean, sigma]
double gaussian(int n, double args[n]){
  double xx = (args[0]-args[1])*(args[0]-args[1]);
  double variance2 = args[2]*args[2]*2.0;
  return exp(-xx/variance2);
}
