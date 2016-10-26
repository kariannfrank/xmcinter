#include <math.h>

//1D Gaussian
//double gaussian(double mean, double stddev, x)
// args = [x, mean, sigma]
double gaussian(int n, double args[n]){
  double xx = (args[0]-args[1])*(args[0]-args[1]);
  double variance2 = args[2]*args[2]*2.0;
  return exp(-xx/variance2);
}


//2D Gaussian
//double gaussian(double meanx, double meany, double stddevx, double stddevy,x,y)
// args = [x, y, meanx, meany, sigmax, sigmay]
double gaussian2d(int n, double args[n]){
  double xx = (args[0]-args[2])*(args[0]-args[2]);
  double variancex2 = args[4]*args[4]*2.0;
  double yy = (args[2]-args[3])*(args[2]-args[3]);
  double variancey2 = args[5]*args[5]*2.0;
  return exp(-xx/variancex2)*exp(-yy/variancey2);
}
