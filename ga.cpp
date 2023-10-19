#include <iostream>
#include <cmath>
#include "ga.h"

using namespace std;

//t is positive

double ga(double t, int toret) {

  double M = .9383;

  double gv = 1.267/pow((1.+t/pow(1.026,2)),2);
  double gs = 0.585/pow((1.+t/pow(1.2,2)),2);

  if(toret == 0) {
    return gv;
  }
  else if(toret == 1) {
    return gs;
  }
}
