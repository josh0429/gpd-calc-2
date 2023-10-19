#include <cmath>
#include <iostream>

using namespace std;

double gpdH(double X,double zeta,double t) {

  double pi = 3.141592653;
  
  double M = 0.9383;
  double m = 1;
  double Mx = 2.;
  double N = 1.;
  
  double t0 = -zeta*zeta*M*M/(1.-zeta);
  double dT = sqrt((t0-t)*(1.-zeta));

  double H = 0.;
  
  for(double kT = 0; kT <= 5.; kT += 0.001) {

    double MMM = X*M*M - (X/(1.-X))*Mx*Mx - 1./(1.-X) * kT*kT;
    double MMMp = (X-zeta)/(1.-zeta)*M*M - (X-zeta)/(1.-zeta)*Mx*Mx - (1-zeta)/(1-X)*kT*kT;
    H += 4*kT*(kT*kT)/(MMM*MMMp);
  }

  return H*0.001;
  
}

int main() {

  double zeta = 0;
  double t = 0;

  double ff = 0.;
  
  for(double X = 0.01; X <=0.99; X += 0.01) {
    cout << X << " " << gpdH(X,zeta,t) << endl;
    ff += gpdH(X,zeta,t);
  }

}
