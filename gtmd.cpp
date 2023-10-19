#include <iostream>
#include <cmath>
#include <fstream>
#include "gpd.h"

using namespace std;

double oamu(double X) {

  double pi = 3.141592653;
  double M = 0.9383;
  double m = 0.420;
  double Mx = 0.604;
  double Ml = 1.018;
  double alpha = 0.210;
  double N = 2.043;

  double f14 = 0;

  double reg = pow(X,-alpha);
  
  for(double kT = 0; kT <= 5; kT += 0.001) {
    double D = X*(1.-X)*M*M - X*Mx*Mx - Ml*Ml - kT*kT;
    f14 += N*pi*kT*kT*kT*reg*pow(1/D,4)*pow(1-X,3);
  }
  return f14*0.001;
  
}

double oamd(double X) {

  double pi = 3.141592653;
  double M = 0.9383;
  double m = 0.275;
  double Mx = 0.913;
  double Ml = 0.860;
  double alpha = 0.0317;
  double N = 1.570;

  double f14 = 0;

  double reg = pow(X,-alpha);

  for(double kT = 0; kT <= 5; kT += 0.001) {
    double D = X*(1.-X)*M*M - X*Mx*Mx - Ml*Ml - kT*kT;
    f14 += N*pi*kT*kT*kT*reg*pow(1/D,4)*pow(1-X,3);
  }
  return f14*0.001;

}

double oamjmu(double X) {

  double pi = 3.141592653;
  double M = 0.9383;
  double m = 0.420;
  double Mx = 0.604;
  double Ml = 1.018;
  double alpha = 0.210;
  double N = 2.043;

  double f14 = 0;

  double reg = pow(X,-alpha);
  for(double kT = 0.01; kT <=5; kT += 0.01) {
    for(double lT = 0.01; lT <=5; lT += 0.01) {

      double D = X*(1.-X)*M*M - X*Mx*Mx - Ml*Ml - kT*kT;
      double num = (D+lT*lT);
      double denom = 2*lT*pow(D,2)*pow(pow(D-lT*lT,2)-4*kT*kT*lT*lT,1.5);
      f14 += N*kT*kT*kT*(num/denom)*pow(1-X,4)*0.01*0.01;
    }
  }
  return f14;
}


double oamjmd(double X) {

  double pi = 3.141592653;
  double M = 0.9383;
  double m = 0.275;
  double Mx = 0.913;
  double Ml = 0.860;
  double alpha = 0.0317;
  double N = 1.570;
  double f14 = 0;

  double reg = pow(X,-alpha);
  for(double kT = 0.01; kT <=5; kT += 0.01) {
    for(double lT = 0.01; lT <=5; lT += 0.01) {

      double D = X*(1.-X)*M*M - X*Mx*Mx - Ml*Ml - kT*kT;
      double num = (D+lT*lT);
      double denom = 2*lT*pow(D,2)*pow(pow(D-lT*lT,2)-4*kT*kT*lT*lT,1.5);
      f14 += N*kT*kT*kT*(num/denom)*pow(1-X,4)*0.01*0.01;
    }
  }
  return f14;
}


int main() {

  ofstream myfile;
  myfile.open("f14.dat");

  double oam_u = 0.;
  double oam_d = 0.;
  double oam_jmu = 0.;
  double oam_jmd = 0.;
  
  for(double X = 0.001; X<0.999; X+= 0.001) {

    double ddxoamu = (oamu(X+0.001) - oamu(X) )/0.001;
    double ddxoamd = (oamd(X+0.001) - oamd(X) )/0.001;
    double gpdH_u = gpdHu(X,0,0);
    double gpdE_u = gpdEu(X,0,0);
    double gpdH_d = gpdHd(X,0,0);
    double gpdE_d = gpdEd(X,0,0);

    double E2Ttilu = 0;
    double E2Ttild = 0;
    
    for(double Y = X; Y <= 0.9999; Y+=0.0001) {

      cout << X << " " << Y << endl;
      
      E2Ttilu += -(gpdHu(Y,0,0) + gpdEu(Y,0,0))/Y + gpdHutil(Y,0,0)/(Y*Y);
      E2Ttild += -(gpdHd(Y,0,0) + gpdEd(Y,0,0))/Y + gpdHdtil(Y,0,0)/(Y*Y);
    }
    
    E2Ttilu = E2Ttilu*0.0001;
    E2Ttild = E2Ttild*0.0001;
    
    E2Ttilu += -gpdHutil(X,0,0)/X;
    E2Ttild += -gpdHdtil(X,0,0)/X;
    
    myfile << X << " " << oamu(X) << " " << oamd(X) << " " << ddxoamu << " " << ddxoamd << " " << gpdH_u << " " << gpdE_u << " " << gpdH_d << " " << gpdE_d << " " << E2Ttilu <<  " " << E2Ttild << endl;
    //    oam_u += oamu(X)*0.001;
    //    oam_d += oamd(X)*0.001;
    //    oam_jmu += oamjmu(X)*0.001;
    //    oam_jmd += oamjmd(X)*0.001;
  }

  //  cout << oam_u << " " << oam_d << endl;
  return 0;
}
