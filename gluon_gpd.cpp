
//Definitely not finished!
#include <iostream>
#include <cmath>
#include <fstream>
#include "gluon_gpd.h"

using namespace std;

double gpdHg(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double Mx = 1.12;
  double Ml = 1.045;
  double alpha = 0.005;
  double N = 1.525 + 0.228;
  double Hg = 0;
  double alphap = 0.275;//- 0.1;
  double p = 0.17;//-0.05;
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-alphap*(pow(1-X,p)*t));
  double betag = 0.;
  double omx2 = (1.-X)*(1.-X);
  double xmz = X-zeta;
  double omx = 1.-X;
  double omxp = (1.-X)/(1.-zeta);
  double xp = (X-zeta)/(1.-zeta);
  for(double kT = 0; kT <= 5; kT+=.001) {
    double kmm = X*(1.-X)*M*M - X*Mx*Mx - (1.-X)*Ml*Ml - kT*kT;
    double kmm2 = kmm*kmm;
    double aaa = xp*omx*M*M-(X-zeta)*Mx*Mx \
      - (1.-zeta)*kT*kT - omxp*omx*dT*dT - omx*Ml*Ml;
    double bbb = -2*kT*dT;
    double Hg1 = -2*pi*N*reg1*reg2*kT*(1.-X)*(1.-X)*(X*xmz*(omx*M-Mx)\
     *(omxp*M-Mx))*(1/(kmm*kmm))*(aaa/(pow(aaa*aaa-omx2*bbb*bbb,1.5)));
    double Hg2 = -2*pi*N*reg1*reg2*kT*(1.-X)*(1.-X)*(1.-zeta+omx2)\
      *kT*kT*(1./(kmm*kmm))*(aaa/(pow(aaa*aaa-omx2*bbb*bbb,1.5)));
    double Hg3 = 2*pi*N*reg1*reg2*kT*(1.-X)*(1.-X)*(1.-zeta+omx2)\
      *kT*dT*omxp*(1./(kmm*kmm))*(bbb/(pow(aaa*aaa-omx2*bbb*bbb,1.5)));
    Hg += Hg1+Hg2+Hg3;
  }
  return Hg*.001 + zeta*zeta/(4*(1.-zeta))*gpdEg(X,zeta,t);
}



double gpdEg(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double Mx = 1.12;
  double Ml = 1.10;
  double alpha = 0.053;
  double N = 3.97 ;
  double Eg = 0;
  double alphap = 0.45;
  double p = -0.2;
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-alphap*(pow(1-X,p)*t));
  double omx2 = (1.-X)*(1.-X);
  double xmz = X-zeta;
  double omx = 1.-X;
  double omxp = (1.-X)/(1.-zeta);
  double xp = (X-zeta)/(1.-zeta);
  for(double kT = 0; kT <= 5; kT+=.001) {
    double mox = X*M*M-Ml*Ml-Mx*Mx*X/(1.-X);
    double moxp = xp*M*M-Ml*Ml-Mx*Mx*xp/(1.-xp);
    double aaa = moxp - kT*kT/(1.-xp) - dT*dT*(1.-xp);
    double bbb = 2*kT*dT;
    double ddd = mox - kT*kT/(1.-X);
    double denom = 1/(ddd*ddd*pow(aaa*aaa-bbb*bbb,1.5));
    Eg += 2*pi*N*kT*(1/(1.-X))*denom*reg1*reg2*(-2*M*(1.-zeta))/(1.-zeta/2.)\
      *(2*kT*kT*(X*(omx*M-Mx)-xp*(1-zeta)*(omxp*M-Mx))-aaa*omxp*X*(omx*M-Mx));
  }
  return Eg*.001;
}

double gpdHgtil(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double Mx = 1.12;
  double Ml = 1.10;
  double alpha = 0.053;
  double N = 1.52;
  double Hgtil = 0;
  for(double kT = 0; kT <= 5; kT+=.001) {
    Hgtil +=2*pi*pow(X,-alpha)*N*kT*(1-X)*(1-X)*((X*X*((1-X)*M-Mx)*((1-X)*M-Mx)\
       +(1-(1-X)*(1-X))*kT*kT)/pow((X*Mx*Mx+(1-X)*Ml*Ml-X*(1-X)*M*M+kT*kT),4)) ;
  }
  return Hgtil*.001;

}

double gpdEgtil(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = 0.9383;
  double Mx = 1.12;
  double Ml = 1.10;
  double alpha = 0.053;
  double N = 1.52;
  double Egtil = 0;
  for (double kT = 0; kT<=5;kT+=0.001){
    double num = 8*M*(pow(1.-X,2)*M-Mx)*kT*kT;
    double denom = pow(X*Mx*Mx + (1.-X)*Ml*Ml - X*(1.-X)*M*M+kT*kT,3);
    Egtil += 2*pi*N*pow(X,-alpha)*kT*pow(1.-X,3)*(num/denom);
  }
  return Egtil*.001;
}
