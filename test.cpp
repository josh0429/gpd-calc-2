//Necessary includes
#include <iostream>
#include <cmath>
#include <fstream>

//Include Header File
#include "gpd.h"

//Necessary for output
using namespace std;

double gpdHu(double X, double zeta, double t) {
  //***************************//
  //      PARAMETERS
  const double pi = 3.14159265358979323846;
  double M = .9383;
  //10.1103/PhysRevC.88.065206 Parameters
  double m = .420;
  double Mx = .604; // M_X
  double Ml = 1.018; // M_\Lambda
  double alpha = .210; // Regge term alpha
  double alphap = 2.448; // Regge term alpha^\prime
  double p = 0.620; // Regge term p
  //  double alphap = 1.814;
  //double p = 0.449;
  double N = 2.043; // Normalization N
  double beta = 0.;
  double a = 1.;
  //***************************//

  //**********************************************************
  //      ODIL VARIABLES
  double Xp = (X-zeta)/(1.-zeta); // X^\prime
  double t0 = -zeta*zeta*M*M/(1.-zeta);
  double dT = sqrt((t0-t)*(1.-zeta)); // \Delta_\perp
  double d = dT * (1. - Xp);
  double MoX = X*Mx*Mx + (1. - X)*Ml*Ml - X*(1. - X)*M*M;
  double MoXp = Xp*Mx*Mx + (1. - Xp)*Ml*Ml - Xp*(1. - Xp)*M*M;
  double S = MoX - MoXp + d*d;
  double R = pow((MoX - MoXp + d*d), 2) + 4*d*d*MoXp;
  //**********************************************************

  //**********************************************************
  //      ODIL INTEGRALS
  double ASarg = (sqrt(R)/ MoX) * (1. + (S / (2*MoXp)));
  double AS = asinh(ASarg);
  double prefac = (S / (R*R));
  double I0 = prefac * ((S / (4*MoXp*(MoXp + d*d))) + (1. / (2*MoX)) + ((MoXp/MoX) \
   - ((MoX + 3*d*d)/(MoXp + d*d))) / S + (3 / (2*sqrt(R)))*AS);
  double I1 = prefac * ((S / (4 * MoXp)) + 1.5 - ((2 * d*d) / S) + ((((4*d*d*MoXp)/S) + S - 3*MoX) / (2 * sqrt(R)))*AS);
  double I2 = prefac * (((d*d*S)/(4*MoXp)) - 1.25*S + 1.5*(3*d*d - MoXp) + (2*d*d*(2*MoXp - d*d)) / S - (MoX / sqrt(R)) * ((4*d*d*MoXp / S) + S - 1.5*MoX) * AS);
  //**********************************************************

  double Hu1 = (m + M*X) * (m + M*Xp) * (MoXp + d*d)*I0;
  double Hu2 = (MoXp + (m + M*X)*(m + M*Xp) - d*d)*I1;
  double Hu = 2*pi*N*(1. - X)*pow((1. - Xp),2)*(1 - zeta / 2) * (Hu1 + Hu2 +I2);
  return Hu + ((zeta*zeta)/(4*(1.-zeta)))*gpdEu(X,zeta,t);
}

double gpdEu(double X, double zeta, double t) {
  //***************************//
  //      PARAMETERS
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .604;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 2.835;
  double p = .969;
  double beta = 0.;
  double a = 1;
  double N = 1.803;
  //***************************//

  //**********************************************************
  //      ODIL VARIABLES
  double Xp = (X-zeta)/(1.-zeta); // X^\prime
  double t0 = -zeta*zeta*M*M/(1.-zeta);
  double dT = sqrt((t0-t)*(1.-zeta)); // \Delta_\perp
  double d = dT * (1. - Xp);
  double MoX = X*Mx*Mx + (1. - X)*Ml*Ml - X*(1. - X)*M*M;
  double MoXp = Xp*Mx*Mx + (1. - Xp)*Ml*Ml - Xp*(1. - Xp)*M*M;
  double S = MoX - MoXp + d*d;
  double R = pow((MoX - MoXp + d*d), 2) + 4*d*d*MoXp;
  //**********************************************************
  
  //**********************************************************
  //      ODIL INTEGRALS
  double ASarg = (sqrt(R)/ MoX) * (1. + (S / (2*MoXp)));
  double AS = asinh(ASarg);
  double prefac = (S / (R*R));
  double I0 = prefac * ((S / (4*MoXp*(MoXp + d*d))) + (1. / (2*MoX)) + ((MoXp/MoX) \
   - ((MoX + 3*d*d)/(MoXp + d*d))) / S + (3 / (2*sqrt(R)))*AS);
  double I1 = prefac * ((S / (4 * MoXp)) + 1.5 - ((2 * d*d) / S) + ((((4*d*d*MoXp)/S) + S - 3*MoX) / (2 * sqrt(R)))*AS);
  //**********************************************************

  double Eu1 = 2*M*(m + M*X)*(MoXp + d*d)*I0;
  double Eu2 = (2*M*(m + M*X) - 4*M*M*(X - Xp))*I1;
  double Eu = 2*pi*N*(1. - zeta/2)*pow((1. - zeta),2)*pow((1. - Xp), 4)*(Eu1 + Eu2);
  
  return Eu;
}

int main(){
    double r = gpdHu(0.4, 0.2, -0.3);
    cout << r << endl;
}