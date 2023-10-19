//************************************************************
//
//   KELLY form factor parametrization
//
//   AUTHOR: Brandon Kriesten
//
//   Calculates kelly parametrization of elastic form factors
//   from data. Dirac and Pauli form factors for proton
//   and neutron from and their errors
//
//   t-momentum transfer
//   must be POSITIVE!!!
//
//************************************************************

//Necessary includes
#include <iostream>
#include <cmath>

//Header File includes
#include "formfactor.h"

//For output
using namespace std;

//
// Kelly Form factors
//
// t MUST BE POSITIVE!!!
//
// ret tells you which part of the form factor you want
// to return
//
double formfactor(double t,int ret) {

  //Parameters used later!
  double a1_ep = -.24;
  double a1_mp = .12;
  double a1_mn = 2.33;
  
  double b1_ep = 10.98;
  double b1_mp = 10.97;
  double b1_mn = 14.72;
  
  double b2_ep = 12.82;
  double b2_mp = 18.86;
  double b2_mn = 24.20;
  
  double b3_ep = .12;
  double b3_mp = 6.55;
  double b3_mn = 84.1;
  
  double ea1ep = .12;
  double ea1mp = .04;
  double ea1mn = 1.4;
  
  double eb1ep = .19;
  double eb1mp = .11;
  double eb1mn = 1.7;

  double eb2ep = 1.10;
  double eb2mp = .28;
  double eb2mn = 9.8;
  
  double eb3ep = 6.8;
  double eb3mp = 1.2;
  double eb3mn = 41.0;

  double aa = 1.7;
  double bb = 3.3;
  double eaa = .04;
  double ebb = .32;

  double amp = .9383;
  double amup = 2.79;
  double amun = -1.91;

  //Tau
  double tau = t/4/amp/amp;

  //GE,GM proton and neutron
  double gep = (1+a1_ep*tau)/(1+b1_ep*tau+b2_ep*tau*tau+b3_ep*tau*tau*tau);
  double gmp = amup*(1+a1_mp*tau)/(1+b1_mp*tau+b2_mp*tau*tau + b3_mp*tau*tau*tau);
  double gmn = amun*(1+a1_mn*tau)/(1+b1_mn*tau+b2_mn*tau*tau + b3_mn*tau*tau*tau);
  double gen = aa*tau/(1+tau*bb)/pow((1+t/.71),2);

  //Errors GEP
  double e1 = tau/(1+b1_ep*tau+b2_ep*tau*tau + b3_ep*tau*tau*tau);
  double e2 = (1+a1_ep*tau)/pow((1+b1_ep*tau+b2_ep*tau*tau + b3_ep*tau*tau*tau),2)*tau;
  double e3 = (1+a1_ep*tau)/pow((1+b1_ep*tau+b2_ep*tau*tau + b3_ep*tau*tau*tau),2)*tau*tau;
  double e4 = (1+a1_ep*tau)/pow((1+b1_ep*tau+b2_ep*tau*tau + b3_ep*tau*tau*tau),2)*tau*tau*tau;
  double egep = sqrt(pow((e1*ea1ep),2)+pow((e2*eb1ep),2)+pow((e3*eb2ep),2)+pow((e4*eb3ep),2));

  //Errors GMP
  double e1m = tau/(1+b1_mp*tau+b2_mp*tau*tau + b3_mp*tau*tau*tau);
  double e2m = (1+a1_mp*tau)/pow((1+b1_mp*tau+b2_mp*tau*tau + b3_mp*tau*tau*tau),2)*tau;
  double e3m = (1+a1_mp*tau)/pow((1+b1_mp*tau+b2_mp*tau*tau + b3_mp*tau*tau*tau),2)*tau*tau;
  double e4m = (1+a1_mp*tau)/pow((1+b1_mp*tau+b2_mp*tau*tau + b3_mp*tau*tau*tau),2)*tau*tau*tau;
  double egmp = sqrt(pow((e1m*ea1mp),2)+pow((e2m*eb1mp),2)+pow((e3m*eb2mp),2)+pow((e4m*eb3mp),2));

  //Errors GMN
  double e1n = tau/(1+b1_mn*tau+b2_mn*tau*tau + b3_mn*tau*tau*tau);
  double e2n = (1+a1_mn*tau)/pow((1+b1_mn*tau+b2_ep*tau*tau + b3_mn*tau*tau*tau),2)*tau;
  double e3n = (1+a1_mn*tau)/pow((1+b1_mn*tau+b2_ep*tau*tau + b3_mn*tau*tau*tau),2)*tau*tau;
  double e4n = (1+a1_mn*tau)/pow((1+b1_mn*tau+b2_ep*tau*tau + b3_mn*tau*tau*tau),2)*tau*tau*tau;
  double egmn = sqrt(pow((e1n*ea1mn),2)+pow((e2n*eb1mn),2)+pow((e3n*eb2mn),2)+pow((e4n*eb3mn),2));

  //Errors GEN
  double e1t = tau/(1+tau*bb)/pow((1+t/.71),2);
  double e2t = aa*tau*tau/pow((1+tau*bb),2)/pow((1+t/.71),2);
  double egen = sqrt(pow(e1t*eaa,2)+pow(e2t*ebb,2));

  //Dirac form factor proton
  double f1p = (tau*gmp+gep)/(1+tau);

  //Pauli form factor proton
  double f2p = (gmp-gep)/(1+tau);

  //Dirac form factor neutron
  double f1n = (tau*gmn+gen)/(1+tau);

  //Pauli form factor neutron
  double f2n = (gmn-gen)/(1+tau);

  //Dipole form factor shape
  double gd = 1/pow((1+t/.71),2);

  //Error Dirac form factor proton
  double ef1p = sqrt(pow((tau*egmp),2)+ pow(egep,2))/(1+tau);

  //Error Pauli form factor proton
  double ef2p = sqrt(pow(egmp,2)+ pow(egep,2))/(1+tau);

  //Error Dirac form factor neutron
  double ef1n = sqrt(pow((tau*egmn),2)+ pow(egen,2))/(1+tau);

  //Error Pauli form factor neutron
  double ef2n = sqrt(pow(egmn,2)+ pow(egen,2))/(1+tau);

  //array with all values
  double form[] = {f1p,f2p,f1n,f2n,ef1p,ef2p,ef1n,ef2n};

  //Control for which value to return
  if(ret == 0) {
    return f1p;
  }
  else if(ret == 1) {
    return f2p;
  }
  else if(ret == 2) {
    return f1n;
  }
  else if (ret == 3) {
    return f2n;
  }

}

