//**********************************************************//
//                                                          //
//    Diquark Calculation of GPDs                           //
// ---------------------------------------------------------//
//                                                          //
//    Author: Brandon Kriesten                              //
//                                                          //
//    Calculates initial scale GPDs in the Diquark model    //
//    using helicity amplitude formalism. Parameters are    //
//    fitted using form factor elastic scattering data.     //
//                                                          //
//    16 total parameters per GPD, each parameter is        //
//    fitted using Cates data, following reference          //
//    Gonzalez,Goldstein,Kathuria,Liuti.                    //
//                                                          //
//    Reference: 10.1103/PhysRevC.88.065206                 //
//                                                          //
//    Input Variables:                                      //
//                                                          //
//    X = longitudinal variable                             //
//    zeta = x_Bj                                           //
//    t = 4-momentum transfer squared                       //
//        always negative                                   //
//                                                          //
//    Output Variables:                                     //
//                                                          //
//    hu = H_u in DGLAP region                              //
//    hd = H_d in DGLAP region                              //
//    eu = E_u in DGLAP region                              //
//    ed = E_d in DGLAP region                              //
//    hutil = H_u tilde in DGLAP region                     //
//    hdtil = H_d tilde in DGLAP region                     //
//    eutil = E_u tilde in DGLAP region                     //
//    edtil = E_d tilde in DGLAP region                     //
//                                                          //
//    hu_plus = H_u upper curve for error band in DGLAP     //
//              region                                      //
//    hd_plus = H_d upper curve for error band in DGLAP     //
//              region                                      //
//    eu_plus = E_u upper curve for error band in DGLAP     //
//              region                                      //
//    ed_plus = E_d upper curve for error band in DGLAP     //
//              region                                      //
//    hutil_plus = H_u tilde upper curve for error band in  //
//              DGLAP region                                //
//    hdtil_plus = H_d tilde upper curve for error band in  //
//              DGLAP region                                //
//    eutil_plus = E_u tilde upper curve for error band in  //
//              DGLAP region                                //
//    edtil_plus = E_d tilde upper curve for error band in  //
//              DGLAP region                                //
//                                                          //
//    contact: btk8bh@virginia.edu                          //
//                                                          //
//**********************************************************//

//Necessary includes
#include <iostream>
#include <cmath>
#include <fstream>

//Include Header File
#include "sea_gpd.h"

//Necessary for output
using namespace std;

//**********************************************************//
//    DGLAP REGION PARAMETERS                               //
//----------------------------------------------------------//
//                                                          //
//    M = Proton Mass                                       //
//    m = Quark Mass                                        //
//    Mx = Diquark Mass                                     //
//    Ml = Diquark Form Factor Mass Cutoff                  //
//    alpha = Regge t-independent exponent                  //
//    alphap = Regge t-depenedent exponent                  //
//    p = Regge t-dependent exponent                        //
//    beta = 0 (this should be fixed if higher moments      //
//           need better fit                                //
//    a = 1 can be fitted later, modifies t-dependent       //
//           regge parameters                               //
//    N = Normalization constant                            //
//    MMp = Mass factor with zeta dependence                //
//    MM = Mass factor zeta-independent                     //
//    Xp = X' which is zeta dependent                       //
//    t0 = minimum t based on zeta                          //
//    dT = Transverse momentum transfer between initial and //
//         final proton states, fourier transform of        //
//         transverse spatial position of quarks.           //
//                                                          //
//**********************************************************//

//Calculation of GPD H for u-quarks fitted to Dirac form factor data
double gpdHub(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double alphap = 0.1;// - 0.06;
  double p = 0.1;// - 0.025;
  double M = .9383;
  double beta = 0;
  double a = 1;
  double alpha = 1.144;
  double Mx = 3.25;
  double m = 0.38;
  double Ml = 1.372;
  double N = 1.206;// + 0.008;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0.; kT < 5.; kT+=.001) {
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    double omx3 = (1-X)*(1-X)*(1-X);
    double zeta2 = zeta*zeta;
    double omz2 = (1-zeta)*(1-zeta);
    double mas = m+M*X;
    double masp = m+M*Xp;
    double kt2 = kT*kT;
    double reg1 = pow(X,-alpha);
    double regp = -alphap*pow(1-X,p);
    double reg2 = pow(X,(regp*t));
    double reg3 = pow(X,-beta*zeta*zeta/(1-zeta)*t);
    
    Hu += -N*2*pi*(1-zeta/2)*reg1*reg2*reg3*(omx3/omz2)*kT\
      *(((mas*masp+kt2)*D1-2*D2*D2)/(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Hu*.001;
}

//Calculation of GPD H for d-quarks fitted to Dirac form factor data
double gpdHdb(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double alphap = 0.125;// - 0.023;
  double p = 0.12;// -0.050;
  double beta = 0;
  double a = 1;
  double m = 0.3;
  double Mx = 2.105;
  double Ml = 1.495;
  double alpha = 1.125;
  double N = 1.230;// + 0.082;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hd = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0.; kT < 5.; kT+=.001) {
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    double omx3 = (1-X)*(1-X)*(1-X);
    double zeta2 = zeta*zeta;
    double omz2 = (1-zeta)*(1-zeta);
    double mas = m+M*X;
    double masp = m+M*Xp;
    double kt2 = kT*kT;
    double reg1 = pow(X,-alpha);
    double regp = -alphap*pow(1-X,p);
    double reg2 = pow(X,(regp*t));
    double reg3 = pow(X,-beta*zeta*zeta/(1-zeta)*t);

    Hd += -N*2*pi*(1-zeta/2)*reg1*reg2*reg3*(omx3/omz2)*kT\
      *(((mas*masp+kt2)*D1-2*D2*D2)/(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Hd*.001;
}

//Calculation of GPD H for d-quarks fitted to Dirac form factor data                                                          
double gpdHsb(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double alphap = 0.;
  double p = 0.;
  double beta = 0;
  double a = 1;
  double m = 0.5;
  double Mx = 2.1;
  double Ml = .7;
  double alpha = 0.6;
  double N = .125;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hsb = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0.; kT < 5.; kT+=.001) {
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    double omx3 = (1-X)*(1-X)*(1-X);
    double zeta2 = zeta*zeta;
    double omz2 = (1-zeta)*(1-zeta);
    double mas = m+M*X;
    double masp = m+M*Xp;
    double kt2 = kT*kT;
    double reg1 = pow(X,-alpha);
    double regp = -alphap*pow(1-X,p);
    double reg2 = pow(X,(regp*t));
    double reg3 = pow(X,-beta*zeta*zeta/(1-zeta)*t);

    Hsb += -N*2*pi*(1-zeta/2)*reg1*reg2*reg3*(omx3/omz2)*kT\
      *(((mas*masp+kt2)*D1-2*D2*D2)/(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Hsb*.001;
}

double gpdHcb(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double alphap = 0.;
  double p = 0.;
  double beta = 0;
  double a = 1;
  double m = 0.25;
  double Mx = 3.00;
  double Ml = 1.35;
  double alpha = 1.3;
  double N = .01;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hcb = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0.; kT < 5.; kT+=.001) {
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    double omx3 = (1-X)*(1-X)*(1-X);
    double zeta2 = zeta*zeta;
    double omz2 = (1-zeta)*(1-zeta);
    double mas = m+M*X;
    double masp = m+M*Xp;
    double kt2 = kT*kT;
    double reg1 = pow(X,-alpha);
    double regp = -alphap*pow(1-X,p);
    double reg2 = pow(X,(regp*t));
    double reg3 = pow(X,-beta*zeta*zeta/(1-zeta)*t);

    Hcb += -N*2*pi*(1-zeta/2)*reg1*reg2*reg3*(omx3/omz2)*kT\
      *(((mas*masp+kt2)*D1-2*D2*D2)/(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Hcb*.001;
}

double gpdEub(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .38;
  double Mx = 3.25;
  double Ml = 1.372;
  double alpha = 1.144;
  double alphap = 0.1;
  double p = 0.1;
  double beta = 0.;
  double a = 1;
  double N = 1.803;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Eu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0; kT < 5; kT+=.001) {
    double omx3 = (1-X)*(1-X)*(1-X);
    double reg1 = pow(X,-alpha);
    double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
    double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    Eu += N*2*pi*(1-zeta/2)*reg1*reg2*reg3*((omx3)/(1-zeta))    \
      *kT*((2*M*((-2*M)*(X-Xp)*kT*kT-(m+M*X)*(1-Xp)*D1))\
           /(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Eu*.001;
}

double gpdEdb(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .3;
  double Mx = 2.105;
  double Ml = 1.495;
  double alpha = 1.125;
  double alphap = 0.125;
  double p = 0.12;
  double beta = 0;
  double a = 1;
  double N = -2.780;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Ed = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  for(double kT = 0; kT < 5; kT+=.001) {
    double omx3 = (1-X)*(1-X)*(1-X);
    double reg1 = pow(X,-alpha);
    double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
    double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);
    double D0 = (1-X)*MM-kT*kT;
    double D1 = (1-Xp)*MMp-kT*kT-dT*dT*(1-Xp)*(1-Xp);
    double D2 = (1-Xp)*kT*dT;
    Ed += N*2*pi*(1-zeta/2)*reg1*reg2*reg3*((omx3)/(1-zeta))    \
      *kT*((2*M*((-2*M)*(X-Xp)*kT*kT-(m+M*X)*(1-Xp)*D1))\
           /(pow((D1*D1-4*D2*D2),1.5)*D0*D0));
  }
  return Ed*.001;
}
