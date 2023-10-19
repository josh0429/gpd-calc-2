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
#include "gpd.h"

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

  double H1 = (m + M*X) * (m + M*Xp) * (MoXp + d*d)*I0;
  double H2 = (MoXp + (m + M*X)*(m + M*Xp) - d*d)*I1;
  double H = 2*pi*N*(1. - X)*pow((1. - Xp), 2)*(1 - zeta / 2) * (H1 + H2 + I2);
  double rexponent = alpha + alphap * pow((1. - Xp), p) * t;
  double regge = pow(X, -rexponent);

  return regge * H + ((zeta*zeta)/(4*(1. - zeta))) * gpdEu(X, zeta, t);
}

double gpdHuplus(double X,double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .603;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 2.448;
  double p = 0.620;
  double beta = 0;
  double a = 1;
  double N = 2.043;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double regp = -alphap*pow(1-X,p);
  double reg2 = pow(X,(regp*t));
  double reg = reg1*reg2;

  double alpha_err = 0.07*alpha;
  double N_err = 0.1*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = .0220;
  double dp = .0170;
  
  double hu_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
		       +pow(reg_err*alpha_err,2)+pow(N_err,2));
  return gpdHu(X,zeta,t) + gpdHu(X,zeta,t)*hu_err/(reg);
  
}

double gpdHuminus(double X,double zeta,double t){

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .603;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 1.814;
  double p = .449;
  //double alphap = 2.448;
  //double p = 0.620;
  double beta = 0;
  double a = 1;
  double N = 2.043;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double regp = -alphap*pow(1-X,p);
  double reg2 = pow(X,(regp*t));
  double reg = reg1*reg2;

  double alpha_err = 0.07*alpha;
  double N_err = 0*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = .0220;
  double dp = .0170;

  double hu_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
		       +pow(reg_err*alpha_err,2)+pow(N_err,2));
  
  return gpdHu(X,zeta,t) - hu_err*gpdHu(X,zeta,t)/(reg);
  
}

//Calculation of GPD H for d-quarks fitted to Dirac form factor data
double gpdHd(double X, double zeta, double t) {
  //***************************//
  //      PARAMETERS
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = 0.860;
  double alpha = .0317;
  double alphap = 2.209;
  double p = 0.658;
  double N = 1.570;
  double beta = 0;
  double a = 1;
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

  double H1 = (m + M*X) * (m + M*Xp) * (MoXp + d*d)*I0;
  double H2 = (MoXp + (m + M*X)*(m + M*Xp) - d*d)*I1;
  double H = 2*pi*N*(1. - X)*pow((1. - Xp), 2)*(1 - zeta / 2) * (H1 + H2 + I2);
  double rexponent = alpha + alphap * pow((1. - Xp), p) * t;
  double regge = pow(X, -rexponent);
  return regge * H + ((zeta*zeta)/(4*(1. - zeta))) * gpdEd(X, zeta, t);
}

double gpdHdplus(double X, double zeta, double t) {

  
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = 0.860;
  double alpha = .0317;
  //double alphap = 1.139;
  //double p = -0.113;
  double alphap = 2.209;
  double p = 0.658;
  double beta = 0;
  double a = 1;
  double N = 1.570;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hd = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));
  
  double reg1 = pow(X,-alpha);
  double regp = -alphap*pow(1-X,p);
  double reg2 = pow(X,(regp*t));
  double reg3 = pow(X,-beta*zeta*zeta/(1-zeta)*t);
  double reg = reg1*reg2*reg3;

  double alpha_err = 0.07*alpha;
  double N_err = 0.1*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = 0.0564;
  double dp = 0.1040;

  double hd_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHd(X,zeta,t) + gpdHd(X,zeta,t)*hd_err/(reg);
    
}

double gpdHdminus(double X, double zeta, double t) {


  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = 0.860;
  double alpha = .0317;
  double alphap = 1.139;
  double p = -0.113;
  //double alphap = 2.209;
  //double p = 0.658;
  double beta = 0;
  double a = 1;
  double N = 1.570;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Hd = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double regp = -alphap*pow(1-X,p);
  double reg2 = pow(X,(regp*t));
  double reg3 = pow(X,-beta*zeta*zeta/(1-zeta)*t);
  double reg = reg1*reg2*reg3;

  double alpha_err = 0.07*alpha;
  double N_err = 0*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = 0.0564 ;
  double dp = 0.1040 ;

  double hd_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHd(X,zeta,t) - gpdHd(X,zeta,t)*hd_err/(reg);

}

//Calculation of GPD E for u-quarks fitted to Pauli form factor data
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

  double rexponent = alpha + alphap * pow((1. - Xp), p) * t;
  double regge = pow(X, -rexponent);
  double E1 = 2*M*(m + M*X)*(MoXp + d*d)*I0;
  double E2 = (2*M*(m + M*X) - 4*M*M*(X - Xp))*I1;
  double E = 2*pi*N*(1. - zeta/2)*pow((1. - zeta),2)*pow((1. - Xp), 4)*(E1 + E2);
  
  return regge*E;
}

double gpdEuplus(double X, double zeta, double t) {

  
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .604;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 2.835;
  double p = .969;
  double beta = 0;
  double a = 1;
  double N = 1.803;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Eu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);

  double reg = reg1*reg2*reg3;

  double alpha_err = 0.07*alpha;
  double N_err = 0.11*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);
  double dap = 0.0509;
  double dp = 0.0307;
  
  double eu_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEu(X,zeta,t) + gpdEu(X,zeta,t)*eu_err/reg;
  
}

double gpdEuminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .420;
  double Mx = .604;
  double Ml = 1.018;
  double alpha = .210;
  double alphap = 2.835;
  double p = .969;
  double beta = 0;
  double a = 1;
  double N = 1.803;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Eu = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);

  double reg = reg1*reg2*reg3;

  double alpha_err = 0.07*alpha;
  double N_err = 0.11*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double dap = 0.0509;
  double dp = 0.0307;
  
  double eu_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEu(X,zeta,t) - gpdEu(X,zeta,t)*eu_err/reg;

}


//Calculation of GPD E for d-quarks fitted to Pauli form factor data
double gpdEd(double X, double zeta, double t) {
  //***************************//
  //      PARAMETERS
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = .860;
  double alpha = .0317;
  double alphap = 1.281;
  double p = 0.726;
  double beta = 0;
  double a = 1;
  double N = -2.780;
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

  double rexponent = alpha + alphap * pow((1. - Xp), p) * t;
  double regge = pow(X, -rexponent);
  double E1 = 2*M*(m + M*X)*(MoXp + d*d)*I0;
  double E2 = (2*M*(m + M*X) - 4*M*M*(X - Xp))*I1;
  double E = 2*pi*N*(1. - zeta/2)*pow((1. - zeta),2)*pow((1. - Xp), 4)*(E1 + E2);
  
  return regge*E;
}

double gpdEdplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = .860;
  double alpha = .0317;
  double alphap = 1.281;
  double p = 0.726;
  double beta = 0;
  double a = 1;
  double N = -2.780;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Ed = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);

  double reg = reg1*reg2*reg3;

  double alpha_err = 0.1*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double dap = 0.0310;
  double dp = 0.0631;

  double ed_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEd(X,zeta,t) + gpdEd(X,zeta,t)*ed_err/reg;
  
}


double gpdEdminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = .275;
  double Mx = .913;
  double Ml = .860;
  double alpha = .0317;
  double alphap = 1.281;
  double p = 0.726;
  double beta = 0;
  double a = 1;
  double N = -2.780;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M-(X/(1-X))*Mx*Mx-Ml*Ml;
  double Ed = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg3 = pow(X,-beta*(pow(zeta,2)/(a-zeta))*t);

  double reg = reg1*reg2*reg3;

  double alpha_err = 0.1*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double dap = 0.0310;
  double dp = 0.0631;

  double ed_err = sqrt(pow(alphap_err*dap,2)+pow(p_err*dp,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEd(X,zeta,t) - gpdEd(X,zeta,t)*ed_err/reg;

}

//Calculation of GPD H_t for u-quarks fitted to axial form factor data
double gpdHutil(double X, double zeta, double t) {
  //***************************//
  //      PARAMETERS
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 1.543;
  double p = .346;
  double N = .0504;
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

  double Ht1 = (m + M*X) * (m + M*Xp) * (MoXp + d*d)*I0;
  double Ht2 = (-MoXp + (m + M*X)*(m + M*Xp) + d*d)*I1;
  double Ht = 2*pi*N*(1. - X)*pow((1. - Xp), 2)*(1 - zeta / 2) * (Ht1 + Ht2 - I2);
  double rexponent = alpha + alphap * pow((1. - Xp), p) * t;
  double regge = pow(X, -rexponent);
  
  return regge * Ht + ((zeta*zeta)/(4*(1. - zeta))) * gpdEutil(X, zeta, t);
}

double gpdHutilplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 1.543;
  double p = .346;
  double N = .0504;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.16*alpha;
  double N_err = 0.16*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double hutil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHutil(X,zeta,t) + gpdHutil(X,zeta,t)*hutil_err/reg;

}

double gpdHutilminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 1.543;
  double p = .346;
  double N = .0504;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.16*alpha;
  double N_err = 0.16*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double hutil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHutil(X,zeta,t) - gpdHutil(X,zeta,t)*hutil_err/reg;

}


//Calculation of GPD H for d-quarks fitted to axial form factor data
double gpdHdtil(double X, double zeta, double t) {
  //***************************//
  //      PARAMETERS
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 1.298;
  double p = .974;
  double N = -.0262;
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

  double Ht1 = (m + M*X) * (m + M*Xp) * (MoXp + d*d)*I0;
  double Ht2 = (-MoXp + (m + M*X)*(m + M*Xp) + d*d)*I1;
  double Ht = 2*pi*N*(1. - X)*pow((1. - Xp), 2)*(1 - zeta / 2) * (Ht1 + Ht2 - I2);
  double rexponent = alpha + alphap * pow((1. - Xp), p) * t;
  double regge = pow(X, -rexponent);
  
  return regge * Ht + ((zeta*zeta)/(4*(1. - zeta))) * gpdEdtil(X, zeta, t);
}

double gpdHdtilplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 1.298;
  double p = .974;
  double N = -.0262;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hdtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.2*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double hdtil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHdtil(X,zeta,t) + gpdHdtil(X,zeta,t)*hdtil_err/reg;
  
}

double gpdHdtilminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 1.298;
  double p = .974;
  double N = -.0262;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Hdtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.2*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double hdtil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdHdtil(X,zeta,t) - gpdHdtil(X,zeta,t)*hdtil_err/reg;

}

//Calculation of GPD E_t for u-quarks fitted to pseudo-scalar form factor data
double gpdEutil(double X, double zeta, double t) {
  //***************************//
  //      PARAMETERS
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 5.130;
  double p = 3.507;
  double N = 1.074;
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

  double rexponent = alpha + alphap * pow((1. - Xp), p) * t;
  double regge = pow(X, -rexponent);
  double E1 = 4*M*(m + M*X)*(MoXp + d*d)*I0;
  double E2 = 4*M*((m + M*X) + 2*(m + M*Xp))*I1;
  double E = 2*pi*N*(1. - zeta/2)*pow((1. - zeta),2)*pow((1. - Xp), 4)*(E1 - E2) / zeta;
  
  return regge*E;
}

double gpdEutilplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 5.130;
  double p = 3.507;
  double N = 1.074;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Eutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.16*alpha;
  double N_err = 0.16*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double eutil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEutil(X,zeta,t) + gpdEutil(X,zeta,t)*eutil_err/reg;
  
}

double gpdEutilminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.624;
  double Mx = .474;
  double Ml = .971;
  double alpha = .219;
  double alphap = 5.130;
  double p = 3.507;
  double N = 1.074;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Eutil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.16*alpha;
  double N_err = 0.16*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double eutil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEutil(X,zeta,t) - gpdEutil(X,zeta,t)*eutil_err/reg;

}



//Calculation of GPD E_t for d-quarks fitted to pseudo-scalar form factor data
double gpdEdtil(double X, double zeta, double t) {
  //***************************//
  //      PARAMETERS
  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 3.385;
  double p = 2.326;
  double N = -0.966;
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

  double rexponent = alpha + alphap * pow((1. - Xp), p) * t;
  double regge = pow(X, -rexponent);
  double E1 = 4*M*(m + M*X)*(MoXp + d*d)*I0;
  double E2 = 4*M*((m + M*X) + 2*(m + M*Xp))*I1;
  double E = 2*pi*N*(1. - zeta/2)*pow((1. - zeta),2)*pow((1. - Xp), 4)*(E1 - E2) / zeta;
  
  return regge*E;
}

double gpdEdtilplus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 3.385;
  double p = 2.326;
  double N = -0.966;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Edtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.2*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double edtil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
                       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEdtil(X,zeta,t) + gpdEdtil(X,zeta,t)*edtil_err/reg;
  
}

double gpdEdtilminus(double X, double zeta, double t) {

  const double pi = 3.14159265358979323846;
  double M = .9383;
  double m = 2.603;
  double Mx = .704;
  double Ml = .878;
  double alpha = .0348;
  double alphap = 3.385;
  double p = 2.326;
  double N = -0.966;
  double MMp = (X-zeta)/(1-zeta)*M*M - ((X-zeta)/(1-X))*Mx*Mx-Ml*Ml;
  double MM = X*M*M - (X/(1-X))*Mx*Mx-Ml*Ml;
  double Edtil = 0;
  double Xp = (X-zeta)/(1-zeta);
  double t0 = -zeta*zeta*M*M/(1-zeta);
  double dT = sqrt((t0-t)*(1-zeta));

  double reg1 = pow(X,-alpha);
  double reg2 = pow(X,-(alphap*pow((1-X),p)*t));
  double reg = reg1*reg2;

  double alpha_err = 0.2*alpha;
  double N_err = 0.15*N;
  double alphap_err = reg*log(X)*pow((1-X),p)*t;
  double p_err = reg*log(X)*log(1-X)*alphap*pow((1-X),p)*t;
  double reg_err = reg*log(X);

  double edtil_err = sqrt(pow(alphap_err*alphap,2)+pow(p_err*p,2)\
		       +pow(reg_err*alpha_err,2)+pow(N_err,2));

  return gpdEdtil(X,zeta,t) - gpdEdtil(X,zeta,t)*edtil_err/reg;

}
