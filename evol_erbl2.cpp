//**************************************************************//
//    Calculates ERBL Region of GPD                             //
//--------------------------------------------------------------//
//                                                              //
//    Input variables:                                          //
//                                                              //
//    X: Longitudinal variable                                  //
//    zeta: Longitudinal momentum transfer ratio                //
//    t: 4-momentum transfer squared                            //
//    sum: DGLAP region integral                                //
//    evolgpdzeta: Value of GPD at X=zeta                       //
//    f1: Value of form factor at t                             //
//    ret: 0 --> Plus symmetry 1 --> Minus symmetry             //
//    ud: 0 --> u-quark 1 --> d-quark                           //
//                                                              //
//    Author: Brandon Kriesten                                  //
//                                                              //
//    contact: btk8bh@virginia.edu                              //
//                                                              //
//**************************************************************//

//Necessary includes
#include <iostream>
#include <cmath>
#include "formfactor.h"
#include "evol_erbl2.h"

//for output
using namespace std;

double Herblevolhd(double X, double zeta, double t, double sum,\
		   double sum2, double sum3, double evolgpdzeta,
		   double f1, double a20, double a30, double a32,
		   double c20, int ret, int ud) {

  //Notice sum has 1/(1-zeta/2) 

  double omz = 1-zeta/2;
  double omz2 = omz*omz;
  double omz3 = omz2*omz;
  double xi = zeta/(2-zeta);
  double xi2 = xi*xi;
  double z2 = zeta*zeta;
  double z3 = z2*zeta;
  double z4 = z3*zeta;
  double z5 = z4*zeta;
  double z6 = z5*zeta;
  
  double Serbl1 = 2*omz*(f1- (sum/omz));
  double Serbl2 = 2*omz2*(a20+4*xi2*c20 - (sum2/omz2));
  double Serbl3 = 2*omz3*(a30+4*xi2*a32 - (sum3/omz3));

  double G = (1/z2)*(Serbl1 - zeta*evolgpdzeta);
  double H = (1/z3)*(Serbl2 - 0.5*z2*evolgpdzeta);
  double J = (1/z4)*(Serbl3 - (1./3.)*z3*evolgpdzeta);
  
  //Minus component is symmetric!
  double A = ((6*G-12*H)-(600/z2)*(J-1.2*H+0.3*G)*(0.1*z2-X*zeta+X*X))/(zeta*(2-600/350))/(0.1*z2-X*zeta+X*X);
  double B = (600/z2)*(J-1.2*H+0.3*G)-(600/350)*zeta*A;
  double C = (1./zeta)*(36*(H-0.6667*G)-1.2*z2*B-1.2*z3*A);
  double D = -0.4*z3*A-0.5*z2*B-0.6667*zeta*C+2*G;
  
  double Hminus = (A*X*X*X*X + B*X*X*X + C*X*X + D*X + evolgpdzeta);

  //plus versus minus
  double aplus;
  if(ud == 0) {
    aplus = 2000;
  }
  else if(ud == 1) {
    aplus = 1000;
  }
  else if(ud == 2) {
    aplus = 2000;
  }
  
  //Plus component is anti-symmetric!
  double c = (2*evolgpdzeta+.5*aplus*zeta*zeta*zeta)/zeta;
  double d = -evolgpdzeta;
  double Hplus = aplus*X*X*X - 1.5*aplus*zeta*X*X + c*X + d;

  //return plus versus minus component
  if(ret == 0){
    return Hplus;
  }
  else if(ret == 1) {
    return Hminus;
  }

}
