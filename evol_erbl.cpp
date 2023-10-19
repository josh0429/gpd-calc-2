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
#include "evol_erbl.h"

//for output
using namespace std;

double Herblevol(double X, double zeta, double t, double sum,\
		 double evolgpdzeta, double f1, int ret, int ud) {

  //Notice sum has 1/(1-zeta/2) 
  double Serbl = (1-(zeta/2))*(f1- (sum/(1-(zeta/2))));

  //Minus component is symmetric!
  double aminus = (6/(zeta*zeta*zeta))*(zeta*evolgpdzeta-2*Serbl);
  double Hminus = aminus*X*X - aminus*zeta*X + evolgpdzeta;

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
