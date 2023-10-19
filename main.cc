#include <iostream>
#include "cff_calc.h"
#include <fstream>

using namespace std;

int main() {

  double M = 0.93828;
  
//  double Q2 = 20.00; // Default 1.92
  double xbj = 0.2; // Default 0.36. Definition of IPPDF is xi = 0 -> zeta = 0
  double zeta = xbj; //zeta is different from xi. See page 12 of Diehl's physics report for why zeta = xbj
  double tmin = -zeta*zeta*M*M/(1-zeta);
//  double t = -1.1; // Default -0.33

  cout << tmin << endl;
  
  double init_q2 = 0.09362;
  double rehpv = 0.;
  double imhpv = 0.;
  double rehuv = 0.;
  double imhuv = 0.;
  double rehdv = 0.;
  double imhdv = 0.;
  double rehp = 0.;
  double imhp = 0.;
  double rehu = 0.;
  double imhu = 0.;
  double rehd = 0.;
  double imhd = 0.;
  double rehg = 0.;
  double imhg = 0.;
  
//  for(int jgpd = 1; jgpd<=4; jgpd++){
//      cff_calc(xbj,t,Q2,jgpd,rehpv,imhpv,rehuv,imhuv,rehdv,imhdv,	\
//	     rehp,imhp,rehu,imhu,rehd,imhd,rehg,imhg);  
//  }

  for(double Q2 = 4.00; Q2<5.00; Q2+=5.00){
    for(double t = 0.3000; t<0.400; t+=1.0000){
      for(int jgpd = 1; jgpd<=4; jgpd++){
        cff_calc(xbj,-t,Q2,jgpd,rehpv,imhpv,rehuv,imhuv,rehdv,imhdv,	\
	      rehp,imhp,rehu,imhu,rehd,imhd,rehg,imhg);  
      }
    }
  }  

//  for(double t = 0.000; t<2.010; t+=0.050){
//    for(int jgpd = 1; jgpd<=4; jgpd++){
//      cff_calc(xbj,-t,Q2,jgpd,rehpv,imhpv,rehuv,imhuv,rehdv,imhdv,	\
//	     rehp,imhp,rehu,imhu,rehd,imhd,rehg,imhg);  
//    }
//  }
  
  return 0;

}
