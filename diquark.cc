#include "diquark.h"
#include "gluon_gpd.h"
#include "sea_gpd.h"
#include "gpd.h"
#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

//All C and C++ STL headers should work fine.
//Additional libs will require tweak of Makefile.

void diquark_sub_ggla_(double *x, double *q2, double *zeta, double *t,
		       double *hu, double *hu_plus,
		       double *eu, double *eu_plus,
		       double *hutil, double *hutil_plus,
		       double *eutil, double *eutil_plus,
		       double *hd, double *hd_plus,
		       double *ed, double *ed_plus,
		       double *hdtil, double *hdtil_plus,
		       double *edtil, double *edtil_plus,
		       double *hg, double *eg, double *htg, double *hub,
		       double *hdb, double *hsb, double *hcb,
		       double *init_q02, double *etg, double *eub, double *edb) {

  if(*init_q02 == 0.09362) {
    
    //VALENCE QUARK PARAMETRIZATION
    *hu = gpdHu(*x,*zeta,*t);
    *hd = gpdHd(*x,*zeta,*t);
    *eu = gpdEu(*x,*zeta,*t);
    *ed = gpdEd(*x,*zeta,*t);
    *hutil = gpdHutil(*x,*zeta,*t);
    *hdtil = gpdHdtil(*x,*zeta,*t);
    *eutil = gpdEutil(*x,*zeta,*t);
    *edtil = gpdEdtil(*x,*zeta,*t);
    *hu_plus = gpdHuplus(*x,*zeta,*t);
    *eu_plus = gpdEuplus(*x,*zeta,*t);
    *hutil_plus = gpdHutilplus(*x,*zeta,*t);
    *eutil_plus = gpdEutilplus(*x,*zeta,*t);
    *hd_plus = gpdHdplus(*x,*zeta,*t);
    *ed_plus = gpdEdplus(*x,*zeta,*t);
    *hdtil_plus = gpdHdtilplus(*x,*zeta,*t);
    *edtil_plus =gpdEdtilplus(*x,*zeta,*t);
    *hg = 0.;
    *eg = 0.;
    *htg = 0.;
    *etg = 0.;
    *hub = 0.;
    *hdb = 0.;
    *hsb = 0.;
    *hcb = 0.;
    *eub = 0.;
    *edb = 0.;
  }
  
  else {

    ifstream inFile("gpd_inter.dat");
    string xs,uvs,dvs,ubs,dbs;
    while(inFile >> xs >> uvs >> dvs >> ubs >> dbs){
      double xd = stod(xs);
      double xtest = *x;
      if(abs(xd-xtest) < 0.000003) {
	*hu = stod(uvs);
	*hu_plus = stod(uvs);
	*eu = stod(uvs);
	*eu_plus = stod(uvs);
	*hutil = stod(uvs);
	*hutil_plus = stod(uvs);
	*eutil = stod(uvs);
	*eutil_plus = stod(uvs);
	*hd = stod(dvs);
	*hd_plus = stod(dvs);
	*ed = stod(dvs);
	*ed_plus = stod(dvs);
	*hdtil = stod(dvs);
	*hdtil_plus = stod(dvs);
	*edtil = stod(dvs);
	*edtil_plus = stod(dvs);
      }
    }
    inFile.close();
    
    // SEA QUARK PARAMETRIZATION
    *hub = gpdHub(*x,*zeta,*t);
    *hdb = gpdHdb(*x,*zeta,*t);
    *hsb = gpdHsb(*x,*zeta,*t);
    *hcb = gpdHcb(*x,*zeta,*t);
    *eub = gpdEub(*x,*zeta,*t);
    *edb = gpdEdb(*x,*zeta,*t);
    
    // GLUON PARAMETRIZATION
    *hg = gpdHg(*x,*zeta,*t); 
    *eg = gpdEg(*x,*zeta,*t);
    *htg = gpdHgtil(*x,*zeta,*t);
    *etg = gpdEgtil(*x,*zeta,*t);
  }
}
