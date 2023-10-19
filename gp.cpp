#include <iostream>
#include <cmath>
#include <fstream>
#include "gp.h"

//t is positive

using namespace std;

double gp(double t, int toret) {

  double mp = .9383;
  double mpi = .140;
  
  int ng = 48;

  double tmm[ng];
  double gu_th[ng];
  double gd_th[ng];
  double emme[ng];
  double emmed[ng];
  double enned[ng];
  double enne[ng];

  double gu, gd;
  
  ifstream inFile("gptheory.dat");
  string line, tms, gu_ths, gd_ths;
  
  int i = 0;

  while(inFile >> tms >> gu_ths >> gd_ths) {
    double tmd = stod(tms);
    double gu_thd = stod(gu_ths);
    double gd_thd = stod(gd_ths);
    
    tmm[i] = tmd;
    gu_th[i] = gu_thd;
    gd_th[i] = gd_thd;
    i++;
  }
  inFile.close();
  
  for(int j=0; j<ng; j++){
    double y1 = gu_th[j];
    double y2 = gu_th[j+1];
    double x1 = tmm[j];
    double x2 = tmm[j+1];
    
    double y1d = gd_th[j];
    double y2d = gd_th[j+1];
    double x1d = tmm[j];
    double x2d = tmm[j+1];

    emme[j] = (y1-y2)/(x1-x2);
    enne[j] = y1 - emme[i]*x1;
    emmed[j] = (y1d-y2d)/(x1d-x2d);
    enned[j] = y1d-emmed[i]*x1d;
    if(t>=tmm[j] and t<=tmm[j+1]){
      gu = emme[j]*t + enne[j];
      gd = emmed[j]*t + enned[j];
    }
  }
  double pole = mpi*mpi+t;
  double g0 = 4*mp*mp/pole*(1-1.7*pole/pow((1+t/2),2));
  double gv = g0;
  double gs = 0;

  if(toret == 0){
    return 0.5*(gs+gv);
  }

  else if(toret == 1) {
    return 0.5*(gs-gv);
  }

}
