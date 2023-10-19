#include <iostream>
#include <cmath>
#include "gluon_ff.h"


double gff_A(double t){

  double alpha_h = 0.58;
  double hmg = 1.13;

  double g_ffa = alpha_h/pow(1+t/(hmg*hmg),2);
  
  return g_ffa;

}

double gff_B(double t) {

  /*  double a0 = -0.42377;
  double a1 = 2.07733;
  double a2 = -2.05108;
 
  double mpi = 0.135;
  double tcut = 4*mpi*mpi;

  double z = (sqrt(tcut+t)-sqrt(tcut))\
    /(sqrt(tcut+t)+sqrt(tcut));

  double gff_B = a0 + a1*z + a2*z*z;
  */
  double alpha_e = 0.0978;
  double emg = -2.5579;

  double gff_B = alpha_e/pow(1+t/(emg*emg),2);

  return gff_B;
  
}

double gff_D(double t) {

  double alpha_d = -10;
  double dmg = 0.48;

  double g_ffd = alpha_d/pow(1+t/(dmg*dmg),2);

  return g_ffd;
  
}
