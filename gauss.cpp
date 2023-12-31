//************************************************************//
//    One dimensional gauss integral                          //
//------------------------------------------------------------//
//                                                            //
//    4 functions GPT,GWT,DINTER,DGAUS1                       //
//                                                            //
//    GPT: Calculates gauss points                            //
//    GWT: Calculates gauss weights                           //
//    DINTER: Interval between points                         //
//    DGAUS1: Calculates Integral                             //
//                                                            //
//    Input variables:                                        //
//                                                            //
//    NPS: Number of gauss points                             //
//    PE: Array of function to integrate values               //
//                                                            //
//    Author: Brandon Kriesten                                //
//                                                            //
//    contact: btk8bh@virginia.edu                            //
//                                                            //
//************************************************************//

//Necessary Includes
#include <iostream>
#include <cmath>
#include "gauss.h"

//For output
using namespace std;

//Calculates gauss points
double gpt(int val) {

  double P[48];

  P[0]=0.9987710072524261e+00;
  P[1]=0.9935301722663508e+00;
  P[2]=0.9841245837228269e+00;
  P[3]=0.9705915925462473e+00;
  P[4]=0.9529877031604309e+00;
  P[5]=0.9313866907065543e+00;
  P[6]=0.9058791367155697e+00;
  P[7]=0.8765720202742479e+00;
  P[8]=0.8435882616243935e+00;
  P[9]=0.8070662040294426e+00;
  P[10]=0.7671590325157403e+0;
  P[11]=0.7240341309238147e+00;
  P[12]=0.6778723796326639e+00;
  P[13]=0.6288673967765136e+00;
  P[14]=0.5772247260839727e+00;
  P[15]=0.5231609747222330e+00;
  P[16]=0.4669029047509584e+00;
  P[17]=0.4086864819907167e+00;
  P[18]=0.3487558862921607e+00;
  P[19]=0.2873624873554556e+00;
  P[20]=0.2247637903946891e+00;
  P[21]=0.1612223560688917e+00;
  P[22]=0.9700469920946270e-01;
  P[23]=0.3238017096286936e-01;

  int J=24;
  int M=24;

  for(int i = M; i < 48; i++){
    J= J-1;
    P[i] = P[J];
    P[J] = -P[J];
  }

  return P[val];
}

//Calculates gauss weights
double gwt(int val) {

  double W[48];

  int L,JP,M,J,I;

  W[0]=0.3153346052305839e-02;
  W[1]=0.7327553901276262e-02;
  W[2]=0.1147723457923454e-01;
  W[3]=0.1557931572294385e-01;
  W[4]=0.1961616045735553e-01;
  W[5]=0.2357076083932438e-01;
  W[6]=0.2742650970835695e-01;
  W[7]=0.3116722783279809e-01;
  W[8]=0.3477722256477044e-01;
  W[9]=0.3824135106583071e-01;
  W[10]=0.4154508294346475e-01;
  W[11]=0.4467456085669428e-01;
  W[12]=0.4761665849249047e-01;
  W[13]=0.5035903555385447e-01;
  W[14]=0.5289018948519367e-01;
  W[15]=0.5519950369998416e-01;
  W[16]=0.5727729210040322e-01;
  W[17]=0.5911483969839564e-01;
  W[18]=0.6070443916589388e-01;
  W[19]=0.6203942315989266e-01;
  W[20]=0.6311419228625403e-01;
  W[21]=0.6392423858464819e-01;
  W[22]=0.6446616443595008e-01;
  W[23]=0.6473769681268392e-01;

  J=24;
  M=24;

  for(int i = M; i < 48; i++){

    J=J-1;
    W[i]=W[J];
    
  }

  return W[val];
  
}

//Interval needed before dgaus1
double dinter(double a, double b, int nps, int val) {

  double df = 0.5*(b-a);
  double ds = 0.5*(b+a);

  double e0;
  
  e0 = df*gpt(val)+ds;

  return e0;
}

//Calculates the integral
double dgaus1(int nps, double a, double b, double* pe) {
  
  double AB;
  double WE[48];
  double DGAUS1;
  double w;
  
  AB = 0.5*(b-a);
  DGAUS1 = 0.;

  for(int i = 0; i<nps; i++ ) {

    w = gwt(i);
    
    WE[i] = w*AB;
    DGAUS1 += WE[i]*pe[i];
  }
  
  return DGAUS1;
  
}
