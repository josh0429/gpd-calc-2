//***************************************************************//
//    Linear Interpolation between points in xbj                 //
//---------------------------------------------------------------//
//                                                               //
//    Input variables:                                           //
//                                                               //
//    xval: Full array of xbj                                    //
//    erbl: Array of cross section points in xbj                 //
//    x: The xbj value you want to calculate at                  //
//                                                               //
//    Output:                                                    //
//                                                               //
//    ret: The value of the linear interpolation at x            //
//                                                               //
//    Author: Brandon Kriesten                                   //
//                                                               //
//    contact: btk8bh@virginia.edu                               //
//                                                               //
//***************************************************************//

//Necessary Includes
#include <iostream>
#include <cmath>
#include "fut.h"

//For output
using namespace std;

double fut(double* xval, double* erbl, double x) {

  //array of slopes and intercepts
  double emme[490];
  double enne[490];

  //loops through xbj values that is given as input
  for(int i=0;i<489;i++) {
    double y1 = erbl[i];
    double y2 = erbl[i+1];
    double x1 = xval[i];
    double x2 = xval[i+1];
    
    //calculates slope
    emme[i] = (y2-y1)/(x2-x1);

    //calculates intercept
    enne[i] = y1 - emme[i]*x1;
  }

  double ret;

  //goes through the arrays and calculates points
  for(int i=0;i<489;i++) {

    //Chooses xbj values around x
    if(x >= xval[i] and x <= xval[i+1]) {

      //linear interpolate between xbj values to give value at x
      ret = emme[i]*x + enne[i];
    }    
  }

  //return value
  return ret;
}
