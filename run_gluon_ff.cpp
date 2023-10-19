#include <iostream>
#include "gluon_ff.h"

using namespace std;

int main() {

  for(double t = 0; t<= 2; t+=0.01){

    cout << t << " " << gff_B(t) << endl;

  }

}
