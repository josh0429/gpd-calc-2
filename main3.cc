//***********************************************************//
//                                                           //
//    DIQUARK CALCULATION -- Compton Form Factors            //
//                                                           //
//    AUTHOR: Brandon Kriesten                               //
//                                                           //
//    REFERENCE: PHYSICAL REVIEW D 84, 034007 (2011)         //
//                                                           //
//    Diquark calculation parametrized by real               //
//    experimental data. Chiral-even GPDs non-singlet        //
//    sector. GPDs are evolved to physical values of         //
//    Q^2.                                                   //
//                                                           //
//    Main file: Calculates Re and Im parts of the CFF       //
//                                                           //
//    Files used: gpd.cpp, ga.cpp, gp.cpp, gluon_gpd.cpp     //
//                evol_erbl.cpp, diquark.cc, formfactor.cpp  //
//                gauss.cpp, fut.cpp                         //
//                                                           //
//    Output files: GPD_H.dat, GPD_E.dat, GPD_Ht.dat         //
//                  GPD_Et.dat, CFF.dat                      //
//                                                           //
//    contact: btk8bh@virginia.edu                           //
//                                                           //
//***********************************************************//


//C++ Libraries
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <sstream>

//**********************************//
//  Necessary Header file includes  //
//**********************************//

//Calculates the GPDs at the initial scale
#include "gpd.h"

//Calculates the axial form factor from data
#include "ga.h"

//Calculates the pseudoscalar form factor from data
#include "gp.h"

//Calculates the gluon GPDs
#include "gluon_gpd.h"

//Calculates the valence ERBL region
#include "evol_erbl.h"
#include "evol_erbl2.h"

//Calculates the gluon/sea ERBL region
#include "evol_erblg.h"

//Wrapper for diquark code to include in evolution
#include "diquark.h"

//Calculates the vector form factors
#include "formfactor.h"

//Gauss integration algorithm
#include "gauss.h"

//Interpolates between xbj values
#include "fut.h"

//Gluon Form Factor
#include "gluon_ff.h"

//ubar and dbar Form Factors
#include "sea_ff.h"

//For output
using namespace std;

//C++ wrapper for Fortran code
extern "C" {
  void struct_ev_gpd_tmd_(double *qsquare, double* zeta, double *t, int *jgpd,\
  int *ierr, int *itmd, int *jt, double *uv0, double *dv0, double *g0,\
  double *ub0, double *db0, double *sb0, double *cb0, double *init_q02);
};

extern "C" {
  double alpha_my_(double *qsquare);
};

//start the madness
int main() {

  //Creates 5 files
  ofstream myfile,myfile2,myfile3,myfile4,myfile5;
  ofstream myfile6,myfile7,myfile8;
  
  //4 of the files hold the calculated GPDs
  myfile.open("GPD_H.dat");
  myfile2.open("GPD_E_q5.dat");
  myfile3.open("GPD_Ht.dat");
  myfile4.open("GPD_Et.dat");

  //Last file contains the flavor separated and proton CFFs
  myfile5.open("CFF.dat");
  
  //Define pi
  const double pi = 3.14159265358979323846;

  //Create all arrays necessary
  
  //xbj values 
  double xval[490];

  //U-valence evolved
  double uv0[490];

  //D-valence evolved
  double dv0[490];

  //Gluons evolved
  double g0[490];

  //U-bar antiquarks
  double ub0[490];

  //D-bar antiquarks
  double db0[490];

  //Strange quarks
  double sb0[490];

  //Charm quarks
  double cb0[490];
  
  //Initializes xbj calculation values
  double xstart;
  double xv;
  double yyy;

  //Proton Mass
  double M = 0.93828;

  //Relic, not needed but too lazy to take it out
  int increment = 1;

  //index used to find crossover point in the 490 x values
  int index =0;

  //Define the initializing kinematics

  //Final Q2 to evolve to
  double final_scale = 5.;
  
  //Returns alpha_s / pi
  double alpha = alpha_my_(&final_scale);
  double mz = 91.1867;
  double mz2 = mz*mz;
  double alpha_s_mz2 = alpha_my_(&mz2)*pi;
  //  cout << "Alpha_s at Mz^2: " << alpha_s_mz2 << endl;
  //zeta = X_bj
  double zeta = 0.0;
  double xi = zeta/(2.-zeta);

  //4-momentum transfer squared
  //  double t = -zeta*zeta*M*M/(1-zeta) - 0.0001;
  double t = 0.;
    
  //Calculate the minimum t value
  double tmin = -zeta*zeta*M*M/(1-zeta);

  //Some formatting output statements
  cout << "Kinematics" << endl;
  cout << endl;
  cout << "Q^2   = " << final_scale << endl;
  cout << "zeta  = " << zeta << endl;
  cout << "t     = " << t << endl;
  cout << "t_min = " << tmin << endl;
  cout << "xi    = " << xi << endl;
  cout << endl;
  cout << "QCD Parameters" << endl;
  cout << endl;
  cout << "alphas/4pi = " << alpha/4 << endl;
  cout << endl;
  //Calculates all of the x values we will use
  //Notice finer spacing for low-x values
  for(int i=0;i<490;i++) {
    if(i <= 290) {
      yyy = log(1.e+4)*(330.-(i)+1.)/330.;
      xv = exp(-yyy);
    }
    else {
      xstart = exp(-log(1.e+4)*41./330.);
      xv = xstart + ((i)-290.)*(1.-xstart)/201.;
    }
    xval[i] = xv;
    if(xval[i] <= zeta) {
      index += 1;
    }
  }

  //Calculates all of the elastic form factors

  //Proton and Neutron Pauli and Dirac Form Factors
  double f1p = formfactor(-t,0);
  double f2p = formfactor(-t,1);
  double f1n = formfactor(-t,2);
  double f2n = formfactor(-t,3);

  //Flavor Separation for Pauli and Dirac Form Factors
  double f1u = 2*f1p + f1n;
  double f1d = 2*f1n + f1p;
  double f2u = 2*f2p + f2n;
  double f2d = 2*f2n + f2p;

  //Axial Form Factor Vector and Scalar Contributions
  double gav = ga(-t,0);
  double gas = ga(-t,1);

  //Flavor Separation for Axial Form Factor
  //Isovector = u-d
  //Isoscalar = u+d
  double gau = .5*(gas+gav);
  double gad = .5*(gas-gav);

  //Pseudo-scalar form factor
  //Needs pion-pole added later
  double gpu = gp(-t,0);
  double gpd = gp(-t,1);

  //Gluon Gravitational Form Factors
  double gluon_a = gff_A(-t);
  double gluon_b = gff_B(-t);
  double gluon_d = gff_D(-t);

  double gluon_F1 = gluon_a + 4*xi*xi*gluon_d;
  double gluon_F2 = gluon_b - 4*xi*xi*gluon_d;

  //Sea Quark Form Factors
  double ubar_F1 = a20_ub_dip(-t,1) + 4*xi*xi*0.5*(c20_upd_zexp(-t) + c20_umd_zexp(-t));
  double dbar_F1 = a20_db_dip(-t,1) + 4*xi*xi*0.5*(c20_upd_zexp(-t) - c20_umd_zexp(-t));

  //Higher order Mellin Moments from the Lattice
  //Will fit to ERBL region soon.
  double a20u = a20_u_dip(-t,1);
  double a20d = a20_d_dip(-t,1);
  double a30u = 0.5*(a30_upd_zexp(-t) + a30_umd_zexp(-t));
  double a30d = 0.5*(a30_upd_zexp(-t) - a30_umd_zexp(-t));
  double c20u = 0.5*(c20_upd_zexp(-t) + c20_umd_zexp(-t));
  double c20d =	0.5*(c20_upd_zexp(-t) - c20_umd_zexp(-t));
  double a32u = 0.0302;
  

  //Sea Quark Form Factor Errors will use later
  double ubar_F1_err = a20_ub_dip(-t,0) - ubar_F1;
  double dbar_F1_err = a20_db_dip(-t,0)	- dbar_F1;

  //Interval for rectangular integration
  double inter;

  //turns error bars on/off (on == 0, off == 1)
  int ierr = 1;

  //not calculating the TMDs
  int itmd = 1;

  //set to 1
  int jt = 1;

  //*********************
  //
  //   JGPD
  //
  //   1 --> H
  //   2 --> E
  //   3 --> Htilde
  //   4 --> Etilde
  //
  //*********************
  int jgpd;

  //Loops over the GPDs to calculate CFFs
  for(jgpd =1; jgpd <=1; jgpd++) { 

    ofstream tempfile;
    tempfile.open("gpd_inter.dat");
    
    //Defines some initial sum variables to calculate the
    //DGLAP area flavor separated
    double sumu = 0;
    double sumd = 0;
    double sumg = 0;
    double sumub = 0;
    double sumdb = 0;
    double sumu2 = 0;
    double sumd2 = 0;
    double sumu3 = 0;
    
    //Initialize Jacobian for low-x and moderate-x
    double xjacobian = 0;
    double xjacobian2 = 0;

    double temp[490];
    double tempd[490];
    double tempg[490];
    double tempub[490];
    double tempdb[490];

    //Calls the evolution code for valence distribution evolves up to scale of gluons
    double init_q02 = 0.09362;
    double qsquare = 0.97;
    cout << "Evolving valence distributions from " << init_q02 << " to " << qsquare << endl;
    double uv0_t[490];
    double dv0_t[490];
    double g0_t[490];
    double ub0_t[490];
    double db0_t[490];
    double sb0_t[490];
    double cb0_t[490];
    struct_ev_gpd_tmd_(&qsquare, &zeta, &t, &jgpd, &ierr, &itmd, &jt, uv0_t, \
		       dv0_t, g0_t, ub0_t, db0_t, sb0_t, cb0_t, &init_q02);

    //Add in section of code that fills file with values of uv0,dv0
    for( int i = 0; i< 490; i+= 1) {
      tempfile << xval[i] << " " << uv0_t[i]/xval[i] << " " << dv0_t[i]/xval[i] << endl;
    }
    tempfile.close();

    init_q02 = 0.97;
    qsquare = final_scale;
    cout << "Evolving valence/glue/sea distributions from " << init_q02 << " to " << qsquare << endl;
    //Calls the evolution code evolves from initial scale of gluons to final Q2
    struct_ev_gpd_tmd_(&qsquare, &zeta, &t, &jgpd, &ierr, &itmd, &jt, uv0,\
                       dv0, g0, ub0, db0, sb0, cb0, &init_q02);
    
    //Calculates integral of evolved gpd in different kinematic regions
    for(int i = 0; i<490; i+=increment) {
      if(xval[i] <= zeta) {
	temp[i] = uv0[i]/xval[i];
	tempd[i] = dv0[i]/xval[i];
	tempg[i] = g0[i];
	tempub[i] = ub0[i];
	tempdb[i] = db0[i];
      }
      else {
	temp[i] = uv0[i]/xval[i];
	tempd[i] = dv0[i]/xval[i];
	tempg[i] = g0[i];
	tempub[i] = ub0[i];
	tempdb[i] = db0[i];
	if(i<=290) {
	  yyy = log(1.e+4)*(330.-(i)+1.)/330.;
	  xv = exp(-yyy);
	  xjacobian = xv*log(1e4)/330;
	  inter = xval[i] - xval[i-increment];	
	  sumu += (uv0[i]/xval[i])*inter;
	  sumd += (dv0[i]/xval[i])*inter;
	  sumg += g0[i]*inter;
	  sumub += ub0[i]*inter;
	  sumdb += db0[i]*inter;
	  sumu2 += uv0[i]*inter;
	  sumd2 += dv0[i]*inter;
	  sumu3 += uv0[i]*inter*xval[i];
	}
	else {
	  xstart = exp(-log(1.e+4)*41./330.);
	  xjacobian2 = (1-xstart)/201.;
	  inter = xval[i] - xval[i-increment];
	  sumu += (uv0[i]/xval[i])*inter;
	  sumd += (dv0[i]/xval[i])*inter;
	  sumg += g0[i]*inter;
	  sumub += ub0[i]*inter;
	  sumdb += db0[i]*inter;
	  sumu2 += uv0[i]*inter;
	  sumd2 += dv0[i]*inter;
	  sumu3 += uv0[i]*inter*xval[i];
	}
      }
    }

    //Starts the integral which we will check at the end
    double erbl_sum = 0.;
    double erbl_sumd = 0.;
    double erbl_sumg = 0.;
    double erbl_sumub = 0.;
    double erbl_sumdb = 0.;

    //division by (1-zeta/2) for jacobian
    erbl_sum += sumu/(1.-zeta/2.);
    erbl_sumd += sumd/(1.-zeta/2.);

    //division by (1-zeta/2)**2 because 2nd mellin moment
    erbl_sumg += sumg/(1.-zeta/2.)/(1.-zeta/2.);

    //Multiply by 2 for symmetry about x = 0
    //Form factor defined for integral dx from -1 to 1
    erbl_sumub += 2*sumub/(1.-zeta/2.)/(1.-zeta/2.);
    erbl_sumdb += 2*sumdb/(1.-zeta/2.)/(1.-zeta/2.);

    //Output to check what we should expect for DGLAP and ERBL
    //region integrals and what total should look like
    if(jgpd == 1) {
      cout << endl;
      cout << "ERBL Uv Area: " << (f1u - (erbl_sum)) << "    ERBL Dv Area: "\
	   << (f1d - (erbl_sumd)) << endl;
      cout << "DGLAP Uv Area: " << sumu/(1-zeta/2) << "    DGLAP Dv Area: "\
	   << sumd/(1-zeta/2) << endl;
      cout << "Total Uv Area: " << f1u << "    Total Dv Area: "\
	   << f1d << endl;
      cout << endl;
      cout << endl;
      cout << "ERBL G Area: " << (gluon_F1 - erbl_sumg)		\
	   << "    ERBL Ub Area: " << ubar_F1 - erbl_sumub\
	   << "    ERBL Db Area: " << dbar_F1 - erbl_sumdb << endl;
      cout << "DGLAP G Area: " << erbl_sumg << "    DGLAP Ub Area: " \
	   << erbl_sumub << "    DGLAP Db Area: " << erbl_sumdb << endl;
      cout << "Total G Area: " << gluon_F1 \
	   << "    Total Ub Area: " << ubar_F1 << "    Total Dbar Area: "\
	   << dbar_F1 << endl;
      cout << endl;
    }
    else if(jgpd == 2) {
      cout << endl;
      cout <<"##################################"<<endl;
      cout << endl;
      cout << "ERBL U Area: " << (f2u - (erbl_sum)) << "    ERBL D Area: "\
	   << (f2d - (erbl_sumd)) << endl;
      cout << "DGLAP U Area: " << sumu/(1-zeta/2) << "    DGLAP D Area: "\
	   << sumd/(1-zeta/2) << endl;
      cout << "Total U Area: " << f2u << "    Total D Area: " << f2d << endl;
      cout << endl;
      cout << "Gluon DGLAP Area: " << erbl_sumg << "    Gluon F2: " << gluon_F2 << endl;
      cout << endl;
    }
    else if(jgpd == 3) {
      cout << endl;
      cout << "#################################"<<endl;
      cout << endl;
      cout << "ERBL U Area: " << (gau - (erbl_sum)) << "    ERBL D Area: "\
	   << (gad - (erbl_sumd)) << endl;
      cout << "DGLAP U Area: " << sumu/(1-zeta/2) << "    DGLAP D Area: "\
	   << sumd/(1-zeta/2) << endl;
      cout << "Total U Area: " << gau << "    Total D Area: " << gad << endl;
      cout << endl;
    }
    else if(jgpd == 4) {
      cout << endl;
      cout << "#################################"<<endl;
      cout << endl;
      cout << "ERBL U Area: " << (gpu - (erbl_sum)) << "    ERBL D Area: "\
	   << (gpd - (erbl_sumd)) << endl;
      cout << "DGLAP U Area: " << sumu/(1-zeta/2) << "    DGLAP D Area: "\
	   << sumd/(1-zeta/2) << endl;
      cout << "Total U Area: " << gpu << "    Total D Area: " << gpd << endl;
      cout << endl;
    }

    //Gives us the GPD evolved at x = zeta
    double gpdzeta = fut(xval,temp,zeta);
    double gpdzetad = fut(xval,tempd,zeta);
    double gpdzetag = fut(xval,tempg,zeta);
    double gpdzetaub = fut(xval,tempub,zeta);
    double gpdzetadb = fut(xval,tempdb,zeta);
    
    //Define the Imaginary part of the CFF as pi*F(zeta,zeta,t)
    double imcff = pi*gpdzeta ;
    double imcffdv = pi*gpdzetad;
    double imcffg = pi*gpdzetag;
    double imcffu = pi*(gpdzeta + gpdzetaub);
    double imcffd = pi*(gpdzetad + gpdzetadb);
    
    //Output what the GPD is at the crossover point
    cout << endl;
    cout << "GPD_u(zeta): " << gpdzeta << endl;
    cout << "GPD_d(zeta): " << gpdzetad << endl;
    cout << "GPD_g(zeta): " << gpdzetag << endl;

    //*****************************
    //
    //  GPD H --> ReH & ImH
    //
    //*****************************
    
    if(jgpd == 1) {
      
      //ERBL GPD flavor separated arrays
      double erbl_huv[490];
      double erbl_hdv[490];
      double erbl_hu[490];
      double erbl_hd[490];
      double erbl_hg[490];

      
      //Symmetry variables for ERBL region
      double Huplus,Huminus,interval,Huq;
      double Hdplus,Hdminus,Hdq;
      double Hgminus;
      
      //*******************************************
      //
      //   Symmetries
      //
      //   q(x) = -qb(-x)
      //
      //   FLAVOR Non-SINGLET valence/gluon
      //   Hminus = Hq(x,zeta,t) - Hqb(x,zeta,t)
      //
      //   FLAVOR SINGLET sea
      //   Hplus = Hq(x,zeta,t) + Hqb(x,zeta,t)
      //
      //*******************************************

      //Goes through xbj and calculates ERBL and DGLAP regions of evolved
      //GPD
      for(int i = 0; i<490; i+=increment) {
	if(xval[i] <= zeta) {
	  interval = xval[i]-xval[i-increment];

	  //Calculates ERBL Region for u-quark
	  Huplus = Herblevol(xval[i],zeta,t,sumu + 2*sumub,gpdzeta + 2*gpdzetaub,f1u,0,0);
	  Huminus = Herblevol(xval[i],zeta,t,sumu,gpdzeta,f1u,1,0);
	  //	  Huminus = Herblevolhd(xval[i],zeta,t,sumu,sumu2,sumu3,\
	  gpdzeta,f1u,a20u,a30u,a32u,c20u,1,0);

	  erbl_sumub += 2*0.5*(Huplus - Huminus)*interval/(1-zeta/2)/(1.-zeta/2);
	  
	  //Total u-quark distribution
	  Huq = 0.5*(Huplus + Huminus);

	  //Total sum should be u-quark Dirac form factor
	  erbl_sum += Huq *interval/(1-zeta/2);

	  //Calculates ERBL Region for d-quark
	  Hdplus = Herblevol(xval[i],zeta,t,sumd + 2*sumdb,gpdzetad+2*gpdzetadb,f1d,0,1);
	  Hdminus = Herblevol(xval[i],zeta,t,sumd,gpdzetad,f1d,1,1);

	  erbl_sumdb += 2*0.5*(Hdplus - Hdminus)*interval/(1-zeta/2)/(1-zeta/2);
	  
	  //Total d-quark distribution
	  Hdq = 0.5*(Hdplus + Hdminus);

	  //Total sum should be d-quark Dirac form factor
	  erbl_sumd += Hdq*interval/(1-zeta/2);

	  //Gluon only have symmetric "minus" component in ERBL region
	  Hgminus = Herblevolg(xval[i],zeta,t,sumg,gpdzetag,gluon_F1,1,2);

	  //Total sum should be gluon Form Factor
	  erbl_sumg += Hgminus*interval/(1.-zeta/2.)/(1.-zeta/2.);
	  
	  //To calculate the Compton form factors for the vector section
	  //We integrate over the minus symmetry component
	  erbl_huv[i] = Huminus;
	  erbl_hu[i] = Huq;
	  erbl_hdv[i] = Hdminus;
	  erbl_hd[i] = Hdq;
	  erbl_hg[i] = Hgminus;
	  
	  //Writes the ERBL region values of the GPD to a file
	  myfile << xval[i] << " " << Huq << " " << Huplus << " " \
		 << Huminus << " " << Hdq << " " << Hdplus << " "\
		 << Hdminus << " " << Hgminus << " " << 0.5*(Huplus-Huminus) << " " \
		 << 0.5*(Hdplus - Hdminus) << " " << sb0[i] << " " << cb0[i] << endl;
	}
	else {
	  
	  //Integration interval
	  interval = xval[i]-xval[i-increment];

	  //Writes the DGLAP region values of the GPD to the file
	  myfile << xval[i] << " " << ( uv0[i]/xval[i]+ ub0[i]) << " " \
		 << uv0[i]/xval[i]+ 2*ub0[i] << " " << uv0[i]/xval[i] << " "	\
		 << dv0[i]/xval[i] + db0[i] << " " << dv0[i]/xval[i] + 2*db0[i]<< " "\
	         << dv0[i]/xval[i] << " " << g0[i] << " " << ub0[i]\
		 << " " << db0[i] <<  " " << sb0[i] << " " << cb0[i] << endl;

	  //DGLAP region added to the integration array
	  erbl_huv[i] = uv0[i]/xval[i];
	  erbl_hu[i] = (uv0[i]/xval[i]) + ub0[i];
	  erbl_hdv[i] = dv0[i]/xval[i];
	  erbl_hg[i] = g0[i];
	}
      }

      //Integration arrays for u and d quark distributions
      double erb[49];
      double erbdv[49];
      double erbg[49];

      double erb2[49];
      double erb2dv[49];
      double erb2g[49];

      double xx[49];
      double xxdv[49];
      double xxg[49];
      
      double xxu[49];
      double erbu[49];
      double erb2u[49];

      double erbd[49];
      double erb2d[49];

      //Fills array for the integration of the u-quark Compton form factors
      for(int i=0;i<49;i++) {

	//Fills xbj array between zeta/2 to 1
	xx[i] = dinter(zeta/2,1,48,i);

	//1/(X-Zeta) portion of the CFF integral
	erb[i] = 0.5*(fut(xval,erbl_huv,xx[i])-fut(xval,erbl_huv,zeta))\
	  /(xx[i]-zeta);

	//1/X portion of the CFF integral
	erb2[i] = 0.5*fut(xval,erbl_huv,xx[i])/xx[i];
      }

      for(int i=0;i<49;i++) {

        //Fills xbj array between zeta/2 to 1
        xxu[i] = dinter(zeta/2,1,48,i);

        //1/(X-Zeta) portion of the CFF integral
        erbu[i] = (fut(xval,erbl_hu,xxu[i])-fut(xval,erbl_hu,zeta))	\
          /(xxu[i]-zeta);

        //1/X portion of the CFF integral
        erb2u[i] = fut(xval,erbl_hu,xxu[i])/xxu[i];
      }
      
      //Fills array for the d-quark CFFs
      for(int i=0;i<49;i++) {

	//Fills xbj between zeta/2 and 1
	xxdv[i] = dinter(zeta/2,1,48,i);

	//1/(X-zeta) portion of teh CFF integral
	erbdv[i] = 0.5*(fut(xval,erbl_hdv,xxdv[i])-fut(xval,erbl_hdv,zeta))\
	  /(xxdv[i]-zeta);

	//1/X portion of the CFF integral
	erb2dv[i] = 0.5*(fut(xval,erbl_hdv,xxdv[i]))/xxdv[i];
      }

      //Fills array for the gluon CFFs
      for(int i=0;i<49;i++) {
	
        //Fills xbj between zeta/2 and 1                                                   
        xxg[i] = dinter(zeta/2,1,48,i);
	
        //1/(X-zeta) portion of teh CFF integral                                           
        erbg[i] = (fut(xval,erbl_hg,xxg[i])-fut(xval,erbl_hg,zeta))\
          /(xxg[i]-zeta);
	
        //1/X portion of the CFF integral                                                  
        erb2g[i] = (fut(xval,erbl_hg,xxg[i]))/xxg[i];
	
      }
      
      //RE CFF_u calculation using gauss integration division by
      //2 for symmetries
      double rehuv = dgaus1(48,zeta/2,1,erb)+0.5*(gpdzeta)*log((1-zeta)	\
        /(zeta/2))+dgaus1(48,zeta/2,1,erb2);
      
      //RE CFF_d calculation using gauss integration division by
      //2 for symmetries
      double rehdv = dgaus1(48,zeta/2,1,erbdv)+0.5*(gpdzetad)*log((1-zeta)\
	/(zeta/2))+dgaus1(48,zeta/2,1,erb2dv);

      double rehd = dgaus1(48,zeta/2,1,erbd) + 0.5*(gpdzetad + gpdzetadb)*log((1-zeta)\
	/(zeta/2)) + dgaus1(48,zeta/2,1,erb2d);
      
      //RE CFF_g calculation using gauss integration division by
      //2 for symmetries
      double rehg = dgaus1(48,zeta/2,1,erbg) + (gpdzetag)*log((1-zeta)  \
	/(zeta/2))+dgaus1(48,zeta/2,1,erb2g);

      double rehu = dgaus1(48,zeta/2,1,erbu) + (gpdzeta + gpdzetaub)*log((1-zeta) \
	/(zeta/2))+dgaus1(48,zeta/2,1,erb2u);
      
      //CFF is multiplied by - alpha_s / 4*pi (minus because of fermion
      //loop in wilson coefficient function
      rehg = -rehg*alpha/4.;
      
      //Imaginary CFF u-quark
      double imhuv = imcff;
      double imhu = imcffu;
      
      //Imaginary CFF d-quark
      double imhdv = imcffdv;
      double imhd = imcffd;

      //Imaginary CFF gluon
      double imhg = -imcffg*alpha/4.;

      //Creates proton CFFs from u and d-quark distributions
      double rehpv = (4./9.)*rehuv + (1./9.)*rehdv;
      double imhpv = (4./9.)*imhuv + (1./9.)*imhdv;
      
      double rehp = (4./9.)*rehu + (1./9.)*rehd;
      double imhp = (4./9.)*imhu + (1./9.)*imhd;
      
      //Rescaling factors for Lattice QCD
      double lattice_rescale_ub = 0.0605407/0.0979;
      double lattice_rescale_db = 0.0718/0.0979;
      
      //Write outputs
      cout << endl;
      cout << "F1_u: " << f1u << " Integral: " << erbl_sum << " F1_d: "\
	   << f1d << " Integral: " << erbl_sumd << endl;
      cout << endl;
      cout << "Re Huv: " << rehuv << "    Im Huv: " << imhuv  << "    Re Hdv: "\
	   << rehdv << "    Im Hdv: " << imhdv << endl;
      cout << "Re Hu: " << rehu << "    Im Hu: " << imhu << "    Re Hd: " \
	   << rehd << "    Im Hd: " << imhd << endl;
      cout << "Re Hp (val): " << rehpv << "    Im Hp (val): " << imhpv << endl;
      cout << "Re Hp (tot): " << rehp << "    Im Hp (tot): " << imhp << endl;

      cout << endl;
      cout << "Gluon Integral: " << erbl_sumg << "   F1_g: " << gluon_F1 << endl;
      cout << "Re Hg: " << rehg << "    Im Hg: " << imhg << endl;

      cout << endl;
      cout << "A20_ub_integral: " << erbl_sumub << "  A20_ub: " << ubar_F1*lattice_rescale_ub << endl;
      cout << "A20_db_integral: " << erbl_sumdb << "  A20_db: " << ubar_F1*lattice_rescale_db << endl;
      cout << endl;
      
      //Writes Compton form factors to the file
      myfile5 << rehp << "," << imhp << ",";
    }

    //**************************
    //
    //    GPD E
    //
    //**************************
    
    else if(jgpd == 2) {

      //ERBL GPD flavor separated arrays
      double erbl_eu[490];
      double erbl_ed[490];
      double erbl_eg[490];

      //Symmetry variables for ERBL region
      double Euplus,Euminus,interval,Euq;
      double Edplus,Edminus,Edq;
      double Egminus;
      //*******************************************
      //                      
      //   Symmetries
      //
      //   q(x) = -qb(-x)
      //
      //   FLAVOR Non-SINGLET valence
      //   Eminus = Eq(x,zeta,t) - Eqb(x,zeta,t)
      //
      //   FLAVOR SINGLET sea
      //   Eplus = Eq(x,zeta,t) + Eqb(x,zeta,t)
      //
      //*******************************************

      //Loops over xbj to calculate the DGLAP and ERBL regions of the
      //GPD to calculate the CFFs
      for(int i = 0; i<490; i+=increment) {
	if(xval[i] <= zeta) {
	  interval = xval[i]-xval[i-increment];

	  //Calculates the u-quark ERBL region for GPD E
	  Euplus = Herblevol(xval[i],zeta,t,sumu,gpdzeta,f2u,0,0);
	  Euminus = Herblevol(xval[i],zeta,t,sumu,gpdzeta,f2u,1,0);

	  //Calculates the GPD E for u-quark distribution
	  Euq = 0.5*(Euplus + Euminus);

	  //Sum should be the u-quark Pauli form factor
	  erbl_sum += Euq *interval/(1-zeta/2);

	  //Calculates the d-quark ERBL region for GPD E
	  Edplus = Herblevol(xval[i],zeta,t,sumd,gpdzetad,f2d,0,1);
	  Edminus = Herblevol(xval[i],zeta,t,sumd,gpdzetad,f2d,1,1);

	  //GPD E for d-quark distribution
	  Edq = 0.5*(Edplus + Edminus);

	  //Sum should be d-quark Pauli form factor
	  erbl_sumd += Edq*interval/(1-zeta/2);

	  //Gluon only have symmetric "minus" component in ERBL region 
          Egminus = Herblevolg(xval[i],zeta,t,sumg,gpdzetag,gluon_F2,1,2);

          //Total sum should be gluon Form Factor         
          erbl_sumg += Egminus*interval/(1.-zeta/2.)/(1.-zeta/2.);
	  
	  //ERBL region added to array
	  erbl_eu[i] = Euminus;
	  erbl_ed[i] = Edminus;
	  erbl_eg[i] = Egminus;

	  //Writes ERBL region to file
	  myfile2 << xval[i] << " " << Euq << " " << " " << Euplus\
		  << " " << Euminus << " " << Edq << " " << Edplus\
	          << " " << Edminus << " " << Egminus << endl;
	}
	else {
	  interval = xval[i]-xval[i-increment];

	  //Writes DGLAP region to file
	  myfile2 << xval[i] << " " << uv0[i]/xval[i] << " "\
		  << uv0[i]/xval[i] << " " << uv0[i]/xval[i] << " " \
		  << dv0[i]/xval[i] << " " << dv0[i]/xval[i] << " " \
	          << dv0[i]/xval[i] << " " << g0[i] << endl;

	  //DGLAP region added to array
	  erbl_eu[i] = uv0[i]/xval[i];
	  erbl_ed[i] = dv0[i]/xval[i];
	  erbl_eg[i] = g0[i];
	}
      }

      //Integration arrays for u and d-quark distributions
      double erb3[49];
      double erb3d[49];
      double erb4[49];
      double erb4d[49];
      double xx2[49];
      double xx2d[49];
      
      //Fills integration arrays for u-quark Compton form factors
      for(int i=0;i<49;i++) {

	//xbj array for gauss integral
	xx2[i] = dinter(zeta/2,1,48,i);

	//1/(X-zeta) array for u-quark filled
	erb3[i] = 0.5*(fut(xval,erbl_eu,xx2[i])-fut(xval,erbl_eu,zeta))\
	  /(xx2[i]-zeta);

	//1/X array for u-quark filled
	erb4[i] = 0.5*fut(xval,erbl_eu,xx2[i])/xx2[i];
      }

      //Fills integration arrays for d-quark Compton form factors
      for(int i=0;i<49;i++) {

	//xbj array for gauss integral
	xx2d[i] = dinter(zeta/2,1,48,i);

	//1/(X-zeta) array for d-quark filled
	erb3d[i] = 0.5*(fut(xval,erbl_ed,xx2d[i])-fut(xval,erbl_ed,zeta))\
	  /(xx2d[i]-zeta);

	//1/X array for d-quark filled
	erb4d[i] = 0.5*(fut(xval,erbl_ed,xx2d[i]))/xx2d[i];
      }

      //Calculates u-quark Re part of the compton form factor
      //1/2 factor is included as part of the symmetry
      double reeu = dgaus1(48,zeta/2,1,erb3)+0.5*(gpdzeta)*log((1-zeta)\
	  /(zeta/2))+dgaus1(48,zeta/2,1,erb4);

      //Calculates d-quark Re part of the compton form factor
      //1/2 factor is included as part of the symmetry
      double reed = dgaus1(48,zeta/2,1,erb3d)+0.5*(gpdzetad)*log((1-zeta)\
	  /(zeta/2))+dgaus1(48,zeta/2,1,erb4d);

      //Imaginary parts of the Compton form factors
      //u-quark
      double imeu = imcff;

      //d-quark
      double imed = imcffd;

      //Creates the Real Compton form factors of the proton
      //Using charges for the u-quark and d-quark
      double reep = (4./9.)*reeu + (1./9.)*reed;
      double imep = (4./9.)*imeu + (1./9.)*imed;

      double cppi2 = (-t/(4*M*M))*((2./3.)*f2u - (1./3.)*f2d)*reep;

      cout << "Cppi_E: " << cppi2 << endl;
      
      //Prints output
      cout << endl;
      cout << "F2_u: " << f2u << " Integral: " << erbl_sum << " F2_d: "\
	   << f2d << " Integral: " << erbl_sumd << endl;
      cout << endl;
      cout << "Re Eu: " << reeu << "    Im Eu: " << imeu << "    Re Ed: "\
	   << reed << "    Im Ed: " << imed << endl;
      cout << "Re Ep: " << reep << "    Im Ep: " << imep << endl;

      cout << endl;
      cout << "Gluon Integral: " << erbl_sumg << "   F1_g: " << gluon_F2 << endl;
      //Writes the Compton form factors to file
      myfile5 << reep << "," << imep << ",";
    }

    //************************
    //
    //    GPD Htil
    //
    //************************

    else if (jgpd ==3) {

      //ERBL GPD flavor separated arrays
      double erbl_hut[490];
      double erbl_hdt[490];
      
      //Defines symmetry arrays
      double Hutplus,Hutminus,interval,Hutq;
      double Hdtplus,Hdtminus,Hdtq;

      //*******************************************        
      //       
      //   Symmetries
      //
      //   Helicity Distributions!
      //
      //   dq(x) = dqb(-x)
      //
      //   FLAVOR Non-SINGLET valence
      //   H_tminus = H_tq(x,zeta,t) - H_tqb(x,zeta,t)
      //
      //   FLAVOR SINGLET sea
      //   H_tplus = H_tq(x,zeta,t) + H_tqb(x,zeta,t)
      //
      //*******************************************

      //Iterates over xbj, separated between ERBL and DGLAP regions
      for(int i = 0; i<490; i+=increment) {

	//ERBL region of X < zeta
        if(xval[i] <= zeta) {
          interval = xval[i]-xval[i-increment];

	  //Plus distribution and minus distribution symmetries for
	  //ERBL u-quark
	  //Definition here is backwards!!!
	  //For the helicity distributions the "plus" distribution is symmetric
	  //and the "minus" distribution is antisymmetric
	  //
	  //I've put it backwards but I pass the anti-symmetric piece - the valence
	  //piece - to the CFF.
          Hutplus = Herblevol(xval[i],zeta,t,sumu,gpdzeta,gau,0,0);
          Hutminus = Herblevol(xval[i],zeta,t,sumu,gpdzeta,gau,1,0);

	  //Total quark distribution which is sum / 2
          Hutq = 0.5*(Hutplus + Hutminus);

	  //Integral which should equal axial u-quark form factor in the end
          erbl_sum += Hutq *interval/(1-zeta/2);

	  //Plus and minus distribution symmetries for ERBL d-quark
          Hdtplus = Herblevol(xval[i],zeta,t,sumd,gpdzetad,gad,0,1);
          Hdtminus = Herblevol(xval[i],zeta,t,sumd,gpdzetad,gad,1,1);

	  //Total d-quark distribution which is sum/2
          Hdtq = 0.5*(Hdtplus + Hdtminus);

	  //Integral which should equal axial d-quark form factor in the end
          erbl_sumd += Hdtq*interval/(1-zeta/2);

	  //Fills ERBL arrays used to integrate later
          erbl_hut[i] = Hutplus;
          erbl_hdt[i] = Hdtplus;

	  //Writes GPDs to file, total quark distribution
	  //separated between u and d
          myfile3 << xval[i] << " " << Hutq << " " << Hutplus << " " \
	          << Hutminus << " " << Hdtq << " " << Hdtplus << " " \
	          << Hdtminus << " " << g0[i] << endl;
        }
	else {

	  //DGLAP region
          interval = xval[i]-xval[i-increment];

	  //put DGLAP GPD in file
          myfile3 << xval[i] << " " << uv0[i]/xval[i] << " "\
	          << uv0[i]/xval[i] << " " << uv0[i]/xval[i] << " "\
		  << dv0[i]/xval[i] << " " << dv0[i]/xval[i] << " "\
	          << dv0[i]/xval[i] << " " << g0[i] << endl;

	  //Fills array with the rest of the GPD 
          erbl_hut[i] = uv0[i]/xval[i];
          erbl_hdt[i] = dv0[i]/xval[i];
        }
      }

      //Integration arrays used for gauss integral algorithm
      double erb5[49];
      double erb5d[49];
      double erb6[49];
      double erb6d[49];
      double xx5[49];
      double xx5d[49];

      //u-quark gauss integral arrays filled
      for(int i=0;i<49;i++) {

	//xbj array filled
        xx5[i] = dinter(zeta/2,1,48,i);

	//1/(X-zeta) integral array filled
        erb5[i] = 0.5*(fut(xval,erbl_hut,xx5[i])-fut(xval,erbl_hut,zeta))/(xx5[i]-zeta);

	//1/X integral array filled
        erb6[i] = 0.5*fut(xval,erbl_hut,xx5[i])/xx5[i];
      }

      //d-quark gauss integral arrays filled
      for(int i=0;i<49;i++) {

	//xbj array filled
        xx5d[i] = dinter(zeta/2,1,48,i);

	//1/(X-zeta) integral array filled for d-quark
        erb5d[i] = 0.5*(fut(xval,erbl_hdt,xx5d[i])-fut(xval,erbl_hdt,zeta))\
	  /(xx5d[i]-zeta);

	//1/X integral array filled for d-quark
        erb6d[i] = 0.5*(fut(xval,erbl_hdt,xx5d[i]))/xx5d[i];
      }

      //ReH_t for u-quark Compton form factor
      double rehut = dgaus1(48,zeta/2,1,erb5)+0.5*(gpdzeta)*log((1-zeta)\
	  /(zeta/2))-dgaus1(48,zeta/2,1,erb6);

      //ReH_t for d-quark Compton form factor
      double rehdt = dgaus1(48,zeta/2,1,erb5d)+0.5*(gpdzetad)*log((1-zeta)\
	  /(zeta/2))-dgaus1(48,zeta/2,1,erb6d);

      //ImH_t u-quark Compton form factor
      double imhut = imcff;

      //ImH_t d-quark Compton form factor
      double imhdt = imcffd;

      //ReH_t for the proton
      double rehpt = (4./9.)*rehut + (1./9.)*rehdt;

      //ImH_t for the proton
      double imhpt = (4./9.)*imhut + (1./9.)*imhdt;

      //Prints output
      cout << endl;
      cout << "GA_u: " << gau << " Integral: " << erbl_sum << " GA_d: "\
	   << gad << " Integral: " << erbl_sumd << endl;
      cout << endl;
      cout << "Re Hu_t: " << rehut << "    Im Hu_t: " << imhut\
	   << "    Re Hd_t: " << rehdt << "    Im Hd_t: " << imhdt << endl;
      cout << "Re Hp_t: " << rehpt << "    Im Hp_t: " << imhpt << endl;

      //Writes Compton form factors to File!
      myfile5 << rehpt << "," << imhpt << "," ;
    }
   
    //***********************
    //
    //    GPD Etil
    //
    //***********************
    
    else if (jgpd == 4) {

      //ERBL GPD flavor separated arrays
      double erbl_eut[490];
      double erbl_edt[490];
      
      //Defines symmetry variables
      double Eutplus,Eutminus,interval,Eutq;
      double Edtplus,Edtminus,Edtq;

      //*******************************************
      //
      //   Symmetries
      //
      //   Helicity Distributions!
      //
      //   dq(x) = dqb(-x)
      //
      //   FLAVOR Non-SINGLET valence
      //   H_tminus = H_tq(x,zeta,t) - H_tqb(x,zeta,t)
      //
      //   FLAVOR SINGLET sea
      //   H_tplus = H_tq(x,zeta,t) + H_tqb(x,zeta,t)
      //
      //*******************************************

      //Iterates over xbj variable between ERBL and DGLAP regions
      for(int i = 0; i<490; i+=increment) {

	//ERBL region
        if(xval[i] <= zeta) {
          interval = xval[i]-xval[i-increment];

	  //Plus and minus distributions called by evol_erbl
          Eutplus = Herblevol(xval[i],zeta,t,sumu,gpdzeta,gpu,0,0);
          Eutminus = Herblevol(xval[i],zeta,t,sumu,gpdzeta,gpu,1,0);

	  //Total quark distribution given by sum/2
          Eutq = 0.5*(Eutplus + Eutminus);

	  //Sum which will integrate to u-quark pseudoscalar form factor
          erbl_sum += Eutq *interval/(1-zeta/2);
	  
	  //Plus and minus distributions for d-quark called by evol_erbl
          Edtplus = Herblevol(xval[i],zeta,t,sumd,gpdzetad,gpd,0,1);
          Edtminus = Herblevol(xval[i],zeta,t,sumd,gpdzetad,gpd,1,1);

	  //Total d-quark distribution given by sum/2
          Edtq = 0.5*(Edtplus + Edtminus);

	  //Integrates to d-quark pseudoscalar form factor 
          erbl_sumd += Edtq*interval/(1-zeta/2);

	  //Fills ERBL array for u and d quark distributions
	  erbl_eut[i] = Eutplus;
          erbl_edt[i] = Edtplus;

	  //Writes ERBL region GPD to the file
          myfile4 << xval[i] << " " << Eutq << " " << Eutplus << " "\
	          << Eutminus << " " << Edtq << " " << Edtplus << " "\
	          << Edtminus << endl;
        }
        else {
          interval = xval[i]-xval[i-increment];

	  //Writes DGLAP region GPDs to file
          myfile4 << xval[i] << " " << uv0[i]/xval[i] << " " \
	          << uv0[i]/xval[i] << " " << uv0[i]/xval[i] << " "\
		  << dv0[i]/xval[i] << " "  << dv0[i]/xval[i] << " "\
                  << dv0[i]/xval[i] << " " << g0[i] << endl;

	  //Fills array with DGLAP GPD
          erbl_eut[i] = uv0[i]/xval[i];
          erbl_edt[i] = dv0[i]/xval[i];
        }
      }

      //Integration arrays for gauss integral declared
      double erb7[49];
      double erb7d[49];
      double erb8[49];
      double erb8d[49];
      double xx7[49];
      double xx7d[49];

      //Gauss integral arrays filled for u-quark
      for(int i=0;i<49;i++) {

	//xbj array filled
        xx7[i] = dinter(zeta/2,1,48,i);

	//1/(x-zeta) array for gauss integral filled for u-quark
        erb7[i] = 0.5*(fut(xval,erbl_eut,xx7[i])-fut(xval,erbl_eut,zeta))/(xx7[i]-zeta);

	//1/X array for gauss integral filled for u-quark
        erb8[i] = 0.5*fut(xval,erbl_eut,xx7[i])/xx7[i];
      }

      //Gauss integral arrays filled for d-quark
      for(int i=0;i<49;i++) {

	//xbj array filled
        xx7d[i] = dinter(zeta/2,1,48,i);

	//1/(X-zeta) array for gauss integral filled for d-quark
        erb7d[i] = 0.5*(fut(xval,erbl_edt,xx7d[i])-fut(xval,erbl_edt,zeta))\
	  /(xx7d[i]-zeta);

	//1/X array for gauss integral filled for d-quark
        erb8d[i] = 0.5*(fut(xval,erbl_edt,xx7d[i]))/xx7d[i];
      }

      //Calculates ReE_t for u-quark using gauss integral routine
      double reeut = dgaus1(48,zeta/2,1,erb7)+0.5*(gpdzeta)\
	*log((1-zeta)/(zeta/2))-dgaus1(48,zeta/2,1,erb8);

      //Calculates ReE_t for d-quark using gauss integral routine
      double reedt = dgaus1(48,zeta/2,1,erb7d)+0.5*(gpdzetad)\
	*log((1-zeta)/(zeta/2))-dgaus1(48,zeta/2,1,erb8d);

      //Im E_t for u-quark distribution
      double imeut = imcff;

      //Im E_t for d-quark distribution
      double imedt = imcffd;

      //ReE_t for the proton using u and d quark charges
      double reept = (4./9.)*reeut + (1./9.)*reedt;

      //ImE_t for the proton using u and d quark charges
      double imept = (4./9.)*imeut + (1./9.)*imedt;\

      //Prints output
      cout << endl;
      cout << "GP_u: " << gpu << " Integral: " << erbl_sum \
	   << " GP_d: " << gpd << " Integral: " << erbl_sumd << endl;
      cout << endl;
      cout << "Re Eu_t: " << reeut << "    Im Eu_t: " << imeut\
	   << "    Re Ed_t: " << reedt << "    Im Ed_t: " << imedt << endl;
      cout << "Re Ep_t: " << reept << "    Im Ep_t: " << imept << endl;

      //Writes Compton form factors to the file
      myfile5 << reept << "," << imept << endl;

    }
  }

  //Endline printed for formatting
  cout << endl;
}

//Madness is over!!!
