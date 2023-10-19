//***********************************************************//
//                                                           //
//    DIQUARK CALCULATION -- Compton Form Factors            //
//                                                           //
//    AUTHOR: Brandon Kriesten                               //
//                                                           //
//    REFERENCE: PHYSICAL REVIEW D 84, 034007 (2011)         //
//               PHYSICAL REVIEW D                           //
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
#include "evol_erbl_new.h"

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
void cff_calc(double xbj, double tt, double Q2, int jgpd,\
	      double& repv, double& impv, double& reuv,\
	      double& imuv, double& redv, double& imdv,\
	      double& rep, double& imp, double& reu,\
	      double& imu, double& red, double& imd,\
	      double &reg, double& img) {

  //Creates 5 files
  ofstream myfile,myfile2,myfile3,myfile4,myfile5;
  ofstream myfile6,myfile7,myfile8;
  
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
  double final_scale = Q2;
  
  //Returns alpha_s / pi
  double alpha = alpha_my_(&final_scale);
  double mz = 91.1867;
  double mz2 = mz*mz;
  double alpha_s_mz2 = alpha_my_(&mz2)*pi;
  cout << "Alpha_s at Mz^2: " << alpha_s_mz2 << endl;

  //zeta = X_bj
  double zeta = xbj;
  double xi = zeta/(2.-zeta);

  //4-momentum transfer squared
  //  double t = -zeta*zeta*M*M/(1-zeta) - 0.0001;
  double t = tt;
    
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

  //DOES THIS NEED A FACTOR OF 4???
  double gluon_F1 = gluon_a + 1*xi*xi*gluon_d;
  double gluon_F2 = gluon_b - 1*xi*xi*gluon_d;

  //  double gluon_F1 = gluon_a + xi*xi*gluon_d;
  //  double gluon_F2 = gluon_b - xi*xi*gluon_d;
  
  //Sea Quark Form Factors
  double ubar_F1 = a20_ub_dip(-t,1) \
    + 4*xi*xi*0.5*(c20_upd_zexp(-t) + c20_umd_zexp(-t));
  double dbar_F1 = a20_db_dip(-t,1) \
    + 4*xi*xi*0.5*(c20_upd_zexp(-t) - c20_umd_zexp(-t));

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

  if(final_scale < 0.97) {

    double init_q02 = 0.09362;
    cout << "Evolving valence distributions from " \
	 << init_q02 << " to " << final_scale << endl;
    
    struct_ev_gpd_tmd_(&final_scale, &zeta, &t, &jgpd,\
		       &ierr, &itmd, &jt, uv0,dv0, g0,\
		       ub0, db0, sb0, cb0, &init_q02);
  }

  else {
    //Calls the evolution code for valence distribution
    //evolves up to scale of gluons
    double init_q02 = 0.09362;
    double qsquare = 0.97;
    cout << "Evolving valence distributions from " \
	 << init_q02 << " to " << qsquare << endl;
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
      tempfile << xval[i] << " " << uv0_t[i]/xval[i] \
	       << " " << dv0_t[i]/xval[i] << " " \
	       << ub0_t[i] << " " << db0_t[i] << endl;
    }
    tempfile.close();

    init_q02 = 0.97;
    qsquare = final_scale;
    cout << "Evolving valence/glue/sea distributions from " \
	 << init_q02 << " to " << qsquare << endl;
    //Calls the evolution code evolves from initial scale of gluons to final Q2
    struct_ev_gpd_tmd_(&qsquare, &zeta, &t, &jgpd, &ierr, &itmd, &jt, uv0, \
		       dv0, g0, ub0, db0, sb0, cb0, &init_q02);
  }
  //Calculates integral of evolved gpd 
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
    cout << "ERBL Uv Area: " << (f1u - (erbl_sum)) << "    ERBL Dv Area: " \
	 << (f1d - (erbl_sumd)) << endl;
    cout << "DGLAP Uv Area: " << sumu/(1-zeta/2) << "    DGLAP Dv Area: " \
	 << sumd/(1-zeta/2) << endl;
    cout << "Total Uv Area: " << f1u << "    Total Dv Area: "	\
	 << f1d << endl;
    cout << endl;
    cout << endl;
    cout << "ERBL G Area: " << (gluon_F1 - erbl_sumg)		\
	 << "    ERBL Ub Area: " << ubar_F1 - erbl_sumub		\
	 << "    ERBL Db Area: " << dbar_F1 - erbl_sumdb << endl;
    cout << "DGLAP G Area: " << erbl_sumg << "    DGLAP Ub Area: "	\
	 << erbl_sumub << "    DGLAP Db Area: " << erbl_sumdb << endl;
    cout << "Total G Area: " << gluon_F1				\
	 << "    Total Ub Area: " << ubar_F1 << "    Total Dbar Area: "	\
	 << dbar_F1 << endl;
    cout << endl;
  }

  else if(jgpd == 2) {
    cout << endl;
    cout <<"##################################"<<endl;
    cout << endl;
    cout << "ERBL U Area: " << (f2u - (erbl_sum)) << "    ERBL D Area: " \
	 << (f2d - (erbl_sumd)) << endl;
    cout << "DGLAP U Area: " << sumu/(1-zeta/2) << "    DGLAP D Area: "	\
	 << sumd/(1-zeta/2) << endl;
    cout << "Total U Area: " << f2u << "    Total D Area: " << f2d << endl;
    cout << endl;
    cout << "Gluon DGLAP Area: " << erbl_sumg << "    Gluon F2: " \
	 << gluon_F2 << endl;
    cout << endl;
  }

  
  else if(jgpd == 3) {
    cout << endl;
    cout << "#################################"<<endl;
    cout << endl;
    cout << "ERBL U Area: " << (gau - (erbl_sum)) << "    ERBL D Area: " \
	 << (gad - (erbl_sumd)) << endl;
    cout << "DGLAP U Area: " << sumu/(1-zeta/2) << "    DGLAP D Area: "	\
	 << sumd/(1-zeta/2) << endl;
    cout << "Total U Area: " << gau << "    Total D Area: " << gad << endl;
    cout << endl;
  }

  
  else if(jgpd == 4) {
    cout << endl;
    cout << "#################################"<<endl;
    cout << endl;
    cout << "ERBL U Area: " << (gpu - (erbl_sum)) << "    ERBL D Area: " \
	 << (gpd - (erbl_sumd)) << endl;
    cout << "DGLAP U Area: " << sumu/(1-zeta/2) << "    DGLAP D Area: "	\
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
  double imcff = pi*gpdzeta;
  double imcffdv = pi*gpdzetad;
  double imcffg = pi*gpdzetag;
  double imcffu = pi*(gpdzeta + 2*gpdzetaub/zeta);
  double imcffd = pi*(gpdzetad + 2*gpdzetadb/zeta);

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
    
    int Qname = (int)(Q2);
    int ttname = (int)(-100 * tt);
    if(ttname < 10) {
      myfile.open("gpd-data/Qsq" + to_string(Qname) + "/GPD_H00" + to_string(ttname) + ".dat");
    }
    else if(ttname < 100) {
      myfile.open("gpd-data/Qsq" + to_string(Qname) + "/GPD_H0" + to_string(ttname) + ".dat");
    }
    else {
      myfile.open("gpd-data/Qsq" + to_string(Qname) + "/GPD_H" + to_string(ttname) + ".dat");
    }
    
//    myfile.open("GPD_H.dat");
    
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

    //SERBL for the plus distributions don't matter because the area of the
    //ERBL region is 0 for the sea quarks.
    
    double Serblhuminus = (1-(zeta/2))*(f1u - (sumu/(1-(zeta/2))));
    double Serblhuplus = (1-(zeta/2))*(f1u - ((sumu+2*sumub)/(1-(zeta/2))));
    double Serblhdminus = (1-(zeta/2))*(f1d - (sumd/(1-(zeta/2))));
    double Serblhdplus = (1-(zeta/2))*(f1d - ((sumd+2*sumdb)/(1-(zeta/2))));

    for(int i = 0; i<491; i+=increment) {
      if(xval[i] <= zeta) {
	interval = xval[i]-xval[i-increment];

	//Calculates ERBL Region for u-quark
	Huplus = Herblevolnew(xval[i],zeta,Serblhuplus,gpdzeta \
			      + (2*gpdzetaub/zeta),0,0);
	Huminus = Herblevolnew(xval[i],zeta,Serblhuminus,gpdzeta,1,0);

	//	Huminus = Herblevolhd(xval[i],zeta,t,sumu\
	,sumu2,sumu3,gpdzeta,f1u,a20u,a30u,a32u,c20u,1,0);
	
	erbl_sumub += 2*0.5*(Huplus-Huminus)*interval/(1-zeta/2)/(1.-zeta/2);
	  
	//Total u-quark distribution
	Huq = 0.5*(Huplus + Huminus);
	
	//Total sum should be u-quark Dirac form factor
	erbl_sum += 0.5*Huminus *interval/(1-zeta/2);
	
	//Calculates ERBL Region for d-quark
	Hdplus = Herblevolnew(xval[i],zeta,Serblhdplus,gpdzetad\
			      +(2*gpdzetadb/zeta),0,1);
	Hdminus = Herblevolnew(xval[i],zeta,Serblhdminus,gpdzetad,1,1);
	
	erbl_sumdb += 2*0.5*(Hdplus-Hdminus)*interval/(1-zeta/2)/(1-zeta/2);
	
	//Total d-quark distribution
	Hdq = 0.5*(Hdplus + Hdminus);
	
	//Total sum should be d-quark Dirac form factor
	erbl_sumd += 0.5*Hdminus*interval/(1-zeta/2);
	
	//Gluon only have symmetric "minus" component in ERBL region
	Hgminus = Herblevolg(xval[i],zeta,t,sumg,gpdzetag,gluon_F1,1,2);
	
	//Total sum should be gluon Form Factor
	erbl_sumg += 0.5*Hgminus*interval/(1.-zeta/2.)/(1.-zeta/2.);
	
	//To calculate the Compton form factors for the vector section
	//We integrate over the minus symmetry component
	erbl_huv[i] = Huminus;
	erbl_hu[i] = Huplus;
	erbl_hdv[i] = Hdminus;
	erbl_hd[i] = Hdplus;
	erbl_hg[i] = Hgminus;
	
	  //Writes the ERBL region values of the GPD to a file
	myfile << xval[i] << " " << Huq << " " << Huplus << " "		\
	       << Huminus << " " << Hdq << " " << Hdplus << " "		\
	       << Hdminus << " " << Hgminus << " " \
	       << 0.5*(Huplus-Huminus) << " "				\
	       << 0.5*(Hdplus - Hdminus) << " " << sb0[i] \
	       << " " << cb0[i] << endl; 
	//	myfile << xval[i] << " " << Huplus << " " << Hdplus << endl;
      }
      else {
	
	//Integration interval
	interval = xval[i]-xval[i-increment];
	
	//Writes the DGLAP region values of the GPD to the file
	myfile << xval[i] << " " << ( uv0[i]/xval[i]+ ub0[i]) << " "	\
	       << uv0[i]/xval[i]+ 2*ub0[i] << " " << uv0[i]/xval[i] << " " \
	       << dv0[i]/xval[i] + db0[i] << " "\
	       << dv0[i]/xval[i] + 2*db0[i]<< " "			\
	       << dv0[i]/xval[i] << " " << g0[i] << " " << ub0[i]	\
	       << " " << db0[i] <<  " " << sb0[i] << " " << cb0[i] << endl;
	/*	myfile << xval[i] << " " << uv0[i]/xval[i]+ 2*ub0[i]/xval[i] << " " \
		<< dv0[i]/xval[i] + 2*db0[i]/xval[i] << endl;*/
	
	//DGLAP region added to the integration array
	  erbl_huv[i] = uv0[i]/xval[i];
	  erbl_hu[i] = (uv0[i]/xval[i]) + 2*ub0[i]/xval[i];
	  erbl_hd[i] = (dv0[i]/xval[i]) + 2*db0[i]/xval[i];
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

    double xxd[49];
    double erbd[49];
    double erb2d[49];
    
    //Fills array for the integration of the u-quark Compton form factors
    for(int i=0;i<49;i++) {
      
      //Fills xbj array between zeta/2 to 1
      xx[i] = dinter(zeta/2,1,48,i);
      
      //1/(X-Zeta) portion of the CFF integral
      erb[i] = (fut(xval,erbl_huv,xx[i])-fut(xval,erbl_huv,zeta))	\
	/(xx[i]-zeta);
      
      //1/X portion of the CFF integral
      erb2[i] = fut(xval,erbl_huv,xx[i])/xx[i];
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
    
    for(int i=0;i<49;i++) {

      //Fills xbj array between zeta/2 to 1
      xxd[i] = dinter(zeta/2,1,48,i);

      //1/(X-Zeta) portion of the CFF integral
      erbd[i] = (fut(xval,erbl_hd,xxu[i])-fut(xval,erbl_hd,zeta))       \
        /(xxd[i]-zeta);

      //1/X portion of the CFF integral
      erb2d[i] = fut(xval,erbl_hd,xxu[i])/xxd[i];
    }
    
    //Fills array for the d-quark CFFs
    for(int i=0;i<49;i++) {
      
      //Fills xbj between zeta/2 and 1
      xxdv[i] = dinter(zeta/2,1,48,i);
      
      //1/(X-zeta) portion of teh CFF integral
      erbdv[i] = (fut(xval,erbl_hdv,xxdv[i])-fut(xval,erbl_hdv,zeta)) \
	/(xxdv[i]-zeta);
      
      //1/X portion of the CFF integral
      erb2dv[i] = (fut(xval,erbl_hdv,xxdv[i]))/xxdv[i];
    }
    
    //Fills array for the gluon CFFs
    for(int i=0;i<49;i++) {
      
      //Fills xbj between zeta/2 and 1
      xxg[i] = dinter(zeta/2,1,48,i);
	
      //1/(X-zeta) portion of teh CFF integral
      erbg[i] = (fut(xval,erbl_hg,xxg[i])-fut(xval,erbl_hg,zeta))	\
	/(xxg[i]-zeta);
      
      //1/X portion of the CFF integral
      erb2g[i] = (fut(xval,erbl_hg,xxg[i]))/xxg[i];
      
    }
    
    //RE CFF valence calculation using gauss integration
    double rehuv = (dgaus1(48,zeta/2,1,erb)+(gpdzeta)*log((1-zeta)	\
       		      /(zeta/2))+dgaus1(48,zeta/2,1,erb2));

    double rehdv = (dgaus1(48,zeta/2,1,erbdv)+(gpdzetad)*log((1-zeta) \
		      /(zeta/2))+dgaus1(48,zeta/2,1,erb2dv));
    
    //RE CFF quark calculation using gauss integration
    double rehu = -1*(dgaus1(48,zeta/2,1,erbu) + (gpdzeta)*log((1-zeta) \
                        /(zeta/2)) + dgaus1(48,zeta/2,1,erb2u));
    
    double rehd = -1*(dgaus1(48,zeta/2,1,erbd) + (gpdzetad)*log((1-zeta) \
			/(zeta/2)) + dgaus1(48,zeta/2,1,erb2d));

    //RE CFF_g calculation using gauss integration
    double rehg = dgaus1(48,zeta/2,1,erbg) + (gpdzetag)*log((1-zeta)	\
			/(zeta/2))+dgaus1(48,zeta/2,1,erb2g);
    rehg = rehg/(2*xi);

    //Imaginary CFF u-quark
    double imhuv = imcff;
    double imhus = imcffu;
    
    //Imaginary CFF d-quark
    double imhdv = imcffdv;
    double imhds = imcffd;

    double imhg = imcffg;
    imhg = imhg;
    
    //Creates proton CFFs from u and d-quark distributions
    double rehpv = (4./9.)*rehuv + (1./9.)*rehdv;
    double imhpv = (4./9.)*imhuv + (1./9.)*imhdv;

    //    cout << rehu << " " << rehd << endl;
    
    double rehp = (4./9.)*rehu + (1./9.)*rehd;
    double imhps = (4./9.)*imhus + (1./9.)*imhds;

    double lattice_rescale_ub = 0.0605/0.0979;
    double lattice_rescale_db = 0.0718/0.0979;

    //    cout << "ReHpv: " << rehpv << " ImHpv: " << imhpv << endl;
    cout << "ReHp: " << rehp << " ImHp: " \
	 << imhps << endl;
    //    cout << "ReHg: " << rehg << " ImHg: " << imhg << endl;
    
    //Writes Compton form factors to the file
    myfile5 << rehp << "," << imhps << ",";
    
    repv = rehpv;
    impv = imhpv;
    reuv = rehuv;
    imuv = imhuv;
    redv = rehdv;
    imdv = imhdv;
    rep = rehp;
    imp = imhps;
    reu = rehu;
    imu = imhus;
    red = rehd;
    imd = imhds;
    reg = rehg;
    img = imhg;

    myfile.close();
    
  }
    //**************************
    //
    //    GPD E
    //
    //**************************
    
  else if(jgpd == 2) {
    myfile2.open("GPD_E.dat");
    //ERBL GPD flavor separated arrays
    double erbl_eu[490];
    double erbl_ed[490];
    double erbl_eg[490];

    double erbl_euplus[490];
    double erbl_edplus[490];
    
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
	erbl_euplus[i] = Euplus;
	erbl_edplus[i] = Edplus;
	
	//Writes ERBL region to file
	myfile2 << xval[i] << " " << Euplus\
		<< " " << Edplus << endl;
      }
      else {
	interval = xval[i]-xval[i-increment];
	
	//Writes DGLAP region to file
	myfile2 << xval[i] << " " << uv0[i]/xval[i]\
		<< " " << dv0[i]/xval[i] << endl;
	
	//DGLAP region added to array
	erbl_eu[i] = uv0[i]/xval[i];
	erbl_ed[i] = dv0[i]/xval[i];
	erbl_euplus[i] = uv0[i]/xval[i];
	erbl_edplus[i] = dv0[i]/xval[i];
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

    double xxuplus[49];
    double xxdplus[49];
    double erbuplus[49];
    double erbuplus2[49];
    double erbdplus[49];
    double erbdplus2[49];
    
    //Fills integration arrays for u-quark Compton form factors
    for(int i=0;i<49;i++) {
      
      //xbj array for gauss integral
      xx2[i] = dinter(zeta/2,1,48,i);
      
      //1/(X-zeta) array for u-quark filled
      erb3[i] = (fut(xval,erbl_eu,xx2[i])-fut(xval,erbl_eu,zeta))	\
	/(xx2[i]-zeta);
      
      //1/X array for u-quark filled
      erb4[i] = fut(xval,erbl_eu,xx2[i])/xx2[i];
    }

    for(int i=0;i<49;i++) {
      xxuplus[i] = dinter(zeta/2,1,48,i);
      erbuplus[i] = (fut(xval,erbl_euplus,xxuplus[i])\
			-fut(xval,erbl_euplus,zeta))/(xx2[i]-zeta);
      erbuplus2[i] = fut(xval,erbl_euplus,xxuplus[i])/xxuplus[i];
    }

    for(int i=0;i<49;i++) {
      xxdplus[i] = dinter(zeta/2,1,48,i);
      erbdplus[i] = (fut(xval,erbl_edplus,xxdplus[i])\
			-fut(xval,erbl_edplus,zeta))/(xx2[i]-zeta);
      erbdplus2[i] = fut(xval,erbl_edplus,xxdplus[i])/xxdplus[i];
    }
    
    
    //Fills integration arrays for d-quark Compton form factors
    for(int i=0;i<49;i++) {
      
      //xbj array for gauss integral
      xx2d[i] = dinter(zeta/2,1,48,i);
      
      //1/(X-zeta) array for d-quark filled
      erb3d[i] = (fut(xval,erbl_ed,xx2d[i])-fut(xval,erbl_ed,zeta))	\
	/(xx2d[i]-zeta);
      
      //1/X array for d-quark filled
      erb4d[i] = (fut(xval,erbl_ed,xx2d[i]))/xx2d[i];
    }
    
    //Calculates u-quark Re part of the compton form factor
    //1/2 factor is included as part of the symmetry
    double reeuv = dgaus1(48,zeta/2,1,erb3)+(gpdzeta)*log((1-zeta)	\
		  /(zeta/2))+dgaus1(48,zeta/2,1,erb4);
    
    //Calculates d-quark Re part of the compton form factor
    //1/2 factor is included as part of the symmetry
    double reedv = dgaus1(48,zeta/2,1,erb3d)+(gpdzetad)*log((1-zeta)\
		   /(zeta/2))+dgaus1(48,zeta/2,1,erb4d);

    //Calculate e sea quark cff contributions
    double reeus = -1*(dgaus1(48,zeta/2,1,erbuplus)+(gpdzeta)*log((1-zeta) \
			       /(zeta/2))+dgaus1(48,zeta/2,1,erbuplus2));

    double reeds = -1*(dgaus1(48,zeta/2,1,erbdplus)+(gpdzeta)*log((1-zeta) \
			       /(zeta/2))+dgaus1(48,zeta/2,1,erbdplus2));
    
    //Imaginary parts of the Compton form factors
    //u-quark
    double imeuv = imcff;
    
    //d-quark
    double imedv = imcffdv;
    
    //Creates the Real Compton form factors of the proton
    //Using charges for the u-quark and d-quark
    double reepv = (4./9.)*reeuv + (1./9.)*reedv;
    double imepv = (4./9.)*imeuv + (1./9.)*imedv;

    double reeps = (4./9.)*reeus + (1./9.)*reeds;
    double imeps = (4./9.)*imeuv + (1./9.)*imedv;

    double reep = reeps;
    double imep = imeps;

    //    cout << reepv << " " << imepv << endl;
    cout << "ReEp: " << reep << " ImEp: " << imep << endl;
    
    myfile5 << reepv << "," << imepv << ",";
        
    repv = reepv;
    impv = imepv;
    reuv = reeuv;
    imuv = imeuv;
    redv = reedv;
    imdv = imedv;
    rep = reep;
    imp = imep;
    reu = reeus;
    imu = imeuv;
    red = reeds;
    imd = imedv;
    reg = 0.;
    img = 0.;

    myfile2.close();
    
  }
    //************************
    //
    //    GPD Htil
    //
    //************************
  
  else if (jgpd ==3) {
    myfile3.open("GPD_Ht.dat");
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
	//For the helicity distributions the "plus" distribution is antisymmetric
	//and the "minus" distribution is symmetric
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
	erbl_hut[i] = Hutminus;
	erbl_hdt[i] = Hdtminus;
	
	//Writes GPDs to file, total quark distribution
	//separated between u and d
	myfile3 << xval[i] << " " << Hutminus << " " << Hdtminus << endl;
      }
      else {
	
	//DGLAP region
	interval = xval[i]-xval[i-increment];
	
	//put DGLAP GPD in file
	myfile3 << xval[i] << " " << uv0[i]/xval[i] << " " << dv0[i]/xval[i] << endl;	
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
      erb5[i] = (fut(xval,erbl_hut,xx5[i])-fut(xval,erbl_hut,zeta))/(xx5[i]-zeta);
      
      //1/X integral array filled
      erb6[i] = fut(xval,erbl_hut,xx5[i])/xx5[i];
    }
    
    //d-quark gauss integral arrays filled
    for(int i=0;i<49;i++) {
      
      //xbj array filled
      xx5d[i] = dinter(zeta/2,1,48,i);
      
      //1/(X-zeta) integral array filled for d-quark
      erb5d[i] = (fut(xval,erbl_hdt,xx5d[i])-fut(xval,erbl_hdt,zeta)) \
	/(xx5d[i]-zeta);
      
      //1/X integral array filled for d-quark
      erb6d[i] = (fut(xval,erbl_hdt,xx5d[i]))/xx5d[i];
    }
    
    //ReH_t for u-quark Compton form factor
    double rehut = (-dgaus1(48,zeta/2,1,erb5)+(gpdzeta)*log((1-zeta)	\
			       /(zeta/2))+dgaus1(48,zeta/2,1,erb6));
    
    //ReH_t for d-quark Compton form factor
    double rehdt = (-dgaus1(48,zeta/2,1,erb5d)+(gpdzetad)*log((1-zeta) \
			       /(zeta/2))+dgaus1(48,zeta/2,1,erb6d));
    
    //ImH_t u-quark Compton form factor
    double imhut = imcff;
    
    //ImH_t d-quark Compton form factor
    double imhdt = imcffdv;
    
    //ReH_t for the proton
    double rehpt = (4./9.)*rehut + (1./9.)*rehdt;
    
    //ImH_t for the proton
    double imhpt = (4./9.)*imhut + (1./9.)*imhdt;

    cout << "ReHtp: " << rehpt << " ImHtp: " << imhpt << endl;
    
    //Writes Compton form factors to File!
    myfile5 << rehpt << "," << imhpt << "," ;
    repv = rehpt;
    impv = imhpt;
    reuv = rehut;
    imuv = imhut;
    redv = rehdt;
    imdv = imhdt;
    rep = rehpt;
    imp = imhpt;
    reu = rehut;
    imu = imhut;
    red = rehdt;
    imd = imhdt;
    reg = 0.;
    img = 0.;
  }


  
  
  //***********************
  //
  //    GPD Etil
  //
  //***********************
  
  else if (jgpd == 4) {
    myfile4.open("GPD_Et.dat");
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
	myfile4 << xval[i] << " " << Eutminus << " " << Edtminus << endl;
      }
      else {
	interval = xval[i]-xval[i-increment];
	
	//Writes DGLAP region GPDs to file
	myfile4 << xval[i] << " " << uv0[i]/xval[i] << " " << dv0[i]/xval[i] << endl;
	
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
      erb7[i] = (fut(xval,erbl_eut,xx7[i])-fut(xval,erbl_eut,zeta))/(xx7[i]-zeta);
      
      //1/X array for gauss integral filled for u-quark
      erb8[i] = fut(xval,erbl_eut,xx7[i])/xx7[i];
    }
    
    //Gauss integral arrays filled for d-quark
    for(int i=0;i<49;i++) {
      
      //xbj array filled
      xx7d[i] = dinter(zeta/2,1,48,i);
      
      //1/(X-zeta) array for gauss integral filled for d-quark
      erb7d[i] = (fut(xval,erbl_edt,xx7d[i])-fut(xval,erbl_edt,zeta)) \
	/(xx7d[i]-zeta);
      
      //1/X array for gauss integral filled for d-quark
      erb8d[i] = (fut(xval,erbl_edt,xx7d[i]))/xx7d[i];
    }
    
    //Calculates ReE_t for u-quark using gauss integral routine
    double reeut = -dgaus1(48,zeta/2,1,erb7)+(gpdzeta)	\
      *log((1-zeta)/(zeta/2))+dgaus1(48,zeta/2,1,erb8);
    
    //Calculates ReE_t for d-quark using gauss integral routine
    double reedt = -dgaus1(48,zeta/2,1,erb7d)+(gpdzetad)	\
      *log((1-zeta)/(zeta/2))+dgaus1(48,zeta/2,1,erb8d);
    
    //Im E_t for u-quark distribution
    double imeut = imcff;
    
    //Im E_t for d-quark distribution
    double imedt = imcffdv;
    
    //ReE_t for the proton using u and d quark charges
    double reept = (4./9.)*reeut + (1./9.)*reedt;
    
    //ImE_t for the proton using u and d quark charges
    double imept = (4./9.)*imeut + (1./9.)*imedt;	\

    cout << "ReEtp: " << 0.5*reept << " ImEtp: " << imept << endl;
    
    myfile5 << reept << "," << imept << endl;

    repv = reept;
    impv = imept;
    reuv = reeut;
    imuv = imeut;
    redv = reedt;
    imdv = imedt;
    rep = reept;
    imp = imept;
    reu = reeut;
    imu = imeut;
    red = reedt;
    imd = imedt;
    reg = 0.;
    img = 0.;

    }
  }

//Madness is over!!!
