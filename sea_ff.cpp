#include <iostream>
#include <cmath>
#include "sea_ff.h"

using namespace std;

double a20_upd_dip(double t, int ierr) {

  double alpha = 0.5768;
  double m = 1.6433;

  double alpha_err = 0.0073;
  double m_err = 0.0503;
  
  double a20;

  if(ierr == 1){
    a20 = alpha / pow(1+t/(m*m),2);
  }

  else if(ierr == 0) {
    a20 = (alpha+alpha_err) / pow(1+t/((m+m_err)*(m+m_err)),2);
  }

  return a20;
  
}


double a20_umd_dip(double t, int ierr) {

  /*  double alpha = 0.2204;
  double m = 1.9657;

  double alpha_err = 0.0074;
  double m_err = 0.1781;
  
  double a20;

  if(ierr == 1){
    a20 = alpha / pow(1+t/(m*m),2);
  }

  else if(ierr == 0) {
    a20 = (alpha+alpha_err) / pow(1+t/((m+m_err)*(m+m_err)),2);
  }
  */

  double M = 0.9383;
  double a20t0 = 0.213;
  double b20t0 = 0.279;
  double gmt0 = a20t0 + b20t0;
  double get0 = a20t0;
  double tau = t/(4*M*M);
  double hmumd = gmt0/((1 + t/1.371)*(1 + t/1.371));
  double heumd = get0/((1 + t/0.704)*(1 + t/0.704));
  double a20 = (tau*hmumd + heumd) / (1+tau);
  return a20;

}

double a20_u_dip(double t, int ierr) {

  double a20_upd = a20_upd_dip(t,1);
  double a20_umd = a20_umd_dip(t,1);
  double a20_upd_err = a20_upd_dip(t,0)-a20_upd;
  double a20_umd_err = a20_umd_dip(t,0)-a20_umd;

  double a20_q_err = 0.5*sqrt(a20_upd_err*a20_upd_err + a20_umd_err*a20_umd_err);

  double a20_u_c = 0.5*(a20_upd_dip(t,1) + a20_umd_dip(t,1));

  double a20u;

  if( ierr == 1) {
    a20u = a20_u_c;
  }
  else if (ierr == 0) {
    a20u = a20_u_c + a20_q_err;
  }

  return a20u;

}

double a20_d_dip(double t, int ierr) {

  double a20_upd = a20_upd_dip(t,1);
  double a20_umd = a20_umd_dip(t,1);
  double a20_upd_err = a20_upd_dip(t,0)-a20_upd;
  double a20_umd_err = a20_umd_dip(t,0)-a20_umd;

  double a20_q_err = 0.5*sqrt(a20_upd_err*a20_upd_err + a20_umd_err*a20_umd_err);

  double a20_d_c = 0.5*(a20_upd_dip(t,1) - a20_umd_dip(t,1));

  double a20d;

  if( ierr == 1) {
    a20d = a20_d_c;
  }
  else if (ierr == 0) {
    a20d = a20_d_c + a20_q_err;
  }

  return a20d;

}

double a20_uv_dip(double t, int ierr) {
  
  double alpha = 0.297;
  double m = 1.187;

  double alpha_err = 0.0001;
  double m_err = 0.0003;

  double a20;

  if(ierr == 1){
    a20 = alpha / pow(1+t/(m*m),2);
  }

  else if(ierr == 0) {
    a20 = (alpha+alpha_err) / pow(1+t/((m+m_err)*(m+m_err)),2);
  }

  return a20;

}

double a20_dv_dip(double t, int ierr) {

  double alpha = 0.1101;
  double m = 1.0704;

  double alpha_err = 0.0002;
  double m_err = 0.0060;

  double a20;

  if(ierr == 1){
    a20 = alpha / pow(1+t/(m*m),2);
  }

  else if(ierr == 0) {
    a20 = (alpha+alpha_err) / pow(1+t/((m+m_err)*(m+m_err)),2);
  }

  return a20;

}

double a20_ub_dip(double t, int ierr) {

  double a20_u = a20_u_dip(t,1);
  double a20_uv = a20_uv_dip(t,1);
  double a20_u_err = a20_u_dip(t,0)-a20_u;
  double a20_uv_err = a20_uv_dip(t,0)-a20_uv;

  double a20_qb_err = sqrt(a20_u_err*a20_u_err + a20_uv_err*a20_uv_err);

  double a20_ub_c = a20_u_dip(t,1) - a20_uv_dip(t,1);

  double a20ub;
  
  if( ierr == 1) {
    a20ub = a20_ub_c;
  }
  else if (ierr == 0) {
    a20ub = a20_ub_c + a20_qb_err;
  }
  return a20ub;

}


double a20_db_dip(double t, int ierr) {

  double a20_d = a20_d_dip(t,1);
  double a20_dv = a20_dv_dip(t,1);
  double a20_d_err = a20_d_dip(t,0)-a20_d;
  double a20_dv_err = a20_dv_dip(t,0)-a20_dv;

  double a20_qb_err = sqrt(a20_d_err*a20_d_err + a20_dv_err*a20_dv_err);

  double a20_db_c = a20_d_dip(t,1) - a20_dv_dip(t,1);

  double a20db;

  if( ierr == 1) {
    a20db = a20_db_c;
  }
  else if (ierr == 0) {
    a20db = a20_db_c + a20_qb_err;
  }
  return a20db;

}

double c20_umd_zexp(double t){

  double mpi = 0.135;
  double tc = 4*mpi*mpi;
  double z = (sqrt(tc+t)-sqrt(tc))/(sqrt(tc+t)+sqrt(tc));

  double a = -0.5517;
  double b = 2.6953;
  double c = -3.3170;

  double c20 = a + b*z + c*z*z;

  return c20;
  
}

double c20_upd_zexp(double t){

  double mpi =	0.135;
  double tc = 4*mpi*mpi;
  double z = (sqrt(tc+t)-sqrt(tc))/(sqrt(tc+t)+sqrt(tc));

  double a = -1.0398;
  double b = 4.1087;
  double c = -4.4119;

  double c20 = a + b*z + c*z*z;

  return c20;

}

double a30_umd_zexp(double t){

  double mpi =	 0.135;
  double tc = 4*mpi*mpi;
  double z = (sqrt(tc+t)-sqrt(tc))/(sqrt(tc+t)+sqrt(tc));

  double a = 0.0907;
  double b = -0.1229;
  double c = 0.1599;

  double a30 = a + b*z + c*z*z;

  return a30;

}

double a30_upd_zexp(double t){

  double mpi =   0.135;
  double tc = 4*mpi*mpi;
  double z = (sqrt(tc+t)-sqrt(tc))/(sqrt(tc+t)+sqrt(tc));

  double a = 0.1488;
  double b = 0.0069;
  double c = -0.2353;

  double a30 = a + b*z + c*z*z;

  return a30;

}


/*
int main(){

  for(double t =0; t<= 2.01; t+=0.01){

    double ubff = a20_ub_dip(t,1);
    double ubfferr = a20_ub_dip(t,0) - a20_ub_dip(t,1);

    double dbff = a20_db_dip(t,1);
    double dbfferr = a20_db_dip(t,0) - a20_db_dip(t,1);

    double upd = a20_upd_dip(t,1);
    double upderr = a20_upd_dip(t,0) - upd;

    double uv = a20_uv_dip(t,1);
    double uverr = a20_uv_dip(t,0) - uv;
    double dv = a20_dv_dip(t,1);
    double dverr = a20_dv_dip(t,0)-dv;

    //    double qb = 0.5*(upd - uv - dv);
    //    double qberr = 0.5*sqrt(upderr*upderr + uverr*uverr + dverr*dverr);
    cout << t << " " << ubff << " " << dbff << endl;
    //    cout << t << " " << a20_umd_dip(t,1) << " " << a20_uv_dip(t,1)-a20_dv_dip(t,1) << endl;
    //cout << t << " " << a20_umd_dip(t,1) << endl;
    //    cout << t << " " << ubff << " " << ubfferr << " " << dbff << " " << dbfferr << endl;
  }
}

*/
