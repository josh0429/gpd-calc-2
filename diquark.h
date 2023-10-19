#ifndef _DIQUARK3_H
#define _DIQUARK3_H

extern "C" {
  void diquark_sub_ggla_temp_fortran_(double *x, double *q2, double *zeta, double *t,
				      double *hu, double *hu_plus, double *hu_minus,
				      double *eu, double *eu_plus, double *eu_minus,
				      double *hutil, double *hutil_plus, double *hutil_minus,
				      double *eutil, double *eutil_plus, double *eutil_minus,
				      double *hd, double *hd_plus, double *hd_minus,
				      double *ed, double *ed_plus, double *ed_minus,
				      double *hdtil, double *hdtil_plus, double *hdtil_minus,
				      double *edtil, double *edtil_plus, double *edtil_minus,
				      double *hg, double *eg, double *htg, double *hub, double *hdb);

  void diquark_sub_ggla_(double *x, double *q2, double *zeta, double *t,
			 double *hu, double *hu_plus,
			 double *eu, double *eu_plus,
			 double *hutil, double *hutil_plus,
			 double *eutil, double *eutil_plus,
			 double *hd, double *hd_plus,
			 double *ed, double *ed_plus,
			 double *hdtil, double *hdtil_plus,
			 double *edtil, double *edtil_plus,
			 double *hg, double *eg, double *htg, double *hub, double* hdb,
			 double *hsb, double *hcb, double *init_q02, double *etg, double *eub, double *edb);
};

#endif
