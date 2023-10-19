F90=gfortran
CXX=g++
# location of the pgplot directory which contains the libpgplot.a library.
#
#PGPLOTLIB=/lv9/users/cdc5m/pgplot
#PGOPTIONS=-L$(PGPLOTLIB) -lpgplot -lX11

# location of the cern libraries, eg. libpawlib.a, ...
#
#CERNLIB=/cern/pro/lib
#CERNLIB=/cern/2006/MacIntel_gcc4/lib/
#CERNOPTIONS=-L$(CERNLIB) -lpdflib804 -lpacklib-shift -lmathlib -lkernlib-shift -lgeant 
#CERNOPTIONS=-L$(CERNLIB) -lpdflib804 -lpacklib -lmathlib -lkernlib -lgeant 

#fflags = -g
# list all the object files which make up the program.  The name "components"
# can be anything.
#*
components=  ga.o gp.o struct_ev_gpd_tmd.o diquark_sub_GGLA_v2.o kelly.o alpha_sub.o main.o  diquark.o gpd.o gluon_gpd.o evol_erbl.o formfactor.o dg.o gauss.o fut.o sea_gpd.o gluon_ff.o sea_ff.o evol_erblg.o evol_erbl2.o cff_calc.o evol_erbl_new.o

#bh_comp.o nucl_bh_gangof8.o fvscal.o fvsum.o nucl_cff.o diquark_BoerM_GGLA.o

#components= breit.o double_d.o ~/simolib/dg.o  
# this is built-in "rule" for the make utility to build ".o" files from ".f" 
# files
#
.f.o:
	$(F90) -c -g -Wall $<

.cc.o:
	$(CXX) -c -g -Wall -O3 -g $<
# this is our "rule" which we have called "all" which depends upon the files in
# $(components).  The action that it takes (listed in the next line) uses 
# those files in $(components) which are out of date.  The "-o inclu" tells f77
# to name the executable file "inclu".
#
all: $(components)
	$(F90) -g $? -o bod -lc -lc++ $(CERNOPTIONS)

clean:
	-rm -rf *.o
	-rm bod
