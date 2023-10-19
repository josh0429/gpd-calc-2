C
C-----------------------------------------------------------------------
C  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  + 
C-----------------------------------------------------------------------
C
C File - Version - Date: diquark_sub_GGL_v2.f - v 1.0 - 4/14/2014
C                   
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C    Calculates GPDs IN THE SPECTATOR -- OR DIQUARK -- MODEL 
C
C    NEW QUANTITATIVE PARAMETRIZATION II
C    PROGRAM ALSO COMPARES WITH OTHER CURRENT CALCULATIONS IN THE LITERATURE
C
C INPUT VARIABLES
C          x       initial quark longitudinal momentum fraction 
C                  (this is X in the literature, different from the symmetric system x)
C          q2      DIS momnetum transfer squared (GeV^2)
C          t       negative momentum transfer squared (GeV^2)
C          zeta    skewness/longitudinal momentum transfer fraction
C
C OUTPUT VARIABLES
***         hu = H_u in DGLAP region
***         eu = E_u in DGLAP region
***         hutil= H_u tilde in DGLAP region
***         eutil = E_u tildein DGLAP region
***
****         
***         hu_plus = H_u upper curve for error band in DGLAP region
***         eu_plus = E_u upper curve for error band in DGLAP region
***         hutil_plus= H_u tilde upper curve for error band in DGLAP region
***         eutil_plus = E_u tilde upper curve for error band in DGLAP region
****
***         hu_minus = H_u lower curve for error band in DGLAP region
***         eu_minus = E_u lower curve for error band in DGLAP region
***         hutil_minus= H_u tilde lower curve for error band in DGLAP region
***         eutil_minus = E_u tilde lower curve for error band in DGLAP region
********
*** SAME SCHEME REPEATED FOR D QUARK
C
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Author: Simonetta Liuti
C
C-----------------------------------------------------------------------
C

      subroutine diquark_sub_GGLA_temp_fortran(x,q2,zeta,t,
     > hu,hu_plus,hu_minus,
     > eu,eu_plus,eu_minus,
     > hutil,hutil_plus,hutil_minus,
     > eutil,eutil_plus,eutil_minus,
     > hd,hd_plus,hd_minus,
     > ed,ed_plus,ed_minus,
     > hdtil,hdtil_plus,hdtil_minus,
     > edtil,edtil_plus,edtil_minus)

      implicit none
      integer*4 n,ik
      parameter(n=8)

      real*8 kperp(n),kmax,kappau,kappad
     
      real*8 amp,amup,pi,conv

      real*8 dgaus1

      real*8 amqu,amqd,amdiq_uh,amdiq_dh,amdiq_ue,amdiq_de

      real*8 amqutil,amqdtil,amdiq_uht,amdiq_dht,amdiq_uet,amdiq_det

      real*8 alam_uh,alam_dh,alam_ue,alam_de

      real*8 alam_uht,alam_dht,alam_uet,alam_det

      real*8 alf_uh,alf_dh,alf_ue,alf_de

      real*8 ah_util,ah_dtil

      real*8 anormu,anormd,anormue,anormde

      real*8 autil,adtil,auetil,adetil

      real*8 bet_uh,bet_dh,bet_ue,bet_de,p1u,p1d,p2u,p2d

      real*8 bet_uhtil,bet_dhtil,bet_uetil,bet_detil,
     > p1_util,p1_dtil,p2_util,p2_dtil

      real*8 dbet_uh,dbet_dh,dbet_ue,dbet_de,dp1u,dp1d,dp2u,dp2d
 
      real*8 dbet_uhtil,dbet_dhtil,dbet_uetil,dbet_detil,
     > dp1_util,dp1_dtil,dp2_util,dp2_dtil

      real*8 x,q2,zeta,t0,xi,t,del,del2,xtil,xtil1

      real*8 emtil0,emtil0t,emtil,emtilt

      real*8 aaa,aaat,bbb,anumer1,anumer2,anume1,anume2,antile1,antile2,
     > antil1,antil2,akntil1,akntil2,akntile1,akntile2
      
      real*8 akfunc1,akfunc2,akfune1,akfune2

      real*8 denom,denomt,fintu(n),finteu(n),fitilu(n),fitileu(n),
     > fintd(n),finted(n),fitild(n),fitiled(n)

      real*8 gu,geu,gutil,geutil,r1u,r2u,r1util,r2util

      real*8 gd,ged,gdtil,gedtil,r1d,r2d,r1dtil,r2dtil

      real*8 er_bet_uh,er_bet_ue,er_bet_uhtil,er_bet_uetil,
     > er_p_uh,er_p_ue,er_p_uhtil,er_p_uetil

      real*8 er_bet_dh,er_bet_de,er_bet_dhtil,er_bet_detil,
     > er_p_dh,er_p_de,er_p_dhtil,er_p_detil

      real*8 dalp_uh,dalp_ue,dalp_uhtil,danormu,danormue,dautil,dauetil,
     >  er_alp_uh,er_alp_ue,er_alp_uhtil,er_alp_uetil

      real*8 dalp_dh,dalp_de,dalp_dhtil,danormd,danormde,dadtil,dadetil,
     >  er_alp_dh,er_alp_de,er_alp_dhtil,er_alp_detil

      real*8 er_uh,er_ue,er_uhtil,er_uetil,er_dh,er_de,er_dhtil,er_detil

      real*8 hu,eu,hd,ed,hutil,eutil,hdtil,edtil

      real*8 hu_plus,hd_plus,eu_plus,ed_plus,hutil_plus,hdtil_plus,
     > eutil_plus,edtil_plus

      real*8 hu_minus,hd_minus,eu_minus,ed_minus,hutil_minus,
     > hdtil_minus,eutil_minus,edtil_minus

** ALL masses in GeV 
      data amp/0.938d0/  ! PROTON MASS
      data amup/2.79d0/  ! PROTON ANOMALOUS MAGNETIC MOMENTS

*** anomalous magnetic moments
      data kappau/1.67/  ! U QUARK ANOMALOUS MAGNETIC MOMENT
      data kappad/-2.03/ ! D QUARK ANOMALOUS MAGNETIC MOMENT

*** regge xi dependent parameter
c      data axi/-0.2/, axit/0.8/
**********************************************************************
*** DGLAP REGION PARAMETERS 
**********************************************************************
***  H_u = G(amqu,amdiq_uh,alf_uh,anormu) R(bet_uh,p1u) , 
***  E_u = G(amqu,amdiq_ue,alf_ue,anormue) R(bet_ue,p2u) 
***  H_u tilde = G(amqu,amdiq_uht,alf_uht,anormuht) R(bet_uht,p1u) 
***  E_u tilde = G(amqu,amdiq_uet,alf_uet,anormuet) R(bet_uht,p2u) 

*** amqu = quark mass
*** amdiq_uh = diquark mass
*** alam_uh = diquark dipole form factor mass 
*** alf_uh = Regge t-independent exponent 
*** anormu = normalization constant 
*** bet_uh = Regge-t-dependent parameter
*** p1u = Regge-t dependent parameter

CCC NOTE THAT: 
CC  1) amdiq_uh=amdiq_ue, alam_uh= alam_ue, alf_uh= alf_ue
CC  2) amdiq_uht=amdiq_uet, alam_uht= alam_uet, alf_uht= alf_uet

*** dbet_uh = error on paramater bet_uh from form factor fit 
*** dp1u = error on parameter p1u from form factor fit 

*** same scheme for H_d, E_d, Htilde_d, Etilde_d
 
***********
*********** H,E *** Osvaldo PRD
***********
*     QUARK MASSES  
      data amqu/0.420d0/
      data amqd/0.275d0/

*     DIQUARK MASSES
      data amdiq_uh/0.603d0/
      data amdiq_dh/0.912d0/
      data amdiq_ue/0.603d0/
      data amdiq_de/0.912d0/

*     DIQUARK FORM FACTOR MASSES
      data alam_uh/1.018d0/
      data alam_dh/0.860D0/
      data alam_ue/1.018d0/
      data alam_de/0.860d0/

*     REGGE PARAMETERS t=0
      data alf_uh/0.210d0/
      data alf_dh/0.0317d0/
      data alf_ue/0.210d0/
      data alf_de/0.0317d0/

*     NORMALIZATIONS
      data anormu/2.043d0/
      data anormd/1.570d0/
      data anormue/1.803d0/
      data anormde/-2.780d0/

******************************************************************************************************
******************************************************************************************************
******************************************************************************************************
*     t-DEPENDENT PARAMETERS: beta (alpha prime in text), p and normalization
**
**
******************************************************************************************************
*** Osvaldo PRD
c      data bet_uh/2.448d0/
c      data dbet_uh/.0845/
* 
c      data bet_dh/2.209d0/
c      data dbet_dh/0.145/
*
c      data bet_ue/2.811d0/
c      data dbet_ue/0.764d0/
*
c      data bet_de/1.362d0/
c      data dbet_de/0.584/
*
c      data p1u/0.620d0/
c      data dp1u/0.0893/
*
c      data p1d/0.658d0/
c      data dp1d/0.370/
*
c      data p2u/0.863d0/
c      data dp2u/0.482/
*
c      data p2d/1.115d0/
c      data dp2u/1.148/

***
*** After Cates DATA 
      data bet_uh/1.814d0/
      data dbet_uh/.0220/
* 
      data bet_dh/1.139d0/
      data dbet_dh/0.0564/
*
      data bet_ue/2.835d0/
      data dbet_ue/0.0509d0/
*
      data bet_de/1.281d0/
      data dbet_de/0.0310/
*
      data p1u/0.449d0/
      data dp1u/0.0170/
*
      data p1d/-0.113d0/
      data dp1d/0.104/
*
      data p2u/0.969d0/
      data dp2u/0.0307/
*
      data p2d/0.726d0/
      data dp2d/0.0631/


***********************************************************
***********
*********** Htilde, Etilde
***********

*     QUARK MASSES
      data amqutil/2.624d0/
      data amqdtil/2.602d0/

*     DIQUARK MASSES
      data amdiq_uht/0.474d0/
      data amdiq_dht/0.703d0/

*     DIQUARK FORM FACTOR MASSES
      data alam_uht/0.971d0/
      data alam_dht/0.878D0/      

*     REGGE PARAMETERS t=0
      data ah_util/0.219d0/
      data ah_dtil/0.0347d0/

*     NORMALIZATIONS THESE ARE DIFFERENT FOR H AND E
      data autil/0.0504d0/
      data adtil/-0.0262d0/
      data auetil/1.074d0/
      data adetil/-0.966d0/


*     t-DEPENDENT PARAMETERS
      data bet_uhtil/1.543d0/
      data dbet_uhtil/0.296/

      data bet_dhtil/1.297d0/
      data dbet_dhtil/0.245/

      data bet_uetil/5.130d0/
      data dbet_uetil/0.101/

      data bet_detil/3.385d0/
      data dbet_detil/0.145d0/

*** Osvaldo PRD
      data p1_util/0.346d0/
      data dp1_util/0.248/

      data p1_dtil/0.974d0/
      data dp1_dtil /0.357/
*** 
      data p2_util/3.507d0/
      data dp2_util/0.054/

      data p2_dtil/2.326d0/
      data dp2_dtil /0.137/

**********
**********
      data kmax/5.0d0/ !kperp max in GeV

*****************************************************************************
*****************************************************************************
*****************************************************************************
*****************************************************************************

      pi = acos(-1.)
      conv=197.3/1000.


*******
*******
******* GPDs in DGLAP region
**
** X > \zeta
**
      xi = zeta/(2.-zeta)

**** asymptotic
      t0 = - amp**2*zeta**2/(1.-zeta)

      del2 = (-t + t0)*(1.-zeta)        ! Delta_T^2 

      del = sqrt (del2)      ! \mid Delta_T \mid

      xtil  = (x-zeta)/(1.-zeta)

      xtil1 = (x-zeta)/(1.-x)

***********
**** GO! 
***********

**************************************************************************************
**************************************************************************************
** u-quark
**************************************************************************************
**************************************************************************************
** k_T independent part of denominator for H_u and E_u   
**********
      emtil0 = x*amp**2 - alam_uh**2 - amdiq_uh**2*x/(1.-x)  ! zeta=0
      emtil = xtil*amp**2 - alam_uh**2 - xtil1*amdiq_uh**2   ! zeta neq 0

**  k_T independent part of denominator for H_u tilde and E_u tilde      
************
      emtil0t = x*amp**2 - alam_uht**2 - amdiq_uht**2*x/(1.-x)  ! zeta=0
      emtilt = xtil*amp**2 - alam_uht**2 - xtil1*amdiq_uht**2   ! zeta neq 0

** Integration in k_perpT = kperp
** (integration in phi angular variable for k_perp prime, done analytically following
**  Gradshtein, 6.11)

      CALL DINTER(kperp,0.d0,kmax,N)
      do 1850 ik=1,n
*** H_u and E_u
       aaa = emtil - kperp(ik)**2*(1.-zeta)/(1.-x)
     >                          - del2*(1.-x)/(1.-zeta)

*** H_u tilde and E_u tilde
       aaat = emtilt - kperp(ik)**2*(1.-zeta)/(1.-x)
     >                          - del2*(1.-x)/(1.-zeta)

       bbb = del*kperp(ik)*2.d0

**** Numerator H_u  
           anumer1 = (amqu+amp*x)*(amqu+amp*xtil) + kperp(ik)*kperp(ik)
           anumer2 = -kperp(ik)*(1.-x)/(1.-zeta)*del

**** Numerator H_u tilde  
        antil1 = (amqutil+amp*x)*(amqutil+amp*xtil)- kperp(ik)*kperp(ik)
        antil2 = kperp(ik)*(1.-x)/(1.-zeta)*del

**** Numerator E_u 
           anume1 =  2.*amp*(-(amqu+amp*x)+(amqu+amp*xtil))*kperp(ik) 
           anume2 = -2.*amp*(amqu+amp*x)*(1.-x)/(1.-zeta)

**** Numerator E_u tilde  
        antile1 = -4.*amp*((amqutil+amp*x)+(amqutil+amp*xtil))*kperp(ik) 
        antile2 = -4.*amp*(amqutil+amp*x)*(1.-x)/(1.-zeta)

          if(aaa.ge.abs(bbb).or.aaa.le.-abs(bbb)) then
*** H
         akfunc1 =  anumer1*2.*pi/(aaa**2-bbb**2)**1.5*abs(aaa)
         akfunc2 =  anumer2*2.*pi/(aaa**2-bbb**2)**1.5*(bbb)

*** E
         akfune1 =  anume1*2.*pi/(aaa**2-bbb**2)**1.5*kperp(ik)*2.d0
         akfune2 =  anume2*2.*pi/(aaa**2-bbb**2)**1.5*abs(aaa)

*** Htilde
         akntil1 =  antil1*2.*pi/(aaat**2-bbb**2)**1.5*abs(aaat)
         akntil2 =  antil2*2.*pi/(aaat**2-bbb**2)**1.5*(bbb)

*** Etilde
         akntile1 =  antile1*2.*pi/(aaat**2-bbb**2)**1.5*kperp(ik)*2.d0
         akntile2 =  antile2*2.*pi/(aaat**2-bbb**2)**1.5*aaat

      else
         endif
       
*****
***** GPD integrand

         denom = (emtil0 - kperp(ik)**2/(1.-x) )**2  
         denomt = (emtil0t - kperp(ik)**2/(1.-x) )**2  

** H_u
           fintu(ik) = kperp(ik)* (akfunc1 + 0.75*akfunc2)/denom  

** E_u
           finteu(ik) = kperp(ik)* (akfune1 + akfune2)/denom

** H_u tilde
           fitilu(ik) = kperp(ik)* (akntil1 + akntil2)/denomt

** E_u tilde
           fitileu(ik) = kperp(ik)* (akntile1 + akntile2)/denomt

 1850   enddo


*** GPD -- G function diquark contribution

           gu = dgaus1(n,0.d0,kmax,fintu)
           geu = dgaus1(n,0.d0,kmax,finteu)
           gutil = dgaus1(n,0.d0,kmax,fitilu)
           geutil = dgaus1(n,0.d0,kmax,fitileu)

*** Regge term and normalization constant

** H_u 
           r1u = anormu*x**(-alf_uh)*x**(bet_uh*(1.-x)**p1u*del2)

** E_u 
           r2u =  anormue*x**(-alf_ue)*x**(bet_ue*(1.-x)**p2u*del2)

** H_u tilde
        r1util = autil*x**(-ah_util)*x**(bet_uhtil*(1.-x)**p1_util*del2)

** E_u tilde
        r2util =auetil*x**(-ah_util)*x**(bet_uetil*(1.-x)**p2_util*del2)


**** 
**** CALCULATE ERROR on Regge term from parameter error
****
******************************************************************
           er_bet_uh = r1u * dlog(x)*(1.-x)**p1u*del2  
           er_bet_ue = r2u * dlog(x)*(1.-x)**p2u*del2  
           er_bet_uhtil = r1util * dlog(x)*(1.-x)**p1_util*del2
           er_bet_uetil = r2util * dlog(x)*(1.-x)**p2_util*del2

           er_p_uh = r1u * dlog(x)* dlog(1.-x) * bet_uh*(1.-x)**p1u*del2  
           er_p_ue = r2u * dlog(x)* dlog(1.-x) * bet_ue*(1.-x)**p2u*del2
           er_p_uhtil = r1util*
     > dlog(x) * dlog(1.-x) * bet_uhtil*(1.-x)**p1_util*del2
           er_p_uetil = r2util*
     > dlog(x) * dlog(1.-x) * bet_uetil*(1.-x)**p2_util*del2
  
*** 12/2013 add in error on x dependent variables alpha and normalization (this is not a real hessian error because we obtained it by minimizing to valence parametrizations)
*** estimated error on alpha parameter
           dalp_uh = 0.07*alf_uh
           dalp_ue = 0.07*alf_ue
           dalp_uhtil = 0.16*ah_util
*** estimated error on normalization
           danormu = 0. !0.10*anormu
           danormue = 0. !0.11*anormue
           dautil = 0. !0.16*autil
           dauetil = 0. !0.16*auetil

           er_alp_uh = r1u * dlog(x)
           er_alp_ue = r2u * dlog(x)
           er_alp_uhtil = r1util * dlog(x)
           er_alp_uetil = r2util * dlog(x)


*** for Compton form factor the errors enter times 2 because the function is squared
           er_uh = sqrt((er_bet_uh*dbet_uh)**2 + (er_p_uh*dp1u)**2 + 
     > (er_alp_uh*dalp_uh)**2 + danormu**2)
*
           er_ue = sqrt((er_bet_ue*dbet_ue)**2 + (er_p_ue*dp2u)**2
     > + (er_alp_ue*dalp_ue)**2 + danormue**2)
*
           er_uhtil = sqrt((er_bet_uhtil*dbet_uhtil)**2 + 
     > (er_p_uh*dp1_util)**2
     > +                    (er_alp_uhtil*dalp_uhtil)**2+ dautil**2)
*
           er_uetil = sqrt((er_bet_uetil*dbet_uetil)**2  
     > +    (er_p_ue*dp2_util)**2
     > +                    (er_alp_uetil*dalp_uhtil)**2+ dauetil**2)
           
************************************************************************************
*** H_u
           
           hu = gu/(1.-x)*r1u
 
           hu_plus = gu/(1.-x)  * (r1u+er_uh)
           hu_minus = gu/(1.-x)  * (r1u-er_uh)

*** E_u
           geu = geu /(1.-x)*(1.-zeta)
           eu =  -geu * r2u 

           eu_plus = -geu * (r2u+er_ue)
           eu_minus = -geu * (r2u-er_ue)

*** Htil_u
           hutil = gutil /(1.-x) *r1util

           hutil_plus = gutil /(1.-x) * (r1util+er_uhtil)
           hutil_minus = gutil /(1.-x) * (r1util-er_uhtil)


*** Etil_u (zeta = 0)
           if(zeta.eq.0.and.t.eq.0.)  then

           eutil = amp*pi*(1.-x)**6*auetil*x**(-ah_util)*
     >  (amp/3./(emtil0t*(1.-x))**3 - 
     >  4.*(amp*x+amqutil)/5./(emtil0t*(1.-x))**4*( (1.-2.*x)*amp**2-
     >   amdiq_uht**2 + alam_uht**2) )  

              else
*** Etil_u (zeta neq 0)  
           eutil = geutil /(1.-x) *r2util *(1.-zeta)/zeta

           eutil_plus = geutil /(1.-x) * (r2util+er_uetil)
           eutil_minus = geutil /(1.-x) * (r2util+er_uetil)
           endif

**********************************************************************************
**********************************************************************************
** d-quark
**********************************************************************************
**********************************************************************************
      

** zeta =0 (and first denominator for \zeta \neq 0)
      emtil0 = x*amp**2 - alam_dh**2 - amdiq_dh**2*x/(1.-x)
      emtil0t = x*amp**2 - alam_dht**2 - amdiq_dht**2*x/(1.-x)

** zeta neq 0
      xtil  = (x-zeta)/(1.-zeta)
      xtil1 = (x-zeta)/(1.-x)
      emtil = xtil*amp**2- alam_dh**2- 
     > xtil1*amdiq_dh**2
      emtilt = xtil*amp**2- alam_dht**2- 
     > xtil1*amdiq_dht**2

* integration in k_perp
      CALL DINTER(kperp,0.d0,kmax,N)
      do 2850 ik=1,n

**      
** integration in phi (angular variable for k_perp prime, done analytically following
**  Gradshtein, 6.11)
*** modify => insert analytical integrations

       aaa = emtil - kperp(ik)**2*(1.-zeta)/(1.-x)
     >                          - del2*(1.-x)/(1.-zeta)
       aaat = emtilt - kperp(ik)**2*(1.-zeta)/(1.-x)
     >                          - del2*(1.-x)/(1.-zeta)

       bbb = del*kperp(ik)*2.d0

**** OPTION GGLA (NEW)
         anumer1 = (amqd+amp*x)*(amqd+amp*xtil) + kperp(ik)*kperp(ik)
         anumer2 = -kperp(ik)*(1.-x)/(1.-zeta)*del

         antil1 = (amqdtil+amp*x)*(amqdtil+amp*xtil)-kperp(ik)*kperp(ik)
         antil2 = kperp(ik)*(1.-x)/(1.-zeta)*del

           anume1 =  2.*amp*
     >              (-(amqd+amp*x)+(amqd+amp*xtil))*kperp(ik) !/del
           antile1 =  -4.*amp*
     >              ((amqdtil+amp*x)+(amqdtil+amp*xtil))*kperp(ik) !/del


           anume2 =  -2.*amp*(amqd+amp*x)*(1.-x)/(1.-zeta)
           antile2 =  -2.*amp*(amqdtil+amp*x)*(1.-x)/(1.-zeta)

          if(aaa.ge.abs(bbb).or.aaa.le.-abs(bbb)) then
*** H
         akfunc1 =  anumer1*2.*pi/(aaa**2-bbb**2)**1.5*abs(aaa)
         akfunc2 =  anumer2*2.*pi/(aaa**2-bbb**2)**1.5*(bbb)
*** E
         akfune1 =  anume1*2.*pi/(aaa**2-bbb**2)**1.5*kperp(ik)*2.d0
         akfune2 =  anume2*2.*pi/(aaa**2-bbb**2)**1.5*abs(aaa)
*** Htilde
         akntil1 =  antil1*2.*pi/(aaat**2-bbb**2)**1.5*abs(aaat)
         akntil2 =  antil2*2.*pi/(aaat**2-bbb**2)**1.5*(bbb)

*** Etilde
         akntile1 =  antile1*2.*pi/(aaat**2-bbb**2)**1.5*kperp(ik)*2.d0
         akntile2 =  antile2*2.*pi/(aaat**2-bbb**2)**1.5*aaat
      else
 
c        akfunc1 =  anumer1*2.*pi/(-aaa**2+bbb**2)**1.5*abs(aaa)
c        akfunc2 =  anumer2*2.*pi/(-aaa**2+bbb**2)**1.5*(bbb)
c
c         akfune1 =  anume1*2.*pi/(-aaa**2+bbb**2)**1.5*kperp(ik)*2.
c         akfune2 =  anume2*2.*pi/(-aaa**2+bbb**2)**1.5*abs(aaa)
c
c        akntil1 =  antil1*2.*pi/(-aaa**2+bbb**2)**1.5*abs(aaa)
c        akntil2 =  antil2*2.*pi/(-aaa**2+bbb**2)**1.5*(bbb)
c
c         akntile1 =  antile1*2.*pi/(-aaa**2+bbb**2)**1.5*kperp(ik)*2.
c         akntile2 =  -anume2*2.*pi/(-aaa**2+bbb**2)**1.5*abs(aaa)


         endif

       
*****
***** GPD
**** Delta dependence, phase space x-dependent factor
         denom = (emtil0 - kperp(ik)**2/(1.-x) )**2
         denomt = (emtil0t - kperp(ik)**2/(1.-x) )**2  

         fintd(ik) = kperp(ik)* (akfunc1 + akfunc2)   
     >  /denom  

           finted(ik) = kperp(ik)* (akfune1 + akfune2)   
     >  /denom

           fitild(ik) = kperp(ik)* (akntil1 + akntil2)   
     >  /denomt

           fitiled(ik) = kperp(ik)* (akntile1 + akntile2)   
     >  /denomt 


 2850   enddo


*** GPD -- diquark contribution

           gd = dgaus1(n,0.d0,kmax,fintd)
           ged = dgaus1(n,0.d0,kmax,finted)

           gdtil = dgaus1(n,0.d0,kmax,fitild)
           gedtil = dgaus1(n,0.d0,kmax,fitiled)

*** Regge term

*** phase space x-dependent factor 
*** (commented out x-zeta dependent phase space factor)
*** back in 4.1.08

           r1d = anormd*x**(-alf_dh)
     > *x**(bet_dh*(1.-x)**p1d*del2)

           r2d = anormde*x**(-alf_dh) !/1.378
     > *x**(bet_de*(1.-x)**p2d*del2)

           r1dtil = adtil*x**(-ah_dtil)
     > *x**(bet_dhtil*(1.-x)**p1_dtil*del2)

           r2dtil = adetil*x**(-ah_dtil)
     > *x**(bet_detil*(1.-x)**p2_dtil*del2)

**** 
**** CALCULATE ERROR on Regge term from parameter error
****
******************************************************************
           er_bet_dh = r1d * dlog(x)*(1.-x)**p1d*del2  
           er_bet_de = r2d * dlog(x)*(1.-x)**p2d*del2  
           er_bet_dhtil = r1dtil * dlog(x)*(1.-x)**p1_dtil*del2
           er_bet_detil = r2dtil * dlog(x)*(1.-x)**p2_dtil*del2

           er_p_dh = r1d * dlog(x)* dlog(1.-x) * bet_dh*(1.-x)**p1d*del2  
           er_p_de = r2d * dlog(x)* dlog(1.-x) * bet_de*(1.-x)**p2d*del2
           er_p_dhtil = r1dtil*
     > dlog(x) * dlog(1.-x) * bet_dhtil*(1.-x)**p1_dtil*del2
           er_p_detil = r2dtil*
     > dlog(x) * dlog(1.-x) * bet_detil*(1.-x)**p2_dtil*del2
  
*** estimated error on alpha parameter
           dalp_dh = 0.1*alf_dh
           dalp_de = 0.1*alf_de
           dalp_dhtil = 0.2*ah_dtil

           er_alp_dh = r1d * dlog(x)
           er_alp_de = r2d * dlog(x)
           er_alp_dhtil = r1dtil * dlog(x)
           er_alp_detil = r2dtil * dlog(x)

*** estimated error on normalization
           danormd = 0.15*anormd
           danormde = 0.15*anormde
           dadtil = 0.15*adtil
           dadetil = 0.15*adetil

           er_dh = sqrt((er_bet_dh*dbet_dh)**2 + (er_p_dh*dp1d)**2 + 
     > (er_alp_dh*dalp_dh)**2 + danormd**2)
*
           er_de = sqrt((er_bet_de*dbet_de)**2 + (er_p_de*dp2d)**2
     > + (er_alp_de*dalp_de)**2 + danormde**2)
*
           er_dhtil = sqrt((er_bet_dhtil*dbet_dhtil)**2 + 
     > (er_p_dh*dp1_dtil)**2
     > +                    (er_alp_dhtil*dalp_dhtil)**2+ dadtil**2)
*
           er_detil = sqrt((er_bet_detil*dbet_detil)**2 + 
     >                     (er_p_de*dp2_dtil)**2
     > +                    (er_alp_detil*dalp_dhtil)**2+ dadetil**2)


************************************************************************************
*** H_d
           gd = gd /(1.-x)
           hd = gd *r1d

           hd_plus = gd * (r1d+er_dh)
           hd_minus = gd * (r1d-er_dh)


*** E_d
           ged= ged /(1.-x)*(1.-zeta)
           ed =  -ged * r2d 

           ed_plus = -ged * (r2d+er_de)
           ed_minus = -ged * (r2d-er_de)

*** Htil_d
           hdtil = gdtil /(1.-x) *r1dtil

           hdtil_plus = gdtil /(1.-x) * (r1dtil+er_dhtil)
           hdtil_minus = gdtil /(1.-x) * (r1dtil-er_dhtil)

*** Etil_d  (zeta = 0)
           if(zeta.eq.0.and.t.eq.0.)  then

           edtil = amp*pi*(1.-x)**6*adetil*x**(-ah_dtil)*
     >  (amp/3./(emtil0t*(1.-x))**3 - 
     >  4.*(amp*x+amqdtil)/5./(emtil0t*(1.-x))**4*( (1.-2.*x)*amp**2-
     >   amdiq_dht**2 + alam_dht**2) )
 
*** Etil_d  (zeta neq 0) 
           else
           edtil = gedtil /(1.-x) *r2dtil *(1.-zeta) /zeta

           edtil_plus = gedtil /(1.-x) * (r2dtil+er_detil)
           edtil_minus = gedtil /(1.-x) * (r2dtil-er_detil)

           endif

c           print 600,x,zeta,q2,zeta,t,hu,r1u,hu_plus,hu_minus
 600	format(1x,10(g10.4,1x))
****************************************************************************************

        return
          end






