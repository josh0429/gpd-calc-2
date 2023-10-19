C     
C-----------------------------------------------------------------------
C  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  + 
C-----------------------------------------------------------------------
C
C File - Version - Date: struct_ev_gpd_tmd.f - v 1.0 - 4/14/2014
C                   
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C     Calculation of PQCD evolution for GPDs and moments of TMDs 
C
C          qsquare  DIS momnetum transfer squared (GeV^2)      
C          zeta     longitudinal momentum transfer fraction
C          t        negative momentum transfer squared (GeV^2)
C          
C          jgpd    GPD type: 
C                  jgpd=1 --> H
C                  jgpd=2 --> E
C                  jgpd=3 --> H tilde
C                  jgpd=4 --> E tilde
CC
*         ierr  flag on error calculation 
**              ierr=1 --> no error
**              ierr=2 --> (GPD) + error 
*
*         itmd  flag on GPD or TMD moment evolution
***             itmd = 1 --> GPD evolution
***             itmd = 2 --> TMD moments evolution
**
***       jt    TMD type for TMD option only
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Author: Simonetta Liuti
C
C-----------------------------------------------------------------------
C
    	subroutine struct_ev_gpd_tmd
     >     (qsquare,zeta,t,jgpd,ierr,itmd,jt,uv0,dv0,g0,ub0,db0,sb0,cb0,
     >     init_q02)

        
**
*
* Calculates Q^2 evolution of parton distributions
* 
**
        implicit none
!
        integer*4 num0,neq,num,nout,kstops,nflavor,jt,itmd,ierr,
     > i,nstep,ntot,jgpd

	parameter(num0=490,ntot=4410)

	real*8 uvi(num0),dvi(num0),ubi(num0),dbi(num0)
     >           ,sbi(num0),cbi(num0),glui(num0)

        real*8 xbj(num0),XMAX(NUM0),q2(num0),q02(num0),X(num0)

        real*8 uv,dv,u_p,d_p,s_p,c_p,b_p,t_p,ub,db,glu

        real*8 uv0(num0),dv0(num0),g0(num0),ub0(num0),db0(num0)
     >           ,sb0(num0),cb0(num0)

        real*8 yprime(ntot),y(ntot),zz(ntot)

        real*8 z,xstart,yyy

*** TMDs
        real*8 tu(9),td(9)

*** GPDs  
        real*8 hu,hu_plus,
     > eu,eu_plus,
     > hutil,hutil_plus,
     > eutil,eutil_plus,
     > hd,hd_plus,
     > ed,ed_plus,
     > hdtil,hdtil_plus,
     > edtil,edtil_plus,hg,eg,htg,hub,hdb,hsb,hcb,
     > init_q02,etg,eub,edb

        real*8 xmass(10)

        real*8 erstep

        real*8 qsquare,t,q002,xlam,zeta

      COMMON/parton/uvi,dvi,ubi,dbi,sbi,cbi,glui
      COMMON/evol/xbj,q2,q02

      COMMON/STASH1/X,XMAX,NEQ,NSTEP
     1  ,ERSTEP,NOUT,KSTOPS
      COMMON/STASH/YPRIME,Y,ZZ
      COMMON/MASS/XMASS,NFLAVOR,NUM,XLAM

***********************************************************************

!				=> put data in common mass
      XMASS(1)=1.4**2
      XMASS(2)=4.5**2
      XMASS(3)=180.**2
!  				=> masses of the charm bottom and top quarks
        IF(qsquare.LT.1.96d0)NFLAVOR=3
        IF(qsquare.GE.1.96d0)NFLAVOR=4
        IF(qsquare.GT.20.25d0)NFLAVOR=5
        IF(qsquare.GT.2.89D4)NFLAVOR=6

c        PRINT *, NFLAVOR
        
! START WITH NF=4 
! MAKE THE CONNECTION WITH LambdaNF=4 AND NF=3,5,6, A LA' MARCIANO
! LATER ON 
!				=> n flavors
c      XLAM=.200d0
c	xlam = 0.300d0       !MRST98,NLO
!				=> XLAM=LAMBDA QCD
      NUM=num0
! 				=> number of x points, it cannot exceed 10000
*****
* now read in the starting values of the quark distributions
* Note that the distributions fill the array
*************************************************************************
*************************************************************************
*** INITIAL SCALE

      if(itmd.eq.1) then
       q002= init_q02
c      q002=0.09362d0
c      q002=0.26d0
c      q002=0.50d0
      else
* Bacchetta
      q002=0.3d0
      endif

       DO I=1,NUM
        IF(I.LE.290)THEN
          YYY=LOG(1.E+4)*(330.-FLOAT(I)+1.)/330.
          Z=EXP(-YYY)
c  note that the small x region is expanded.
c  the range of bjorken x (here called Z) is 1.e-4 to 1.
        ELSE
          XSTART=EXP(-LOG(1.E4)*41./330.)
          Z=XSTART+(FLOAT(I)-290.)*(1.-XSTART)/201.
        END IF

*****
***** define Lambda a la Marciano (LO)
*****
      xbj(i)=z
      q2(i)=qsquare
      if(q2(i).le.xmass(1)) then 
         xlam = 0.247d0         !!! Ahmad et al, PRD 2005
      else
      if(q2(i).le.xmass(2)) then 
         xlam = 0.215d0
        else
            xlam=0.165d0  
       endif
       endif


**** INSERT INITIAL PARAMETRIZATION
      if(z.ge.zeta) then

         go to (1200,1300)itmd
**************
************** GPD PARAMETRIZATION
 1200  continue   
       call diquark_sub_GGLA(z,qsquare,zeta,t,
     > hu,hu_plus,
     > eu,eu_plus,
     > hutil,hutil_plus,
     > eutil,eutil_plus,
     > hd,hd_plus,
     > ed,ed_plus,
     > hdtil,hdtil_plus,
     > edtil,edtil_plus,
     > hg,eg,htg,hub,hdb,hsb,hcb,init_q02,etg,eub,edb)
         go to 1600
1300    continue
     
 1600   continue

      q02(i) = q002
      X(i)=LOG(q02(i))


********* for NLO evolution convert into V and T linear combinations
********* OF PARTON DISTRIBUTIONS
********* V and T are NS

      if(itmd.eq.1) then
      if(ierr.eq.1) then

***************************************************************************************

      go to (110,120,130,140)jgpd
 110  uv= z*hu                  !(hu-hu_diq_c)
      dv= z*hd !(hd-hd_diq_c)
      ub= z*hub
      db= z*hdb
      u_p = uv + 2.*ub
      d_p = dv + 2.*db
      s_p= z*hsb
	c_p=z*hcb
	b_p=0.d0
	t_p=0.d0
        glu = hg
       go to 200

 120  uv=z* eu ! eu_diq_c 
      dv=z* ed  !ed_diq_c 
      UB=z*eub
      db=z*edb
      u_p = uv + 2.*ub
      d_p = dv + 2.*db
      s_p=0.
	c_p=0.d0
	b_p=0.d0
	t_p=0.d0
        glu = eg
        go to 200

c 130  u_p=z* hutil  
c      d_p=z* hdtil
**** Construct F_1 
 130  uv=z* hutil/(0.404*(1.+1.47*z**1.41))
      dv=-z* hdtil/(0.274*(1.+2.65*z**1.25))
      UB=0.D0
      db=0.d0
      u_p = uv + 2.*ub
      d_p = dv + 2.*db
      s_p=0.
	c_p=0.d0
	b_p=0.d0
	t_p=0.d0
        glu = htg
        go to 200

 140  uv=z* eutil/(0.404*(1.+1.47*z**1.41))
      dv=-z* edtil/(0.274*(1.+2.65*z**1.25))
      UB=0.D0
      db=0.d0
      u_p = uv + 2.*ub
      d_p = dv + 2.*db
      s_p=0.
	c_p=0.d0
	b_p=0.d0
	t_p=0.d0
        glu=etg
 200    continue

        else

      go to (111,121,131,141)jgpd
 111  uv=z* hu_plus
      dv=z* hd_plus
       UB=0.D0
      db=0.d0
      u_p = uv + 2.*ub
      d_p = dv + 2.*db
      s_p=0.
	c_p=0.d0
	b_p=0.d0
	t_p=0.d0
        glu=0.d0
        go to 201

 121  uv=z* eu_plus
      dv=z* ed_plus
      UB=0.D0
      db=0.d0
      u_p = uv + 2.*ub
      d_p = dv + 2.*db
      s_p=0.
	c_p=0.d0
	b_p=0.d0
	t_p=0.d0
        glu=0.d0
        go to 201

c 130  u_p=z* hutil  
c      d_p=z* hdtil 
 131   uv=z* hutil_plus/(0.404*(1.+1.47*z**1.41))
       dv=-z* hdtil_plus/(0.274*(1.+2.65*z**1.25))
      UB=0.D0
      db=0.d0
       u_p = uv + 2.*ub
       d_p = dv + 2.*db
       s_p=0.
	c_p=0.d0
	b_p=0.d0
	t_p=0.d0
        glu=0.d0
        go to 201
 141   uv=z* eutil_plus/(0.404*(1.+1.47*z**1.41))
       dv=-z* edtil_plus/(0.274*(1.+2.65*z**1.25))
      UB=0.D0
      db=0.d0
       u_p = uv + 2.*ub
       d_p = dv + 2.*db
       s_p=0.
	c_p=0.d0
	b_p=0.d0
	t_p=0.d0
        glu=0.d0
 201    continue

        endif
        endif

*** V_u and V_d
	y(i) = uv
	y(i+num) = dv
*** T3
      Y(i+2*NUM)= u_p - d_p    
*** T8
      Y(i+3*NUM)= u_p + d_p - 2.*s_p
*** T15
      Y(i+4*NUM)= u_p + d_p + s_p! - 3.*c_p
*** T24
      Y(i+5*NUM)= u_p + d_p + s_p! + c_p !- 4.*b_p
*** T35
      Y(i+6*NUM)= u_p + d_p + s_p! + c_p !+ b_p - 5.*t_p

******* Singlet part Sigma
*** Sigma
      Y(i+8*NUM)= u_p + d_p + s_p! + c_p !+ b_p + t_p
***** K=7 is gluons
      Y(i+7*NUM) = glu

      XMAX(i) = LOG(Q2(i))

      else
	y(i) = 0.d0
	y(i+num) = 0.d0
         endif 
      END DO

*************************************************************************
*************************************************************************
*************************************************************************
      NEQ=NTOT
C num is the number of x points on the grid, neq is the total number
C of equations of which there are the number of X points times the total
c number of  sets of distributions (up valence, down valence, anti up,
c anti down, gluons, charm, strange, top and bottom)
      NSTEP=128
!					==> number of Q^2 points
c      NOUT=1
      NOUT=0
c      ERSTEP=.00001
      ERSTEP=.1
!
      CALL SPACEOUT
!
******* Convert back to parton distributions from linear combinations
*******
c      KIP=NUM

      if(jgpd.eq.1.or.jgpd.eq.2) then
      DO I=1,num

	uvi(i) = y(i)
	dvi(I) = y(i+num)
	glui(i) = y(i+7*num)

	ubi(i) =(y(i+2*num)+y(i+3*num)/3.+y(i+4*num)/6.+
     >         y(i+5*num)/10.+y(i+6*num)/15.+y(i+8*num)/3.-2.*y(i))/4.

	dbi(i) =(-y(i+2*num)+y(i+3*num)/3.+y(i+4*num)/6.+
     >    y(i+5*num)/10.+y(i+6*num)/15.+y(i+8*num)/3.-2.*y(i+num))/4.

	sbi(i) =-y(i+3*num)/6.+y(i+4*num)/24.+
     >       y(i+5*num)/40.+y(i+6*num)/60.+y(i+8*num)/12.
        
        cbi(i) = -y(i+4*num)/8.+y(i+5*num)/40.+y(i+6*num)/60.+
     >      y(i+8*num)/12.
        
        uv0(i) = uvi(i)
        dv0(i) = dvi(i)
        g0(i) = glui(i)
        ub0(i) = ubi(i)
        db0(i) = dbi(i)
        sb0(i) = sbi(i)
        cb0(i) = cbi(i)
      END DO

        else

      DO I=1,num

	uvi(i) = y(i)*0.404*(1.+1.47*xbj(i)**1.41)
	dvi(I) = - y(i+num)*0.274*(1.+2.65*xbj(i)**1.25)
	glui(i) = y(i+7*num)

	ubi(i) =(y(i+2*num)+y(i+3*num)/3.+y(i+4*num)/6.+
     >         y(i+5*num)/10.+y(i+6*num)/15.+y(i+8*num)/3.-2.*y(i))/4.

	dbi(i) =(-y(i+2*num)+y(i+3*num)/3.+y(i+4*num)/6.+
     >    y(i+5*num)/10.+y(i+6*num)/15.+y(i+8*num)/3.-2.*y(i+num))/4.

	sbi(i) =-y(i+3*num)/6.+y(i+4*num)/24.+
     >      y(i+5*num)/40.+y(i+6*num)/60.+y(i+8*num)/12.

        uv0(i) = uvi(i)
        dv0(i) = dvi(i)
        g0(i) = glui(i)

      END DO


        endif

c      CLOSE(UNIT=5)
c      CLOSE(UNIT=7)
	RETURN
      END
**************************************************************************
**************************************************************************
**************************************************************************
      SUBROUTINE FAROUT
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/STASH1/X,XMAX,NEQ,NSTEP
     1   ,ERSTEP,NOUT,KSTOPS
      COMMON/STASH/YPRIME(4410),Y(4410),ZZ(4410)
      COMMON/MASS/XMASS(10),NFLAVOR,NUM,XLAM
c      user supplied code to test and/or output values
c       set kstops .ne. 0 to end integration
c
C this routine is called every nout steps
C       ****************************************************************
C       ****************************************************************
c      WRITE(10,'(1x,g9.3)')exp(X)*XLAM*XLAM
c      WRITE(6,'(1x,g9.3)')exp(X)*XLAM*XLAM
      KIP=NUM
      DO I=1,KIP
	val = y(i) + y(i+num)
	glu = y(i+2*num)
	sea = y(i+3*num)+y(i+4*num)
c      WRITE(8,'(1x,7(g9.3,1x))')ZZ(I),(Y(I+(K-1)*NUM),K=1,6)
c      WRITE(6,'(1x,7(g9.3,1x))')ZZ(I), (Y(I+(K-1)*NUM),K=1,6)
      END DO
	end
************************************************************************
************************************************************************
      SUBROUTINE CHANGES(I)
C
C       user supplied code to compute yprime(k)=f(k,y(l))
C
C       ****************************************************************
C       ****************************************************************

      implicit none

      integer*4 num,i,j,nlo,nstep,nflavor,k,kk,neq,nout,ntot,kstops,
     > num0

      parameter (num0=490,ntot=4410)

        real*8 xbj(num0),XMAX(NUM0),q2(num0),q02(num0),X(num0)

        real*8 yprime(ntot),y(ntot),zz(ntot)

        real*8 xstart,yy,xpoint,ypoint,ypp,yppp,yyp,xjacobian

        real*8 xmass(10)

        real*8 qsq,xlam,zeta

        real*8 zeta1

        real*8 alpha_my,alp0,alp02

        real*8 pi

        real*8 cf,ca,tr,tf

*** NLO
        real*8 zeta3

        real*8 alx,alx2,alx1,alx12,alxp

        real*8 pqq_p,pqq_m,pqqv,pqqv_1,pqqbarv
 
        real*8 pgq_p,pgq_m,pgq1

        real*8 pqg_p,pqg_m,pa,pb,pc,pqq1s1,pqq1s,pqq1s0,pqg1

        real*8 pgg_p1,pgg_p2,pgg_m,pgg1,pgg1_1,pgg1_0

***
*** ENDPOINTS
      real*8 extra,extra2,es2

      real*8 erstep


      COMMON/STASH1/X,XMAX,NEQ,NSTEP
     1   ,ERSTEP,NOUT,KSTOPS
      COMMON/STASH/YPRIME,Y,ZZ
      common /evol/xbj,q2,q02

      COMMON/MASS/XMASS,NFLAVOR,NUM,XLAM
        common/gpdkin/zeta1

c      dimension FLAG(10)

*******************************************************************************
      PI=4.*ATAN(1.)

      QSQ=EXP(X(I))

	ca=3.
	cf=4./3.
	tr=0.5
	tf=tr*float(nflavor)

	zeta3 = 1.202057d0

        num=num0
        zeta=zeta1
*************************
*************************
****** ALPHA

C  10/2009 insert my own code for alpha
C
c          alp0 = 0.125
        alp0=ALPHA_MY(qsq)/2.
c        print *, ALPHA_MY(qsq)
c        print *, alp0*2.
c        alp0 = 0.212
c          print *, i,xbj(i),qsq,alp0*2.
          alp02 = alp0*alp0

****
      DO K=1,9
         KK=K-1
        YPRIME(I+KK*NUM)=0.d0
      END DO
C
C INITIALISE DERIVATIVES TO ZERO
C
C derivatives have two parts one integral and one delta fn
C FIRST DO INTEGRALS
C
C to avoid 0/0 in plus dist. integrals we will omit the first point
C and then include it below in the non plus pieces
C **********************************************
      XSTART=EXP(-LOG(1.E4)*41./330.)
        IF(I.LE.290)THEN
          YPP=LOG(1.E+4)*(330.-FLOAT(I)+1.)/330.
          XPOINT=EXP(-YPP)
        ELSE
          XPOINT=XSTART+(FLOAT(I)-290.)*(1.-XSTART)/201.
        END IF

        DO J=I+1,NUM
          IF(J.LE.290)THEN
            YPPP=LOG(1.E4)*(330.-FLOAT(J)+1.)/330.
            YPOINT=EXP(-YPPP)
            XJACOBIAN=YPOINT*LOG(1.E4)/330.
C the variable transform is ypoint to j
          ELSE
            YPOINT=XSTART+(FLOAT(J)-290.)*(1.-XSTART)/201.
            XJACOBIAN=(1.-XSTART)/201.
          END IF

!
          YY=XPOINT/YPOINT
          YYP = (XPOINT-ZETA)/(YPOINT-ZETA)
          IF(YY.GT.1.)THEN
            PRINT 123,I,J
  123       FORMAT(2I5)
          ELSE
          END IF

C
C first do NS quark distributions
C ********************** ************************************
	  alx =log(yy)
	  alx2 =alx*alx
	  alx1 =log(1.-yy)
	  alx12 =alx1*alx1
	  alxp = alx*alx1

	  pqq_p = (1.+yy*yy)     !times 1-z
	  pqq_m = (1.+yy*yy)/(1.+yy)
**
* separate out "+" terms 
* this is for both non singlet and singlet
*** 
*** this is not divided by (1-yy) yet

	  pqqv = ( 
     1      cf*ca*(67./18.-pi*pi/6.)
     1      + cf*tf*(-10./9.) )*pqq_p

* add in non "+" terms
* non singlet
***
         pqqv_1 =    
     2   cf*cf*((-2.*alxp-1.5*alx)*pqq_p/(1.-yy)
     1  -(1.5+3.5*yy)*alx -0.5*(1.+yy)*alx*alx-         
     2      5.*(1.-yy) ) 
     2  + cf*ca*((0.5*alx*alx+11./6.*alx)*pqq_p/(1.-yy)
     1  + (1.+yy)*alx+20./3.*(1.-yy)) + 
     2         cf*tf*((-2./3.*alx)*pqq_p/(1.-yy)-4./3.*(1.-yy) )  

	  pqqbarv=cf*(cf-ca/2.)*
     1	        (pqq_m*2.*es2(yy)
     1    +  2.*(1.+yy)*alx+4.*(1.-yy) )      


C
C u valence
C ********************** ************************************
*
* Notice: I do this only for the uv and dv distributions. I will take 
* care of the rest later (for F_L).
* 

          nlo=0.
          if(nlo.eq.0) then
             YPRIME(I) =  alp0* cf*( (1.+YY*YYP)*Y(J)-(1.+yyp/yy)*Y(I) )
     3    /(1.-YY) *XPOINT/YPOINT/YPOINT 
     4   *XJACOBIAN + YPRIME(I)
             
          else
          YPRIME(I) = ( alp0 * cf*( (1.+YY*YYP)*Y(J)-(1.+yyp/yy)*Y(I) )
     3    /(1.-YY) 
     3   +  alp02 * ( ( pqqv*Y(J) - 
     1   2.*(cf*ca*(67./18.-pi*pi/6.)-cf*tf*10./9.)*Y(I) )
     1   /(1.-yy) +  
     2   (pqqv_1 - pqqbarv)*Y(J) ) 
     2  ) *XPOINT/YPOINT/YPOINT 
     4   *XJACOBIAN + YPRIME(I)
       endif
C
C d valence
C ********************** ************************************
          if(nlo.eq.0) then
          YPRIME(I+NUM) =  alp0* cf*( (1.+YY*YYP)*Y(J+NUM)-
     >   (1.+yyp/yy)*Y(I+NUM) )
     3    /(1.-YY) *XPOINT/YPOINT/YPOINT 
     4   *XJACOBIAN + YPRIME(I+NUM)

          else
          YPRIME(I+NUM) = ( alp0 * cf*( (1.+YY*YYP)*Y(J+NUM)-
     >                      (1.+yyp/yy)*Y(I+NUM) )/(1.-YY) 
     3   +  alp02 * ( ( pqqv*Y(J+NUM) - 
     1   2.*(cf*ca*(67./18.-pi*pi/6.)-cf*tf*10./9.)*Y(I+NUM) )
     1   /(1.-yy) +  
     2   (pqqv_1 - pqqbarv)*Y(J+NUM) ) 
     2  ) *XPOINT/YPOINT/YPOINT 
     4   *XJACOBIAN + YPRIME(I+NUM)
          endif


C
C T3,T8,T15,T24,T35 (K=2,6)
C ********************** ************************************
          DO K=2,6

             if(nlo.eq.0)then

          YPRIME(K*NUM+I) = alp0 * cf* 

     1   ((1.+YY*YY)*Y(J+K*NUM)-2.*Y(I+K*NUM))/(1.-YY)

     1    *XPOINT/YPOINT/YPOINT*XJACOBIAN 

     >  + YPRIME(I+K*NUM)

          else
          YPRIME(K*NUM+I)=(alp0*cf*

     1   ((1.+YY*YY)*Y(J+K*NUM)-2.*Y(I+K*NUM))/(1.-YY)

     2   + alp02*( (pqqv*Y(J+K*NUM) - 

     1	  2.*(cf*ca*(67./18.-pi*pi/6.)-cf*tf*10./9.)*Y(I+K*NUM) )

     3   /(1.-yy)   

     2   + (pqqv_1 + pqqbarv)*Y(J+K*NUM) 

     3   ) )*XPOINT/YPOINT/YPOINT

     4   *XJACOBIAN+YPRIME(I+K*NUM)
          endif

        END DO
C ******************************   *******************************
C now for gluon distributions (K=7) and SIGMA  K=8
C ******************              ******************************
	pgq_p = 1.+(1.-yy)**2 !times yy
	pgq_m = -1.-(1.+yy)**2 !times yy

	pgg_p1= yy    ! times 1-yy
	pgg_p2= (1.-yy)*(1.+yy*yy)     !times yy
	pgg_m= (-yy**2-(1.+yy)**2*(1.+yy**2))/(1.+yy) !times yy

	pgq1=cf*cf * ( (-2.5-3.5*yy+(2.+3.5*yy)*alx-
     1     (1.-0.5*yy)*alx*alx-2.*yy*alx1) *xpoint/ypoint**2 -
     1     (3.*alx1+alx1*alx1)*pgq_p /ypoint ) +
     2     cf*ca*( (28./9.+65./18.*yy+44./9.*yy*yy-
     3          (12.+5.*yy+8./3.*yy*yy)*alx+
     4          (4.+yy)*alx*alx+2.*yy*alx1) *xpoint/ypoint**2 +
     5          es2(yy)*pgq_m/ypoint  +
     6      ( 0.5-2.*alxp+0.5*alx*alx+11./3.*
     7      alx1+ alx1*alx1-pi*pi/6.)*pgq_p/ypoint )+
     8     cf*tf*( -4./3.*yy *xpoint/ypoint**2
     9     -(20./9.+4./3.*alx1)*pgq_p/ypoint ) 

        pgg1=
     5   (-ca*tf*20./9.+ 
     1   ca*ca*(67./9.-pi*pi/3.))*pgg_p1

	pgg1_1 = cf*tf*(-16.+8.*yy+20./3.*yy*yy-(6.+10.*yy)*alx
     4        -(2.+2.*yy)*alx*alx) +
     1     ca*tf*(2.-2.*yy+26./9.*yy*yy-4./3.*(1.+yy)*alx) +
     1     ca*ca*(13.5*(1.-yy)+67./9.*yy*yy-
     8         (25.-11.*yy+44.*yy*yy)/3.*alx+
     1	     4.*(1.+yy)*alx*alx -4.*alxp+alx*alx*pgg_p1/(1.-yy))  

	pgg1_0 = cf*tf*4./3. 
     1       + ca*tf*(-20./9.*pgg_p2 -26./9.) +
     1          ca*ca*( -67./9. + 2.*pgg_m*es2(yy)+ 
     2       (67./9.-4.*alxp+alx*alx-pi*pi/3.)*pgg_p2)
 
***
*** Pgq
        if(nlo.eq.0) then

         YPRIME(I+7*NUM)=alp0*cf/YPOINT*(1.+(1.-YY)*(1.-YYP))*Y(8*NUM+J)
     1	   *(1./(1.-zeta/2.))*XJACOBIAN + YPRIME(I+7*NUM)       

         else

         YPRIME(I+7*NUM) = (alp0*cf/YPOINT*(1.+(1.-YY)**2)*Y(8*NUM+J)

     2       + ALP02 * pgq1 * Y(8*NUM+J) )  

     3    *XJACOBIAN +YPRIME(I+7*NUM)

            endif
***
*** Pgg
            if(nlo.eq.0) then
         YPRIME(I+7*NUM)=ALP0*(XPOINT/YPOINT/YPOINT*(3.*(YY+YYP*YYP/YY)
     1 * Y(J+7*NUM)- 3.*(1.+YYP/YY)*Y(I+7*NUM))
     1  /(1.-YY) + 6.*(1.-YYP)/YPOINT*(1.+YY*YYP)
     2  *Y(7*NUM+J) 
     4   ) *XJACOBIAN+YPRIME(I+7*NUM)
          else
       	  YPRIME(I+7*NUM)=(ALP0*(XPOINT/YPOINT/YPOINT*(6.*YY*Y(J+7*NUM)-
     1   6.*Y(I+7*NUM))
     1  /(1.-YY) + 6.*(1.-YY)/YPOINT*(1.+YY*YY)
     2  *Y(7*NUM+J)    
     2  + ALP02 *
     3  (xpoint/ypoint/ypoint*(
     1  (pgg1*Y(J+7*NUM)- Y(I+7*NUM)*
     3   (ca*ca*(67./9.-pi*pi/3.)-ca*tf*20./9.) )/(1.-yy) 
     1   + pgg1_1*Y(J+7*NUM) ) +pgg1_0*Y(J+7*NUM)/ypoint))
     4   *XJACOBIAN) +  YPRIME(I+7*NUM)
          endif
C *****************************   ***************
C now SIGMA  K=8
C*****************************************************
	  pqg_p=yy*yy+(1.-yy)**2
	  pqg_m=yy*yy+(1.+yy)**2

	  pa=cf*cf*(-(1.5*alx+2.*alxp)*pqq_p/(1.-yy)
     1                -1. + yy + (0.5-1.5*yy)*alx
     1       - 0.5*(1.+yy)*alx*alx + 2.*pqq_m*es2(yy) )
          pb=cf*ca*( 14./3.*(1.-yy) - pqq_m*es2(yy)+
     1 (11./6.*alx+0.5*alx*alx)*pqq_p/(1.-yy) ) 
          pc=cf*tf*(-16./3.+ 40./3.*yy-2./3.*alx*pqq_p/(1.-yy) +
     1                (10.*yy+16./3.*yy*yy+2.)*alx
     1            -112./9.*yy*yy-2.*(1.+yy)*alx*alx )

	  pqq1s1 = pa + pb + pc

         pqq1s= (
     1   cf*ca*(67./18.-pi*pi/6)-cf*tf*10./9. )*pqq_p 

	  pqq1s0 = cf*tf*40./9.

	  pqg1=( cf*tf*(4.-9.*yy-(1.-4.*yy)*alx-(1.-2.*yy)*
     1         alx*alx+4.*alx1+
     2        (2.*(alx1-alx)*(alx1-alx)-4.*(alx1-alx)-2./3.*pi*pi+
     1      10.)*pqg_p )+
     3      ca*tf*(182./9.+14./9.*yy+(136.*yy-38.)/3.
     3      *alx-4.*alx1-(2.+8*yy)*alx*alx+2.*pqg_m*es2(yy)
     4     +(-alx*alx+44./3.*alx-2.*alx1*alx1+4.*alx1+
     4      pi*pi/3.-218./9.)*pqg_p) )*xpoint/ypoint/ypoint
     5     +ca*tf*40./9./ypoint          
***
*** Pqq

c	  pqq1=pqq1s

          if(nlo.eq.0) then
          YPRIME(I+8*NUM)=alp0*cf*XPOINT/YPOINT/YPOINT*((1.+YY*YYP)*

     1   Y(J+8*NUM) - (1.+YYP/YY)*Y(I+8*NUM))/(1.-YY)

     4   *XJACOBIAN + YPRIME(I+8*NUM)

          else
          YPRIME(I+8*NUM)=( alp0*cf*XPOINT/YPOINT/YPOINT*
     > ( (1.+YY*YY)*Y(J+8*NUM) - 2.*Y(I+8*NUM) )/(1.-YY)

     2   + ALP02*( (pqq1s*Y(J+8*NUM)-

     1   2.*(ca*cf*(67./18.-pi*pi/6.)-cf*tf*10./9.)*Y(I+8*NUM))

     >   /(1.-yy) 

     1    + pqq1s1*Y(J+8*NUM)*xpoint/ypoint/ypoint 

     >    + pqq1s0*Y(J+8*NUM)/ypoint    

     4   ) ) *XJACOBIAN + YPRIME(I+8*NUM)

          endif
***
*** Pqg
          if(nlo.eq.0) then
           YPRIME(I+8*NUM) = alp0*(XPOINT/YPOINT**2

     2   *Y(J+7*NUM)*(.5-YY*((1. + YYP/YY)/2.-YYP) ) *2.*float(nflavor) 

     2  )*(1.-zeta/2.)*XJACOBIAN + YPRIME(I+8*NUM)

              else

              YPRIME(I+8*NUM) = (alp0 * XPOINT/YPOINT**2

     2   *Y(J+7*NUM)*(.5-YY*(1.-YY) ) *2.*float(nflavor) 

     3   + ALP02* 

     3   Y(J+7*NUM)*pqg1 ) *XJACOBIAN + YPRIME(I+8*NUM)  
              endif

        END DO
****************
****************
****************
C ******************************************
C
C NOW ADD DELTA FN PIECES
C ******************    *****************

        IF(I.EQ.NUM)THEN
          EXTRA=0.
          EXTRA2=0.
        ELSE
          EXTRA=LOG(((1.-XPOINT)**2)/(1.-zeta))
          extra2=EXTRA*EXTRA
        END IF

        DO KK=1,2
           if(nlo.eq.0) then
           K=KK-1
           YPRIME(I+K*NUM)=YPRIME(I+K*NUM) + Y(I+K*NUM)* 
     1          cf*( alp0*(1.5 + EXTRA))
        else

           YPRIME(I+K*NUM) = YPRIME(I+K*NUM) +Y(I+K*NUM)*
     1      (cf*( alp0*(1.5 + EXTRA))
     1     + alp02 *

     2       ( cf*cf*(3./8.-pi*pi/2.+6.*zeta3) +

     3        cf*ca*(17./24. +11.*pi*pi/18.-3.*zeta3)-

     4        cf*tf*(1./6.+2.*pi*pi/9.)

     5 + 1.*(ca*cf*(67./18.-pi*pi/6.)-cf*tf*10./9.)*extra) )
              endif
        END DO

        DO K=2,6
           if(nlo.eq.0) then
          YPRIME(I+K*NUM)=YPRIME(I+K*NUM) +Y(I+K*NUM)*
     1      cf*( alp0*(1.5 + EXTRA))        

          else

           YPRIME(I+K*NUM) = YPRIME(I+K*NUM) +Y(I+K*NUM)*
     1      (2.*(1.+4./3.*EXTRA)* alp0
     1     + alp02 *

     2       ( cf*cf*(3./8.-pi*pi/2.+6.*zeta3) +

     3        cf*ca*(17./24. +11.*pi*pi/18.-3.*zeta3)-

     4        cf*tf*(1./6.+2.*pi*pi/9.)  

     5 + 2.*(ca*cf*(67./18.-pi*pi/6.)-cf*tf*10./9.)*extra) )

             endif
       END DO

       if(nlo.eq.0) then
        YPRIME(I+7*NUM) = YPRIME(I+7*NUM) +Y(I+7*NUM)*alp0* 
     1         (-(NFLAVOR/6.-11./4.) + 1.5*EXTRA)*2.

        else
        YPRIME(I+7*NUM)=YPRIME(I+7*NUM) +Y(I+7*NUM)*(alp0* 
     1   (-(NFLAVOR/6.-11./4.)+1.5*EXTRA)*2.
     1     + alp02 * 
     2       ( ca*ca*(8./3.+3.*zeta3)-cf*tf-4./3.*ca*tf 
     3      + (ca*ca*(67./9.-pi*pi/3.)-ca*tf*20./9.)*extra) )
        endif

        if(nlo.eq.0)then
          YPRIME(I+8*NUM)=YPRIME(I+8*NUM) +Y(I+8*NUM)
     1        *cf*( alp0*(1.5 + EXTRA))

          else

         YPRIME(I+8*NUM)= YPRIME(I+8*NUM) + Y(I+8*NUM)

     1  * (alp0*2.*( 1.+4./3.*EXTRA) 

     1     + alp02 *

     2       (cf*cf* (3./8.-pi*pi/2.+6.*zeta3)+

     3        cf*ca* (17./24. +11.*pi*pi/18.-3.*zeta3)-

     4        cf*tf* (1./6.+2.*pi*pi/9.)  

     5     + 2.*(ca*cf*(67./18.-pi*pi/6.)-cf*tf*10./9.)*extra) )

      endif


C*****************************************
C
C    NOW ADD THE FIRST BIT TO THE NON PLUS INTEGRALS
C **************************** ********************
c        IF(I.LE.290)THEN
c          XJACOBIAN=XPOINT*LOG(1.E4)/330.
c        ELSE
c          XJACOBIAN=(1.-XSTART)/201.
c        END IF

c       yprime(I) = yprime(i) -
c     >   y(i)*fac0*alp02 *(8./3.) /2./pi/xpoint
c     > *xjacobian 
c       yprime(I+NUM) = yprime(i+NUM) - 
c     >   y(i+NUM)*fac0* (8./3.) *alp02/2./pi/xpoint 
c     > *xjacobian 

c        YPRIME(I+7*NUM)=(4./3. / XPOINT*
c     1   Y(8*NUM+I) + ALPHA(I)/2./xpoint* (
c     1  (-cf*cf*6.+cf*ca*(218./18.-pi*pi/6.)
c     1   -cf*tf*32./9.
c     1 + (-5.*cf**2+17./3.*cf*ca-4./3.*cf*tf)*log(0.98)
c     1 + (-cf**2+cf*ca)*log(0.98)**2)*Y(8*NUM+I) + 
c     1  (ca*tf*20./9.-ca*ca*(67./9.-pi*pi/3.))*Y(7*NUM+I)
c     2  ))*XJACOBIAN+ YPRIME(I+7*NUM)
        
c        YPRIME(8*NUM+I)=(.5*2.*float(nflavor)/
c     2  XPOINT * Y(I+7*NUM) 
c     1  + ALPHA(I)/2./xpoint* ( (cf*cf*4.-cf*ca*(67./9.-pi*pi/3.)
c     1 + cf*tf*20./9.)*Y(I+8*NUM) 
c     1 +(cf*tf*(5.-2./3.*pi*pi)+ca*tf*(-2.+pi*pi/3.)+
c     1  (cf-ca)*tf*log(0.98)**2)
c     1  *Y(I+7*NUM) ))
c     1  *XJACOBIAN+YPRIME(I+8*NUM)

c        DO K=1,9
c        DO K=3,9
c          KK=K-1
c          YPRIME(I+KK*NUM)=YPRIME(I+KK*NUM)   !*ALPHA/2.
c        END DO

c        aa=- y(i)*fac0*alp02 *(8./3.) /2./pi/xpoint
c     > *xjacobian 
c        bb= y(I)*cf*( alp0*(1.5 + 2.*EXTRA) )/ (2.*pi)
c        cc= y(i)*cf*
c     2  fac0*alp02*(pi*pi/3. + 0.5 - extra2)
c     1  / (2.*pi)  
c        dd= y(i)*alp02/(2.*pi)**2 *
c     1       ( cf*cf*(3./8.-pi*pi/2.+6.*zeta3)+
c     3        cf*ca*(17./24. +11.*pi*pi/18.-3.*zeta3)-
c     4        cf*tf*(1./6.+2.*pi*pi/9.)  
c     5       +2.*(ca*cf*(67./18.-pi*pi/6.)-cf*tf*10./9.)*extra) 

c        tot=aa+bb+cc+dd

c        write(6,'(1x,8(1x,g10.4),1x,i3)')xpoint,x(i),qsq,yprime(i),
c     >  bb,cc,dd,yprime(i)-tot,i
c        write(8,'(1x,8(1x,g10.4))')xpoint,x(i),qsq,yprime(i),bb,cc,dd
c     > ,yprime(i)-tot

c      END DO

      RETURN
      END

***************
***************
	double precision function es2(x)
	implicit real *8 (a-h,o-z)

	pi=acos(-1.)
	xm=-x
	es2 = -2.*dil(xm)+0.5*(log(x))**2-2.*log(x)*log(1.+x)
     >        -pi*pi/6.
	return
	end
***************
*************** Dilogarithm
***************
	double precision function dil(x)
	implicit real *8 (a-h,o-z)
        dimension y(8),f1(8)
        call dinter(y,0.d0,x,8)
        do i=1,8
           f1(i)=log(1.d0-y(i))/y(i)
           enddo
        dil= -dgaus1(8,0.d0,x,f1)
	return
	end

      SUBROUTINE SPACEOUT
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/STASH1/X,XMAX,NEQ,NSTEP
     1    ,ERSTEP,NOUT,KSTOPS
      COMMON/STASH/YPRIME(4410),Y(4410),ZZ(4410)
      DIMENSION YPBUF(4410,4),YLAST(4410),YCHK(4410)
      DIMENSION D0(4410),D1(4410),D2(4410),D3(4410)
	DIMENSION DX(4410),XMAX(490),X(490)
C
C       subroutine numerically integrates first order simultaneous
C       differential equations:
C                       yprime(k)=f(k,y(l),x)
C
C       before calling this routine, user must fill the following
C       variables in common block /stash/:
C       y(k),x : fill with initial conditions
C       xmax = end of integration interval
C       neq = number of differential equations
C       nstep = number of integration steps
C       erstep = error/step : corrector repeats up to 20 times or
C                until successive values agree within erstep
C       nout : subroutine farout is called every nout steps
C              nout=0 means no calls to farout
C
C       user must also supply code in subroutines changes and farout
C       to compute yprime(k) and to test/output values
C       setting kstops .ne. 0 in farour flags for return
C
c
c        
c     NUM is the number of x steps, right now set to 4410 which divided
c     by 9 is 490. So this is going x by x in the distribution
        
      NUM= NEQ/9
      KSTOPS=0
      KSTP=0
      KBUF=1
      
c     Iterates x by x in the distribution      
      do I=1,NUM

c     Iterates over q2 steps. Xmax is the log of the final q2
c     X(i) is log of the initial q2. Nstep is the number of
c     steps in q2 you are iterating over.
c         
c     Creates an interval in q2 array and transposes that array
c     over the total number of quark distributions.
      DX(I)=(XMAX(I)-X(I))/NSTEP
	enddo
	do I=1,NUM
	   do L=1,8
	   DX(I+L*NUM)=DX(I)
	   enddo
	   enddo

c     Iterates over total number of points in the dx array which is
c     now 4410 or NTOT.
c
c     Setting all the derivatives to 0 in the beginning.
      DO 100 L=1,NEQ
        D0(L)=0.
        D1(L)=0.
        D2(L)=0.
        D3(L)=0.
        DO 50 K=1,4
 50        YPBUF(L,K)=0.
  100 CONTINUE
      DO 195 I=1,NUM
      CALL CHANGES(I)
 195	 ENDDO
	 DO 200 L=1,NEQ
        YPBUF(L,KBUF)=YPRIME(L)
  200 D0(L)=YPBUF(L,KBUF)
C
  250 CONTINUE
      IF(NOUT.EQ.0) GO TO 300
c      IOUT=NOUT*(KSTP/NOUT)-KSTP
c      IF((IOUT.NE.0) .AND. (KSTP.NE.NSTEP)) GO TO 300
c      CALL FAROUT
  300 IF(KSTOPS.NE.0) RETURN
      IF(KSTP.GE.NSTEP) RETURN
C
      DO 350 L=1,NEQ
        YLAST(L)=Y(L)
        Y(L)=Y(L)+DX(L)*D0(L)
        IF(KSTP.GE.1) Y(L)=Y(L)+(1./2.)*DX(L)*D1(L)
        IF(KSTP.GE.2) Y(L)=Y(L)+(5./12.)*DX(L)*D2(L)
        IF(KSTP.GE.3) Y(L)=Y(L)+(3./8.)*DX(L)*D3(L)
  350 CONTINUE
C
      KSTP=KSTP+1
      KBUF=MOD(KBUF,4)+1
      KB1=MOD(KBUF+2,4)+1
      KB2=MOD(KB1+2,4)+1
      KB3=MOD(KB2+2,4)+1
      KERR=0
C

C 400  CONTINUE

      DO 450 L=1,NUM
       X(L)=X(L)+DX(L)
      CALL CHANGES(L)
 450	ENDDO

      DO 500  L=1,NEQ
       YPBUF(L,KBUF)=YPRIME(L)
        D0(L)=YPBUF(L,KBUF)
        D1(L)=YPBUF(L,KBUF)-YPBUF(L,KB1)
        D2(L)=YPBUF(L,KBUF)-2.*YPBUF(L,KB1)
     1                           +YPBUF(L,KB2)
        D3(L)=YPBUF(L,KBUF)-3.*YPBUF(L,KB1)
     1                           +3.*YPBUF(L,KB2)-YPBUF(L,KB3)
        YCHK(L)=Y(L)
        Y(L)=YLAST(L)+DX(L)*(D0(L)-D1(L)/2.)
        IF(KSTP.GE.2) Y(L)=Y(L)-DX(L)*D2(L)/12.
        IF(KSTP.GE.3) Y(L)=Y(L)-DX(L)*D3(L)/24.
  500 CONTINUE
C
      DO 600 L=1,NEQ
        ERROR=ABS(Y(L)-YCHK(L))
        IF(ERROR.LT.ERSTEP) GO TO 600
        KERR=KERR+1
C        IF(KERR.LT.20) GO TO 400
C                    print 901,l,x,error
  600 CONTINUE
C
      GO TO 250
C
C 901   format(' convergence error for y(',i2,') at x = ',e15.5,/
C     1' after 20 corrections, error = ',e15.5)
C
      END

