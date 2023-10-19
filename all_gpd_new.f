C
C-----------------------------------------------------------------------
C  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  + 
C-----------------------------------------------------------------------
C
C File - Version - Date: all_gpd_new.f - v 1.0 - 4/14/2014
C                   
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Author: Simonetta Liuti
C
C-----------------------------------------------------------------------
C
C Calculation of  the Compton form factors 
C
C This code provides also the GPDs in both the ERBL and DGLAP regions (these are not part of the output but they can be
C extracted from the code, as indicated)   
C
C INPUT
C          t       negative momentum transfer squared (GeV^2)
C          zeta    (SKEWNESS PARAMETER) connected to xi
C          scale   DIS momentum transfer squared (GeV^2)
C          jgpd    GPD type: 
C                  jgpd=1 --> H
C                  jgpd=2 --> E
C                  jgpd=3 --> H tilde
C                  jgpd=4 --> E tilde
C
C
C Careful: definitions of pseudoscalar form factors gpu and gpd  vary depending on whether you have a pion pole (in pi+) or not
C K+ has been calculated without pion pole 
C
***
*** OUTPUT IS CFFs OBTAINED FROM QUARK-DIQUARK BASED MODELs.
***
*** -- THE CHIRAL EVEN GPDs IN THE DGLAP REGION ARE CALCULATED AT A LOW Q^2 SCALE AND THEY
***    ARE EVOLVED TO THE DSIRED Q^2 IN THE SUBROUTINE struct_ev_gpd_tmd
*** -- NEXT THE ERBL PART IS CALCULATED Q^2 BY Q^2 FOLLOWING: 
***
C    J.~O.~Gonzalez-Hernandez, S.~Liuti, G.~R.~Goldstein and K.~Kathuria,
C   ``Interpretation of the Flavor Dependence of Nucleon Form Factors in a Generalized Parton Distribution Model,''
C    Phys.\ Rev.\ C {\bf 88}, 065206 (2013) [arXiv:1206.1876 [hep-ph]]. (PARAMETRIZATION II)
C 
C    G.~R.~Goldstein, J.~O.~Hernandez and S.~Liuti,
C  ``Flexible Parametrization of Generalized Parton Distributions from Deeply Virtual Compton Scattering Observables,''
C Phys.\ Rev.\ D {\bf 84}, 034007 (2011)[arXiv:1012.3776 [hep-ph]]. (PARAMETRIZATION I)
*** 
*** IN THE FINAL STEP THE COMPTON FORM FACTORS ARE CALCULATED 
***
      subroutine all_gpd_new(T,ZETA,SCALE,jgpd,HUIM,HDIM,HURE,
     > HDRE,DHUIM,DHDIM,DHURE,DHDRE)

      implicit none

      integer*4 i,num,ng,jgpd,ierr

      real*8 pi

      real*8 scale

      real*8 zeta,t,tm

      parameter(num=490,ng=48)
  
      real*8 q02(num),xbj(num),q2(num)

      real*8 uvi(num),dvi(num),ubi(num),dbi(num),sbi(num),glui(num)

      real*8 xjacobian,xstart,yyy,z

      real*8 bu,bd,f1u,f1d,f2u,f2d,gau,gad,areau,aread,
     > areau_m,aread_m

c      real*8 uvint,dvint

      real*8 huzz,hdzz,xuz

      real*8 asu,bsu,csu,asd,bsd,csd,aau,bau,cau,dau,aad,bad,cad,dad

      real*8 f1pp,ef1p,f1np,ef1n,f2p,ef2p,f2n,ef2n

      real*8 gpv,gps,gpu,gpd,gav,gas

      real*8 hu_erbl(num),hd_erbl(num),hus_erbl(num),hua_erbl(num),
     >  hds_erbl(num),hda_erbl(num),hubar_erbl(num),hdbar_erbl(num)

      real*8 hut(num),hdt(num)

      real*8 az,bz,s

      real*8 f(ng),xx(ng),fut,fdt,dgaus1,
     > fintu,hintu,fintd,hintd,fu1,fd1,fintu_tot,fintd_tot

      real*8 huim(2),hdim(2),hure(2),hdre(2)

      real*8 dhuim,dhdim,dhure,dhdre


      external fut,fdt

      data aau/2000./,aad/1000./ 

*	common /pdf/ npt,ngroup,nset
      common /parton/uvi,dvi,ubi,dbi,sbi,glui
      common /evol/xbj,q2,q02
      common/gpds/hu_erbl,hd_erbl,hus_erbl,hua_erbl,
     >     hds_erbl,hda_erbl,hut,hdt
  
************************************************************************
************************************************************************ 
******** GO!!!
************************************************************************
************************************************************************ 

**** positive t 
        tm=-t

      	pi = dacos(-1.d0)

************************************************************************

*** struct_ev puts NUM values of DGLAP distributions for each jgpd in COMMON parton 

***** LOOP OVER ERROR EVALUATION (CAN BE COMMMENTED OUT SETTING IERR=1)
*** ierr=1 central value
*** ierr=2 upper value 

        do 2000 ierr=1,2


*** DGLAP CALCULATION EVOLVED TO Q^2 
*** set itmd=1 (for GPD =1, for TMD moment =2)
*** set jt = 1 (type of TMD, NA here)  

           call struct_ev_gpd_tmd(scale,zeta,t,jgpd,ierr,1,1)

**************************** Check: u-quark,d-quark Baryon number
c      uvint=0.d0
c      DO I=1,NUM
c            if(xbj(i).lt.zeta)then
c            uvi(i)=0.d0
c            dvi(i)=0.d0
c            else
c        IF(I.LE.290)THEN
c             xjacobian = xbj(i)*LOG(1.E4)/330.
c             else
c          XSTART=EXP(-LOG(1.E4)*41./330.)
c            XJACOBIAN=(1.-XSTART)/201.
c            endif
c            uvint = uvint + uvi(i)/xbj(i)*xjacobian
c            dvint = dvint + dvi(i)/xbj(i)*xjacobian
c         endif
c               enddo
c         write(6,'(a,2x,2(g10.4,2x),3x,I3)')'Baryon Numbers ',
c     > uvint,dvint,jgpd
********************************************************************************************
******* ERBL
**
***** CALCULATION OF MISSING AREA AT SCALE Q^2 -> AREA IN THE ERBL REGION = bu,bd

        if(zeta.ne.0.) then

        bu=0.d0
        bd=0.d0
        DO I=1,NUM
        IF(I.LE.290)THEN
          YYY=LOG(1.E+4)*(330.-FLOAT(I)+1.)/330.
          Z=EXP(-YYY)
          XJACOBIAN=Z*LOG(1.E4)/330.
c  note that the small x region is expanded.
c  the range of bjorken x (here called Z) is 1.e-4 to 1.
        ELSE
          XSTART=EXP(-LOG(1.E4)*41./330.)
          Z=XSTART+(FLOAT(I)-290.)*(1.-XSTART)/201.
          XJACOBIAN=(1.-XSTART)/201.
        END IF
        
        if(z.gt.zeta) then

        bu = bu + uvi(i)/z*xjacobian
        bd = bd + dvi(i)/z*xjacobian

        else
           endif
          
        enddo

c         print *,t,tm,zeta,bu 
*** total area

*** Kelly's electromagnetic form factor parametrization
        call kelly(tm,f1pp,ef1p,f1np,ef1n,f2p,ef2p,f2n,ef2n)

*** axial form factor
        call ga(tm,gav,gas)

*** pseudo-scalar form factor (pion pole part needs to be fixed)
        call gp(tm,gpv,gps,gpu,gpd)

*** MISSING AREAS 

        go to (11,12,13,14),jgpd
 11     f1u = 2.*f1pp+f1np
        f1d = 2.*f1np+f1pp
        areau = f1u - bu/(1.-zeta/2.)
        aread = f1d - bd/(1.-zeta/2.)
        go to 20
 12     f2u = 2.*f2p+f2n
        f2d = 2.*f2n+f2p
        areau = f2u - bu/(1.-zeta/2.)
        aread = f2d - bd/(1.-zeta/2.)
        go to 20
 13     gau = (gav+gas)/2.
        gad = (-gav+gas)/2.
        areau = gau - bu/(1.-zeta/2.)
        aread = gad - bd/(1.-zeta/2.)
        go to 20
 14     continue

        areau = gpu - bu/(1.-zeta/2.)
        aread = gpd - bd/(1.-zeta/2.)

**** Notice factor 2. This is because areau is calculated for the "q" component, while what obeys
**** the symmetry property is the "+" component 
 20     areau_m = areau *(1.-zeta/2.)*2.
        aread_m = aread *(1.-zeta/2.)*2.

*** DETERMINE CROSSOVER POINT, also for Im part of CFF, huzz, hdzz

        do i=1,num-1
        if(uvi(i+1).ne.0.and.uvi(i).eq.0.d0)then

         huzz = uvi(i+1)/xbj(i+1)
         xuz = xbj(i+1)  !notice that we do not need to define xdz because it is the same
         hdzz = dvi(i+1)/xbj(i+1)

         else
            endif
            enddo

            do 300 i=1,num

            if(xbj(i).lt.xuz) then

*************** SYMMETRIC PART (FLAVOR NON SINGLET)
*** ==> u quark

         asu = 6.*(huzz*xuz-areau_m)/xuz**3
         bsu = -asu*xuz
         csu = huzz 
         hus_erbl(i) = asu*xbj(i)**2 + bsu*xbj(i) + csu

*** ==> d quark
         asd = 6.*(hdzz*xuz-aread_m)/xuz**3
         bsd = -asd*xuz
         csd = hdzz 
         hds_erbl(i) = asd*xbj(i)**2 + bsd*xbj(i) + csd

*************** ANTISYMMETRIC PART (FLAVOR SINGLET)

         bau = -aau*1.5*xuz
         cau = (2.*huzz + 0.5*aau*xuz**3)/xuz
         dau = - huzz
         hua_erbl(i) = aau*xbj(i)**3+bau*xbj(i)**2 + cau*xbj(i) + dau

*** ==> d quark
         bad = -aad*1.5*xuz
         cad = (2.*hdzz + 0.5*aad*xuz**3)/xuz
         dad = - hdzz
         hda_erbl(i) = aad*xbj(i)**3+bad*xbj(i)**2 + cad*xbj(i) + dad
*************************************************************************
**** Quark Distributions = (H^+ + H^-)/2
         hu_erbl(i) = (hus_erbl(i)+hua_erbl(i))/2.
         hd_erbl(i) = (hds_erbl(i)+hda_erbl(i))/2.
**** Anti-Quark Distributions = (H^+ - H^-)/2
         hubar_erbl(i) = (hus_erbl(i)-hua_erbl(i))/2.
         hdbar_erbl(i) = (hds_erbl(i)-hda_erbl(i))/2.

            else
          hu_erbl(i)=0.d0
          hd_erbl(i)=0.d0
          hua_erbl(i)=0.d0
          hda_erbl(i)=0.d0
          hus_erbl(i)=0.d0
          hds_erbl(i)=0.d0
           endif

**** QUARK Distribution in ERBL + DGLAP regions

           hut(i) = hu_erbl(i) + uvi(i)/xbj(i)
           hdt(i) = hd_erbl(i) + dvi(i)/xbj(i)
           
 300    continue

****** at this point hut and hdt vectors are filled for a number "num" of  X points covering both the 
****** ERBL and DGLAP regions, for each kinematics given by: zeta, Q^2 and t, and for each value of 
****** jg=H,E,Htilde 

***********************************************************************************
********
******** Now Compton Form Factors
********
***********************************************************************************
*** Real part of CFF: CAUCHY PV

***** Old method
*** standard done at eps=0.05
c      epsabs =1.0d0
c      epsrel=1.0d0
c      limit=4
c      eps=0.1d0
c      s=zeta
c      az=zeta/2.
c      bz=1.d0
c        fintu = dcauch(fut,az,bz,s,eps) 
c        fintd = dcauch(fdt,az,bz,s,eps) 
*********************
c      call dqawce(fut,az,bz,s,epsabs,epsrel,limit,result,abserr,neval,
c     *   ier,alist,blist,rlist,elist,iord,last)
c
c      fintu=result

***** Uses fut and fdt external functions that interpolate linearly between the points given by
***** struct_ev 
***** NOTICE: fut and fdt need the whole set of num points in X before being called
*****
      s=zeta
      az=zeta/2.
      bz=1.d0
      call dinter(xx,az,bz,ng)
      do i=1,ng
         f(i)= (fut(xx(i))-fut(s))/(xx(i)-s)
         enddo
         hintu = dgaus1(ng,az,bz,f)
         fintu = hintu + fut(s)*dlog((1.-zeta)/az) 
      do i=1,ng
         f(i)= (fdt(xx(i))-fdt(s))/(xx(i)-s)
         enddo
         hintd = dgaus1(ng,az,bz,f)
         fintd = hintd + fdt(s)*dlog((bz-s)/az) 

*********************
***   ADDITIONAL (NON SINGULAR) PART OF INTEGRAL FROM 1/X term

      fu1=0.d0
      fd1=0.d0
        DO I=1,NUM
        IF(I.LE.290)THEN
          YYY=LOG(1.E+4)*(330.-FLOAT(I)+1.)/330.
          Z=EXP(-YYY)
          XJACOBIAN=Z*LOG(1.E4)/330.
c  note that the small x region is expanded.
c  the range of bjorken x (here called Z) is 1.e-4 to 1.
        ELSE
          XSTART=EXP(-LOG(1.E4)*41./330.)
          Z=XSTART+(FLOAT(I)-290.)*(1.-XSTART)/201.
          XJACOBIAN=(1.-XSTART)/201.
        END IF

                if(z.gt.az) then

        fu1 = fu1 + (hus_erbl(i)+uvi(i)/xbj(i))/z*xjacobian
        fd1 = fd1 + (hds_erbl(i)+dvi(i)/xbj(i))/z*xjacobian
 
       else
           fu1=0.d0
           fd1=0.d0
           endif

        enddo
********************************************************************************************
********** FINAL Real CFFs
***** H and E are defined with C^+ = 1/(X-zeta) + 1/X
***** Htilde and Etilde are defined with C^-=1/(X-zeta) - 1/X
        if(jgpd.eq.3.or.jgpd.eq.4) then
      fintu_tot = fintu-fu1 
      fintd_tot = fintd-fd1
      else
      fintu_tot = fintu +fu1 
      fintd_tot = fintd +fd1

      endif

       else

          huzz=uvi(1)/xbj(1)
          hdzz=dvi(1)/xbj(1)
          fintu_tot=0.d0
          fintd_tot=0.d0

          endif

          huim(ierr) = pi*huzz
          hdim(ierr) = pi*hdzz

          hure(ierr) = fintu_tot
          hdre(ierr) = fintd_tot
2000     continue

          dhuim = abs(huim(2)-huim(1))
          dhdim = abs(hdim(2)-hdim(1))
          dhure = abs(hure(2)-hure(1))
          dhdre = abs(hdre(2)-hdre(1))
c          if(jgpd.eq.3) 
c        print *, huim(1),huim(2),dhuim,hdim(1),hdim(2),dhdim
         return
         end

********************************************************************************
********************************************************************************
**********************************************************************************************************
**********************************************************************************************************
**********************************************************************************************************
      double precision function fut(x)
C  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  + 
C-----------------------------------------------------------------------
C
C File - Version - Date: fut.f - v 1.0 - 4/14/2014
C                   
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C     Calculates the PV of the symmetric (-) d quark GPD 
C
C          x = X = initial quarks longitudinal momentum fraction
C          
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Author: Simonetta Liuti
C
C-----------------------------------------------------------------------
C
      implicit none 

      integer*4 num,i
      parameter(num=490)

      real*8 x
      
      real*8 uvi(num),dvi(num),ubi(num),dbi(num),sbi(num),glui(num)

      real*8 q02(num),xbj(num),q2(num)

      real*8 hu_erbl(num),hd_erbl(num),hus_erbl(num),hua_erbl(num),
     >  hds_erbl(num),hda_erbl(num),hut(num),hdt(num)

      real*8 hus(num)

      real*8 y1,y2,x1,x2,emme(num),enne(num)

      common /parton/uvi,dvi,ubi,dbi,sbi,glui
      common /evol/xbj,q2,q02
      common/gpds/hu_erbl,hd_erbl,hus_erbl,hua_erbl,
     >  hds_erbl,hda_erbl,hut,hdt

**** Linear interpolation of the function in the bins
**** given above: n1=490 bins. This generates my function f(x)
         do i=1,num
            hus(i) = hus_erbl(i)+uvi(i)/xbj(i) 
            enddo

         do i=1,num-1
            y1 = hus(i) 
            y2 = hus(i+1)
            x1 = xbj(i)
            x2 = xbj(i+1)

            emme(i) = (y1-y2)/(x1-x2)
            enne(i) = y1 - emme(i)*x1
            enddo
 
          do i=1,num-1
            if(x.ge.xbj(i).and.x.le.xbj(i+1)) then

c         f = -(emme(i)*x + enne(i))/x/(x-zeta)*(2.*x-zeta)
*** For dispersion relations
          fut = (emme(i)*x + enne(i))  !/(x-zeta1)

            else
            endif

            enddo

      return
      end
***********************************************************************************************************************
**********************************************************************************************************
      double precision function fdt(x)
C  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  + 
C-----------------------------------------------------------------------
C
C File - Version - Date: fdt.f - v 1.0 - 4/14/2014
C                   
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C     Calculates the PV of the symmetric (-) d quark GPD 
C
C          x = X = initial quarks longitudinal momentum fraction
C          
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Author: Simonetta Liuti
C
C-----------------------------------------------------------------------
C
      implicit none 

      integer*4 num,i
      parameter(num=490)

      real*8 x
      
      real*8 uvi(num),dvi(num),ubi(num),dbi(num),sbi(num),glui(num)

      real*8 q02(num),xbj(num),q2(num)

      real*8 hu_erbl(num),hd_erbl(num),hus_erbl(num),hua_erbl(num),
     >  hds_erbl(num),hda_erbl(num),hut(num),hdt(num)

      real*8 hds(num)

      real*8 y1,y2,x1,x2,emme(num),enne(num)

      common /parton/uvi,dvi,ubi,dbi,sbi,glui
      common /evol/xbj,q2,q02
      common/gpds/hu_erbl,hd_erbl,hus_erbl,hua_erbl,
     >  hds_erbl,hda_erbl,hut,hdt

**** Linear interpolation of the function in the bins
**** given above: n1=490 bins. This generates my function f(x)
         do i=1,num
            hds(i) = hds_erbl(i)+dvi(i)/xbj(i) 
            enddo

         do i=1,num-1
            y1 = hds(i) 
            y2 = hds(i+1)
            x1 = xbj(i)
            x2 = xbj(i+1)

            emme(i) = (y1-y2)/(x1-x2)
            enne(i) = y1 - emme(i)*x1
            enddo
 
          do i=1,num-1
            if(x.ge.xbj(i).and.x.le.xbj(i+1)) then

c         f = -(emme(i)*x + enne(i))/x/(x-zeta)*(2.*x-zeta)
*** For dispersion relations
          fdt = (emme(i)*x + enne(i))  !/(x-zeta1)

            else
            endif

            enddo

      return
      end




