C
C-----------------------------------------------------------------------
C  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  +  + 
C-----------------------------------------------------------------------
C
C File - Version - Date: all_gpd_main.f - v 1.0 - 8/11/2018
C                   
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
C
C Author: Simonetta Liuti
C
C-----------------------------------------------------------------------
C
C Main file that calls all_gpd_new.f
C


CCCCC HERE FOLLOWS THE DESCRIPTION OF VARIABLES IN all_gpd_new.f
CCCCC
C INPUT
C          t       negative momentum transfer squared (GeV^2)
C          zeta    longitudinal momentum transfer fraction
C          scale   DIS momnetum transfer squared (GeV^2)
C          jgpd    GPD type: 
C                  jgpd=1 --> H
C                  jgpd=2 --> E
C                  jgpd=3 --> H tilde
C                  jgpd=4 --> E tilde
C
***
*** OUTPUT: CFFs 
***
*** -- THE CHIRAL EVEN GPDs IN THE DGLAP REGION ARE CALCULATED AT A LOW Q^2 SCALE AND THEY
***    ARE EVOLVED TO THE DESIRED Q^2 IN THE SUBROUTINE struct_ev_gpd_tmd
*** -- THE ERBL PART IS CALCULATED FOLLOWING: 
***
C    J.~O.~Gonzalez-Hernandez, S.~Liuti, G.~R.~Goldstein and K.~Kathuria,
C   ``Interpretation of the Flavor Dependence of Nucleon Form Factors in a Generalized Parton Distribution Model,''
C    Phys.\ Rev.\ C {\bf 88}, 065206 (2013) [arXiv:1206.1876 [hep-ph]]. (PARAMETRIZATION II)
C 
C    G.~R.~Goldstein, J.~O.~Hernandez and S.~Liuti,
C  ``Flexible Parametrization of Generalized Parton Distributions from Deeply Virtual Compton Scattering Observables,''
C Phys.\ Rev.\ D {\bf 84}, 034007 (2011)[arXiv:1012.3776 [hep-ph]]. (PARAMETRIZATION I)
*** 
***
      program all_gpd_main

      implicit none

      integer*4 nkin,jgpd,jg,it
      parameter(nkin=4)

      real*8 pi

      real*8 scale

      real*8 zeta,t,dt,tm,xm

      real*8 q2loop(nkin),tloop(nkin),xloop(nkin)

      real*8 huim(2),hdim(2),hure(2),hdre(2)
! two values because of the error evaluation

      real*8 huim1(4),hdim1(4),hure1(4),hdre1(4)

      real*8 dhuim,dhdim,dhure,dhdre


*** possible kinematic choices
      data tloop/0.35, 0.1, 0.438, 1.133/
      data xloop/0.01, 0.2, 0.266, 0.262/
      data q2loop/2.000, 2.00, 2.469, 2.436/

*****
      PI = DACOS(-1.D+00)

        scale = q2loop(1)
        t = -tloop(1)
        zeta = xloop(1)

        write (6,'(//,3(a,f9.3,/))')
     > 't (GeV)',t,'x_Bj ',zeta,'Q^2 ',scale

        do jg=1,4

          if(jg.eq.1) write(6,'(8x,5(a,7x))')
     > 't (GeV)','Im H_u','Re H_u','Im H_d','Re H_d'
          if(jg.eq.2) write(6,'(8x,4(a,7x))')
     > 't (GeV)','Im E_u','Re E_u','Im E_d','Re E_d'
          if(jg.eq.3) write(6,'(8x,4(a,7x))')
     > 't (GeV)','Im Ht_u','Re Ht_u','Im Ht_d','Re Ht_d'
          if(jg.eq.4) write(6,'(8x,4(a,7x))')
     >  't (GeV)','Im Et_u','Re Et_u','Im Et_d','Re Et_d'
 
*** check in one kinematics
         jgpd=jg 
         tm = 0.35
         xm = 0.3 ! positive t
         dt = 0.05
        do it=1,5
           zeta = xm
           t = -tm
        call all_gpd_new(T,ZETA,SCALE,jgpd,HUIM,HDIM,HURE,
     > HDRE,DHUIM,DHDIM,DHURE,DHDRE)

        huim1(jg)=huim(1)
        hdim1(jg)=hdim(1)
        hure1(jg)=hure(1)
        hdre1(jg)=hdre(1)

        write(6,'(i3,3x,5(g11.4,2x))')
     > jg,t,2./3.*huim1(jg)-1./3*hdim1(jg),2./3.*hure1(jg)
     >-1./3.*hdre1(jg),0,0
C  jg,t,2./3.*huim1(jg)-1./3.*hdim1(jg),2./3.*hure1(jg)-1./3.*hdre1(jg)
        xm = xm + dt
        enddo

        enddo

        stop
        end

