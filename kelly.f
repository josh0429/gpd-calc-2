      subroutine kelly(t,f1p,ef1p,f1n,ef1n,f2p,ef2p,f2n,ef2n)

*** t is positive and in GeV^2

      implicit real *8(a-h,o-z)

!      parameter(nd=179)

*** parameters
      data a1_ep/-0.24/,a1_mp/0.12/,a1_mn/2.33/
      data b1_ep/10.98/,b1_mp/10.97/,b1_mn/14.72/
      data b2_ep/12.82/,b2_mp/18.86/,b2_mn/24.20/
      data b3_ep/21.97/,b3_mp/6.55/,b3_mn/84.1/

      data ea1ep/0.12/,ea1mp/0.04/,ea1mn/1.4/
      data eb1ep/0.19/,eb1mp/0.11/,eb1mn/1.7/
      data eb2ep/1.10/,eb2mp/0.28/,eb2mn/9.8/
      data eb3ep/6.80/,eb3mp/1.20/,eb3mn/41./
      
      data aa/1.7/,bb/3.30/
      data eaa/0.04/,ebb/0.32/
** Walcher: aa=1.73 bb=4.59
      data amp/0.938/
      data amup/2.79/,amun/-1.91/

c       do t=0.,11.,0.2
c      write(6,*)'t value in GeV2 '
c      read (5,*)t

          tau=t/4./amp/amp

          gep= (1.+a1_ep*tau)/(1.+b1_ep*tau+b2_ep*tau**2+b3_ep*tau**3)

          gmp= amup*
     > (1.+a1_mp*tau)/(1.+b1_mp*tau+b2_mp*tau**2+b3_mp*tau**3)

          gmn= amun*
     > (1.+a1_mn*tau)/(1.+b1_mn*tau+b2_mn*tau**2+b3_mn*tau**3)

          gen = aa*tau/(1.+tau*bb)/(1.+t/0.71)**2


**** errors on gmp,gep,gmn
          e1 = tau/(1.+b1_ep*tau+b2_ep*tau**2+b3_ep*tau**3)
          e2 =(1.+a1_ep*tau)/(1.+b1_ep*tau+b2_ep*tau**2+b3_ep*tau**3)**2
     >        *tau
          e3 =(1.+a1_ep*tau)/(1.+b1_ep*tau+b2_ep*tau**2+b3_ep*tau**3)**2
     >        *tau**2
          e4 =(1.+a1_ep*tau)/(1.+b1_ep*tau+b2_ep*tau**2+b3_ep*tau**3)**2
     >        *tau**3
          egep = sqrt((e1*ea1ep)**2+(e2*eb1ep)**2+(e3*eb2ep)**2
     >            +(e4*eb3ep)**2)
***
          e1 = tau/(1.+b1_mp*tau+b2_mp*tau**2+b3_mp*tau**3)
          e2 =(1.+a1_mp*tau)/(1.+b1_mp*tau+b2_mp*tau**2+b3_mp*tau**3)**2
     >        *tau
          e3 =(1.+a1_mp*tau)/(1.+b1_mp*tau+b2_mp*tau**2+b3_mp*tau**3)**2
     >        *tau**2
          e4 =(1.+a1_mp*tau)/(1.+b1_mp*tau+b2_mp*tau**2+b3_mp*tau**3)**2
     >        *tau**3
          egmp = sqrt((e1*ea1mp)**2+(e2*eb1mp)**2+(e3*eb2mp)**2
     >            +(e4*eb3mp)**2)
***
          e1 = tau/(1.+b1_mn*tau+b2_mn*tau**2+b3_mn*tau**3)
          e2 =(1.+a1_mn*tau)/(1.+b1_mn*tau+b2_mn*tau**2+b3_mn*tau**3)**2
     >        *tau
          e3 =(1.+a1_mn*tau)/(1.+b1_mn*tau+b2_mn*tau**2+b3_mn*tau**3)**2
     >        *tau**2
          e4 =(1.+a1_mn*tau)/(1.+b1_mn*tau+b2_mn*tau**2+b3_mn*tau**3)**2
     >        *tau**3
          egmn = sqrt((e1*ea1mn)**2+(e2*eb1mn)**2+(e3*eb2mn)**2
     >            +(e4*eb3mn)**2)

*** error on Gen
          e1 = tau/(1.+tau*bb)/(1.+t/0.71)**2
          e2 = aa*tau*tau/(1.+tau*bb)**2/(1.+t/0.71)**2
          egen = sqrt((e1*eaa)**2+(e2*ebb)**2)

***** 
***** DIRAC AND PAULI
          f1p = (tau*gmp+gep)/(1.+tau)
          f2p = (gmp-gep)/(1.+tau)  !/(amup-1.)

          f1n = (tau*gmn+gen)/(1.+tau)
          f2n = (gmn-gen)/(1.+tau)  !/(amun-1.)
c          f2n = (gmn-gen)/(1.+tau)/(amun)


          gd = 1./(1.+t/0.71)**2
**** error          
          ef1p = sqrt((tau*egmp)**2+ egep**2)/(1.+tau) 
          ef2p = sqrt(egmp**2+ egep**2)/(1.+tau) 
          ef1n = sqrt((tau*egmn)**2+ egen**2)/(1.+tau) 
          ef2n = sqrt(egmn**2+ egen**2)/(1.+tau) 

c      write(6,'(6(1x,(g10.4,1x)))') !fintu,fintd,fintu0,fintd0,
c     > xju,dint_u,xjd,dint_d
c     > t,gep,gmp,gmn,gen
c     > t,f1p,f2p,f1n,f2n,sqrt(t)*f2p/f1p

c      write(8,'(6(1x,(g10.4,1x)))') !fintu,fintd,fintu0,fintd0,
c     > xju,dint_u,xjd,dint_d
c     > t,gep/gd,gmp/gd,gmn/gd,gen/gd
c     > t,f1p,f2p,f1n,f2n,sqrt(t)*f2p/f1p
c          enddo
          
c 600     format(1x,10(g10.4,1x))

          return
            end

