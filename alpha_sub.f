	DOUBLE PRECISION FUNCTION ALPHA_MY(SCALE)

****
**** CALCULATES PI/ALPHA using W. Marciano's prescription
**** RETURNS ALPHA_MY=ALPHA/PI
**** SCALE=Q^2 in GeV^2
**** 
	implicit none

	integer*4 nflavor,imass

	real*8 pi,xlam,xlam3,xlam5

	real*8 xmass(3),xmas2(3),scale

	real*8 x3,x5,x6,add3,add5,add6

	real*8 beta0,beta1

	real*8 x

	real*8 alphainv

	data xmass/1.4d0, 4.5d0, 170./

      PI=4.*ATAN(1.)

**** REFERENCE SCALE XLAM(NFLAVOR=4)
****
       XLAM=.215d0 !Nf=4
c	XLAM = 0.176
c	XLAM = 0.158
	do imass=1,3
	xmas2(imass)=xmass(imass)*xmass(imass)
	enddo

	if(scale.le.xmas2(1)) imass=1
	if(scale.gt.xmas2(1).and.scale.le.xmas2(2)) imass=2
	if(scale.gt.xmas2(2).and.scale.le.xmas2(3)) imass=3
	if(scale.gt.xmas2(3)) imass=4

	go to (110,120,130,140),imass

 110	NFLAVOR=3
	xlam3 = xlam*(xmass(1)/xlam)**(2./27.)
c     >         *(dlog((xmass(1)/xlam)**2))**(107./2025.)

*** additive procedure
c	x3 = dlog(xmas2(1)/xlam/xlam)
c        add3 =   25./12.*x3 /(1. - 462./625.*dlog(X3)/X3 )
c     1  - 27./12.*x3 /(1. - 576./729.*dlog(X3)/X3 )
c        add3 =   25./12.*x3 - 27./12.*x3 
c        X=DLOG(scale/XLAM/XLAM)

	beta0 = 11.d0 - float(nflavor)*2./3.
	beta1 = 51.d0 - float(nflavor)*19./3.
c	beta1 = 102.d0 - float(nflavor)*38./3.

        X=DLOG(scale/XLAM3/XLAM3)

	ALPHAINV  = (beta0/4.)*x!/(1. - 2.*beta1/beta0**2*dlog(X)/X )
c     > + add3

c	ALPHAINV = (beta0*x/4.)*((1 - ((beta1*dlog(x))/(beta0**2*x)))**-1)
	
c	print 600, scale,sqrt(scale),xlam3,x,pi/alphainv,
c     > 4./beta0/x,beta0,beta1,add3

	go to 100


 120	NFLAVOR=4

	beta0 = 11.d0 - float(nflavor)*2./3.
	beta1 = 51.d0 - float(nflavor)*19./3.
c	beta1 =	102.d0 - float(nflavor)*38./3.

        X=DLOG(scale/XLAM/XLAM)
c	ALPHAINV = (4*(1/(beta0*X)) - (8*beta1/beta0)*DLOG(X)/(beta0*X)**2)**-1.
	alphainv = (beta0/4.)*x!/(1. - 2.*beta1/beta0**2*dlog(X)/X )
c	ALPHAINV = (beta0*x/4.)*((1 - ((beta1*dlog(x))/(beta0**2*x)))**-1)

c	print 600, scale,sqrt(scale),xlam,x,pi/alphainv,4./beta0/x

	go to 100

 130	NFLAVOR=5
	xlam5 = xlam*(xmass(2)/xlam)**(-2./23.)
c     >         *(dlog((xmass(2)/xlam)**2))**(-963./13225.)

** additive procedure
c	x5 = dlog(xmas2(2)/xlam/xlam)
c        add5 =  25./12.*x5 /(1. - 462./625.*dlog(X5)/X5 )
c     1  - 23./12.*x5 /(1. - 348./529.*dlog(X5)/X5 )
c        add5 =  25./12.*x5 - 23./12.*x5 
c        X=DLOG(scale/XLAM/XLAM)

	beta0 = 11.d0 - float(nflavor)*2./3.
	beta1 = 51.d0 - float(nflavor)*19./3.
c	beta1 =	102.d0 - float(nflavor)*38./3.

        x=DLOG(scale/XLAM5/XLAM5)
c	ALPHAINV = (4*(1/(beta0*X)) - (8*beta1/beta0)*DLOG(X)/(beta0*X)**2)**-1.
       ALPHAINV  = (beta0/4.)*x!/(1. - 2.*beta1/beta0**2*dlog(X)/X )
c	ALPHAINV = (beta0*x/4.)*((1 - ((beta1*dlog(x))/(beta0**2*x)))**-1)
c     >  + add5

c	print 600, sqrt(scale),xlam5,x,dacos(-1.d0)/alphainv

	go to 100

 140	NFLAVOR=6
	beta0 = 11.d0 - float(nflavor)*2./3.
	beta1 = 51.d0 - float(nflavor)*19./3.

	x6 = dlog(xmass(3)/xlam/xlam)
        add6 =  23./12.*x6 /(1. - 348./529.*dlog(X6)/X6 )
     1  - 21./12.*x6 /(1. - 234./441.*dlog(X6)/X6 )

	beta0 = 11.d0 - float(nflavor)*2./3.
	beta1 = 51.d0 - float(nflavor)*19./3.

        X=DLOG(scale/XLAM/XLAM)
	alphainv= (beta0/4.)*x! /(1. - 2.*beta1/beta0**2*dlog(X)/X )
c     >  + add6

 100	continue


*****
*****
***** alpha/pi

        ALPHA_MY=1./ALPHAINV

c	print *,scale,x,xlam**2,beta0,alphainv,alpha_my
****
****
****
c 600	format(1x,10(g10.4,1x))
	return 
	end
