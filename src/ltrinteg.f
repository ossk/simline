c     **************************************************************
c     Subroutines to carry out integration of the intensity 
c     convolved with the beam function for a given central beam
c     displacement
c     The background intensity is removed here.
c     **************************************************************
      subroutine intbeam(p,iout,ubg,iconv,nnp,nff,deltap,sigbeam)
      implicit none
      integer nnp,nff

      include 'fsizes.inc'
      real p(nstretch),iout(nstretch,-nf:nf),iconv(-nf:nf)
      real deltap,sigbeam,ubg

      real zz(nstretch),uf(nstretch),iff(nstretch),skip,pskip
      real plow,pup,sigsqu,dintstepwise,specerr
      integer i,j,imin,imax,iskip1,iskip2

      real negexp,negint
      common /accgauss/ negexp,negint
      save /accgauss/

c     Criterion for narrow intervals (linear integration): 1%
      data skip /0.02/

c     **************************************************************
c     Compute the integration limits in p
c     **************************************************************
      sigsqu=sigbeam**2
      plow=max(0.0,deltap-negint*sigbeam)
      pup=deltap+negint*sigbeam
      do 200 i=2,nnp
      if (p(i).gt.plow) then
      imin=i-1
      goto 250
      endif
200   continue
      imin=nnp-1
250   do 280 i=imin,nnp
      if (p(i).ge.pup) then
      imax=i
      goto 290
      endif
280   continue
      imax=nnp
c     **************************************************************
c     Check for pencil beams - extend the integration range
c     **************************************************************
290   if (imax-imin.lt.2) then
       if (imin.gt.1) then
        imin=imin-1
        goto 290
       endif
       if (imax.lt.nnp) then
        imax=imax+1
        goto 290
       endif
      endif
c     **************************************************************
c     Set up the scales for the numeric integration routine
c     **************************************************************
      uf(imin)=2.0*p(imin)*exp(-(p(imin)-deltap)**2/sigsqu)
     %         *specerr(p(imin)*deltap/sigsqu)
      do 300 i=imin+1,imax
      zz(i)=p(i)-p(i-1)
      uf(i)=2.0*p(i)*exp(-(p(i)-deltap)**2/sigsqu)
     %         *specerr(p(i)*deltap/sigsqu)
300   continue
c     **************************************************************
c     In case of extreme steps, the cubic spline integration may
c     go bananas. At the edges, it is typically still fine, but in
c     the range it may be necessary to replace it by trapezium rule.
c     This is opposite to uzint.
c     **************************************************************
c     Find the small-step region
      pskip=skip*2.0*sigbeam
      iskip1=imax
      do 340 i=imin+1,imax
      if (zz(i).lt.pskip) then
      iskip1=i-1
      goto 360
      endif
340   continue
360   if (iskip1.eq.imax) then
c     no linear integration range needed
      iskip2=imin-1
      else
      iskip2=imin
      do 350 i=imax,iskip1+1,-1
      if (zz(i).lt.pskip) then
      iskip2=i
      goto 370
      endif
350   continue
370   continue
      endif
c     check outcome
      if (iskip1.lt.imin+3) iskip1=imin
      if (iskip2.gt.imax-3) iskip2=imax
c     **************************************************************
c     Loop for the frequencies to integrate the intensities
c     **************************************************************
      do 100 j=-nff,nff
      do 150 i=imin,imax
      iff(i)=uf(i)*(iout(i,j)-ubg)
150   continue
      iconv(j)=dintstepwise(zz,iff,imin,imax,iskip1,iskip2)
100   continue
      return
      end

c     ***************************************************************
c     Function to integrate stepwise through cubic splines and 
c     trapezium rule depending on given limits
c     ***************************************************************
      real function dintstepwise(x,y,n1,n2,nc1,nc2)
      implicit none
      real x(*),y(*),z
      real dintcub,dintlin
      integer n1,n2,nc1,nc2

c     integration logic
      if (nc1.ge.n2) then
c       no steps: only cubic spline
        z=dintcub(x,y,n1,n2)
      else 
        if (nc1.gt.n1) then
c         start with spline
          z=dintcub(x,y,n1,nc1)+dintlin(x,y,nc1,nc2)
        else 
c         start with trapezium
          z=dintlin(x,y,n1,nc2)
        endif
        if (nc2.lt.n2) then
c         end with spline
          z=z+dintcub(x,y,nc2,n2)
        endif
      endif

      dintstepwise=z
      return
      end

c     ***************************************************************
c     Function for approximative values of the integral
c     integrate(exp(-2*g*(1-cos(x)))dx) from 0 to pi
c     depending on the value of g
c     ***************************************************************
      real function specerr(g)
      implicit none
      real g
      real spi,pi,a,c

c     Sqrt(pi)/2,pi
      data spi,pi/ 0.886227,3.141593/
c     Coefficients of the power fit
      data a,c/ 0.96921,-0.55247/

c     Three different approximations g<0.3, 0.3<g<4, 4<g
      if (g.ge.4.0) then
      specerr=spi/sqrt(g)
      else if (g.ge.0.3027) then
      specerr= a*g**c
      else
      specerr=pi*(1.0-2.0*g+3.0*g**2-3.33333*g**3+2.91667*g**4
     %   -2.1*g**5)
      endif
      return
      end

c     **************************************************************
c     Subroutines to carry out integration of the intensity for
c     the energy density (J=2*pi/c*int(int(I(ny,z)*phi(ny,z),dny),dz)
c     **************************************************************
c     First routine - integration over the frequencies for a given 
c     z value and radius
c     Carried out within the radiative tranfer loop
c     The number of frequency points has to be restricted to <= 512
c     **************************************************************
      subroutine unyint(iplus,uny,lin,vs,sig,dny)
      implicit none
      integer lin

      include 'fsizes.inc'
      real iplus(-nf:nf,ntra),uny(ntra)
      real vs,sig,phi,quintgauss,dny
      real ff(2*nf),prof(-nf:nf)
      integer j,l,ii,j1,j2,nj

      real negexp,negint
      common /accgauss/ negexp,negint
      save /accgauss/

c     Frequency range
      j1=int((vs-sig*negint)/dny)-1
      j2=int((vs+sig*negint)/dny)+1
      nj=j2-j1+1
c     initialize profile parameters
      do 100 j=j1,j2
       prof(j)=phi(j*dny,vs,sig)
100   continue
c     Treat the different lines
      do 200 l=1,lin
      ii=0
      do 250 j=j1,j2
      ii=ii+1
       ff(ii)=iplus(j,l)*prof(j)
250   continue
      uny(l)=quintgauss(ff,dny,nj)
200   continue
      return
      end

c     **************************************************************
c     Second routine - integration over z for a given radius
c     Attention - a prefactor 4*pi/c is shifted to the Einstein 
c     coefficients to use the same energy density units as Kruegel
c     **************************************************************
      subroutine uzint(uny,z,nz,ulev,nn,lin)
      implicit none
      integer nn,lin

      include 'fsizes.inc'
      real uny(ns,2*nstretch,ntra),z(ns,2*nstretch),ulev(ns,ntra)
      integer nz(ns)
      real zz(2*nstretch),uz(2*nstretch)
      real pre,dpre,skip,dintcub,dintlin
      integer i,j,k,ninz,kskip,kend

c     Prefactor from the translation of a phi- into a z-integration
c     Normally pre=4*pi/c, here pre=1
      data pre /1.0/
c     Criterion for narrow intervals (linear integration): 1%
      data skip /0.02/

c     Most inner point
      do 400 j=1,lin
      ulev(1,j)=pre*uny(1,1,j)
400   continue
c     **************************************************************
c     Radius loop - set up the z scale for a given r
c     **************************************************************
      do 100 i=2,nn
      ninz=nz(i)
      dpre=pre/(z(i,ninz)-z(i,1))
      do 120 k=1,ninz
      zz(k)=dpre*z(i,k)
120   continue
c     Make the field of differences (separately for simpler indices)
      do 130 k=ninz,2,-1
      zz(k)=zz(k)-zz(k-1)
130   continue
      zz(1)=0.0
c     **************************************************************
c     Due to the concentration of the points at the intervall ends
c     for large radii, numerical instabilities in the cubic integration
c     may appear in extreme cases. In these regions of dense points
c     the linear integration will be taken.
c     **************************************************************
c     Find the region separation
      kskip=1
      do 140 k=2,ninz
      if (zz(k).gt.skip) then
      kskip=k-1
      goto 160
      endif
140   continue
160   kend=ninz-kskip+1
c     ***************************************************************
c     Level loop - set up the function vector to be integrated  
c     Cubic integration is always executed
c     ***************************************************************
      do 500 j=1,lin
      do 530 k=1,ninz
       uz(k)=uny(i,k,j)
530   continue
      ulev(i,j)=dintcub(zz,uz,kskip,kend)
c     Linear integration part
      if (kskip.gt.1) ulev(i,j)=ulev(i,j)
     %   +dintlin(zz,uz,1,kskip)+dintlin(zz,uz,kend,ninz)
500   continue
c     ***************************************************************
c     End of radius loop
c     ***************************************************************
100   continue
      return
      end

c     *************************************************************
c     Integration routine for equidistant smooth tabulated functions
c     I have found that the most simple trapezoidal rule is 
c     in general sufficient. Third order rules change the
c     weighting at the ends which plays no role here (->0).
c     More sophisticated algorithms take more time than they win.
c     *************************************************************
      real function quintgauss(y,dx,nn)
      implicit none
      real y(*),dx
      real*8 z
      integer nn,i

      z=0.5*(y(1)+y(nn))
      do 200 i=2,nn-1
      z=z+y(i)
200   continue
      quintgauss=z*dx
      return
      end

C     ***********************************************************
C     DINTCUB IS A CUBIC SPLINE INTEGRATION PROGRAM
C     X has to be the field of differences, not of abscissa values 
C     ***********************************************************
      REAL FUNCTION DINTCUB(X,Y,N1,N2)
      IMPLICIT NONE
      REAL X(*),Y(*)
      REAL D(4095)
      REAL*8 A1,A2,A3,A4,Z,CONT,H2
      INTEGER N1,N2,I

C     SPLINE INTERPOLATION
      CALL SPLINE(X,Y,N1,N2,D)
C     COEFFICIENTS OF CUBIC FUNCTIONS
      Z=0.0D0
      DO 200 I=N1+1,N2
      H2=X(I)
        A1=((D(I-1)+D(I))*H2-2.0*(Y(I)-Y(I-1)))/H2**3
        A2=(D(I)-D(I-1)-3.0*A1*H2**2)/(2.0*H2)
        A3=D(I-1)
        A4=Y(I-1)
C     STEPWISE ANALYTIC INTEGRATION
      CONT=A1*0.25D0*H2**4+A2/3.0D0*H2**3+A3*0.5D0*H2**2+A4*H2
      Z=Z+CONT
200   CONTINUE
      DINTCUB=Z
      RETURN
      END

C     ***********************************************************
C     DINTLIN IS A LINEAR INTEGRATION PROGRAM
C     ***********************************************************
      REAL FUNCTION DINTLIN(X,Y,N1,N2)
      IMPLICIT NONE
      REAL X(*),Y(*)
      REAL*8 Z,A
      INTEGER N1,N2,I

      Z=0.0d0
      DO 200 I=N1+1,N2
      A=0.5*(Y(I)+Y(I-1))
      Z=Z+A*X(I)
200   CONTINUE
      DINTLIN=Z
      RETURN
      END

C     ***********************************************************
C     SPLINE INTERPOLATION
C     MODIFIED FROM NUMERICAL RECIPIES
C     X(N),Y(N) are the input parameters characterizing
C     the function to be considered
C     X(N) is the field of differences between grid points and
C     not the field of the grid points themselves
C     Y2(N) is the output giving the first derivative of the
C     function at the points X(N)
C     ***********************************************************
      SUBROUTINE SPLINE(X,Y,N1,N2,Y2)
      IMPLICIT NONE
      REAL X(*),Y(*),Y2(*)
      REAL U(4095),SIG,P,YN
      INTEGER I,N1,N2,K

        Y2(N1)=0.
        U(N1)=0.
      DO 11 I=N1+1,N2-1
        SIG=X(I)/(X(I+1)+X(I))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/X(I+1)-(Y(I)-Y(I-1))
     *      /X(I))/(X(I+1)+X(I))-SIG*U(I-1))/P
11    CONTINUE
      Y2(N2)=0.0
      DO 12 K=N2-1,N1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      YN=Y2(N2-1)
      DO 13 I=N1,N2-1
        Y2(I)=(Y(I+1)-Y(I))/X(I+1)-X(I+1)*(2.0*Y2(I)+Y2(I+1))/6.0
13    CONTINUE
      Y2(N2)=(Y(N2)-Y(N2-1))/X(N2)+YN*X(N2)/6.0
      RETURN
      END
