c     **********************************************************
c     Subroutine to estimate the level populations from the
c     velocity profile in Sobolev approximation
c     An appropriate first guess has to be given !
c     **********************************************************
      subroutine sobolev(r,vs,vt,rc,t,dh,dco,uext,dd,nn,lin,lev)
      implicit none
      integer nn,lev,lin

      include 'fsizes.inc'
      real r(ns),vs(ns),vt(ns),t(ns),dh(ns),dco(ns),rc(ns)
      real*8 dd(ns,nlev),dlocal(nlev),zerolev
      real uext(ntra),uu(ntra),abu(ntra),abu1(ntra)
      real uout(ntra)
      real vh(ns),vr,dvdr,rcorr,ed
      real cap,dcap,dcap1,s,ds,ds1,beta,abbeta,datau,rerr
      real*8 conv,neglev
      integer itermax,iter,i,l,levtmp,redlev
      logical hii
      common /acclev/ conv,neglev
      save /acclev/
        data zerolev /1.0d-14/

c     **********************************************************
c     Internal formats
c     **********************************************************
11    format(/,' The Sobolev approximation did not reach ',
     %  'convergence at radial point: ',I3)
12    format(' Maximum number of iterations in Sobolev ',
     %  'approximation: ',I3)
c     **********************************************************
c     Initializations of constant quantities
c     **********************************************************
      itermax=1
c     Initialize points for the angular integration
      call gauleg
c     dimensionless velocities
      do 150 i=2,nn
       vh(i)=vt(i)/r(nn)
150   continue
c     Total external intensities from background
      do 160 l=1,lin
       uout(l)=uext(l)
160   continue
c     Save and restore central vs in case of an HII region
      hii = (vt(1).ge.1.0)
      ed=vs(1)
      vs(1)=vs(2)

c     **********************************************************
c     Outer loop - radial points
c     **********************************************************
      do 100 i=2,nn
c     Fixed quantities
      rcorr=rc(i)
      if (hii) call ufromhii(r(i),r(2),uext,uout,lin)
      do 300 l=1,lev
       dlocal(l)=dd(i,l)
300   continue

c     Velocity gradients are only approximate (not analytical)
      if (i.eq.nn) then
      dvdr=(vs(nn)-vs(nn-1))/(r(nn)-r(nn-1))
      vr=vs(nn)/r(nn)
      else
      dvdr=(vs(i+1)-vs(i-1))/(r(i+1)-r(i-1))
      vr=vs(i)/r(i)
      endif
c     **********************************************************
c     To stay at least a bit physical, the velocity should not
c     fall below a certain limit in an LVG regime
c     Here, we use vturb/rmax as this limit
c     **********************************************************
      if(abs(dvdr).lt.vh(i)) dvdr=vh(i)

c     **********************************************************
c     Iteration to find the local solution
c     **********************************************************
      iter=0
      redlev=0
      levtmp=lev
400   continue

c     **********************************************************
c     Loop for the levels
c     Compute kappa and S, reduce kappa for turbulent clumping,
c     integrate finally the energy densities
c     **********************************************************
      do 450 l=1,lin
      call sobkapsrc(dco(i),dlocal,cap,dcap,dcap1,s,ds,ds1,l)
      call sobeffkap(cap,datau,rcorr)

      call gaussob(vr,dvdr,cap,beta,abbeta)
      uu(l)  =(1.0-beta)*s + beta*uout(l)
      abu(l) =(1.0-beta)*ds + (abbeta*datau*dcap)*(uout(l)-s)
      abu1(l)=(1.0-beta)*ds1 + (abbeta*datau*dcap1)*(uout(l)-s)
450   continue
c     **********************************************************
c     Set up the matrix and solve the system 
c     Look for convergence
c     **********************************************************
      call sobmatrix(t(i),dh(i),uu,abu,abu1,dlocal,lin,levtmp,rerr)
      iter=iter+1
      if (iter.ge.100) then
      write(6,11) i
      goto 500
      endif
c     ************************************************************
c     Reduce the level number temporarily if necessary
c     Neglect level if unpopulated three subsequent times
c     ************************************************************
      if (dlocal(levtmp).le.zerolev) then
      redlev=redlev+1
      if (redlev.ge.3) then
      levtmp=levtmp-1
      redlev=0
      endif
      else
      redlev=0
      endif
      if (rerr.ge.conv) goto 400
c     ************************************************************
c     End of the iteration, assign the local values
c     to the global field of level populations
c     ************************************************************
500   if (iter.gt.itermax) itermax=iter
      do 550 l=1,lev
      dd(i,l)=dlocal(l)
550   continue
c     ************************************************************
100   continue
      write(6,12) itermax
c     ************************************************************
c     It makes no sense to compute Sobolev for the most inner point
c     The values from the second point are taken there
c     ************************************************************
      do 200 l=1,lev
       dd(1,l)=dd(2,l)
200   continue
      vs(1)=ed
      return
      end

c     *************************************************************
c     Subroutine to carry out the mu-integration for the
c     Sobolev approximation, beta and dbeta/dkappa are treated
c     The points and weightings from gauleg are used.
c     The symmetry of the escape probability at 0 is expoited.
c     *************************************************************
      subroutine gaussob(vr,dvdr,kappa,y1,y2)
      implicit none
      real vr,dvdr,kappa,y1,y2

      integer m,mmax
      parameter(mmax=50)
      real x(mmax),w(mmax),x1,x2
      integer i

      common /gausspt/m,x,w
      save /gausspt/

      y1=0.0
      y2=0.0 
      do  500  i = 1, m
      call escape(vr,dvdr,kappa, x(i), x1, x2)
      y1=y1+w(i)*x1
      y2=y2+w(i)*x2
500   continue
c     The division by two is already done by taking the half range
      return
      end

c     ***************************************************************
c     Subroutine to compute the escape probability and its 
c     derivative from a given point into a certain direction
c     beta = (1-exp(-kappa/q)) / (kappa/q)
c     q = abs( mu^2 * dvdr + (1-mu^2) * v/r )
c     vr and dvdr have signs, q has none
c     abbeta = dbeta/dkappa
c     ***************************************************************
      subroutine escape(vr,dvdr,kappa,mu,beta,abbeta)
      implicit none
      real vr,dvdr,kappa,mu,beta,abbeta
      real q,tau,exptau

      q = abs(mu**2*(dvdr-vr) + vr)
      tau=kappa/q

c     ***************************************************************
c     avoid the accuracy problems for small tau
c     ***************************************************************
      if (tau.lt.0.01) then
      beta = 1.0 - tau*0.5
      abbeta = -0.5/q
      else
      exptau = exp(-tau)
      beta = (1.0-exptau) / tau
      abbeta = (exptau*(1.0+tau)-1.0) / (q*tau**2)
      endif
c     ***************************************************************
      return
      end

c     ************************************************************
c     Compute the set of points for the Gaussian integration
c     using a given accuracy
c     This subroutine is taken from "Numerical recipies", p. 125
c     and modified to use the fixed interval [-1,1] and only
c     half the set of grid points (>0)
c     ************************************************************
      SUBROUTINE gauleg
      implicit none
      integer m,mmax

      parameter(mmax=50)
      real x(mmax),w(mmax)
      integer n,i,j
      real*8 p1,p2,p3,pp,z,z1,EPS
      PARAMETER (EPS=3.d-14)

      real deltav,epsz
      common /acctrans/ deltav,epsz
      common /gausspt/ m,x,w
      save /acctrans/,/gausspt/

c     Approximate the number of points required for the integration
      m=min(int(2/epsz),mmax)
      n=2*m

c     The Gauss-Legendre loop from "Numerical recipies"
       do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if (abs(z-z1).gt.EPS) goto 1
        x(m-i+1)=z
        w(m-i+1)=2.d0/((1.d0-z*z)*pp*pp)
12     continue

      return
      end      

c     *************************************************************
c     Subroutine to compute the derivatives of the absorption 
c     and the source function from the level populations
c     at one radial point
c     *************************************************************
      subroutine sobkapsrc(dco,dd,cap,dcap,dcap1,s,ds,ds1,l)
      implicit none
      integer l

      include 'fsizes.inc'
      real*8 dd(nlev)
      real*8 bij,bji,b1,b2,rat,mrat,hep,zerolev
        real lgth,cap,dcap,dcap1
      real dco,s,ds,ds1
      real aij,pc
      integer lfrom(ntra),lto(ntra)

      common /lintable/ lfrom, lto
      save /lintable/

c     The length translation to pc is done via kappa
      data pc /3.0856e18/
        data zerolev /1.0d-14/

      b1=bij(l)
      b2=bji(l)
      lgth=pc*dco
      mrat=zerolev*b2
      rat=dd(lto(l))*b2-dd(lfrom(l))*b1
      if (abs(rat).lt.mrat) rat=mrat
      cap=lgth*rat
      dcap=lgth*b2
      dcap1=-lgth*b1
      hep=aij(l)/rat 
      s=hep*dd(lfrom(l))
      ds=-s*(b2/rat)
      ds1=hep + s*(b1/rat)
      return
      end

c     *************************************************************
c     Subroutine to get effective absorption coefficients for
c     a clumpy medium determined by a turbulent correlation length
c     tau=cap(l)*r(clump)/sigma(therm)=cap(l)*r(corr)/sigma(turb)
c     (The source function is not influenced by clumpiness.)
c     *************************************************************
c     In contrast to effkap, the routine is unified with atau
c     ************************************************************* 
      subroutine sobeffkap(cap,dcap,rc)
      implicit none
      real cap,dcap,rc

      real tau,ttau,dtt,mintau,deltat,atau
      real tabatau(57),tabdtau(57)
      integer i,it1,it2

c     *************************************************************
c     Since the reduction function is well behaving but represents
c     a nasty integral, the logarithmic values are tabulated.
c     The table goes from -4 to 10 with a step size of 0.25.
c     *************************************************************
      data deltat,mintau,it1,it2 /0.25, 0.01850, 17, 56/
      data tabatau / -0.00323227, -0.00414831, -0.00532323, -0.00682972,
     % -0.00876059, -0.0112341, -0.0144006, -0.018451, -0.0236262,
     % -0.0302296, -0.0386404, -0.0493291, -0.0628737, -0.0799749, 
     % -0.101467, -0.128322, -0.161636, -0.202594, -0.252403, -0.312191, 
     % -0.382874, -0.465019, -0.558747, -0.663717, -0.779211, -0.904295, 
     % -1.03798, -1.17932, -1.32748, -1.48173, -1.64142, -1.806, -1.975, 
     % -2.148, -2.32464, -2.5046, -2.6876, -2.87339, -3.06177, -3.25252, 
     % -3.44549, -3.6405, -3.83743, -4.03615, -4.23654, -4.43849,
     % -4.64192, -4.84673, -5.05285, -5.26021, -5.46873, -5.67836,
     % -5.88904, -6.10071, -6.31333, -6.52685, -6.74123 /
c     *************************************************************
c     Table of the differentiated reduction function (linear)
c     *************************************************************
      data tabdtau / -0.175607, -0.175277, -0.174854, -0.174313,
     % -0.173621, -0.172739, -0.171615, -0.170185, -0.168374,-0.166086,
     % -0.163209, -0.159613, -0.155152, -0.149669, -0.14301, -0.135042, 
     % -0.125685, -0.114946, -0.102962, -0.0900296, -0.0766157,
     % -0.0633147, -0.0507607, -0.0395066, -0.0299147, -0.0221102,
     % -0.0160085, -0.0113919, -0.00798999, -0.00553579, -0.00379561, 
     % -0.00257921, -0.00173904, -0.00116459, -0.000775221,-0.000513288, 
     % -0.000338241, -0.000221938, -0.000145064, -0.0000944844, 
     % -0.0000613441, -0.0000397114, -0.0000256383, -0.0000165116, 
     % -0.0000106095, -6.80268e-6, -4.35319e-6, -2.78059e-6,-1.77305e-6, 
     % -1.12878e-6, -7.1753e-7, -4.55467e-7, -2.88731e-7, -1.82803e-7, 
     % -1.156e-7, -7.30203e-8, -4.60751e-8 /

      tau=cap*rc
c     *************************************************************
c     Everything beyond the small-tau end of the table is ignored
c     *************************************************************
      if (tau.lt.mintau) then
      dcap=1.0
      return
      endif
c     *************************************************************
c     Read the tables
c     *************************************************************
      ttau=alog(tau)
      dtt=ttau/deltat+it1
      i=min(int(dtt),it2)
      dtt=dtt-i
      ttau=tabatau(i+1)*dtt+tabatau(i)*(1.0-dtt)
      atau=exp(ttau)
      ttau=tabdtau(i+1)+dtt+tabdtau(i)*(1.0-dtt)
c     *************************************************************
c     Assign to the cappa values
c     *************************************************************
      cap=cap*atau
      dcap=atau+tau*ttau
      return
      end

c     *************************************************************
c     Computation of the new level populations from the local
c     radiation field in the Sobolev approximation
c     The level populations are normalized to 1 here!
c     *************************************************************
c     At first the matrix of transitional balance equations
c     is set up. Then it is translated to differences of level
c     populations. Finally it is solved by lusolver
c     ************************************************************* 
      subroutine sobmatrix(T,sdens,uu,abu,abu1,vek,lin,lev,rerr)
      implicit none
      integer lev,lin
      real T,sdens,rerr

      include 'fsizes.inc'
      real uu(ntra),abu(ntra),abu1(ntra)
      real*8 vek(nlev),v(nlev)
      real*8 am(nlev,nlev),sumcol(nlev),urel,dens,rate, zerolev
      real*8 bij,bji
      real aij,cij
      integer j,k,kup,kdo
      integer lfrom(ntra),lto(ntra)

        common /lintable/ lfrom, lto
        save /lintable/
        data zerolev /1.0d-14/

      dens=sdens
c     ***********************************************************
c     Compute the original matrix
c     Same approach like in matixup
c     ***********************************************************
c     Initialize the matrix with the constructive collisional 
c     transitions ( col(j,j)=0 has to be fulfilled )
c     *************************************************************
      do j=1,lev
      sumcol(j)=0d0
      do k=1,lev
        am(k,j)=dens*cij(T,j,k)
        sumcol(j)=sumcol(j)+am(k,j)
      enddo
      enddo
      do j=1,lev
        am(j,j)=am(j,j)-sumcol(j)
      enddo
c     *************************************************************
c     Set up the different lines of the matrix  with radiative
c     transitions. Here, I go through the transitions, not levels.
c     *************************************************************
      do j=1,lin
        kup=lfrom(j)
        if (kup.le.lev) then
        kdo=lto(j)
        rate=aij(j)+bij(j)*uu(j)
        am(kdo,kup)=am(kdo,kup)+rate
        am(kup,kup)=am(kup,kup)-rate
        rate=bji(j)*uu(j)
        am(kup,kdo)=am(kup,kdo)+rate
        am(kdo,kdo)=am(kdo,kdo)-rate
        endif
      enddo
c     ***********************************************************
c     The right hand side
c     ***********************************************************
      do j=1,lev-1
        v(j)=0d0
        do k=1,lev
          v(j)=v(j)-vek(k)*am(j,k)
        enddo
      enddo
c     Last line
      v(lev)=1d0
      do k=1,lev 
        v(lev)=v(lev)-vek(k)
      enddo
c     ***********************************************************
c     Add the derivatives to the matrix
c     ***********************************************************
      do j=1,lin
        kup=lfrom(j)
        if (kup.le.lev) then
        kdo=lto(j)
        rate=abu1(j)*(bij(j)*vek(kup)-bji(j)*vek(kdo))
        am(kdo,kup)=am(kdo,kup)+rate
        am(kup,kup)=am(kup,kup)-rate
        rate=abu(j)*(bji(j)*vek(kdo)-bij(j)*vek(kup))
        am(kup,kdo)=am(kup,kdo)+rate
        am(kdo,kdo)=am(kdo,kdo)-rate
        endif
      enddo
c     *************************************************************
c     The last line. It is not linear independent from the 
c     others when taking the transition equations.
c     Here, the normalization condition has to be used.
c     *************************************************************
      do j = 1,lev
        am(lev,j) = 1d0 
      enddo
c     ***********************************************************
c     The solver
c     ***********************************************************
      call lusolver(am,lev,v)
c     ***********************************************************
c     Avoid negative populations which may occur during the
c     Newton-Raphson fit , compute the maximum change
c     ***********************************************************
      do j=1,lev
        vek(j)=vek(j)+v(j)
        if (vek(j).lt.zerolev) vek(j)=zerolev
      enddo
      rerr=0.0d0
      do j=1,lev-2
        urel=abs(v(j))/vek(j)
        if (urel.gt.rerr) rerr=urel
      enddo
      return
      end

