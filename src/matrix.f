c     ***********************************************************
c     Caller routine for ulevmic and dlevmic 
c     Depending on the number of lines the lambda iteration
c     is either done for the energy densities or the levels
c     ***********************************************************
      subroutine levmic(t,dh,ulev,dd,nn,lin,lev,istat)
      implicit none
      integer nn,lin,lev,istat

      include 'fsizes.inc'
      real t(ns),dh(ns),ulev(ns,ntra)
      real*8 dd(ns,nlev)

      if (lin.gt.lev) then
        call dlevmic(t,dh,ulev,dd,nn,lin,lev,istat)
      else
        call ulevmic(t,dh,ulev,dd,nn,lin,lev,istat)
      endif
      return
      end

c     ***********************************************************
c     Routine to compute the level populations and the new
c     emmissivity/absorptivity from the radiation density
c     at each radial grid point
c     In this version, the accelerated lambda iteration is
c     applied to the energy densities ulev
c     ***********************************************************
      subroutine ulevmic(t,dh,ulev,dd,nn,lin,lev,istat)
      implicit none
      integer nn,lin,lev,istat

      include 'fsizes.inc'
      real t(ns),dh(ns),ulev(ns,ntra)
      real ualt(3,ns,ntra)
      real*8 conv,neglev
      real*8 dd(ns,nlev)

      real ulin(ntra),urel
      real*8 vek(nlev),udens
      integer i,l,is,il

      common /acclev/ conv,neglev
      save /acclev/
      save ualt

c     ************************************************************
c     Call lambda iteration if 3 steps have been made
c     ************************************************************
      istat=istat+1
      if (istat.eq.4) call lambda(ualt,ulev,nn,lin)
c     ************************************************************
c     Compute new populations - initialize
c     ************************************************************
      do 400 i=1,nn
      do 450 l=1,lin
       ulin(l)=ulev(i,l)
450   continue
      call matrixup(T(i),dh(i),ulin,vek,lin,lev)
      do 480 l=1,lev
       dd(i,l)=vek(l)
480   continue
400   continue
c     ************************************************************
c     Compute maximum change in this step
c     The most outer transition(closure determined) is not counted
c     ************************************************************
      if (istat.gt.1) then
      udens=0.0d0
      is=nn
      il=1
      do 110 i=1,nn
      do 110 l=1,lin
      urel=abs((ulev(i,l)-ualt(istat-1,i,l))/ulev(i,l))
      if (urel.gt.udens) then
      udens=urel
      is=i
      il=l
      endif
110   continue
      else
      udens=1.0d0
      is=nn
      il=1
      endif
c     ************************************************************
c     Look for convergence
c     The removal of superfluid levels is now done externally
c     ************************************************************
      if (istat.eq.4) then
      if (udens.lt.conv) then
      istat=-1
      else
      istat=1
      call removelevels(dd,nn,lin,lev,istat)
      endif
      endif
c     ************************************************************
c     Store old values for accelerated lambda-iteration
c     ************************************************************
      if (istat.gt.0) then
      do 100 i=1,nn
      do 100 l=1,lin
       ualt(istat,i,l)=ulev(i,l)
100   continue
      endif
C     Control Output (Kruegel)
17    format (' Maximum change in this step: ',1PE9.2,' shell: ',
     %  I3,' transition: ',I3)
      write(6,17) udens,is,il
      return
      end

c     ***********************************************************
c     Routine to compute the level populations
c     at each radial grid point
c     In this version, the accelerated lambda iteration is
c     applied to the level populations
c     ***********************************************************
      subroutine dlevmic(t,dh,ulev,dd,nn,lin,lev,istat)
      implicit none
      integer nn,lin,lev,istat

      include 'fsizes.inc'
      real t(ns),dh(ns),ulev(ns,ntra)
      real*8 djkalt(3,ns,nlev)
      real*8 conv,neglev
      real*8 dd(ns,nlev)

      real ulin(ntra)
      real*8 vek(nlev),ddens,drel
      integer i,l,is,il

      common /acclev/ conv,neglev
      save /acclev/
      save djkalt

c     ************************************************************
c     Store old values in case of a new grid
c     ************************************************************
      if (istat.eq.0) then
        do i=1,nn
        do l=1,lev
          djkalt(istat,i,l)=dd(i,l)
        enddo
        enddo
        istat=1
      endif
c     ************************************************************
c     Compute new populations - initialize
c     ************************************************************
      do 400 i=1,nn
      do 450 l=1,lin
       ulin(l)=ulev(i,l)
450   continue
      call matrixup(T(i),dh(i),ulin,vek,lin,lev)
      do 480 l=1,lev
       dd(i,l)=vek(l)
480   continue
400   continue
c     ************************************************************
c     Call lambda iteration if 3 steps have been made
c     ************************************************************
      istat=istat+1
      if (istat.eq.4) call dlambda(djkalt,dd,nn,lev)
c     ************************************************************
c     Compute maximum change in this step
c     The most outer transition(closure determined) is not counted
c     ************************************************************
      ddens=0.0d0
      do 110 i=1,nn
      do 110 l=1,lev-1
      drel=abs((dd(i,l)-djkalt(istat-1,i,l))/dd(i,l))
      if (drel.gt.ddens) then
      ddens=drel
      is=i
      il=l
      endif
110   continue
c     ************************************************************
c     Look for convergence
c     The removal of superfluid levels is now done externally
c     ************************************************************
      if (istat.eq.4) then
      if (ddens.lt.conv) then
      istat=-1
      else
      istat=1
      call removelevels(dd,nn,lin,lev,istat)
      endif
      endif
c     ************************************************************
c     Store old values for accelerated lambda-iteration
c     ************************************************************
      if (istat.gt.0) then
      do 100 i=1,nn
      do 100 l=1,lev
       djkalt(istat,i,l)=dd(i,l)
100   continue
      endif

C     Control Output (Kruegel)
17    format (' Maximum change in this step: ',1PE9.2,' shell: ',
     %  I3,' level: ',I3)
      write(6,17) ddens,is,il
      return
      end

c     ************************************************************
c     Throw out superfluid levels
c     ************************************************************
      subroutine removelevels(dd,nn,lin,lev,istat)
      implicit none
      integer nn,lin,lev,istat

      include 'fsizes.inc'
      real*8 dd(ns,nlev),negvek(nlev),zerolev
      real*8 conv,neglev
      integer i,l,olev

      common /acclev/ conv,neglev
      save /acclev/
      data zerolev /1.0d-14/

      do 150 l=1,lev
       negvek(l)=zerolev
150   continue

c     Loop of radii
      do i=1,nn
        do l=1,lev
         if (dd(i,l).gt.negvek(l)) negvek(l)=dd(i,l)
        enddo
      enddo
c     Now check all levels
      olev=lev
      l=lev
300   if (negvek(l).gt.neglev) then
      lev=l
      else
      l=l-1
      goto 300
      endif

      if (lev.ne.olev) then
       call linfromlev(lin,lev)
       write(6,'('' Superfluid high levels removed: '',I3)') olev-lev
       istat=0
      endif
      return
      end

c     *************************************************************
c     Subroutine to remove superfluid lines and rearrange fields
c     *************************************************************
      subroutine linfromlev(lin,lev)
      implicit none
      integer lin,lev

      include 'fsizes.inc'
      real fr(ntra),eina(ntra),einb(ntra),gewratio(ntra)
      integer i,l,lfrom(ntra),lto(ntra)

        common /lintable/ lfrom, lto
        common /coeffein/ fr,eina,einb,gewratio
      save /coeffein/,/lintable/

c     Remove superfluid lines
      l=1
400   if (lfrom(l).gt.lev) then
       do i=l+1,lin
         lfrom(i-1)=lfrom(i)
         lto(i-1)=lto(i)
         fr(i-1)=fr(i)
         eina(i-1)=eina(i)
         einb(i-1)=einb(i)
         gewratio(i-1)=gewratio(i)
       enddo
       lin=lin-1
       l=l-1
      endif
      l=l+1
      if (l.le.lin) goto 400

      return
      end

c     *************************************************************
c     Subroutine for lambda acceleration
c     according to Auer (1987), in Kalkofen (Ed.) "Numerical
c     Methods in Radiative Transfer", p. 101
c     Second order residual minimization
c     First version with real input fields
c     *************************************************************
      subroutine lambda(dalt,dd,nn,lin)
      implicit none
      integer nn,lin

      include 'fsizes.inc'
      real dalt(3,ns,ntra), dd(ns,ntra), dneu(ns)
      real*8 aa,bb,cc,a1,b1,b2,c1,c2,d0,d1,d2,wt,ddif,drel
      integer i,l

c     *************************************************************
c     The selection of the vectorization is arbitrary.
c     Either all levels at one point or all points at one level
c     may be taken. Tests have favoured the second choice.
c     *************************************************************
      do 100 l=1,lin
      a1=0d0
      b1=0d0
      b2=0d0
      c1=0d0
      c2=0d0
      wt=0d0
      do 200 i=1,nn
      wt=wt+dble(dd(i,l)**2)
      d0=dble(dd(i,l)-dalt(3,i,l))
      d1=dble(dd(i,l)-2.0*dalt(3,i,l)+dalt(2,i,l))
      d2=dble(dd(i,l)-dalt(3,i,l)-dalt(2,i,l)+dalt(1,i,l))
      a1=a1+d1*d1
      b1=b1+d1*d2
      b2=b2+d2*d2
      c1=c1+d0*d1
      c2=c2+d0*d2
200   continue
c     ************************************************************
c     Compute the resulting changes and apply them
c     ************************************************************
      wt=1d0/wt
      a1=a1*wt
      b1=b1*wt
      b2=b2*wt
      c1=c1*wt
      c2=c2*wt
c     ************************************************************
c     If the changes are below the machine accuracy - skip 
c     ************************************************************
      cc=b2*a1-b1*b1
      if (abs(cc).gt.1d-38) then
      aa=(b2*c1-b1*c2)/cc
      bb=(a1*c2-b1*c1)/cc
c     ************************************************************
c     Possible loop to prevent overshooting
c     ************************************************************
250   cc=1d0-aa-bb
      ddif=0.6d0
      do 300 i=1,nn
      dneu(i)=cc*dd(i,l)+aa*dalt(3,i,l)+bb*dalt(2,i,l)
      drel=dneu(i)/dd(i,l)
      if (drel.lt.ddif) ddif=drel
300   continue
      if (ddif.lt.0.5d0) then
      drel=2.5d0*(1d0-ddif)
      aa=aa/drel
      bb=bb/drel
      goto 250
      endif
      do 310 i=1,nn
       dd(i,l)=dneu(i)
310   continue
      endif
100   continue
      return
      end

c     *************************************************************
c     Subroutine for lambda acceleration
c     according to Auer (1987), in Kalkofen (Ed.) "Numerical
c     Methods in Radiative Transfer", p. 101
c     Second order residual minimization
c     Second version with real*8 input fields
c     *************************************************************
      subroutine dlambda(dalt,dd,nn,lin)
      implicit none
      integer nn,lin

      include 'fsizes.inc'
      real*8 dalt(3,ns,nlev), dd(ns,nlev), dneu(ns)
      real*8 aa,bb,cc,a1,b1,b2,c1,c2,d0,d1,d2,wt,ddif,drel
      integer i,l

c     *************************************************************
c     The selection of the vectorization is arbitrary.
c     Either all levels at one point or all points at one level
c     may be taken. Tests have favoured the second choice.
c     *************************************************************
      do 100 l=1,lin
      a1=0d0
      b1=0d0
      b2=0d0
      c1=0d0
      c2=0d0
      wt=0d0
      do 200 i=1,nn
      wt=wt+dd(i,l)**2
      d0=dd(i,l)-dalt(3,i,l)
      d1=dd(i,l)-2.0*dalt(3,i,l)+dalt(2,i,l)
      d2=dd(i,l)-dalt(3,i,l)-dalt(2,i,l)+dalt(1,i,l)
      a1=a1+d1*d1
      b1=b1+d1*d2
      b2=b2+d2*d2
      c1=c1+d0*d1
      c2=c2+d0*d2
200   continue
c     ************************************************************
c     Compute the resulting changes and apply them
c     ************************************************************
      wt=1d0/wt
      a1=a1*wt
      b1=b1*wt
      b2=b2*wt
      c1=c1*wt
      c2=c2*wt
c     ************************************************************
c     If the changes are below the machine accuracy - skip 
c     ************************************************************
      cc=b2*a1-b1*b1
      if (abs(cc).gt.1d-38) then
      aa=(b2*c1-b1*c2)/cc
      bb=(a1*c2-b1*c1)/cc
c     ************************************************************
c     Possible loop to prevent overshooting
c     ************************************************************
250   cc=1d0-aa-bb
      ddif=0.6d0
      do 300 i=1,nn
      dneu(i)=cc*dd(i,l)+aa*dalt(3,i,l)+bb*dalt(2,i,l)
      drel=dneu(i)/dd(i,l)
      if (drel.lt.ddif) ddif=drel
300   continue
      if (ddif.lt.0.5d0) then
      drel=2.5d0*(1d0-ddif)
      aa=aa/drel
      bb=bb/drel
      goto 250
      endif
      do 310 i=1,nn
       dd(i,l)=dneu(i)
310   continue
      endif
100   continue
      return
      end

c     *************************************************************
c     Subroutine to compute the absorption coefficient kappa and 
c     the source function sc from the level populations for all 
c     radii and levels
c     The frequency dependence phi(ny)*ny0 is multiplied outside
c     *************************************************************
      subroutine kapsrc(dco,dd,cap,src,nn,lin)
      implicit none
      integer nn,lin

      include 'fsizes.inc'
      real*8 dd(ns,nlev),rat,mrat,zerolev
      real dco(ns),cap(ns,ntra),src(ns,ntra)
      real*8 bij,bji
      real aij,pc
      integer i,l,l1,l2
      integer lfrom(ntra),lto(ntra)

        common /lintable/ lfrom, lto
      save /lintable/
      data zerolev /1.0d-14/

c     The prefactor h/(4*pi) is shifted towards the Einstein Bs 
c     i.e. the intensity units are divided by h/(4pi)
c     Here, the length translation to pc is done
      data pc /3.0856e18/

      do  400  l = 1, lin
      l1=lfrom(l)
      l2=lto(l)
      mrat=zerolev*bji(l)
      do  400  i = 1, nn
       rat=dd(i,l2)*bji(l)-dd(i,l1)*bij(l)
       if (abs(rat).lt.mrat) rat=mrat
       cap(i,l)=pc*dco(i)*rat
       src(i,l)=dd(i,l1)*aij(l)/rat
400   continue
      return
      end

c     *************************************************************
c     Subroutine to get effective absorption coefficients for
c     a clumpy medium determined by a turbulent correlation length
c     tau=cap(l)*r(clump)/sigma(therm)=cap(l)*r(corr)/sigma(turb)
c     (The source function is not influenced by clumpiness.)
c     *************************************************************
      subroutine effkap(cap,rc,nn,lin)
      implicit none
      integer nn,lin

      include 'fsizes.inc'
      real cap(ns,ntra),rc(ns)
      real tau,atau
      integer i,l

      do 100 i=1,nn
      do 100 l=1,lin
      tau=cap(i,l)*rc(i)
      cap(i,l)=cap(i,l)*atau(tau)
100   continue
      return
      end

c     *************************************************************
c     Function giving the effective kappa reduction depening on
c     the normalized clump optical depth: A(tau)/tau
c     *************************************************************
      real function atau(tau)
      implicit none
      real tau,ttau,dtt,mintau,deltat
      real tabatau(57)
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
c     Everything beyond the small-tau end of the table is ignored
c     *************************************************************
      if (tau.lt.mintau) then
      atau=1.0
      return
      endif
c     *************************************************************
c     Read the table
c     *************************************************************
      ttau=alog(tau)
      dtt=ttau/deltat+it1
      i=min(int(dtt),it2)
      dtt=dtt-i
      ttau=tabatau(i+1)*dtt+tabatau(i)*(1.0-dtt)
      atau=exp(ttau)
      return
      end

c     *************************************************************
c     Make the excitation in an HII region smooth from the viewpoint 
c     of other routines which do not need to know about it
c     *************************************************************
      subroutine hiismth(dd,lev)
      implicit none
      integer lev

      include 'fsizes.inc'
      real*8 dd(ns,nlev)
      integer l

      do 100 l=1,lev
       dd(1,l)=dd(2,l)
100   continue
      return
      end

c     *************************************************************
c     Compute kappa and S for the HII region - field adaptor
c     *************************************************************
      subroutine hiikap(cap,src,lin)
      implicit none
      integer lin

      include 'fsizes.inc'
      real cap(ns,ntra),src(ns,ntra)
      real kp(ntra),s(ntra)
      integer l

      common /hiidata/ kp,s
      save /hiidata/

      do 100 l=1,lin
      cap(1,l)=kp(l)
      src(1,l)=s(l)
100   continue
      return
      end

c     *************************************************************
c     Computation of the new level populations from the given
c     radiation field
c     The level populations are normalized to 1 here!
c     *************************************************************
c     At first the matrix of transitional balance equations
c     is set up
c     In a second step it is solved by lusolver
c     ************************************************************* 
      subroutine matrixup(T,sdens,uu,vek,lin,lev)
      implicit none
      integer lev,lin
      real T,sdens

      include 'fsizes.inc'
      real uu(ntra)
      real*8 vek(nlev)
      real*8 am(nlev,nlev),sumcol(nlev),dens,rate,zerolev
      real*8 bij,bji
      real aij,cij
      integer j,k,kup,kdo
      integer lfrom(ntra),lto(ntra)

        common /lintable/ lfrom, lto
      save /lintable/
      data zerolev /1.0d-14/

      dens=sdens
c     ************************************************************
c     Set up the matrix
c     ************************************************************
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
        kdo=lto(j)
        rate=aij(j)+bij(j)*uu(j)
        am(kdo,kup)=am(kdo,kup)+rate
        am(kup,kup)=am(kup,kup)-rate
        rate=bji(j)*uu(j)
        am(kup,kdo)=am(kup,kdo)+rate
        am(kdo,kdo)=am(kdo,kdo)-rate
      enddo
c     *************************************************************
c     The last line. It is not linear independent from the 
c     others when taking the transition equations.
c     Here, the normalization condition has to be used.
c     *************************************************************
      do  250  j = 1,lev
       am(lev,j) = 1d0 
250   continue
c     *************************************************************
c     The right hand sight of the equation system
c     ************************************************************
      do 258 j=1,lev-1
       vek(j)=0d0
258   continue
      vek(lev)=1d0
c     ***********************************************************
c     The solver
c     ***********************************************************
      call lusolver(am,lev,vek)
c     ***********************************************************
c     Under extreme cases it may happen that some elements which
c     are many orders of magnitude smaller than the others
c     become negative due to limited real*8 accuracy. 
c     Values which are smaller than one digit of 1 are fixed. 
c     ***********************************************************
      do 222 j=1,lev
      if (vek(j).lt.zerolev) vek(j)=zerolev
222   continue
      return
      end

c     **************************************************************
c     The matrix equation for the number densities of the level
c     populations is solved by LU decomposition with improvement.
c     Global caller routine for the Numerical recipies subroutines
c     **************************************************************
      SUBROUTINE lusolver(a,n,b)
      implicit none
      integer n
      include 'fsizes.inc'

      INTEGER i,j,indx(nlev),count
      REAL*8 a(nlev,nlev),astor(nlev,nlev)
      REAL*8 b(nlev),bstor(nlev),mxchg
      real*8 conv,neglev

      common /acclev/ conv,neglev
      save /acclev/
c     **************************************************************
c     First save a and b for later improvement
c     **************************************************************
      do i=1,n
        do j=1,n
          astor(i,j)=a(i,j)
        enddo
        bstor(i)=b(i)
      enddo
c     **************************************************************
c     Call the Numerical recipies functions
c     **************************************************************
      call ludcmp(a,indx,n)
      call lubksb(a,indx,b,n)
c     **************************************************************
c     Improvement only executed for large matrices
c     Iterative improvement up to convergence limit
c     **************************************************************
      if (n.lt.15) return
      count=0
100   continue
        call mprove(astor,a,indx,bstor,b,n,mxchg)
        count=count+1
      if ((mxchg.gt.conv).and.(count.lt.5)) goto 100

      return
      end

c     **************************************************************
c     Here, we use the routines ludcmp, lubksb, and mprove from the
c     'Numerical recipies in FORTRAN'
c     They were slightly modified to reduce the number of 
c     parameters for use with unique field sizes
c     **************************************************************
c     Third part - imrovement according to mprove
c     **************************************************************
      SUBROUTINE mprove(a,alud,indx,b,x,n,mxchg)
      implicit none
      include 'fsizes.inc'
      INTEGER n,indx(nlev),i,j
      REAL*8 a(nlev,nlev),alud(nlev,nlev),b(nlev),x(nlev),r(nlev)
      real*8 sdp,mxchg

      do 12 i=1,n
        sdp=-b(i)
        do 11 j=1,n
          sdp=sdp+a(i,j)*x(j)
11      continue
        r(i)=sdp
12    continue
      call lubksb(alud,indx,r,n)
      mxchg=0.0d0
      do 13 i=1,n
        mxchg=max(mxchg,abs(r(i)))
        x(i)=x(i)-r(i)
13    continue
      return
      END

c     **************************************************************
c     First part - matrix decomposition according
c     **************************************************************
c     a is the matrix (will be destroyed)
c     n is the really used dimension of the matrix 
c     nlev is the physical dimension of all vectors
c     index is a vector of permutations 
c      - will be used by the second part
c     **************************************************************
      SUBROUTINE ludcmp(a,indx,n)
      implicit none
      include 'fsizes.inc'
      INTEGER n,indx(nlev),i,imax,j,k
      REAL*8 a(nlev,nlev),TINY,aamax,dum,sum,vv(nlev)
      PARAMETER (TINY=1.0e-20)

      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) then
        write(6,*) 'singular matrix in ludcmp'
        stop 1
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END

c     **************************************************************
c     Second part - solution of the two triagular matrices
c     by forward and back substitutaion
c     **************************************************************
c     b is the vector with inhomogeneous part 
c      - will be returned as the solution vector
c     n is the used dimension of the vector
c     **************************************************************
      SUBROUTINE lubksb(a,indx,b,n)
      implicit none
      include 'fsizes.inc'
      INTEGER n,indx(nlev),i,ii,j,ll
      REAL*8 a(nlev,nlev),b(nlev),sum

      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END

