c     ****************************************************************
c     Routine to set up the initial radial grid with the values for 
c     density, temperature, turbulent, and systematic velocity
c     ****************************************************************
      subroutine initrad(r,t,dh,dco,vs,vt,rc,nn,kerr)
      implicit none
      integer nn,kerr

c     output parameter
      include 'fsizes.inc'
      real r(ns),t(ns),dh(ns),dco(ns),vt(ns),vs(ns),rc(ns)
      real epsv,epsn,dvgauss,hyster,minspac,maxspac
      
c     input parameter (from common)
      integer ishell,fixed(0:ns)
      real rs(0:ns),parm(ns,6,2)

c     internal variables
      real ldelta,dhtmp,exv,dv,v0,v1,rdelta,dhrad
      real thick,small,minedge,toosmall
      integer i,j,jj,na,n0,np,k,nk
      integer new(ns),add(ns),place(ns)
      integer cnt,sumadd,i0,ind

      common /model/ rs,fixed,parm,ishell
      common /accgrid/ epsn,epsv,dvgauss,hyster,minspac,minedge
      save /model/,/accgrid/

c     The limit for the smallest layer will not be used very often,
c     so that it is not treated as a normal numerical parameter
      data small / 0.02/
      data toosmall /1.0e-6/

      ldelta=hyster*log(epsn)
      rdelta=epsn-1.0
c     ************************************************************
c     Exception handling for zero width shells - remove
c     ************************************************************
14    format(/,' Zero width shell ',i3,' removed')
c     Determine the maximum shell size
      maxspac=0.0
      do 300 j=1,ishell
       thick=rs(j)-rs(j-1)
       if (thick.gt.maxspac) maxspac=thick
300   continue
      minspac=toosmall*maxspac
      j=0
      cnt=0
c     loop to remove shells
320   continue
        j=j+1
        thick=rs(j)-rs(j-1)
        if ((thick.lt.minspac).and.(fixed(j).eq.0)) then
            write(6,14) j-1+cnt
            cnt=cnt+1
            do 330 k=j,ishell-1
             rs(k)=rs(k+1)
             fixed(k)=fixed(k+1)
             do 330 nk=1,6
             parm(k,nk,1)=parm(k+1,nk,1)
             parm(k,nk,2)=parm(k+1,nk,2)
330         continue
            ishell=ishell-1
            j=j-1
        endif
      if (j.lt.ishell) goto 320
c     ************************************************************
c     At first, a minimum spacing will be chosen which is used
c     to remove possible jumps in the parameters
c     ************************************************************
11    format(/,' Discontinuity detected between shells',i3,' and',i3)
c     Determine the minimum shell size
      minspac=1e30
      do 310 j=1,ishell
      thick=rs(j)-rs(j-1)
      if (thick.lt.minspac) minspac=thick
310   continue
      minspac=small*minspac
      minedge=rs(2)
c     ************************************************************
c     The grid is adjusted to the velocity and density structure 
c     of the cloud. The density at neighbouring points may deviate 
c     by less than a factor epsn^hyster. The velocity may deviate 
c     by less than epsv.
c     The parameters within the core have to be constant
c     ************************************************************
c     The core is treated separately
      kerr=0
      nn=2
      r(1)=0.0
      call values(r(1),dh(1),dco(1),t(1),vt(1),rc(1),vs(1))
c     No discontinuity can occur at the most inner shell radius
      r(2)=rs(1)
      call values(r(2),dh(2),dco(2),t(2),vt(2),rc(2),vs(2))
c     ************************************************************
c     Outer loop - treat all shells
c     ************************************************************
      do 500 j=2,ishell
c     Look for discontinuities at fixed shell boundaries
      if (fixed(j-1).gt.0) then
      np=nn+1
      r(np)=rs(j-1)
      call values(r(np),dh(np),dco(np),t(np),vt(np),rc(np),vs(np))
      if ((abs(log(dh(np)/dh(nn))).gt.ldelta).or.
     %      (abs(vs(np)/vt(np)-vs(nn)/vt(nn)).gt.epsv)) then
      write(6,11) j-2,j-1
      nn=np
      else
      r(nn)=r(np)
      dh(nn)=dh(np)
      dco(nn)=dco(np)
      t(nn)=t(np)
      vt(nn)=vt(np)
      rc(nn)=rc(np)
      vs(nn)=vs(np)
      endif
      endif
c     ************************************************************
c     The first grid is set up on the basis of the velocity
c     structure - linear steps in velocity difference
c     ************************************************************
c     Since we compare the quotient of systematic and turbulent 
c     velocity, the turbulent velocity should vary slower than the 
c     systematic velocity to obtain meaningful spacings.
c     ************************************************************
      n0=nn
      np=n0+1
      call values(rs(j),dh(np),dco(np),t(np),vt(np),rc(np),vs(np))
      v0=vs(nn)/vt(nn)
      v1=vs(np)/vt(np)
      dv=abs(v1-v0)
      na=int(dv/epsv)+1
c     Check whether the number of radial points is sufficient
      if (n0+na.ge.ns) then
      kerr=1 
      return
      endif
c     Assign the new radii
      if (na.gt.1) then
        if (v1*v0.le.0.0) then
        dv=(rs(j)-rs(j-1))/na
        do i=1,na-1
          nn=nn+1
          r(nn)=rs(j-1)+i*dv
        enddo
        else
        exv=alog(rs(j)/rs(j-1))/abs(alog(v1/v0))
        dv=dv/(na*v0)
        do i=1,na-1
          nn=nn+1
          r(nn)=rs(j-1)*(1.0+i*dv)**exv
        enddo
        endif
      endif
      nn=nn+1
      if (fixed(j).eq.0) then
       r(nn)=rs(j)
      else
       r(nn)=rs(j)-0.5*minspac
      endif
c     ************************************************************
c     Compute the values of the other quantities at these points
c     ************************************************************
      do 200 i=n0+1,nn
      call values(r(i),dh(i),dco(i),t(i),vt(i),rc(i),vs(i))
200   continue
c     ************************************************************
c     Second loop - it is tested whether I need additional points
c     due to the density structure, and whether additional points
c     are required to get a reasonable radial coverage and these
c     points are added logarithmically
c     ************************************************************
230   cnt=0
      sumadd=0
      do 400 i=n0,nn-1
      place(i)=i+sumadd
c     The radial criterion
      dhrad=(r(i+1)-r(i))/max(r(i+1),minedge)
c     The density criterion
      dhtmp=abs(log(dh(i+1)/dh(i)))
      if ((dhtmp.gt.ldelta).or.(dhrad.gt.rdelta)) then
       cnt=cnt+1
       new(cnt)=i
       add(cnt)=int(max(dhtmp/ldelta,dhrad/rdelta))
       sumadd=sumadd+add(cnt)
      endif
400   continue
      place(nn)=nn+sumadd
c     No or too many additional points
      if (cnt.eq.0) goto 500
      if (nn+sumadd.ge.ns) then
      kerr=1 
      return
      endif
c     ************************************************************
c     Add additional points in those intervals 
c     ************************************************************
c     First make holes in the field - backward scan
      do 420 i=nn,n0,-1
      i0=place(i)
      r(i0)=r(i)
      dh(i0)=dh(i)
      dco(i0)=dco(i)
      t(i0)=t(i)
      vt(i0)=vt(i)
      rc(i0)=rc(i)
      vs(i0)=vs(i)
420   continue
c     Now fill the holes
      do 450 jj=1,cnt
      i0=place(new(jj))
      exv=(r(i0+add(jj)+1)/r(i0))**(1.0/(add(jj)+1))
      do 480 ind=1,add(jj)
c     Assign the new values
      i=i0+ind
      r(i)=r(i0)*exv**ind
      call values(r(i),dh(i),dco(i),t(i),vt(i),rc(i),vs(i))
480   continue
450   continue
      nn=nn+sumadd
c     ************************************************************
c     End of the loop for all shells - Correct the last point
c     ************************************************************
500   continue
      r(nn)=rs(ishell)
      call values(r(nn),dh(nn),dco(nn),t(nn),vt(nn),rc(nn),vs(nn))
      return
      end
 
c     ************************************************************
c     Subroutine to add and remove additional points in the radial
c     grid depending on the variation of the level populations
c     ************************************************************
c     Additionally required points are always added and the
c     accelerated lambda iteration cycle is restarted,
c     Points which are not necessary are removed only in the first
c     two steps of an accelerated lambda iteration cycle
c     ************************************************************
      subroutine adaptgrid(r,vs,vt,rc,dco,dh,t,dd,nn,lev,istat,kerr)
      implicit none
      integer nn,lev,istat,kerr

      include 'fsizes.inc'
      real r(ns),vs(ns),vt(ns),dco(ns),dh(ns),t(ns),rc(ns)
      real*8 dd(ns,nlev)
      integer nold,cntadd

c     ************************************************************
c     First remove superfluid points to spare grid points
c     ************************************************************
      if (istat.lt.3)
     %     call removepoints(r,vs,vt,rc,dco,dh,t,dd,nn,lev,istat)
      nold=nn
c     Now add necessary points
      call addpoints(r,vs,vt,rc,dco,dh,t,dd,nn,lev,istat,kerr)
      if (kerr.ne.0) return
      cntadd=nn-nold
c     Again remove redundant points appearing at inset edges
      if (cntadd.ne.0)
     %     call removepoints(r,vs,vt,rc,dco,dh,t,dd,nn,lev,istat)
      return
      end

c     *************************************************************
c     First subroutine - add necessary points
c     *************************************************************
      subroutine addpoints(r,vs,vt,rc,dco,dh,t,dd,nn,lev,istat,kerr)
      implicit none
      integer nn,lev,istat,kerr

      include 'fsizes.inc'
      real r(ns),vs(ns),vt(ns),dco(ns),dh(ns),t(ns),rc(ns)
      real*8 dd(ns,nlev)
      real*8 d1,diff
      real epsn,epsv,dvgauss,hyster,minspac,minedge
      real*8 conv,neglev

      integer i,l,j,i0,cnt, loop
      integer place(ns),new(ns)

      common /accgrid/ epsn,epsv,dvgauss,hyster,minspac,minedge
      common /acclev/ conv,neglev
      save /accgrid/, /acclev/

      kerr=0
c     *************************************************************
c     First radius scan
c     Find the maximum difference factors, look for the need of
c     additional points (only included for spaces>minspac)
c     The last level is not considered since it is often determined
c     by the cut of the equation system.
c     Values below neglev are ignored
c     *************************************************************
      loop=0
1000  cnt=0
      do 100 i=1,nn-1
      if (dd(i,1)*0.0d0.ne.0.0d0) then
       kerr=-1
       return
      endif
      diff=0d0
      do 150 l=1,lev-1
      if ((dd(i+1,l).gt.neglev).and.(dd(i,l).gt.neglev)) then
        d1=dco(i+1)*dd(i+1,l)/(dco(i)*dd(i,l))
        d1=max(d1,1d0/d1)
        if (d1.gt.diff) diff=d1
      endif  
150   continue
c     Compute the numbers for the field translation
      place(i)=i+cnt
      if ((diff.gt.epsn).and.(r(i+1)-r(i).gt.minspac)) then
      cnt=cnt+1
      new(cnt)=i
      endif
100   continue
      place(nn)=nn+cnt
      loop=loop+1
c     ************************************************************
c     Anything/nothing to be done
c     ************************************************************
      if (cnt.eq.0) return
c     Check for field overflow
      write(6,'('' Additional grid points needed: '',I4)') cnt
      if (nn+cnt.gt.ns) then
      kerr=1 
      return
      endif
c     ************************************************************
c     Include the additional grid points
c     Due to the power law behaviour of most quantities, the
c     additional points are added as geometric means
c     ************************************************************
c     First make holes in the field - backward scan
      do 420 i=nn,1,-1
      i0=place(i)
      r(i0)=r(i)
      dh(i0)=dh(i)
      dco(i0)=dco(i)
      t(i0)=t(i)
      vt(i0)=vt(i)
      rc(i0)=rc(i)
      vs(i0)=vs(i)
      do 420 l=1,lev
      dd(i0,l)=dd(i,l)
420   continue
c     Now fill the holes
      do 450 j=1,cnt
      i0=place(new(j))
      i=i0+1
      if (i0.eq.1) then
      r(i)=0.5*(r(i0+2)+r(i0))
      else
      r(i)=sqrt(r(i0+2)*r(i0))
      endif
      call values(r(i),dh(i),dco(i),t(i),vt(i),rc(i),vs(i))
      do 455 l=1,lev
      dd(i,l)=dsqrt(dd(i0+2,l)*dd(i0,l))
455   continue
450   continue
      nn=nn+cnt
      istat=0
c     *************************************************************
c     End of the loop 
c     *************************************************************
      goto 1000
      end

c     ************************************************************
c     Second part - find and remove superfluid grid points
c     ************************************************************
      subroutine removepoints(r,vs,vt,rc,dco,dh,t,dd,nn,lev,istat)
      implicit none
      integer nn,lev,istat

      include 'fsizes.inc'
      real r(ns),vs(ns),vt(ns),dco(ns),dh(ns),t(ns),rc(ns)
      real*8 dd(ns,nlev)
      real*8 den,denu,denl,d1,diff
      real epsn,epsv,dvgauss,hyster,minspac,minedge
      real*8 conv,neglev
      real dv,epshyst,epsmon,dhrad,rdelta
      integer i,l,i0,ilocat,cnt
      integer place(ns)

      integer ishell,fixed(0:ns)
      real rs(0:ns),parm(ns,6,2)
      common /model/ rs,fixed,parm,ishell

      common /accgrid/ epsn,epsv,dvgauss,hyster,minspac,minedge
      common /acclev/ conv,neglev
      save /model/,/accgrid/,/acclev/

c     ************************************************************
c     find the maximum difference factors
c     ************************************************************
      epshyst=epsn**hyster
      rdelta=epsn-1.0
      epsmon=1.0+conv
100   cnt=0
      place(1)=1
      i=2
650   continue
c     Points which must not be removed due to fixed shell edges
      i0=ilocat(rs,ishell+1,r(i))-1
      if ((abs((r(i)-rs(i0))/rs(i0)).lt.1e-6).and.(fixed(i0).gt.0))
     %    goto 680
c     Points which must not be removed due to large velocity gradients
c     In case of an HII region this point must not be reached
      dv=abs(vs(i+1)/vt(i+1)-vs(i-1)/vt(i-1))
      if (dv.gt.epsv) goto 680
c     Points which must not be removed due to keep the geometric
c     resolution fine enough
      dhrad=(r(i+1)-r(i-1))/max(r(i+1),minedge)
      if (dhrad.gt.rdelta) goto 680
c     Scan the levels
      diff=0d0
      do 250 l=1,lev-1
      denl=dco(i-1)*dd(i-1,l)
      den=dco(i)*dd(i,l)
      denu=dco(i+1)*dd(i+1,l)
c     Points which must not be removed due to non-monotony in dd*dh
      if (epsmon*abs(denu-denl).lt.abs(denu-den)+abs(den-denl)) 
     %     goto 680
      d1=denu/denl
      d1=max(d1,1d0/d1)
      if (d1.gt.diff) diff=d1
250   continue
      if (diff.lt.epshyst) then
      place(i)=i-cnt
      cnt=cnt+1
      i=i+1
      endif
c     End of the loop and jump address when nothing may happen
680   place(i)=i-cnt
      i=i+1
      if (i.lt.nn) goto 650
      place(nn)=nn-cnt
c     **************************************************************
c     Exit in case of no superfluid points
c     **************************************************************
      if (cnt.eq.0) return
      write(6,'('' Superfluid grid points removed: '',I4)') cnt
      if ((cnt.lt.0.1*nn).and.(istat.gt.1)) return
c     **************************************************************
c     Remove the superfluid points using the table created
c     The points are overwritten
c     **************************************************************
      do 720 i=1,nn
      i0=place(i)
      r(i0)=r(i)
      dh(i0)=dh(i)
      dco(i0)=dco(i)
      t(i0)=t(i)
      vt(i0)=vt(i)
      rc(i0)=rc(i)
      vs(i0)=vs(i)
      do 720 l=1,lev
      dd(i0,l)=dd(i,l)
720   continue
      nn=nn-cnt
      istat=0
c     *************************************************************
c     End of the loop 
c     *************************************************************
      goto 100
      end

c     ************************************************************
c     Subroutine to assign the physical parameters to a given 
c     radial point. (Here, the reduction of the velocities 
c     relative to c is performed.)
c     ************************************************************
      subroutine values(rx,pp1,pp3,pp2,pp4,pp5,pp6)
      implicit none

      include 'fsizes.inc'
      integer is,fixed(0:ns)
      real rx,rs(0:ns),parm(ns,6,2),pp(6)
      real pp1,pp2,pp3,pp4,pp5,pp6

      integer ind,ilocat,j
      real vtherm,cst,clicht
      common /model/ rs,fixed,parm,is
      save /model/

c     Speed of light in km/s
      data clicht/ 2.99792e5 /
c     Translation FWHM-sigma cst=1/(2*sqrt(ln 2))
      data cst /0.600561/

c     find the shell to a given radius
c     rx has to be positive !
      ind=ilocat(rs,is+1,rx)
c     Compute the values
c     It has to be guaranteed that the exponents in the core are 0!
      do 300 j=1,6
      if(abs(parm(ind,j,2)).lt.0.01) then
        pp(j)=parm(ind,j,1)
      else
        pp(j)=parm(ind,j,1)*(dble(rx/rs(ind-1)))**dble(parm(ind,j,2))
      endif
300   continue
c     ************************************************************
c     Adapt field treatment to separate variables
c     ************************************************************
      pp1=pp(1)
      pp2=pp(2)
c     Compute molecule density from depletion ratios
      pp3=pp(3)*pp1
c     ************************************************************
c     Special treatment for HII region, detect by negative rc
c     ************************************************************
      if (pp(5).lt.0.0) then
      pp4=max(pp(4),1.0)
      pp5=0.0
      pp6=pp(6)
      else
c     ************************************************************
c     The turbulent velocity is input as FWHM. The thermal
c     velocity has to be added to obtain the complete random part.
c     ************************************************************
      pp4=sqrt(vtherm(pp2)+(cst*pp(4))**2)/clicht
c     Get the clump size from the velocity correlation length
      pp5=pp(5)/pp4
      pp6=pp(6)/clicht
      endif
      return
      end

c     **************************************************************
c     Routine to check for line overlap
c     At the moment, only a warning is produced. Later a regular
c     overlap treatment has to be added
c     **************************************************************
      subroutine overlap(vmax,lin)
      implicit none
      real vmax,deltav,cfreq,fr,vhig,vlow
      integer lin,i,j

      include 'fsizes.inc'
      real molmass,energy(nlev),stgew(nlev)
      character*10 levnam(nlev)
      integer lfrom(ntra),lto(ntra)
        common /molecule/ molmass, energy, stgew, levnam
      common /lintable/ lfrom, lto
      save /molecule/, /lintable/


11    format(/,' Line overlap detected between:')
12    format(' transition ',i3,' from ',A10,' to ',A10,' at ',
     %           1PE12.5,' Hz')
      deltav=2.0*vmax
      do i=1,lin-1
       fr=cfreq(i)
       vlow=fr*(1.0-deltav)
       vhig=fr*(1.0+deltav)
       do j=i+1,lin
        fr=cfreq(j)
        if ((fr.lt.vhig).and.(fr.gt.vlow)) then
         write(6,11)
         write(6,12) i,levnam(lfrom(i)),levnam(lto(i)),cfreq(i)
         write(6,12) j,levnam(lfrom(j)),levnam(lto(j)),fr
        endif
       enddo
      enddo
      return
      end

c     **************************************************************
c     Routine to determine additional geometry parameters required
c     to calculate the frequency scale in the radiative transfer
c     computations
c     **************************************************************
      subroutine initfreq(vs,vt,nn,vmax,dny,kerr)
      implicit none

c     input parameter
      include 'fsizes.inc'
      real vt(ns),vs(ns)
      integer nn

c     output parameter
      real dny(ns),vmax
      integer kerr
      
c     internal variables
      real negexp,negint,dvgauss,epsn,epsv,hyster,minspac,minedge
      real dvmin
      integer i,nff
      common /accgrid/ epsn,epsv,dvgauss,hyster,minspac,minedge
      common /accgauss/ negexp,negint
      save /accgrid/,/accgauss/

c     ***************************************************************
c     Find the minimum line width for all radii outside of r(i)
c     Find the extreme systematic velocities
c     ***************************************************************
      kerr=0
      dvmin=1e30
      vmax=0.0
      do 100 i=nn,2,-1
      dvmin=min(dvmin,vt(i))
      dny(i)=dvmin*dvgauss
      vmax=max(vmax,abs(vs(i))+vt(i)*negint)
100   continue
      dny(1)=dny(2)
c     **************************************************************
c     Extreme velocities - check whether nf is sufficient
c     **************************************************************
      nff=int(vmax/dny(1))+1
      if (nff.gt.nf) kerr=1
      return
      end

c     ************************************************************
c     Subroutine to set up the first guess for the level
c     populations at the different radial points
c     Superfluid levels are only removed in cases 1 and 3
c     ************************************************************
      subroutine initn(inum,r,vs,vt,rc,t,dh,dco,uext,dd,nn,lin,lev)
      implicit none
      integer inum,nn,lin,lev

      include 'fsizes.inc'
      real r(ns),vs(ns),vt(ns),t(ns),dh(ns),dco(ns),rc(ns)
      real*8 dd(ns,nlev)
      real uext(ntra),utot(ntra)

      real*8 dnn(nlev),sumdn,probtherm
      integer i,l,istat

      if (abs(inum).eq.1) then
c     ************************************************************
c     Thermal population
c     ************************************************************
      do 100 i=1,nn
      sumdn=0.0d0
      do 150 l=1,lev
      dnn(l)=probtherm(t(i),l)
      sumdn=sumdn+dnn(l)
150   continue
      do 180 l=1,lev
       dd(i,l)=dnn(l)/sumdn
180   continue
100   continue
      call removelevels(dd,nn,lin,lev,istat)
      else
c     ************************************************************
c     Use the same routines as in the radiative transfer
c     for computing the populations without radiation
c     (only the cosmic background is taken)
c     ************************************************************
c     Take central HII region into account (if present)
      if (vt(1).gt.1.0) then
      do 250 i=2,nn
       call ufromhii(r(i),r(2),uext,utot,lin)
       call matrixup(t(i),dh(i),utot,dnn,lin,lev)
       do 260 l=1,lev
        dd(i,l)=dnn(l)
260    continue
250   continue
      do 270 l=1,lev
       dd(1,l)=dd(2,l)
270   continue
      else
c     No HII region, only external radiation
      do 200 i=1,nn
       call matrixup(t(i),dh(i),uext,dnn,lin,lev)
       do 220 l=1,lev
        dd(i,l)=dnn(l)
220    continue
200   continue
      endif
      if (abs(inum).eq.3) then
c     ************************************************************
c     Sobolev program for the most sophisticated initial guess
c     The solution without radiation is used as first guess.
c     ************************************************************
      call sobolev(r,vs,vt,rc,t,dh,dco,uext,dd,nn,lin,lev)
      call removelevels(dd,nn,lin,lev,istat)
      endif
      endif
      return
      end

c     ************************************************************
c     Subroutine to estimate the energy density at a given point
c     ignoring all radiative transfer but including the HII region
c     ************************************************************
      subroutine ufromhii(r,rc,uext,uout,lin)
      implicit none
      integer lin

      include 'fsizes.inc'
      real uext(ntra),uout(ntra),s(ntra),kp(ntra)
      real r,rc,tau,muc,dil
      integer l

      common /hiidata/ kp,s
      save /hiidata/

      muc=sqrt(1.0-(rc/r)**2)
      do 100 l=1,lin
      tau=2.0*kp(l)*r
      if (tau.lt.0.01) then
      dil=1.0
      else
      dil=muc+(1.0-muc)*(1.0-exp(-tau))/tau
      endif
      uout(l)=0.5*((uext(l)+s(l))+(uext(l)-s(l))*dil)
100   continue
      return
      end

c     ***********************************************************
c     Subroutine to read precomputed levels from a file as input
c     for another radiative transfer problem
c     ***********************************************************
      subroutine filedens(r,t,dh,dco,vs,vt,rc,dd,nrad,lev,kerr)
      implicit none

      include 'fsizes.inc'
      real r(ns),t(ns),dh(ns),dco(ns),vt(ns),vs(ns),rc(ns)
      real*8 dd(ns,nlev),rat,zerolev
      integer nrad,lev,kerr

      real rf(ns),tb,outrad,dh1,dh2,dco1,t1,vt1,vt2,rc1,vs1,vs2
      integer i,j,l,i0,nn,np,hlev,ilocat
      character*250 mol,name

c     parameters in common
      integer ishell,fixed(0:ns)
      real rs(0:ns),parm(ns,6,2)
      real epsv,epsn,dvgauss,hyster,minspac,minedge,thick,small,ldelta
      common /model/ rs,fixed,parm,ishell
      common /accgrid/ epsn,epsv,dvgauss,hyster,minspac,minedge
      save /model/,/accgrid/
c     The limit for the smallest layer will not be used very often,
c     so that it is not treated as a normal numerical parameter
      data small / 0.02/
      data zerolev / 1.0d-14/

11    format(/,' Discontinuity detected between shells',i3,' and',i3)
c     ************************************************************
c     First read the data table
c     ************************************************************
      call readdens(mol,tb,r,vs,vt,rc,dco,dd,nrad,hlev,name,kerr)
      if (kerr.ne.0) return

c     ************************************************************
c     Now we have to check the radial spacings. This is done
c     in the same way as by initrad - parts of initrad are
c     repeated here.
c     ************************************************************
c     At first, the minimum spacing (minspac) has to be defined
c     ************************************************************
      minspac=1e30
      do j=1,ishell
        thick=rs(j)-rs(j-1)
        if (thick.lt.minspac) minspac=thick
      enddo
      minspac=small*minspac
      ldelta=hyster*log(epsn)
c     ************************************************************
c     Now look for fixed radii which have to be matched
c     ************************************************************
c     The core is treated separately
c     No discontinuity can occur at the most inner shell radius
      nn=2
      rf(1)=0.0
      rf(2)=rs(1)
c     Look for discontinuities at fixed shell boundaries
      do j=2,ishell-1
        if (fixed(j).gt.0) then
        nn=nn+1
        np=nn+1
        if (np.ge.ns) then
          kerr=1
          return
        endif
        rf(nn)=rs(j)-0.5*minspac
        call values(rf(nn),dh1,dco1,t1,vt1,rc1,vs1)
        rf(np)=rs(j)
        call values(rf(np),dh2,dco1,t1,vt2,rc1,vs2)
        if ((abs(log(dh2/dh1)).gt.ldelta).or.
     %        (abs(vs2/vt2-vs1/vt1).gt.epsv)) then
          write(6,11) j-1,j
          nn=np
        else
          rf(nn)=rf(np)
        endif
        endif
      enddo
      nn=nn+1
      rf(nn)=rs(ishell)
c     ************************************************************
c     Now, we throw away points outside of rs(ishell)
c     I assume that the file is written by the program so that the
c     smallest radius is 0 -> no additional check for underflow
c     ************************************************************
      nrad=nrad+1
      outrad=rs(ishell)+0.25*minspac
600     nrad=nrad-1
      if (r(nrad).gt.outrad) goto 600
c     ************************************************************
c     Look for entries in the table for the fixed radii
c     Create them if necessary
c     Special treatment for radii beyond the table
c     ************************************************************
      outrad=r(nrad)-0.25*minspac
      do j=1,nn
       i0=ilocat(r,nrad,rf(j))
       if (rf(j).gt.outrad) then
         if (abs(rf(j)-r(i0+1)).gt.0.25*minspac) then
           nrad=nrad+1
           if (nrad.gt.ns) then
             kerr=1
             return
           endif
           i0=i0+1
           r(nrad)=rf(j)
           do l=1,hlev
             dd(nrad,l)=dd(i0,l)
           enddo
         endif
       else
         if (abs(rf(j)-r(i0)).gt.0.25*minspac) then
           nrad=nrad+1
           if (nrad.gt.ns) then
             kerr=1
             return
           endif
           i0=i0+1
           do i=nrad,i0+1,-1
            r(i)=r(i-1)
            do l=1,hlev
              dd(i,l)=dd(i-1,l)
            enddo
           enddo
           r(i0)=rf(j)
           rat=(r(i0)-r(i0-1))/(r(i0+1)-r(i0-1))
           do l=1,hlev
             dd(i0,l)=(1.0d0-rat)*dd(i0-1,l)+rat*dd(i0+1,l)
           enddo
         endif
       endif
      enddo
c     ************************************************************
c     Radial table complete - fill with physical parameter values
c     Fill the levels from hlev+1 to lev with 1d-14
c     ************************************************************
      do i=1,nrad
        call values(r(i),dh(i),dco(i),t(i),vt(i),rc(i),vs(i))
        do l=hlev+1,lev
          dd(i,l)=zerolev
        enddo
      enddo

      return
      end

