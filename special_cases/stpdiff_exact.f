c     ***********************************************************
c     Subroutine to carry out radiative transport computation
c     It determines the energy density at all given radial points
c     ***********************************************************
      subroutine stpmic(r,v,dlb,dvx,vx,caps,ss,uext,ulev,nn,lin,ke)
      implicit none
      integer nn,lin,ke

      include 'fsizes.inc'
      real r(ns),v(ns),dlb(ns),dvx(ns),vx
      real caps(ns,ntra),ss(ns,ntra)
      real uext(ntra),ulev(ns,ntra)

      real uny(ntra)
      real cap(-nf:nf,ntra),iplus(-nf:nf,ntra)
      real cp(ntra),sca(ntra)
      real cpo(ntra),sco(ntra)
      real p2(nstretch)
      real zz(ns,2*nstretch),ufull(ns,2*nstretch,ntra)
      integer nz(ns),indr(nstretch)

      real y,vs,dv,vv,yo,vso,dvo,vo,rr,dr,dy2,dny,fak
      integer i,j,l,j1,j2,it
      integer nnp,idown,kk,nff,itstp,iback

c     Always allocate ufull statically to avoid stack problems
      save ufull

c     *******************************************************
c     Initialize a field of displacement parameters depending
c     on the grid of radial points so that the z integration
c     has a sufficient efficiency
c     *******************************************************
      call makps(r,nz,nn,p2,indr,nnp,ke)
c     write(6,'('' Number of rays in radiative transfer: '',i3)')nnp
      if (ke.ne.0) return

c     ********************************************************
c     If the last ray is touching - separated treatment
c     ********************************************************
      if (indr(nnp).eq.2*nn) then
      dny=dvx(nn)
      nff=int(vx/dny)+1
      do 200 l=1,lin
      do 200 j=-nff,nff
200   iplus(j,l) = uext(l)
      call unyint(iplus,uny,lin,0.0,dlb(nn),dny)
      zz(nn,nnp)=0.0
      do 250 l=1,lin
      ufull(nn,nnp,l)=uny(l)
250   continue
      nnp=nnp-1
      endif

c     ********************************************************
c     This is the beginning of the p (displacement parameter) 
c     loop. 
c     ********************************************************
      do 400 i=1,nnp
      if (indr(i).gt.nn) then
      idown=indr(i)-nn+1
      else
      idown=indr(i)
      endif

      dny=dvx(idown-1)
      nff=int(vx/dny)+1

c     ********************************************************
c     Beginning of the z-loop (integration along a ray)
c     The organization of the p and z loops has to guaranty
c     that growing z values for the same r are treated
c     in a monotonically increasing order
c     ********************************************************
c     Special treatment for the most outer left point
c     ********************************************************
c     Geometry parameters
      kk=nn
      rr = r(kk)
      yo = -sqrt(rr**2-p2(i))
      vo = v(kk)
      vso = vo*yo /rr
      dvo = dlb(kk)
c     Background for all frequencies
      do 500 l=1,lin
      do 500 j=-nff,nff
500   iplus(j,l) = uext(l)
c     Initialize kappa and source function fields
      do 510 l=1,lin
      cpo(l)=caps(kk,l)
510   sco(l)=ss(kk,l)
      call leftkp(cpo,cap,j1,j2,lin,vso,dvo,dny)
c     Frequency integration for all lines
      call unyint(iplus,uny,lin,vso,dvo,dny)
c     **********************************************************
c     Store the results into large fields
c     In contrast to their production in lines with constant p
c     they are stored in fields lines with constant r.
c     The first index characterizes the radius r. 
c     The second index is z.
c     **********************************************************
      zz(kk,i)=yo
      do 550 l=1,lin
      ufull(kk,i,l)=uny(l)
550   continue

c     ********************************************************
c     Now integration - treat all points through the cloud
c     Split into two parts - positive and negative z values
c     ********************************************************
c     First part - negative z
c     ********************************************************
      do  1600  kk = nn-1,idown,-1
c     Initialize geometry parameters
      dr= r(kk)-rr
      rr= r(kk)
      y = -sqrt(rr**2-p2(i))
      vv= v(kk)
      vs= vv*y/rr
      dv= dlb(kk)
      do 610 l=1,lin
      cp(l)=caps(kk,l)
610   sca(l)=ss(kk,l)
c     *******************************************************
c     Radiative transfer - Multiple steps if the velocity
c     gradient is too large
c     *******************************************************
      call firststep(yo,vo,vso,dvo,cpo,sco,
     %                 y,vv,vs,dr,dv,cp,sca,p2(i),dy2,lin,itstp)
      call transfer(cpo,sca,iplus,cap,sco,j1,j2,lin,vso,dvo,dny,dy2)
      if (itstp.gt.0) then
      do 640 it=1,itstp
      call nextstep(vso,dvo,cpo,sca,dy2,lin)
      call transfer(cpo,sca,iplus,cap,sco,j1,j2,lin,vso,dvo,dny,dy2)
640   continue
      endif
c     ******************************************************
c     Frequency integration for all lines
c     Store the results into large fields
c     ******************************************************
      call unyint(iplus,uny,lin,vso,dvo,dny)
      zz(kk,i)=y
      do 650 l=1,lin
      ufull(kk,i,l)=uny(l)
650   continue
1600  continue

c     ********************************************************
c     Special treatment for zero z
c     ********************************************************
      kk=idown-1
      dv=dlb(kk)
      vv=v(kk)
c     ********************************************************
c     Treat the central HII region if there is one
c     ********************************************************
      if (dv.ge.1.0) then
      dy2 = -2.0*yo
      do 910 l=1,lin
      cp(l)=caps(kk,l)
910   sca(l)=ss(kk,l)
      call ctransfer(cp,sca,iplus,nff,lin,dy2)
c     ********************************************************
c     Treatment of the point at the other side of the HII region
c     by a 0 length transfer in the next transfer step
c     ********************************************************
      yo = -yo
      vo = v(kk+1)
      vso = vo*yo /rr
      dvo = dlb(kk+1)

      else
c     ********************************************************
c     Normal radiative transfer when no HII region is touched
c     ********************************************************
c     Geometry parameters
      y = 0.0
      vs = 0.0
      if (indr(i).gt.nn) then
c     The zero point is given itself
      dr = r(kk)-rr
      rr = r(kk)
      do 810 l=1,lin
      cp(l)=caps(kk,l)
810   sca(l)=ss(kk,l)
      else
c     interpolate the values at z=0
      rr=sqrt(p2(i))
      dr=rr-r(idown)
      fak=(rr-r(kk))/(r(idown)-r(kk))
      dv=fak*(dlb(idown)-dv)+dv
      vv=fak*(v(idown)-vv)+vv
      do 820 l=1,lin
      cp(l)=fak*(caps(idown,l)-caps(kk,l))+caps(kk,l)
      sca(l)=(fak*(ss(idown,l)*caps(idown,l)
     %    -ss(kk,l)*caps(kk,l))+ss(kk,l)*caps(kk,l))/cp(l)
820   continue
      endif
c     *******************************************************
c     Radiative transfer
c     *******************************************************
c     Extra treatment for p=0 since the velocity has
c     to be constant there with a jump at z=0
c     *******************************************************
      if (i.eq.1) then
      vs=vso
      call firststep(yo,vo,vso,dvo,cpo,sco,
     %                 y,vv,vs,dr,dv,cp,sca,p2(i),dy2,lin,itstp)
      call transfer(cpo,sca,iplus,cap,sco,j1,j2,lin,vso,dvo,dny,dy2)
c     New kappa field in case of a velocity jump
      if (vso.ne.0.0) then
      vso=-vso
      call leftkp(cpo,cap,j1,j2,lin,vso,dvo,dny)
      endif
      else
c     *******************************************************
c     Multiple steps if the velocity gradient is too large
c     *******************************************************
      call firststep(yo,vo,vso,dvo,cpo,sco,
     %                 y,vv,vs,dr,dv,cp,sca,p2(i),dy2,lin,itstp)
      call transfer(cpo,sca,iplus,cap,sco,j1,j2,lin,vso,dvo,dny,dy2)
      if (itstp.gt.0) then
      do 840 it=1,itstp
      call nextstep(vso,dvo,cpo,sca,dy2,lin)
      call transfer(cpo,sca,iplus,cap,sco,j1,j2,lin,vso,dvo,dny,dy2)
840   continue
      endif
      endif
      endif

c     ******************************************************
c     This radiative transfer point is only used in the
c     intensity integration if it falls on a given r-shell	
c     ********************************************************
c     To enable uniform treatment assign the same fields as 
c     for normal points to the central point in an HII region
c     ********************************************************
      if (indr(i).gt.nn) then
c     Frequency integration for all lines
      call unyint(iplus,uny,lin,vso,dvo,dny)
c     Store the results into large fields
      zz(kk,i)=0.0
      do 850 l=1,lin
      ufull(kk,i,l)=uny(l)
850   continue
      endif

c     ********************************************************
c     Second part - positive z
c     ********************************************************
      do  1700  kk = idown,nn
c     Geometry parameters
      dr=r(kk)-rr
      rr=r(kk)
      y = -zz(kk,i)
      vv= v(kk)
      vs= vv*y/r(kk)
      dv= dlb(kk)
      do 710 l=1,lin
      cp(l)=caps(kk,l)
710   sca(l)=ss(kk,l)
c     *******************************************************
c     Radiative transfer - Multiple steps if the velocity
c     gradient is too large
c     *******************************************************
      call firststep(yo,vo,vso,dvo,cpo,sco,
     %                 y,vv,vs,dr,dv,cp,sca,p2(i),dy2,lin,itstp)
      call transfer(cpo,sca,iplus,cap,sco,j1,j2,lin,vso,dvo,dny,dy2)
      if (itstp.gt.0) then
      do 740 it=1,itstp
      call nextstep(vso,dvo,cpo,sca,dy2,lin)
      call transfer(cpo,sca,iplus,cap,sco,j1,j2,lin,vso,dvo,dny,dy2)
740   continue
      endif
c     **********************************************************
c     Frequency integration for all lines,
c     Store the results into large fields
c     In the part of positive z, the values have to be stored
c     in backward order into the r lines in the fields.
c     **********************************************************
      call unyint(iplus,uny,lin,vso,dvo,dny)
      iback=nz(kk)+1-i
      zz(kk,iback)=y
      do 750 l=1,lin
      ufull(kk,iback,l)=uny(l)
750   continue
1700  continue
c     **************************************************************
c     End of the z loop
c     **************************************************************
400   continue

c     **************************************************************
c     Calculate integral over mean intensity times profile function
c     **************************************************************
      call uzint(ufull,zz,nz,ulev,nn,lin)
      return
      end

c     **************************************************************
c     Subroutine to compute the physical parameters at the end of 
c     the first radiative transfer step within one interval.
c     It checks how many steps have to be done there.
c     **************************************************************
      subroutine firststep(yo,vo,vso,dvo,cpo,sco,y,v,vs,dr,dv,cp,sc,
     %                       p,z,lev,istp)
      implicit none
      integer lev,istp

      include 'fsizes.inc'
      real cpo(ntra),sco(ntra),cp(ntra),sc(ntra)
      real y,yo,v,vo,vs,vso,dv,dvo,z,p,dr
      real deltav,epsz

      real ddcp(ntra),ddsc(ntra),ddv,ddt,ddz
      real aa,bb,cadd,cpre,vzgeo,pp
      real fi,rv,rvo,dd,dz
      integer l
      logical lin

c     Internal variables are transmitted to nextstep for all
c     later steps
      common /step/ ddcp,ddsc, ddv,ddt,ddz,lin
      common /nonlin/ bb,cadd,cpre,vzgeo,pp
      common /acctrans/ deltav,epsz
      save /step/,/nonlin/,/acctrans/

c     *******************************************************
c     First call - compute the number of steps
c     *******************************************************
      ddv=vs-vso
c     *******************************************************
c     Single step transfer
c     *******************************************************
      if (abs(ddv).le.deltav*dv) then
c     New step size
      z=0.5*(y-yo)
c     Assign values for the recent point
      do 610 l=1,lev
610   cpo(l)=cp(l)
      yo=y
      vo=v
      vso=vs
      dvo=dv
      istp=0
      return
      endif
c     *******************************************************
c     Multiple steps - look whether linear approx is possible
c     *******************************************************
      ddz=y-yo
      if (p*v*vo.eq.0.0) goto 1000
      aa=(v-vo)/dr
      rvo=vso/vo
      rv=vs/v
      if (y.gt.0.0) then
      dd=(vs-vso)-(y-yo)*(vs/y+rv**2*(aa-vs/y))
      else
      dd=(vs-vso)-(y-yo)*(vso/yo+rvo**2*(aa-vso/yo))
      endif
c     ******************************************************
c     Full treatment of the geometry for strong nonlinearity
c     ******************************************************
      if (abs(0.5*dd).gt.deltav*dv) then
c     Globally necessary constants
      bb=0.5*(v+vo)
      fi=rv-rvo
      istp=int(abs(bb*fi)/(deltav*dv))+1
      ddv=bb*fi/istp
      ddt=(dv-dvo)/ddz
      do 820 l=1,lev
      ddcp(l)=(cp(l)-cpo(l))/ddz
820   ddsc(l)=(sc(l)-sco(l))/ddz
c     Constants for the velocity correction
      pp=p
      cadd=(v-bb)*rvo
      cpre=(v-bb)*(rv+rvo)/ddz
c     *****************************
c     Compute the first step
c     *****************************
c     Step size
      vzgeo=bb*rvo+ddv
      ddz=sign(sqrt(pp*vzgeo**2/(bb**2-vzgeo**2)),yo+y)
      dz=ddz-yo
      z=0.5*dz
c     Correct the velocity at that point
      cadd=cadd-cpre*dz
      vso=vzgeo-cadd
c     Assign the linearily interpolated values
      do 840 l=1,lev
      cpo(l)=cpo(l)+ddcp(l)*dz
840   sc(l)=sco(l)+ddsc(l)*dz
      dvo=dvo+ddt*dz
c     Values which are not changed/needed in the following steps
      yo=y
      vo=v
      istp=istp-1
      lin=.false.
      return
      endif
c     ******************************************************
c     If the nonlinearity is small take the linear approx
c     Linear approx is also taken for p=0 or v,vo=0
c     ******************************************************
1000  istp=int(abs(ddv)/(deltav*dv))+1
c     Constants for the next steps
      ddv=(vs-vso)/istp
      ddt=(dv-dvo)/istp
      ddz=ddz/istp
      do 620 l=1,lev
      ddcp(l)=(cp(l)-cpo(l))/istp
620   ddsc(l)=(sc(l)-sco(l))/istp
c     *****************************
c     Compute the first step
c     *****************************
c     New step size
      z=0.5*ddz
c     Assign the new values
      do 640 l=1,lev
      cpo(l)=cpo(l)+ddcp(l)
640   sc(l)=sco(l)+ddsc(l)
      vso=vso+ddv
      dvo=dvo+ddt
c     Values which are not changed/needed in the following steps
      yo=y
      vo=v
      istp=istp-1
c     Detect this path in nextstep
      lin=.true.
      return
      end

c     *********************************************************
c     Subroutine to compute the physical parameters at the
c     end of the following radiative transfer steps
c     *********************************************************
      subroutine nextstep(vso,dvo,cpo,sco,z,lev)
      implicit none
      integer lev

      include 'fsizes.inc'
      real cpo(ntra),sco(ntra)
      real vso,dvo,z
      real y,dz
      integer l

c     Internal variables computed in firststep 
      real ddcp(ntra),ddsc(ntra),ddv,ddt,ddz
      real bb,cadd,cpre,vzgeo,pp
      logical lin
      common /step/ ddcp,ddsc, ddv,ddt,ddz,lin
      common /nonlin/ bb,cadd,cpre,vzgeo,pp
      save /step/,/nonlin/

c     *********************************************************
c     If linear approximation is possible
c     *********************************************************
      if (lin) then
c     Update the physical parameters with firststep variables
      z=0.5*ddz
      do 640 l=1,lev
      cpo(l)=cpo(l)+ddcp(l)
640   sco(l)=sco(l)+ddsc(l)
      vso=vso+ddv
      dvo=dvo+ddt
c     *********************************************************
c     Nonlinear - compute the possible step size
c     *********************************************************
      else
c     Step size
      vzgeo=vzgeo+ddv
      y=sign(sqrt(pp*vzgeo**2/(bb**2-vzgeo**2)),ddz)
      dz=y-ddz
      ddz=y
      z=0.5*dz
c     Correct the velocity at that point
      cadd=cadd-cpre*dz
      vso=vzgeo-cadd
c     Assign the linearily interpolated values
      do 840 l=1,lev
      cpo(l)=cpo(l)+ddcp(l)*dz
840   sco(l)=sco(l)+ddsc(l)*dz
      dvo=dvo+ddt*dz
      endif
      return
      end

c     ***********************************************************
c     Initialization of the kappa field on the left cloud edge 
c     ***********************************************************
      subroutine leftkp(cp,cappa,j1,j2,lev,vs,dv,dny)
      implicit none
      integer lev,j1,j2

      include 'fsizes.inc'
      real cp(ntra),cappa(-nf:nf,ntra)
      real vs,dv,dny
      real negexp,negint
      real phi,prof
      integer j,l

      common /accgauss/ negexp,negint
      save /accgauss/

c     Frequency range
      j1=int((vs-dv*negexp)/dny)-1
      j2=int((vs+dv*negexp)/dny)+1
c     Frequency dependent quantities
      do 500 j=j1,j2
      prof = phi(j*dny,vs,dv)
      do 500 l=1,lev
      cappa(j,l)  = cp(l) * prof
500   continue
      return
      end

c     ***********************************************************
c     Subroutine for radiative transfer between two points
c     ***********************************************************
      subroutine transfer(cp,sc,iplus,cap,ss,j1,j2,lev,vs,dv,dny,z)
      implicit none
      integer lev,j1,j2

      include 'fsizes.inc'
      real cp(ntra),sc(ntra),ss(ntra)
      real cap(-nf:nf,ntra)
      real iplus(-nf:nf,ntra), pp(-nf:nf)
      real vs,dv,dny,z
      real negexp,negint,explin,softzone
      real phi,prof,tau1,tau2,tau,tauh1,tauh2,st1,st2
      real g1,ex1,ex2
      real ihlp1,ihlp2
      integer jn1,jn2,jj1,jj2,j,l

      common /accstep/ explin,softzone
      common /accgauss/ negexp,negint
      save /accstep/,/accgauss/

c     *******************************************************
c     Determine indices of the frequency range
c     *******************************************************
      jn1=int((vs-dv*negexp)/dny)-1
      jn2=int((vs+dv*negexp)/dny)+1
c     Find the overlapping interval
      jj1=max(j1,jn1)
      jj2=min(j2,jn2)
c     *******************************************************
c     Treat the non-overlapping left end
c     *******************************************************
      if (jn1.lt.jj1) then
      do 510 j=jn1,jj1-1
      prof = phi(j*dny,vs,dv)
      do 510 l=1,lev
      cap(j,l) = cp(l)*prof
      tau = z*cap(j,l)
      iplus(j,l)=(iplus(j,l)+sc(l)*tau)/(1.0+tau)
510   continue
      else if (j1.lt.jj1) then
      do 520 j=j1,jj1-1
      do 520 l=1,lev
      tau = z*cap(j,l)
      iplus(j,l)=(iplus(j,l)+ss(l)*tau)/(1.0+tau)
520   continue
      endif
c     ******************************************************
c     Treat the main, fully overlapping part
c     Only in this part, it is controlled whether the lines 
c     are optically thick, i.e. it is assumed that they are
c     thin within the non-overlapping wings
c     ******************************************************
      prof=2.0*z*phi(0.0,0.0,dv)
      do 500 j=jj1,jj2
500   pp(j) = phi(j*dny,vs,dv)
      do 600 l=1,lev
c     ******************************************************
c     When even the centre of the line is optically 
c     thin the time consuming if in the inner loop shall
c     be avoided
c     ******************************************************
c     Starting from version 2.0, I prevent exponential maser
c     amplification for stability reasons. This is unphysical
c     in situations with real masers - no final solution yet.
c     Uncomment the old lines below for maser amplification.
c     ******************************************************
c     if (abs(cp(l)*prof).lt.explin) then
      if (cp(l)*prof.lt.explin) then
      do 650 j=jj1,jj2
      tau1 = z*cap(j,l)
      cap(j,l) = cp(l)*pp(j)
      tau2 = z*cap(j,l)
      tauh1 = 0.5*tau1
      st1 = ss(l)*tau1
      tauh2 = 0.5*tau2
      st2 = sc(l)*tau2
      iplus(j,l) = (iplus(j,l)+st1+st2+st1*tauh1+st2*tauh2
     %	 +st1*tauh2*0.33333+st2*tauh1*1.66667)
     %   /(1.0+tau1+tau2+tau1*tauh1+tau2*tauh2+tau1*tau2)
650   continue
      else
c     ******************************************************
c     In each step there is independently asked whether
c     the line is optically thick
c     It has to be guaranteed that none of the cappa values
c     is 0 in this case. (Should be always fulfilled for
c     the overlapping part.)
c     ******************************************************
      do 680 j=jj1,jj2
      tau1 = z*cap(j,l)
      cap(j,l) = cp(l)*pp(j)
      tau2 = z*cap(j,l) 
      tau = tau1+tau2
c     rel = abs(tau)
c     if (rel.lt.explin) then
      if (tau.lt.explin) then
c     Thin case - linear approx.
      tauh1 = 0.5*tau1
      st1 = ss(l)*tau1
      tauh2 = 0.5*tau2
      st2 = sc(l)*tau2
      iplus(j,l) = (iplus(j,l)+st1+st2+st1*tauh1+st2*tauh2
     %	 +st1*tauh2*0.33333+st2*tauh1*1.66667)
     %   /(1.0+tau1+tau2+tau1*tauh1+tau2*tauh2+tau1*tau2)
c     else if (rel.lt.softzone) then
      else if (tau.lt.softzone) then
c     Small transition regime - superposition
      tauh1 = 0.5*tau1
      st1 = ss(l)*tau1
      tauh2 = 0.5*tau2
      st2 = sc(l)*tau2
      ihlp1 = (iplus(j,l)+st1+st2+st1*tauh1+st2*tauh2
     %	 +st1*tauh2*0.33333+st2*tauh1*1.66667)
     %   /(1.0+tau1+tau2+tau1*tauh1+tau2*tauh2+tau1*tau2)
      ex1=exp(-0.33333*tau1)
      ex2=exp(-0.33333*tau2)
      g1=0.5*(ss(l)-sc(l))
      ihlp2=sc(l)+ex2*ex2*(g1+ex2*ex1*
     %              (g1+ex1*ex1*(iplus(j,l)-ss(l))))
c     iplus(j,l)=(ihlp1*(softzone-rel)+ihlp2*(rel-explin))/
c     %              (softzone-explin)
      iplus(j,l)=(ihlp1*(softzone-tau)+ihlp2*(tau-explin))/
     %              (softzone-explin)
      else
c     Optically thick case - exponential regime
      ex1=exp(-0.33333*tau1)
      ex2=exp(-0.33333*tau2)
      g1=0.5*(ss(l)-sc(l))
      iplus(j,l)=sc(l)+ex2*ex2*(g1+ex2*ex1*
     %              (g1+ex1*ex1*(iplus(j,l)-ss(l))))
      endif
680   continue
      endif
600   continue
c     *******************************************************
c     Treat the non-overlapping right end
c     *******************************************************
      if (jn2.gt.jj2) then
      do 560 j=jj2+1,jn2
      prof = phi(j*dny,vs,dv)
      do 560 l=1,lev
      cap(j,l) = cp(l)*prof 
      tau = z*cap(j,l)
      iplus(j,l)=(iplus(j,l)+sc(l)*tau)/(1.0+tau)
560   continue
      else if (j2.gt.jj2) then
      do 570 j=jj2+1,j2
      do 570 l=1,lev
      tau = z*cap(j,l)
      iplus(j,l)=(iplus(j,l)+ss(l)*tau)/(1.0+tau)
570   continue
      endif
c     *******************************************************
c     New parameters
c     *******************************************************
      j1=jn1
      j2=jn2
      do 580 l=1,lev
580   ss(l)=sc(l)
      return
      end

c     ***********************************************************
c     Subroutine for the radiative transfer in the continuum
c     ***********************************************************
      subroutine ctransfer(cp,sc,iplus,nff,lev,z)
      implicit none
      integer lev,nff

      include 'fsizes.inc'
      real cp(ntra),sc(ntra)
      real iplus(-nf:nf,ntra)
      real z,tau,tauq,xx
      real explin,softzone
      integer j,l

      common /accstep/ explin,softzone
      save /accstep/

c     ******************************************************
c     Treat all levels and frequencies
c     ******************************************************
      do 600 l=1,lev
      tau=cp(l)*z
      if (tau.lt.explin) then
c     ******************************************************
c     Linear approximation for optically thin ranges
c     ******************************************************
      do 650 j=-nff,nff
      tauq = 0.5*tau**2
      iplus(j,l) = (iplus(j,l)+sc(l)*(tau+tauq))/(1.0+tau+tauq)
650   continue
      else
c     ******************************************************
c     Exponential function for thick ranges
c     ******************************************************
      do 680 j=-nff,nff
      xx=exp(-tau)
      iplus(j,l)=xx*(iplus(j,l)-sc(l)) + sc(l)
680   continue
      endif
600   continue
      return
      end

c     ************************************************************
c     Subroutine to compute the efficient scaling in the 
c     displacement parameter p to cover the z scale for all
c     radial points
c     The output is p**2 instead of p!
c     indr contains the number of crossings with radii for all p
c     ************************************************************
c     Here, this is a purely geometrical approach requiring that
c     the relative z distance between two points is below epsz.
c     In a more sophisticated later version, this should be 
c     replaced by adapting the p scale to the intensity structure.
c     ************************************************************
      subroutine makps(r,nz,nn,p,indr,kk,kerr)
      implicit none
      integer nn,kk,kerr

      include 'fsizes.inc'
      real r(ns),p(nstretch)
      real deltav,epsz
      integer nz(ns),indr(nstretch)
      real relz,r2,numerrz
      integer i
      logical odd

      common /acctrans/ deltav,epsz
      save /acctrans/

c     Prefactor to avoid very close double points due to numerical
c     round-off errors
      data numerrz/ 1.02/

      kerr=0
      kk=1
      p(kk)=0.0
      indr(1)=1+nn
      nz(1)=1
c     ***************************************************************
c     Beginning of the radial loop - test for p field overflow
c     ***************************************************************
      do 100 i=2,nn
      odd=.false.
200   if (kk.eq.nstretch) then
      kerr=1
      return
      endif
      r2=r(i)**2
      relz=sqrt(1.0-p(kk)/r2)
c     Include a p line which crosses the circle with r(i)
      if (relz.gt.numerrz*epsz) then
      kk=kk+1
      p(kk)=r2*(1.0-(relz-epsz)**2)
      indr(kk)=i
      odd=.false.
      goto 200
c     Include a p line which touches the circle with r(i)
      else if (2.0*relz.gt.epsz) then
      kk=kk+1
      p(kk)=r2
      indr(kk)=i+nn
      odd=.true.
      endif
c     ***************************************************************
c     End of r loop - compute the number of crossings for each radius
c     ***************************************************************
      if (odd) then
      nz(i)=2*kk-1
      else
      nz(i)=2*kk
      endif
100   continue
      return
      end
