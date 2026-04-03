c     ***********************************************************
c     Subroutine to carry out radiative transport computation
c     for one transition
c     This is a reduced version of stpmic for only one line
c     ***********************************************************
      subroutine stpline(r,v,dlb,caps,ss,nn,uext,p,indr,iout,
     %                     tautot,nnp,nff,dny)
      implicit none
      integer nn,nnp,nff

      include 'fsizes.inc'
      real r(ns),v(ns),dlb(ns)
      real caps(ns),ss(ns)
      real uext,dny
      real iout(nstretch,-nf:nf),tautot(nstretch,-nf:nf),p(nstretch)
      integer indr(nstretch)

      real cap(-nf:nf),iplus(-nf:nf),tauaccu(-nf:nf)
      real cpo,sco,cp,sca

      real y,vs,vv,dv,yo,vso,vo,dvo,rr,dr,dy2,p2,fak
      integer i,j,j1,j2,it
      integer idown,kk,itstp,nnk

c     *******************************************************
c     Here, the field of displacement parameters has to be
c     initialized outside (psforline called by line)
c     ********************************************************
c     If the last ray is touching - separated treatment
c     ********************************************************
      if (indr(nnp).eq.2*nn) then
      do 200 j=-nff,nff
      tautot(nnp,j) = 0.0
       iout(nnp,j) = uext
200   continue
      nnk=nnp-1
      else
      nnk=nnp
      endif

c     ********************************************************
c     This is the beginning of the p (displacement parameter) 
c     loop. 
c     ********************************************************
      do 400 i=1,nnk
      if (indr(i).gt.nn) then
      idown=indr(i)-nn+1
      else
      idown=indr(i)
      endif

      p2=p(i)*p(i)

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
      yo = -sqrt(rr**2-p2)
      vo = v(kk)
      vso = vo*yo /rr
      dvo = dlb(kk)
c     Background for all frequencies
      do 500 j=-nff,nf
      tauaccu(j)= 0.0
      iplus(j) = uext
500   continue
c     Initialize kappa and source function fields
      cpo=caps(kk)
      sco=ss(kk)
      call leftkpl(cpo,cap,j1,j2,vso,dvo,dny)

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
      y = -sqrt(rr**2-p2)
      vv = v(kk)
      vs= vv*y/rr
      dv = dlb(kk)
      cp=caps(kk)
      sca=ss(kk)
c     *******************************************************
c     Radiative transfer - Multiple steps if the velocity
c     gradient is too large
c     *******************************************************
      call firstlstp(yo,vo,vso,dvo,cpo,sco,
     %                 y,vv,vs,dr,dv,cp,sca,p(i),dy2,itstp)
      call transfl(cpo,sca,iplus,tauaccu,cap,sco,j1,j2,vso,dvo,dny,dy2)
      if (itstp.gt.0) then
      do 640 it=1,itstp
      call nextlstp(vso,dvo,cpo,sca,dy2)
      call transfl(cpo,sca,iplus,tauaccu,cap,sco,j1,j2,vso,dvo,dny,dy2)
640   continue
      endif
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
      cp=caps(kk)
      sca=ss(kk)
      call ctransfl(cp,sca,iplus,tauaccu,nff,dy2)
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
      vs= 0.0
      if (indr(i).gt.nn) then
c     The zero z point is directly treated
      dr= r(kk)-rr
      rr= r(kk)
      cp=caps(kk)
      sca=ss(kk)
      else
c     interpolate the values at z=0
      rr=p(i)
      dr=rr-r(idown)
      fak=(rr-r(kk))/(r(idown)-r(kk))
      dv=fak*(dlb(idown)-dv)+dv
      vv=fak*(v(idown)-vv)+vv
      cp=fak*(caps(idown)-caps(kk))+caps(kk)
      sca=(fak*(ss(idown)*caps(idown)
     %   -ss(kk)*caps(kk))+ss(kk)*caps(kk))/cp
      endif
c     *******************************************************
c     Radiative transfer
c     *******************************************************
c     Extra treatment for p=0 since the velocity has
c     to be constant there with a jump at z=0
c     *******************************************************
      if (p(i).eq.0.0) then
      vs=vso
      call firstlstp(yo,vo,vso,dvo,cpo,sco,
     %                 y,vv,vs,dr,dv,cp,sca,p(i),dy2,itstp)
      call transfl(cpo,sca,iplus,tauaccu,cap,sco,j1,j2,vso,dvo,dny,dy2)
c     New kappa field in case of a velocity jump
      if (vso.ne.0.0) then
      vso=-vso
      call leftkpl(cpo,cap,j1,j2,vso,dvo,dny)
      endif
      else
c     *******************************************************
c     Multiple steps if the velocity gradient is too large
c     *******************************************************
      call firstlstp(yo,vo,vso,dvo,cpo,sco,
     %                 y,vv,vs,dr,dv,cp,sca,p(i),dy2,itstp)
      call transfl(cpo,sca,iplus,tauaccu,cap,sco,j1,j2,vso,dvo,dny,dy2)
      if (itstp.gt.0) then
      do 840 it=1,itstp
      call nextlstp(vso,dvo,cpo,sca,dy2)
      call transfl(cpo,sca,iplus,tauaccu,cap,sco,j1,j2,vso,dvo,dny,dy2)
840   continue
      endif
      endif
      endif

c     ********************************************************
c     Second part - positive z
c     ********************************************************
      do  1700  kk = idown,nn
c     Geometry parameters
      dr=r(kk)-rr
      rr=r(kk)
      y = sqrt(rr**2-p2)
      vv= v(kk)
      vs= vv*y/rr
      dv=dlb(kk)
      cp=caps(kk)
      sca=ss(kk)
c     *******************************************************
c     Radiative transfer - Multiple steps if the velocity
c     gradient is too large
c     *******************************************************
      call firstlstp(yo,vo,vso,dvo,cpo,sco,
     %                 y,vv,vs,dr,dv,cp,sca,p(i),dy2,itstp)
      call transfl(cpo,sca,iplus,tauaccu,cap,sco,j1,j2,vso,dvo,dny,dy2)
      if (itstp.gt.0) then
      do 740 it=1,itstp
      call nextlstp(vso,dvo,cpo,sca,dy2)
      call transfl(cpo,sca,iplus,tauaccu,cap,sco,j1,j2,vso,dvo,dny,dy2)
740   continue
      endif
1700  continue

c     **********************************************************
c     End of the z loop
c     Store the results into a large field
c     **********************************************************
      do 420 j=-nff,nff
      tautot(i,j) = tauaccu(j)
      iout(i,j) = iplus(j)
420   continue
400   continue
c     **********************************************************
      return
      end

c     **********************************************************
c     Subroutine to estimate the central optical depth at the
c     line centre, symmetry is expoited.
c     **********************************************************
      subroutine taucentral(r,v,dlb,caps,nn,tau)
      implicit none
      integer nn

      include 'fsizes.inc'
      real r(ns),v(ns),dlb(ns),caps(ns),tau

      real*8 kappa
      real phi,z,vs,vso,y,yo,dv,dvo,cpo,ddv,ddd,ddcp
      real deltav,epsz
      integer i,kk,itstp
      common /acctrans/ deltav,epsz
      save /acctrans/

c     *******************************************************
c     Initialize
c     *******************************************************
      kappa=0d0
      yo = -r(nn)
      vso = -v(nn)
      dvo = dlb(nn)
      cpo = caps(nn)*phi(0.0,vso,dvo)
c     *******************************************************
c     Scan for negative z
c     *******************************************************
      do 600 kk = nn-1,2,-1
      y = -r(kk)
      vs= -v(kk)
      dv = dlb(kk)
      ddv=vs-vso
      if (abs(ddv).le.deltav*dv) then
c     *******************************************************
c     One step
c     *******************************************************
      z=(y-yo)
      kappa=kappa+dble(cpo*z)
      cpo = caps(kk)*phi(0.0,vs,dv)
      kappa=kappa+dble(cpo*z)
      else
c     *******************************************************
c     Multiple steps
c     *******************************************************
      itstp=int(abs(ddv/(deltav*dv)))+1
      z=2.0*(y-yo)/itstp
      kappa = kappa+dble(0.5*cpo*z)
      ddv=ddv/itstp
      ddd=(dv-dvo)/itstp
      ddcp=(caps(kk)-caps(kk+1))/itstp
      do 500 i=1,itstp-1
      cpo = (caps(kk+1)+i*ddcp)*phi(0.0,vso+i*ddv,dvo+i*ddd)
      kappa = kappa+dble(cpo*z)
500   continue
      cpo = caps(kk)*phi(0.0,vs,dv)
      kappa = kappa+dble(0.5*cpo*z)
      endif
c     Values for next step
      yo=y
      vso=vs
      dvo=dv
600   continue
c     *******************************************************
c     Zero z - consider HII region
c     Otherwise the velocities have to be constant
c     *******************************************************
      kk=1
      z=r(2)
      dv=dlb(kk)
      if (dv.ge.1.0) then
      cpo = 2.0*caps(kk)
      kappa = kappa+dble(cpo*z)
      else
      vs= -v(kk)
      kappa=kappa+dble(cpo*z)
      cpo = caps(kk)*phi(0.0,vs,dv)
      kappa = kappa+dble(cpo*z)
      endif
c     *******************************************************
      tau=kappa
      return
      end

c     **************************************************************
c     Subroutine to compute the physical parameters at the end of 
c     the first radiative transfer step within one interval.
c     It checks how many steps have to be done there.
c     This is a reduced version of firststep
c     **************************************************************
      subroutine firstlstp(yo,vo,vso,dvo,cpo,sco,y,v,vs,dr,dv,cp,sc,
     %                       p,z,istp)
      implicit none
      integer istp

      real cpo,sco,cp,sc
      real y,yo,v,vo,vs,vso,dv,dvo,z,p,dr
      real deltav,epsz
      real ddcp,ddsc,ddv,ddt,ddz
      real aa,bb,cadd,cpre,vzgeo,pp
      real fi,rv,rvo,dd,dz
      logical lin

c     Internal variables are transmitted to nextlstp for all
c     later steps
      common /stepl/ ddcp,ddsc,ddv,ddt,ddz,lin
      common /nonlinl/ bb,cadd,cpre,vzgeo,pp
      common /acctrans/ deltav,epsz
      save /stepl/,/nonlinl/,/acctrans/

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
      cpo=cp
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
      ddcp=(cp-cpo)/ddz
      ddsc=(sc-sco)/ddz
c     Constants for the velocity correction
      pp=p
      cadd=(v-bb)*rvo
      cpre=(v-bb)*(rv+rvo)/ddz
c     *****************************
c     Compute the first step
c     *****************************
c     Step size
      vzgeo=bb*rvo+ddv
      ddz=sign(pp*sqrt(vzgeo**2/(bb**2-vzgeo**2)),yo+y)
      dz=ddz-yo
      z=0.5*dz
c     Correct the velocity at that point
      cadd=cadd-cpre*dz
      vso=vzgeo-cadd
c     Assign the linearily interpolated values
      cpo=cpo+ddcp*dz
      sc=sco+ddsc*dz
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
      ddcp=(cp-cpo)/istp
      ddsc=(sc-sco)/istp
c     ****************************
c     Compute the first step
c     ****************************
      z=0.5*ddz
      cpo=cpo+ddcp
      sc=sco+ddsc
      vso=vso+ddv
      dvo=dvo+ddt
c     Values which are not changed in further steps
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
c     This is a reduced version of nextstep
c     *********************************************************
      subroutine nextlstp(vso,dvo,cpo,sco,z)
      implicit none
      real cpo,sco
      real vso,dvo,z
      real y,dz

c     Internal variables computed in firststep 
      real ddcp,ddsc,ddv,ddt,ddz
      real bb,cadd,cpre,vzgeo,pp
      logical lin
      common /stepl/ ddcp,ddsc,ddv,ddt,ddz,lin
      common /nonlinl/ bb,cadd,cpre,vzgeo,pp
      save /stepl/,/nonlinl/

c     *********************************************************
c     If linear approximation is possible
c     *********************************************************
      if (lin) then
c     Update the physical parameters with firststep variables
      z=0.5*ddz
      cpo=cpo+ddcp
      sco=sco+ddsc
      vso=vso+ddv
      dvo=dvo+ddt
c     *********************************************************
c     Nonlinear - compute the possible step size
c     *********************************************************
      else
c     Step size
      vzgeo=vzgeo+ddv
      y=sign(pp*sqrt(vzgeo**2/(bb**2-vzgeo**2)),ddz)
      dz=y-ddz
      ddz=y
      z=0.5*dz
c     Correct the velocity at that point
      cadd=cadd-cpre*dz
      vso=vzgeo-cadd
c     Assign the linearily interpolated values
      cpo=cpo+ddcp*dz
      sco=sco+ddsc*dz
      dvo=dvo+ddt*dz
      endif
      return
      end

c     ***********************************************************
c     Initialization of the kappa field on the left cloud edge
c     This is a reduced version of leftkp for only one line
c     ***********************************************************
      subroutine leftkpl(cp,cappa,j1,j2,vs,dv,dny)
      implicit none
      integer j1,j2

      include 'fsizes.inc'
      real cp, cappa(-nf:nf)
      real vs,dv,dny
      real negexp,negint
      real phi,prof
      integer j

      common /accgauss/ negexp,negint
      save /accgauss/

c     Frequency range
      j1=int((vs-dv*negexp)/dny)-1
      j2=int((vs+dv*negexp)/dny)+1
c     Frequency dependent quantities
      do 500 j=j1,j2
      prof = phi(j*dny,vs,dv)
      cappa(j)  = cp * prof
500   continue
      return
      end

c     ***********************************************************
c     Subroutine for radiative transfer between two points
c     This is a reduced version of transfer for only one line
c     But additional accumulator for tau added
c     ***********************************************************
      subroutine transfl(cp,sc,iplus,tauaccu,cap,ss,j1,j2,vs,dv,dny,z)
      implicit none
      integer j1,j2

      include 'fsizes.inc'
      real cp,sc,ss
      real cap(-nf:nf)
      real iplus(-nf:nf),tauaccu(-nf:nf)
      real vs,dv,dny,z
      real negexp,negint,explin,softzone
      real phi,prof,tau1,tau2,tau,tauq
      real g1,ex1,ex2
      real ihlp1,ihlp2,rel
      integer jn1,jn2,jj1,jj2,j

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
      cap(j) = cp*prof 
      tau = z*cap(j)
      tauaccu(j)=tauaccu(j)+tau
      iplus(j)=(iplus(j)+sc*tau)/(1.0+tau)
510   continue
      else if (j1.lt.jj1) then
      do 520 j=j1,jj1-1
      tau = z*cap(j)
      tauaccu(j)=tauaccu(j)+tau
      iplus(j)=(iplus(j)+ss*tau)/(1.0+tau)
520   continue
      endif
c     ******************************************************
c     Treat the main, fully overlapping part
c     Only in this part, it is controlled whether the lines 
c     are optically thick, i.e. it is assumed that they are
c     thin within the non-overlapping wings
c     ******************************************************
c     When even the centre of the line is optically 
c     thin the time consuming if in the inner loop shall
c     be avoided
c     ******************************************************
      if (abs(2.0*z*cp*phi(0.0,0.0,dv)).lt.explin) then
      do 500 j=jj1,jj2
      prof = phi(j*dny,vs,dv)
      tau1 = z*cap(j)
      cap(j) = cp*prof 
      tau2 = z*cap(j)
      tauq = tau1*tau2
      tau = tau1+tau2
      tauaccu(j) = tauaccu(j)+tau
      iplus(j) = (iplus(j)+sc*tau2+ss*tau1
     %   +0.66667*tauq*(2.0*sc+ss))/(1.0+tau+2.0*tauq)
500   continue
      else
c     ******************************************************
c     In each step there is independently asked whether
c     the line is optically thick
c     It has to be guaranteed that none of the cappa values
c     is 0 in this case. (Should be always fulfilled for
c     the overlapping part.)
c     ******************************************************
      do 1500 j=jj1,jj2
      prof = phi(j*dny,vs,dv)
      tau1 = z*cap(j)
      cap(j) = cp*prof 
      tau2 = z*cap(j)
      tau=tau1+tau2
      rel = abs(tau)
      tauaccu(j) = tauaccu(j)+tau
      if (rel.lt.explin) then
c     Thin case - linear approx.
      tauq = tau1*tau2
      iplus(j) = (iplus(j)+sc*tau2+ss*tau1
     %   +0.66667*tauq*(2.0*sc+ss))/(1.0+tau+2.0*tauq)
      else if (rel.lt.softzone) then
c     Small transition regime - superposition
      tauq = tau1*tau2
      ihlp1 = (iplus(j)+sc*tau2+ss*tau1
     %  +0.66667*tauq*(2.0*sc+ss))/(1.0+tau+2.0*tauq)
      ex1=exp(-0.33333*tau1)
      ex2=exp(-0.33333*tau2)
      g1=0.5*(ss-sc)
      ihlp2=sc+ex2*ex2*(g1+ex2*ex1*(g1+ex1*ex1*(iplus(j)-ss)))
      iplus(j)=(ihlp1*(softzone-rel)+ihlp2*(rel-explin))/
     %              (softzone-explin)
      else
c     Optically thick case - exponential regime
      ex1=exp(-0.33333*tau1)
      ex2=exp(-0.33333*tau2)
      g1=0.5*(ss-sc)
      iplus(j)=sc+ex2*ex2*(g1+ex2*ex1*(g1+ex1*ex1*(iplus(j)-ss)))
      endif
1500  continue
      endif
c     *******************************************************
c     Treat the non-overlapping right end
c     *******************************************************
      if (jn2.gt.jj2) then
      do 560 j=jj2+1,jn2
      prof = phi(j*dny,vs,dv)
      cap(j) = cp*prof
      tau = z*cap(j) 
      tauaccu(j)=tauaccu(j)+tau
      iplus(j)=(iplus(j)+sc*tau)/(1.0+tau)
560   continue
      else if (j2.gt.jj2) then
      do 570 j=jj2+1,j2
      tau = z*cap(j)
      tauaccu(j)=tauaccu(j)+tau
      iplus(j)=(iplus(j)+ss*tau)/(1.0+tau)
570   continue
      endif
c     *******************************************************
c     New parameters
c     *******************************************************
      j1=jn1
      j2=jn2
      ss=sc
      return
      end

c     ***********************************************************
c     Subroutine for the radiative transfer in the continuum
c     This is a reduced version of ctransfer for only one line
c     But additional accumulator for tau added
c     ***********************************************************
      subroutine ctransfl(cp,sc,iplus,tauaccu,nff,z)
      implicit none
      integer nff

      include 'fsizes.inc'
      real cp,sc
      real iplus(-nf:nf),tauaccu(-nf:nf)
      real z,tau,tauq,xx
      real explin,softzone
      integer j

      common /accstep/ explin,softzone
      save /accstep/

c     ******************************************************
c     Treat all frequencies
c     ******************************************************
      tau=cp*z
      if (tau.lt.explin) then
c     ******************************************************
c     Linear approximation for optically thin ranges
c     ******************************************************
      tauq = 0.5*tau**2
      do 650 j=-nff,nff
      tauaccu(j) = tauaccu(j)+tau
      iplus(j) = (iplus(j)+sc*(tau+tauq))/(1.0+tau+tauq)
650   continue
      else
c     ******************************************************
c     Exponential function for thick ranges
c     ******************************************************
      xx=exp(-tau)
      do 680 j=-nff,nff
      tauaccu(j) = tauaccu(j)+tau
      iplus(j)=xx*(iplus(j)-sc) + sc
680   continue
      endif
      return
      end

c     ************************************************************
c     Subroutine to compute the efficient scaling in the 
c     displacement parameter p to have sufficiently dense points.
c     In contrast to makps, the scale will be always denser than
c     the radial point scale. No points are omitted.
c     indr contains the number of crossings with radii for all p
c     ************************************************************
      subroutine psforline(r,vs,dv,nn,p,indr,kk,sig,dmap,doff,nb,kerr)
      implicit none
      integer nn,kk,nb,kerr

      include 'fsizes.inc'
      real r(ns),vs(ns),dv(ns),p(nstretch),sig,doff,dmap
      real vx(ns),r2,relp,relz,dp,deps,velp,vlast,numerrz
      real deltav,epsz,negexp,negint
      real plow,pup
      integer indr(nstretch)
      integer i,nadd,j,imin,imax

      common /acctrans/ deltav,epsz
      common /accgauss/ negexp,negint
      save /acctrans/,/accgauss/

c     Prefactor to avoid very close double points due to numerical
c     round-off errors
      data numerrz/ 1.02/

c     ***************************************************************
c     Restrict the range of p values needed for the chosen map
c     ***************************************************************
      plow=max(0.0,doff-negint*sig)
      if (nb.gt.0) then
      pup=min(r(nn),doff+nb*dmap+negint*sig)
      else
      pup=r(nn)
      nb=max(int((r(nn)-doff)/dmap)+3,1)
      endif
c     Find radial indices according to pup and plow
      do 400 i=1,nn
      if (r(i).gt.plow) then
      imin=i-1
      goto 410
      endif
400   continue
      imin=nn-1
410   do 450 i=imin,nn
      if(r(i).gt.pup) then
      imax=i
      goto 460
      endif
450   continue
      imax=nn
460   continue
c     ***************************************************************
c     Compute the changes in maximum velocities for tangential rays
c     Any present HII region may use only one point.
c     ***************************************************************
      do 300 i=1,nn-1
      vx(i)=0.0
      do 300 j=i+1,nn
      vx(i)=max(vx(i),abs(vs(j)*sqrt(1.0-(r(i)/r(j))**2)))
300   continue
      vx(nn)=0.0
c     ***************************************************************
c     First values
c     ***************************************************************
      kerr=0
      deps=epsz*sig
      kk=1
      p(kk)=r(imin)
      indr(1)=imin+nn
c     ***************************************************************
c     Radial loop for all shells 
c     ***************************************************************
      do 100 i=imin+1,imax
      r2=r(i)**2
      vlast=vx(i-1)
500   if (kk.eq.nstretch) then
       kerr=1
       return
      endif
      relp=r(i)-p(kk)
      dp=deps
c     First criterion spatial density
      relz=sqrt(1.0-p(kk)**2/r2)
      if (relz.gt.numerrz*epsz) then
          dp=min(dp,sqrt(r2*(1.0-(relz-epsz)**2))-p(kk))
      endif
c     Second criterion - frequency/projected velocity density
      velp=vx(i)-vlast
      if (abs(velp).gt.dv(i-1)) then 
        nadd=int(abs(velp)/dv(i-1))+1
        dp=min(dp,relp/nadd)
      endif
      if (dp*numerrz.lt.relp) then
        kk=kk+1
        p(kk)=p(kk-1)+dp
        vlast=vlast+velp*dp/relp
        indr(kk)=i
        goto 500
      else
c     Always include a p line with p=r
c     Test for p field overflow
        if (kk.eq.nstretch) then
         kerr=1
         return
        endif
        kk=kk+1
        p(kk)=r(i)
        indr(kk)=i+nn
      endif
100   continue
      return
      end

