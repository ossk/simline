c     ***********************************************************
c     Program to compute the line profiles for a given cloud
c     In a first step the level populations are determined
c     self-consistently. In a second step the line profiles
c     are computed using the parameters of the telescope.
c     ***********************************************************
      program linetrans

      implicit none
      integer lev,nn

      include 'fsizes.inc'
      real r(ns),dco(ns),vt(ns),vs(ns),rc(ns),tb
      logical aexit,part1,tauwrite
      real*8 dd(ns,nlev)
      integer inum,kerr
      character*250 mol,names
      character yn

c     ************************************************************
c     Formats for the input controlling the general program flow
c     ************************************************************
12    format(//,
     %  ' Do you want to store the results in a file ? [y/n]: ',$)
21    format(/,' Do you want to compute a new model for the cloud ',
     %  ' ? [y/n]: ',$)
23    format(/,' Do you want to compute another line ? [y/n]: ',$)
24    format(/,' Do you want to treat another model ? [y/n]: ',$)
59    format(//,' Stop on all errors requested. Program halted.')
41    format(A1)
c     ************************************************************
c     Loops of computation
c     ************************************************************
      call present(1)
      call predefine(aexit,part1,tauwrite)
      inum=0
110   kerr=0
      write(6,21)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) then
       call micro(inum,mol,tb,r,vs,vt,rc,dco,dd,nn,lev,names,
     %              aexit,kerr)
       if (kerr.ne.0) goto 250
       write(6,12)
       read (5,41) yn
       if ((yn.eq.'y').or.(yn.eq.'Y')) 
     %      call writedens(mol,tb,r,vs,vt,rc,dco,dd,nn,lev,names)
      else
       call readdens(mol,tb,r,vs,vt,rc,dco,dd,nn,lev,names,kerr)
       if (kerr.ne.0) goto 250
      endif
        if (part1) goto 260
200   call line(inum,mol,tb,r,vs,vt,rc,dco,dd,nn,lev,names,
     %        tauwrite,kerr)
      if (kerr.ne.0) goto 250
      write(6,23)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) goto 200
240   write(6,24)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) goto 110
260   call present(2)
      stop
250   if (aexit) then
       write(6,59)
       stop 1
      endif
      goto 240
      end

c     ***********************************************************
c     Subroutine to compute the level populations for one kind
c     of molecules in a dense cloud without internal heating
c     ***********************************************************
c     Statistical turbulence approximation
c     ***********************************************************
      subroutine micro(inum,mol,tb,r,vs,vt,rc,dco,dd,nn,lev,name,
     %                   aexit,kerr)
      implicit none
      integer lev,nn,inum,kerr

      include 'fsizes.inc'
      real r(ns),t(ns),dh(ns),dco(ns),vt(ns),vs(ns),rc(ns),dny(ns)
      real caps(ns,ntra),srcs(ns,ntra),tb
      real*8 dd(ns,nlev)
      character*250 mol,name

      real ulev(ns,ntra),uext(ntra)
      real vmax
      integer istat,all,lin,olin
      logical aexit,hii

c     ************************************************************
c     Formats for the error handling
c     ************************************************************
11    format(/,' Number of radial points insufficient to treat the',
     %  ' problem with the',/,' required accuracy !')
12    format(/,' Number of frequency points insufficient to treat',
     %  ' the problem with the',/,' required accuracy !')
13    format(/,' Number of ray points insufficient to treat the',
     %  ' problem with the',/,' required accuracy !')
14    format(/,' Convergence is not reached !')
15    format(/,' Error reading the molecule file !')
16    format(/,' Numerical singularity occured for the given',
     %  ' parameters !')
17    format(/,' Kinetic temperatures below the cosmic background',
     %  ' temperature are not allowed !')
c     ************************************************************
c     Control output formats
c     ************************************************************
21    format(/,' *** Model iteration is running ***',/)
22    format(' Initial number of radial grid points: ',i3)
23    format(/,' Lambda cycle: ',i3,' Levels: ',i3,' Lines: ',
     %  i3,' Radial grid points: ',i3)
24    format(/,' Total number of iterations used: ',i4)
25    format(' Molecule can be treated with ',I3,' levels and ',
     %          I3, ' lines.')
c     ************************************************************
c     Read the input data of the program
c     An HII region is indicated by vt(0)>1.0
c     ************************************************************
      mol=' '
1100  call readspecies(mol,lin,lev,kerr)
      if (kerr.ne.0) then
        write(6,15)
        return
      endif
      write(6,25) lev,lin
      olin=lin

      call readphys(tb,name,kerr)
      if (kerr.eq.3) then
        write(6,17)
        goto 1000
      else if (kerr.ne.0) then
        return
      endif
      call readnum(inum)
      if (inum.eq.0) then
       kerr=1
       return
      endif

c     ************************************************************
c     Initialize all fields
c     ************************************************************
c     In case of file input for the initial guess 
c     the file determines the initial radial spacings
c     ************************************************************
      if (inum.eq.4) then
       call filedens(r,t,dh,dco,vs,vt,rc,dd,nn,lev,kerr)
       if (kerr.ne.0) goto 1000
      else
      call initrad(r,t,dh,dco,vs,vt,rc,nn,kerr)
      if (kerr.ne.0) then
        write(6,11)
        goto 1000
      endif
      endif
      write(6,22) nn

c     Exceptional treatment in case of central HII region
      hii = ((vt(1).ge.1.0).or.(rc(1).lt.0.0))
      if (hii) call hiiinit(vs(1),vt(1),lin)

      call uback(uext,tb,lin)
      if (inum.ne.4) 
     %     call initn(inum,r,vs,vt,rc,t,dh,dco,uext,dd,nn,lin,lev)
      write(6,21)

c     **********************************************************
c     Cycle for iterations
c     Check the initial values with negative inum
c     **********************************************************
      all=0
      istat=0
      if (inum.lt.0) goto 100
400   continue
c     Change the radial grid to follow the level density changes
      call adaptgrid(r,vs,vt,rc,dco,dh,t,dd,nn,lev,istat,kerr)
       if (kerr.gt.0) then
       write(6,11)
       goto 1000
       endif
       if (kerr.lt.0) then
       write(6,16)
       goto 1000
       endif
      write(6,23) istat,lev,lin,nn

c     Set up the frequency scale for all radial grid points
      call initfreq(vs,vt,nn,vmax,dny,kerr)
       if (kerr.ne.0) then
       write(6,12)
       goto 1000
       endif
      call overlap(vmax,lin)
c     Initialize cappa and S from the populations, effective values
      call kapsrc(dco,dd,caps,srcs,nn,lin)
      call effkap(caps,rc,nn,lin)
c     HII region for most inner point (if present)
      if (hii) call hiikap(caps,srcs,lin)
c     recompute external field if transition numbering changed
      if (lin.ne.olin) then
        call uback(uext,tb,lin)
        olin=lin
      endif
c     Calculate radiative transfer
      call stpmic(r,vs,vt,dny,vmax,caps,srcs,uext,ulev,nn,lin,kerr)
       if (kerr.ne.0) then
       write(6,13)
       goto 1000
       endif
c     Calculate level population of molecule
      call levmic(t,dh,ulev,dd,nn,lin,lev,istat)
      if (hii) call hiismth(dd,lev)
      all=all+1
      if ((istat.ge.0).and.(all.lt.500)) goto 400
c     **********************************************************
c     Write out the results
c     **********************************************************
100   write(6,24) all
      if (all.ge.500) write(6,14)
      return
1000  if (.not.aexit) goto 1100
      return
      end

C     ***********************************************************
C     Line profile phi(ny) * ny0
c     ! in the observers frame !
c     ***********************************************************
c     For Falgarone intermittency profiles add:
c     Falgarone parameters f1=1+a, f2=a/3.3, f3=1/3.3^2
c     data f1,f2,f3 /1.3, 0.090909, 0.091827/
c     phi=pis/(sig*f1)*(exp(-nys)+f2*exp(-nys*f3))
c     ***********************************************************
      real function phi(ny,vs,sig)
      implicit none
      real ny,nys,vs,sig
      real pis

c     1/sqrt(pi)
      data pis /0.5641896/

c     The cut-off for large negative exponents has to be treated 
c     in the calling routine
      nys=((ny-vs)/sig)**2
      phi=pis/sig*exp(-nys)
      return
      end

c     ***********************************************************
c     Subroutine to compute the observed line profiles
c     Microturbulent approximation 
c     ***********************************************************
      subroutine line(inum,mol,tb,r,vs,vt,rc,dco,dd,nn,lev,name,
     %    tauwrite,kerr)
      implicit none
      integer lev,nn,inum

      include 'fsizes.inc'
      real r(ns),dco(ns),vt(ns),vs(ns),rc(ns)
      real caps(ns,ntra),srcs(ns,ntra),tb
      real cap(ns),sc(ns)
      real*8 dd(ns,nlev)
      real p(nstretch),pn(nstretch)
      real iout(nstretch,-nf:nf),tautot(nstretch,-nf:nf)
      real iconv(-nf:nf)
      integer indr(nstretch)
      character*250 mol,name
      character*22 mesg
      logical tauwrite

      real uext(ntra),dhlp(ns)
      real vmax,deltap,dp,dpp,factor,tau0,nu0
      real sigbeam,dny,ubg,tbeam,cfreq,dv,fwhm,dmap,doff,dist
      real spi,pnorm
      real ibeam(ns,-nf:nf),taubeam(ns,-nf:nf)
      integer lin,nff,kerr,ilin,hlev,nbeam,nnp
      integer i,j,il1,il2

c     Normalization factor spi=1/sqrt(pi)
      data spi/ 0.5641896/
c     Translation from I to T: tbeam = h*c^2/(8*pi*k) [K*cm^2]
      data tbeam / 1.716215e9/
c     Empty message at the beginning
      data mesg /'                      '/
      
c     ************************************************************
c     Formats for the error handling
c     ************************************************************
12    format(/,' Number of frequency points insufficient to treat',
     %  ' the problem with the',/,' required accuracy !')
13    format(/,' Number of ray points insufficient to treat the',
     %  ' problem with the',/,' required accuracy !')
15    format(/,' Error reading the molecule file !')
c     ************************************************************
c     Control output formats
c     ************************************************************
21    format(/,' *** Radiative transfer is running ***',/)
22    format(' Number of beam offsets used: ',i4)
23    format(' Number of rays in radiative transfer: ',i4)
24    format(' Optical depth at cloud and line centre: ',1PE11.4)
33    format(/,' The transition: ',A22,' will be computed.')
34    format(' Computing transition: ',A22, ' at ', 1PE11.4, ' Hz')
c     ************************************************************
c     Initialize all fields
c     ************************************************************ 
      call readspecies(mol,lin,hlev,kerr)
        if (kerr.ne.0) then
        write(6,15)
        return
        endif
      if (hlev.ne.lev) call linfromlev(lin,lev)
c     ************************************************************
c     Set up fixed parameters
c     ************************************************************
      call readobs(fwhm,dmap,doff,nbeam,dist,dv,ilin,lin)
      if (ilin.lt.0) return
      if (ilin.gt.0) then
       call translabel(ilin,mesg)
       write(6,33) mesg
      endif
      call translat(fwhm,dmap,doff,dist,dv,sigbeam,dp,dpp,dny)
      call readnum(inum)
      if (inum.eq.0) return
c     **********************************************************
c     Compute physical quantities needed for radiative transfer
c     **********************************************************
c     Set up the frequency scale for all radial grid points
      call initfreq(vs,vt,nn,vmax,dhlp,kerr)
      nff=int(vmax/dny)+1
       if (nff.gt.nf) then
       kerr=1
       write(6,12)
       return
       endif
c     Initialize cappa and S from the populations, effective values
      call kapsrc(dco,dd,caps,srcs,nn,lin)
      call effkap(caps,rc,nn,lin)
      if (vt(1).ge.1.0) then
       call hiiinit(vs(1),vt(1),lin)
       call hiikap(caps,srcs,lin)
      endif
c     Initialize a field of displacement parameters
      p(1)=dpp 
      call psforline(r,vs,dhlp,nn,p,indr,nnp,sigbeam,dp,dpp,nbeam,
     %                                                          kerr)
      write(6,22) nbeam
      write(6,23) nnp
      if (nbeam.gt.ns) then
       kerr=1
       write(6,13)
       return
      endif
      if (kerr.ne.0) then
       write(6,13)
       return
      endif
c     ************************************************************
c     Loop in case of all transitions
c     ************************************************************
      if (ilin.eq.0) then
       il1=1
       il2=lin
      else
       il1=ilin
       il2=ilin
      endif
c     ***********************************************************
      do ilin=il1,il2
      write(6,21)
      nu0=cfreq(ilin)
      call translabel(ilin,mesg)
      write(6,34) mesg, nu0
c     Background radiation
      call uback(uext,tb,ilin)
      ubg=uext(ilin)
c     Reduce the caps and srcs fields to the transition considered
      do i=1,nn
        cap(i)=caps(i,ilin)
        sc(i)=srcs(i,ilin)
      enddo
c     Calculate real radiative transfer
      call stpline(r,vs,vt,cap,sc,nn,ubg,p,indr,iout,tautot,nnp,nff,dny)
c     **********************************************************
c     Beam integration
c     **********************************************************
c     Renormalize displacements to have the beam area=1
      pnorm=sigbeam/spi
      do i=1,nnp
        pn(i)=p(i)/pnorm
      enddo
      factor=tbeam/nu0**2
c     Loop for the integration over the beams at different 
c     displacements, reduction, translation into beam
c     temperature
      do i=1,nbeam
       deltap=(dpp+(i-1)*dp)/pnorm
       call intbeam(pn,iout,ubg,iconv,nnp,nff,deltap,spi)
       do j=-nff,nff
        ibeam(i,j)=iconv(j)*factor
       enddo
      enddo
c     **********************************************************
c     Double effort in case of writing tau files but allows
c     for complete code reuse
c     **********************************************************
      if (tauwrite) then
        do i=1,nbeam
         deltap=(dpp+(i-1)*dp)/pnorm
         call intbeam(pn,tautot,0.0,iconv,nnp,nff,deltap,spi)
         do j=-nff,nff
          taubeam(i,j)=iconv(j)
         enddo
        enddo
      endif
c     **********************************************************
c     For a better cloud control the optical thickness at
c     cloud centre and line centre was computed in a seperate run
c     This is obsolete with the implememntation of the tau profiles
c     taucentral is to be removed in future versions
c     **********************************************************
c     call taucentral(r,vs,vt,cap,nn,tau0)
c     write(6,24) tau0
      tau0=tautot(1,0)
      write(6,24) tau0
c     **********************************************************
c     Write out the results
c     **********************************************************
      call writeline(ibeam,nbeam,nff,dmap,doff,dv,tau0,nu0,.false.,
     %   mesg,name)
      if (tauwrite) then
        call writeline(taubeam,nbeam,nff,dmap,doff,dv,tau0,nu0,.true.,
     %   mesg,name)
      endif
c     **********************************************************
c     End of possible loop
c     **********************************************************
      enddo
      return
      end

c     ***************************************************************
c     Short subroutine translating from observable to physical
c     quantities
c     ***************************************************************
      subroutine translat(fwhm,dmap,doff,dist,dv,sigbeam,dsig,dpp,dny)
      implicit none
      real fwhm,dmap,doff,dist,dv,sigbeam,dsig,dpp,dny

      real clicht,astor,ftosg

c     Speed of light in km/s (required to obtain relative steps)
      data clicht/ 2.99792e5 /
c     astor translates arcsec to rad 
c     astor = pi/(180*60*60)
c     ftosg translates FWHM to sigma
c     ftosg = 2*sqrt(ln2)
      data astor,ftosg / 4.848137e-6,1.66511/

      dny=dv/clicht
      sigbeam=astor*dist*fwhm/ftosg
      dsig=astor*dist*dmap
      dpp=astor*dist*doff
      return
      end

c     ***********************************************************
c     Message screens (Indroduction, credits)
c     ***********************************************************
      subroutine present(i)
      implicit none
      integer i

c     ************************************************************
c     Formats
c     ************************************************************
11    format(/////////////////////////////)
12    format(/,' ------------------------------------------------',
     %  '--------------------')
21    format(//,24X,'SimLine - Version 2.18')
22    format(/,5X,' Radiative transfer in molecular lines through ',
     %  'turbulent media')
23    format(/,26X,' Volker Ossenkopf-Okada')
24    format(20X,' University Observatory Jena /')
32    format(14X,' Cologne Observatory for Sub-mm Astronomy')
25    format(//,27X,'March 22, 2026')
27    format(/,' Credits:',/,' --------')
28    format(' A precursor of this program was written ',
     %  'by E. Kruegel.')
29    format(' The accelerated lambda method was published by ',
     %  'L. Auer.')
30    format(' The treatment of the turbulence correlation length',
     %  ' roughly',/,' follows a method from H.M. Martin, D.B. Sanders',
     %  ' & R.E. Hills.')
31    format(' The LU decomposition, spline interpolation, and ',
     %  'Gaussian',/,' quadrature were ',
     %  'taken from "Numerical Recipies".')
c     ************************************************************
c     Output
c     ************************************************************
      if (i.eq.1) then
      write(6,11)
      write(6,12)
      write(6,21)
      write(6,22)
      write(6,23)
      write(6,24)
      write(6,32)
      write(6,25)
      write(6,12)
      else
      write(6,12)
      write(6,27)
      write(6,28)
      write(6,29)
      write(6,30)
      write(6,31)
      write(6,12)
      endif
      return
      end
