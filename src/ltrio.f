c     ***************************************************************
c     Input of physical parameters - either from an input file
c     or directly from the terminal
c     ***************************************************************
      subroutine readphys(tb,lfile,kerr)
      implicit none
      integer ishell,kerr

      include 'fsizes.inc'
      real rs(0:ns),parm(ns,6,2),ed,et,tb,rat1,tl1,tl2
      integer i,j,l,fixed(0:ns),ipshells
      integer precomm,prestruc,ioerror
      logical oldformat,overwr,fileinp,numdef,anames,allstop
      character*250 lfile,dummyline
      character*50 quant(6)
      character*90 qlin(6)
      character*40 sprm(2)
      character yn

      common /model/ rs,fixed,parm,ishell
      common /options/ precomm,prestruc,overwr,fileinp,numdef,
     %                   anames,allstop
      save /model/,/options/

11    format(//,' *** Input of physical parameters ***',/)
12    format(' Do you want to use an input file ? [y/n]: ',$)
13    format(' Name of the parameter file: ',$)
c     Formats for interactive parameter input
14    format(' File containing the molecular data: ',$)
15    format(' Number of shells for the cloud representation: ',$)
19    format(' Use shells as power law regions or as data table ?',
     %          ' [0/1]: ',$)
16    format(' Radius of the inner core [pc]: ',$)
18    format(' Background radiation temperature [K]: ',$)
21    format(' Use a central HII region ? [y/n]: ',$)
22    format(' Electron density in the HII region [e/cm^3]: ',$)
23    format(' Electron temperature in the HII region [K]: ',$)
c     Repeated parameters within the shells
17    format(/,' Input of power laws for the density, temperature,',
     %  ' and velocity distributions')
24    format(/,' Input of the data table')
25    format(' Radial grid point [pc]: ',$)
26    format(' r [pc], n_H2 [cm^-3], T [K], X [mol/H2], v_turb',
     %  ' [km/s], l_corr [pc] , v [km/s]')
29    format(/,I2,'. shell')
28    format(' Outer radius [pc]: ',$)
c     Input formats
41    format(A1)
42    format(A250)
43    format(' ')
44    format(I4)
45    format(1PG12.5)
46    format(1P7G13.5)
47    format(': ',$)
c     Special control output formats
80    format(I2,'. shell - Outer radius [pc]: ')
81    format(' Core electron density [e/cm^3] and temperature [K]:')

c     Initialize strings for variable format statements
      quant(1)='(/,'' Hydrogen density [H2/cm^3]'',$)'
      quant(2)='(/,'' Temperature [K]'',$)'
      quant(3)='(/,'' Relative molecular abundance [X/H2]'',$)'
      quant(4)='(/,'' FWHM of turbulent velocity [km/s]'',$)'
      quant(5)='(/,'' Turbulence correlation length [pc]'',$)'
      quant(6)='(/,'' Systematic radial velocity [km/s]'',$)'
      sprm(1)='('' Inner value: '',$)'
      sprm(2)='('' Radius exponent: '',$)'
      qlin(1)='(I2,''. shell - H2 density [H2/cm^3] - inner value,'
     %  //' exponent'')'
      qlin(2)='(I2,''. shell - Temperature [K] - inner value,'
     %  //' exponent'')'
      qlin(3)='(I2,''. shell - Molecule abundance [X/H2] - inner '
     %  //'value, exponent'')'
      qlin(4)='(I2,''. shell - Turbulent velocity [km/s] - inner '
     %  //'value, exponent'')'
      qlin(5)='(I2,''. shell - Correlation length [pc] - inner value,'
     %  //' exponent'')'
      qlin(6)='(I2,''. shell - Radial velocity [km/s] - inner value,'
     %  //' exponent'')'
      
c     First selection
      write(6,11)
      kerr=0
100   if (.not.fileinp) then
      write(6,12)
      read (5,41) yn
      if ((yn.ne.'y').and.(yn.ne.'Y')) goto 310
      endif
c     *************************************************************
c     Use an input file
c     *************************************************************
      oldformat=.false.
      write(6,13)
      read (5,42) lfile
      open(14,err=120,file=lfile,status='old')
c     Global parameters
c     read (14,42,err=110,end=110) dummyline
c     read (14,42,err=110,end=110) mol
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) tb
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) ishell
      if (ishell.gt.ns-1) ishell=ns-1
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) ipshells
      if (ipshells.eq.0) oldformat=.true.
      if (oldformat) then
        read (14,42,err=110,end=110) dummyline
        read (14,*,err=110,end=110) rs(1)
      endif
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) ed,et
c     *************************************************************
c     Treat the various shells and the six quantities within them
c     *************************************************************
      if (oldformat) then
      do i=1,ishell
        read (14,42,err=110,end=110) dummyline
        read (14,*,err=110,end=110) rs(i+1)
        do j=1,6
          read (14,42,err=110,end=110) dummyline
          read (14,*,err=110,end=110) parm(i+1,j,1),parm(i+1,j,2)
        enddo
      enddo
c     *************************************************************
c     New format data table
c     Check for wrong order
c     *************************************************************
      else
      read (14,42,err=110,end=110) dummyline
      do i=1,ishell
        read (14,*,err=110,end=110) rs(i),(parm(i,j,1),j=1,6)
      enddo
      if (rs(1).gt.rs(ishell)) then
        do i=1,ishell/2
          l=ishell+1-i
          rat1=rs(i)
          rs(i)=rs(l)
          rs(l)=rat1
          do j=1,6
            rat1=parm(i,j,1)
            parm(i,j,1)=parm(l,j,1)
            parm(l,j,1)=rat1
          enddo
        enddo
      endif
c     End of special treatment for new table format 
      endif
      close(14)
      goto 500
c     *************************************************************
c     Use interactive input
c     *************************************************************
310   oldformat=.false.
      write(6,43)
c     write(6,14)
c     read (5,42,err=180) mol
      write(6,18)
      read (5,*,err=180) tb
      write(6,15)
      read (5,*,err=180) ishell
      if (ishell.gt.ns-1) ishell=ns-1
      write(6,19)
      read (5,*,err=180) ipshells
      if (ipshells.eq.0) oldformat=.true.
      if (oldformat) then
        write(6,16)
        read (5,*,err=180) rs(1)
      endif
c     Question for HII region
      write(6,21)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) then
        write(6,22)
        read (5,*,err=180) ed
        write(6,23)
        read (5,*,err=180) et
      else
        ed=0.0
        et=0.0
      endif
c     *************************************************************
c     Treat the various shells and the six quantities within them
c     *************************************************************
      if (oldformat) then
      write(6,17)
      do i=1,ishell
        write(6,29) i
c     Shell radii
        write(6,28)
        read (5,*,err=180) rs(i+1)
        do j=1,6
          write(6,quant(j))
          write(6,43)
          do l=1,2
            write(6,sprm(l))
            if ((i.gt.1).and.(l.eq.1)) then
c     Continue smoothly
              parm(i+1,j,1)=parm(i,j,1)*(rs(i)/rs(i-1))**parm(i,j,2)
              write(6,45) parm(i+1,j,1)
            else
c     Read other parameters
              read (5,*,err=180) parm(i+1,j,l)
            endif
          enddo
        enddo
      enddo
c     *************************************************************
c     New format data table
c     *************************************************************
      else
      write(6,24)
      do i=1,ishell
        write(6,29) i
        write(6,25)
        read (5,*,err=180) rs(i)
        do  j=1,6
          write(6,quant(j))
          write(6,47)
          read (5,*,err=180) parm(i,j,1)
        enddo
      enddo
c     *************************************************************
c     Correct for wrong order
c     *************************************************************
      if ((rs(1).gt.rs(ishell)).and.(rs(2).ge.rs(ishell))) then
      do i=1,ishell/2
        l=ishell+1-i
        rat1=rs(i)
        rs(i)=rs(l)
        rs(l)=rat1
        do j=1,6
          rat1=parm(i,j,1)
          parm(i,j,1)=parm(l,j,1)
          parm(l,j,1)=rat1
        enddo
      enddo
      endif
      endif
c     *************************************************************
c     In case of interactive input store the parameters for
c     the next time in a normal input file
c     *************************************************************
72    format(//,' Do you want to store these parameters for later ',
     %  'use ? [y/n]: ',$)
      write(6,72)
      read (5,41) yn
      if ((yn.ne.'y').and.(yn.ne.'Y')) goto 500
350   write(6,13)
      read (5,42) lfile
      open(14,err=130,file=lfile,status='new')
250   continue
c     write(14,14)
c     write(14,43)
c     write(14,42) mol
      write(14,18)
      write(14,43)
      write(14,45) tb
      write(14,15)
      write(14,43)
      write(14,44) ishell
      write(14,19)
      write(14,43)
      write(14,44) ipshells
      if (oldformat) then
        write(14,16)
        write(14,43)
        write(14,46) rs(1)
      endif
      write(14,81)
      write(14,46) ed,et
c     Write out shells in old format
      if (oldformat) then
      do i=1,ishell
        write(14,80) i
        write(14,46) rs(i+1)
        do j=1,6
          write(14,qlin(j)) i
          write(14,46) parm(i+1,j,1),parm(i+1,j,2)
        enddo
      enddo
c     Shell table in new format
      else
      write(14,26)
      do i=1,ishell
        write(14,46) rs(i), (parm(i,j,1),j=1,6)
      enddo
      endif
      close(14)
c     *************************************************************
c     Final treatment:
c     Fill the core with constant parameters to allow a simplified 
c     treatment, create exponents from the data table
c     *************************************************************
500   rs(0)=0.0
      if (oldformat) then
      ishell=ishell+1
      do i=0,ishell
       fixed(i)=1
      enddo
c     Check for too low temperatures
      do i=2,ishell
       tl1=parm(i,2,1)
       tl2=tl1*(rs(i)/rs(i-1))**parm(i,2,2)
       if ((tl1.lt.tb).or.(tl2.lt.tb)) then
         kerr=3
         return
       endif
      enddo
      else
c     New format - compute exponents
      fixed(0)=1
      fixed(1)=1
      fixed(ishell)=1
c     Do not allow for zero radii as input
      if (rs(1).eq.0.0) rs(1)=0.1e-3*rs(2)
c     Loop over all shells for exponents
      do i=ishell,2,-1
c       Check for too low temperatures
        if (parm(i,2,1).lt.tb) then
          kerr=3
          return
        endif
c       Compute exponents, reassign
c       Switch to definitions used in old format with shift of
c       radius relative to parameter (inner parameter vs. outer radius)
        rat1=alog(rs(i)/rs(i-1))
        fixed(i)=0
        do j=1,6
          if (1.0d0*parm(i,j,1)*parm(i-1,j,1).le.0.0d0) then
            if (parm(i,j,1)-parm(i-1,j,1).ne.0.0) fixed(i)=1
            parm(i,j,2)=0.0
            parm(i,j,1)=parm(i-1,j,1)
          else
            parm(i,j,2)=alog(parm(i,j,1)/parm(i-1,j,1))/rat1
            parm(i,j,1)=parm(i-1,j,1)
          endif
        enddo
      enddo
      endif
c     ************************************************************
c     Core
c     ************************************************************
      do j=1,6
        parm(1,j,1)=parm(2,j,1)
        parm(1,j,2)=0.0
      enddo
c     Assign electron density and temperature to systematic and
c     turbulent velocities in case of an HII region.
      if ((ed.gt.0.0).and.(et.gt.0.0)) then
        parm(1,6,1)=ed
        parm(1,4,1)=et
        parm(1,5,1)=-1.0
        if (et.lt.tb) kerr=3
      endif
      return

c     *************************************************************
c     Error handling
c     *************************************************************
120   goto (100) ioerror(lfile,1)
      kerr=1
      return
110   goto (100) ioerror(lfile,2)
      kerr=1
      return
130   goto (250,350) ioerror(lfile,3)+1
      goto 500
58    format(/,' A numerical error occured - ',
     %  'Please correct your data input!')
180   write(6,58)
      goto 310
      end

c     ***************************************************************
c     Input of numerical parameters - either from an input file
c     or directly from the terminal
c     ***************************************************************
      subroutine readnum(inum)
      implicit none
      integer inum
      real epsn,epsv,dvgauss,hyster,minspac,minedge
      real negexp,negint,explin,softzone,deltav,epsz
      real*8 conv,neglev

      integer precomm,prestruc,ioerror
      logical overwr,fileinp,numdef,anames,allstop
      character*250 lfile,dummyline
      character yn

      common /accgrid/ epsn,epsv,dvgauss,hyster,minspac,minedge
      common /acctrans/ deltav,epsz
      common /acclev/ conv,neglev
      common /accgauss/ negexp,negint
      common /accstep/ explin,softzone
      common /options/ precomm,prestruc,overwr,fileinp,numdef,
     %                   anames,allstop
      save /accgrid/,/acctrans/,/acclev/,/accgauss/,/accstep/,/options/

11    format(//,' *** Input of numerical parameters ***',/)
12    format(' Do you want to use built-in default values ? [y/n]: ',$)
13    format(' Do you want to use an input file ? [y/n]: ',$)
14    format(' Name of the parameter file: ',$)
15    format(' Do you want to use the same numerical parameters ',
     %  'as in the previous',/,' computation ? [y/n]: ',$)
22    format(' Cut-off for Gaussians in the line emission [sigma]: ',$)
23    format(' Cut-off for Gaussians in freq. integration [sigma]: ',$)
24    format(' Resolution in scanning the Gaussians [sigma]: ',$)
31    format(' Maximum change factor for the level densities between',
     %  ' grid points: ',$)
32    format(' Hysteresis in the radial grid adjustment [0..1]: '
     %  ,$)
33    format(' Maximum change of the radial velocity between grid',
     %  ' points [sigma]: ',$)
34    format(' Maximum relative distance of points on the',
     %  ' spatial integration scale: ',$)
35    format(' Maximum velocity shift between radiative transfer',
     %  ' points [sigma]: ',$)
36    format(' Initial populations [1=thermal,2=no radiation,',
     %  '3=Sobolev,4=file input]: ',$)
37    format(' Neglection threshold for the level populations: ',$)
38    format(' Relative accuracy as convergence criterion: ',$)
c     Input formats
41    format(A1)
42    format(A250)
43    format(' ')
44    format(I4)
45    format(1PG12.5)

c     *************************************************************
c     Skip terminal interaction and take defaults for -numdefaults
c     *************************************************************
      if (numdef) goto 200
      write(6,11)
c     *************************************************************
c     Use the old values if called multiple times
c     *************************************************************
      if (inum.ne.0) then
      write(6,15)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) return
      endif
100   write(6,12)
      read (5,41) yn
      if ((yn.ne.'y').and.(yn.ne.'Y')) goto 410
c     *************************************************************
c     The default values
c     *************************************************************
200   negexp=3.2
      negint=4.2
      dvgauss=0.55
      epsn=1.3
      hyster=0.2
      epsv=2.5
      epsz=0.2
      deltav=dvgauss
      inum=2
      neglev=1d-8
      conv=2d-6
      goto 500
c     *************************************************************
c     Use an input file
c     *************************************************************
410   if (.not.fileinp) then
      write(6,13)
      read (5,41) yn
      if ((yn.ne.'y').and.(yn.ne.'Y')) goto 310
      endif
      write(6,14)
      read (5,42) lfile
      open(14,err=120,file=lfile,status='old')
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) negexp
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) negint
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) dvgauss
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) epsn
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) hyster
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) epsv
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) epsz
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) deltav
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) inum
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) neglev
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) conv
      goto 500
c     *************************************************************
c     Use interactive input
c     *************************************************************
310   write(6,43)
      write(6,22)
      read (5,*,err=180) negexp
      write(6,23)
      read (5,*,err=180) negint
      write(6,24)
      read (5,*,err=180) dvgauss
      write(6,43)
      write(6,31)
      read (5,*,err=180) epsn
      write(6,32)
      read (5,*,err=180) hyster
      write(6,33)
      read (5,*,err=180) epsv
      write(6,34)
      read (5,*,err=180) epsz
      write(6,35)
      read (5,*,err=180) deltav
      write(6,43)
      write(6,36)
      read (5,*,err=180) inum
      write(6,37)
      read (5,*,err=180) neglev
      write(6,38)
      read (5,*,err=180) conv
c     *************************************************************
c     In case of interactive input store the parameters for
c     the next time in a normal input file
c     *************************************************************
72    format(/,' Do you want to store these parameters for later ',
     %  'use ? [y/n]: ',$)
      write(6,72)
      read (5,41) yn
      if ((yn.ne.'y').and.(yn.ne.'Y')) goto 500
350   write(6,14)
      read (5,42) lfile
      open(14,err=130,file=lfile,status='new')
250   write(14,22)
      write(14,43)
      write(14,45) negexp
      write(14,23)
      write(14,43)
      write(14,45) negint
      write(14,24)
      write(14,43)
      write(14,45) dvgauss
      write(14,31)
      write(14,43)
      write(14,45) epsn
      write(14,32)
      write(14,43)
      write(14,45) hyster
      write(14,33)
      write(14,43)
      write(14,45) epsv
      write(14,34)
      write(14,43)
      write(14,45) epsz
      write(14,35)
      write(14,43)
      write(14,45) deltav
      write(14,36)
      write(14,43)
      write(14,44) inum
      write(14,37)
      write(14,43)
      write(14,45) neglev
      write(14,38)
      write(14,43)
      write(14,45) conv
      close(14)
c     *************************************************************
c     Both cases:
c     Check the parameters for being in a reasonable range
c     Use default values if not.
c     *************************************************************
500   if (negexp.le.0.0) negexp=3.2
      if (negint.le.0.0) negint=4.2
      if (dvgauss.le.0.0) dvgauss=0.55
      if (epsn.le.1.0) epsn=1.3
      if ((hyster.lt.0.0).or.(hyster.ge.1.0)) hyster=0.2
      if (epsv.le.0.0) epsv=2.5
      if (epsz.le.0.0) epsz=0.2
      if (deltav.le.0.0) deltav=dvgauss
      if (inum.eq.0) inum=2
c     *************************************************************
c     Threshold for the transition from the linearily interpolated
c     thin approximation to the stepwise constant thick treatment
c     in radiative transfer. This is too dangerous as independent
c     parameter. An appropriate treatment is coupling to the
c     maximum density step
c     *************************************************************
      if (epsn.lt.1.03) then
      softzone=0.06
      else if (epsn.lt.1.15) then
      softzone=0.07
      else if (epsn.lt.1.4) then
      softzone=0.15
      else
      softzone=0.21
      endif
      explin=0.995*softzone 
c     The hysteresis factor in the program is 1 for zero hysteresis
      hyster=1.0-hyster
      return
c     *************************************************************
c     Error handling
c     *************************************************************
120   goto (100) ioerror(lfile,1)
      inum=0
      return
110   goto (100) ioerror(lfile,2)
      inum=0
      return
130   goto (250,350) ioerror(lfile,3)+1
      goto 500
58    format(/,' A numerical error occured - ',
     %  'Please correct your data input!')
180   write(6,58)
      goto 310
      end

c     ***************************************************************
c     Input of observational antenna parameters - either from an
c     input file or directly from the terminal
c     ***************************************************************
      subroutine readobs(fwhm,dmap,doff,nmap,dist,dny,ilin,lin)
      implicit none
      integer ilin,lin,nmap,iuplev,ilolev
      real fwhm,dist,dmap,doff,dny
      character*250 lfile,dummyline
      character yn
      integer precomm,prestruc,ioerror,theline
      logical overwr,fileinp,numdef,anames,allstop

      common /options/ precomm,prestruc,overwr,fileinp,numdef,
     %                   anames,allstop
      save /options/

11    format(//,' *** Input of observational parameters ***',/)
12    format(' Do you want to use an input file ? [y/n]: ',$)
13    format(' Name of the parameter file: ',$)
c     Formats for interactive parameter input
14    format(' Upper level of the transition to be observed: ',$)
21    format(' Lower level of the transition to be observed: ',$)
15    format(' Frequency resolution [km/s]: ',$)
16    format(' Beam width [FWHP in ``]: ',$)
17    format(' Step size for mapping [``]: ',$)
18    format(' Central offset for the first point [``]: ',$)
19    format(' Number of map points [0=map whole cloud]: ',$)
20    format(' Distance of the cloud [pc]: ',$)
c     Input formats
41    format(A1)
42    format(A250)
43    format(' ')
44    format(I4)
45    format(1PG12.5)
      
c     First selection
      write(6,11)
100   if (.not.fileinp) then
      write(6,12)
      read (5,41) yn
      if ((yn.ne.'y').and.(yn.ne.'Y')) goto 310
      endif
c     *************************************************************
c     Use an input file
c     *************************************************************
      write(6,13)
      read (5,42) lfile
      open(14,err=120,file=lfile,status='old')
c     Global parameters
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) iuplev
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) ilolev
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) dny
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) fwhm
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) dmap
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) doff
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) nmap
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) dist
      close(14)
c     Check transition
      ilin=theline(iuplev,ilolev,lin)
      if (ilin.lt.0) goto 100
      goto 500
c     *************************************************************
c     Use interactive input
c     *************************************************************
310   write(6,43)
      write(6,14)
      read (5,*,err=180) iuplev
      write(6,21)
      read (5,*,err=180) ilolev
c     Check transition immediately in case of interactive input
      ilin=theline(iuplev,ilolev,lin)
      if (ilin.lt.0) goto 310
      write(6,15)
      read (5,*,err=180) dny
      write(6,16)
      read (5,*,err=180) fwhm
      write(6,17) 
      read (5,*,err=180) dmap
      write(6,18) 
      read (5,*,err=180) doff
      write(6,19) 
      read (5,*,err=180) nmap
      write(6,20) 
      read (5,*,err=180) dist
c     *************************************************************
c     In case of interactive input store the parameters for
c     the next time in a normal input file
c     *************************************************************
72    format(//,' Do you want to store these parameters for later ',
     %  'use ? [y/n]: ',$)
      write(6,72)
      read (5,41) yn
      if ((yn.ne.'y').and.(yn.ne.'Y')) goto 500
350   write(6,13)
      read (5,42) lfile
      open(14,err=130,file=lfile,status='new')
250   write(14,14)
      write(14,43)
      write(14,44) iuplev
      write(14,21)
      write(14,43)
      write(14,44) ilolev
      write(14,15)
      write(14,43)
      write(14,45) dny
      write(14,16)
      write(14,43)
      write(14,45) fwhm
      write(14,17)
      write(14,43)
      write(14,45) dmap
      write(14,18)
      write(14,43)
      write(14,45) doff
      write(14,19)
      write(14,43)
      write(14,44) nmap
      write(14,20)
      write(14,43)
      write(14,45) dist
      close(14)
c     *************************************************************
c     Both cases:
c     check for errors
c     *************************************************************
32    format(/,' Only finite beam widths are possible !',/)
500   continue
      if (fwhm.le.0.0) then
       write(6,32)
       goto 100
      endif
      return
c     *************************************************************
c     Error handling
c     *************************************************************
120   goto (100) ioerror(lfile,1)
      ilin=-1
      return
110   goto (100) ioerror(lfile,2)
      ilin=-1
      return
130   goto (250,350) ioerror(lfile,3)+1
      goto 500
58    format(/,' A numerical error occured - ',
     %  'Please correct your data input!')
180   write(6,58)
      goto 310
      end

c     ***********************************************************
c     Function to find the corresponding line for levels
c     ***********************************************************
      integer function theline(lup,ldo,lin)
      implicit none
      integer lup,ldo, lin

      include 'fsizes.inc'
      integer lfrom(ntra),lto(ntra),i

      common /lintable/ lfrom, lto
      save /lintable/

31    format(/,' There is no excited transition from level ',I3,
     %   ' to level ',I3,'!',/)
32    format(/,' All excited transitions will be computed.')

c     ***********************************************************
c     default, zero-zero = all lines
c     ***********************************************************
      theline=-1
      if ((lup.eq.0).and.(ldo.eq.0)) then
       theline=0
       write(6,32)
       return
      endif
c     ***********************************************************
c     Check all lines
c     ***********************************************************
      do i=1,lin
        if ((lup.eq.lfrom(i)).and.(ldo.eq.lto(i))) theline=i
      enddo
        if (theline.eq.-1) write(6,31) lup,ldo

      return
      end

c     ***********************************************************
c     Subroutine to get the transition label 
c     ***********************************************************
      subroutine translabel(ilin,mesg)
      implicit none
      integer ilin
      character*22 mesg

      include 'fsizes.inc'
      real molmass,energy(nlev),stgew(nlev)
      integer lfrom(ntra),lto(ntra),endu,endo
      character*10 levnam(nlev),levnu,levno

      common /molecule/ molmass, energy, stgew, levnam
      common /lintable/ lfrom, lto
      save /molecule/, /lintable/

c     ***********************************************************
c     Messages on the result
c     ***********************************************************
      if (ilin.eq.0) then
        mesg=' '
      else
        levnu=levnam(lfrom(ilin))
        call clearstring(levnu,endu,.true.)
        levno=levnam(lto(ilin))
        call clearstring(levno,endo,.true.)
        mesg=levnu(1:endu)//'--'//levno(1:endo)
      endif
      return
      end

c     ***********************************************************
c     Subroutine to write out the results in a form that can be
c     used to compute the observed profiles of certain lines
c     The energy density ulev is only written in a preliminary 
c     state for test purposes
c     Here, the maximum level number is hard wired by a format !
c     ***********************************************************
      subroutine writedens(mol,tb,r,vs,vt,rc,dco,dd,nn,lev,lfile)
      implicit none
      integer lev,nn

      include 'fsizes.inc'
      real r(ns),dco(ns),vt(ns),vs(ns),rc(ns),tb
      real*8 dd(ns,nlev)
c     real ulev(ns,ntra)
      integer i,l,ioerror,lend

      character*250 mol,lfile,leer,newname
      character*40 fmt46,fmt47,fmt43

      integer precomm,prestruc
      logical overwr,fileinp,numdef,anames,allstop
      common /options/ precomm,prestruc,overwr,fileinp,numdef,
     %                   anames,allstop
      save /options/

13    format(' Name of the file for the level populations: ',$)
c     Input formats
42    format(A250)
c     Output lines
44    format(2i4,1P,E11.4)
45    format(i4,1P,5E14.6)
c     Create format codes dynamically fitting to the number of levels
      write(fmt46,'("(i4,1P,",i3,"E14.7)")') nlev+1
      write(fmt47,'("(i4,1P,",i3,"E11.4)")') nlev
      leer=' '
c     *************************************************************
c     Construct or get file name
c     *************************************************************
      if (.not.anames) lfile=leer
350   if (lfile.eq.leer) then
       write(6,13)
       read (5,42) lfile
      else
       call clearstring(lfile,lend,.false.)
       newname=lfile(1:lend)//'.n'
       lfile=newname
      endif
      open(14,err=130,file=lfile,status='new')
250   continue
c     *************************************************************
c     Real output
c     *************************************************************
      call clearstring(mol,lend,.true.)
      write(fmt43,'("(A",i0,")")'), lend
      write(14,fmt43) mol(1:lend)
      write(14,44) lev,nn,tb
      write(14,45) (i,r(i),dco(i),vs(i),vt(i),rc(i),i=1,nn)
      do 100 i=1,nn
      write(14,fmt46) i,(dd(i,l), l=1,lev)
100   continue
c     do 200 i=1,nn
c     write(14,fmt47) i,(ulev(i,l), l=1,lin)
c 200   continue
      close(14)
      return

c     *************************************************************
c     Error handling
c     *************************************************************
130   goto (250,260) ioerror(lfile,3)+1
      return
260   lfile=leer
      goto 350
      end

c     ***********************************************************
c     Subroutine to read the results of a previous calculation
c     giving the level populations at various radii
c     ***********************************************************
      subroutine readdens(mol,tb,r,vs,vt,rc,dco,dd,nn,lev,lfile,kerr)
      implicit none
      integer lev,nn,kerr

      include 'fsizes.inc'
      real r(ns),dco(ns),vt(ns),vs(ns),rc(ns),tb
      real*8 dd(ns,nlev)
      integer i,l,idum,ioerror

      character*250 mol,lfile

13    format(' Name of the file giving the level populations: ',$)
c     Input formats
42    format(A250)

350   write(6,13)
      read (5,42) lfile
      open(14,err=120,file=lfile,status='old')
c     *************************************************************
c     Read everything
c     *************************************************************
      read(14,42,err=110,end=110) mol
      read(14,*,err=110,end=110) lev,nn,tb
      read(14,*,err=110,end=110) 
     %      (idum,r(i),dco(i),vs(i),vt(i),rc(i),i=1,nn)
      do 100 i=1,nn
      read(14,*,err=110,end=110) idum,(dd(i,l), l=1,lev)
100   continue
      close(14)
      return

c     *************************************************************
c     Error handling
c     *************************************************************
120   goto (350) ioerror(lfile,1)
      kerr=1
      return
110   goto (350) ioerror(lfile,2)
      kerr=1
      return
      end

c     ***********************************************************
c     Subroutine to write out the line profiles in a form that
c     can be read by idl, greg or gnuplot
c     ***********************************************************
      subroutine writeline(ibeam,nbeam,nff,dtrans,doff,dny,tau,
     %                                            nu0,istau,mesg,lfile)
      implicit none
      integer nbeam,nff

      include 'fsizes.inc'
      real ibeam(ns,-nf:nf)
      real dtrans,doff,dny,tau,nu0
      integer istruc,icomm,i,j,ioerror
      integer precomm,prestruc
      logical overwr,fileinp,numdef,anames,allstop,istau

c     Variables for FITS file write
      real daten(ns*(2*nf+1))
      real v0(0:3),dv(0:3),convert
      integer nax(3),kerr,lend,mend,usedlen
      character*8 units(0:3)
      character*40 comm(0:3)
      character*4 stdext(3),useext

      character*250 lfile,newname,leer
      character*210 info
      character*10 optdepth
      character*12 freqstring
      character*22 mesg
      character*26 transname
      character*2 comml(2)

      common /options/ precomm,prestruc,overwr,fileinp,numdef,
     %                   anames,allstop
      save /options/

c     Translate from arcsec to degree
      data convert / 0.00027777/
      data comml / '##','!!'/
      data stdext / 'fits','scan','line' /

11    format(/,' Writing line intensity profiles to output file')
12    format(/,' Writing optical depth profiles to output file')
13    format(/,' Name of the output file: ',$)
14    format(' File structure [FITS = 0/lines = 1/compact = 2]: ',$)
15    format(' Comment line format [none = 0/UNIX = 1/VMS = 2]: ',$) 
c     Input formats
42    format(A250)
c     Output lines
81    format(A2,' Offset of first entry: ',E10.3)
82    format(A2,' Frequency of first entry: ',E10.3)
91    format(A2,' Offset step size: ',E10.3)
92    format(A2,' Frequency step size: ',E10.3)
85    format(A2,' Position-velocity map produced by SimLine '
     %    //'(V.Ossenkopf-Okada 1999-2026)')
83    format(A2,' Transition: ',A22)
84    format(A2,' Frequency: ',A20)
86    format(A2,' Number of offset points: ',I4)
87    format(A2,' Number of frequency points: ',I4)
89    format(A2,' Central optical depth: ',E10.3)
88    format(A2,' ')
45    format(1P,E9.2,E11.3,E12.4)
46    format(1P,223E11.4)

c     Opening
      if (istau) then
        write(6,12)
        call clearstring(mesg,usedlen,.true.)
        transname=mesg(1:usedlen)//'-tau'
      else
        write(6,11)
        transname=mesg
      endif
      leer=' '
c     *************************************************************
c     The format of the output files may be selected on compile 
c     time, by command line switch or interactively 
c     *************************************************************
      if (prestruc.ge.0) then
      istruc=prestruc
      else
      write(6,14)
      read (5,*) istruc
      if ((istruc.gt.2).or.(istruc.lt.0)) istruc=0
      endif
      useext=stdext(istruc+1)
c     *************************************************************
c     Construct or get file name
c     *************************************************************
      if (.not.anames) lfile=leer
      call clearstring(transname,mend,.true.)
350   if (lfile.eq.leer) then
        write(6,13)
        read (5,42) lfile
      else
       call clearstring(lfile,lend,.false.)
       newname=lfile(1:lend)//'.'//transname(1:mend)//'.'//useext
       lfile=newname
      endif
c     *************************************************************
c     Separated treatment for ASCII files and Fits files
c     *************************************************************
      if (istruc.eq.0) then
      write(optdepth,'(E10.3)') tau
      write(freqstring,'(E12.5)') nu0
c     Make strings that always fill the 80 character lines in FITS
      info= 'Position-velocity map produced by SimLine'//
     %     ' (V.Ossenkopf-Okada 1997-2026)'//
     %     ' Transition: '//mesg//
     %     '       Rest frequency: '//freqstring//
     %     '   Central optical depth: '//optdepth//
     %     '                '
      nax(1)=nbeam
      nax(2)=2*nff+1
      if (istau) then
        units(0)=''
        comm(0)='Optical depth'
      else
        units(0)='K'
        comm(0)='Beam temperature'
      endif
      units(1)='RA---TAN'
      comm(1)='in degrees'
      units(2)='VELOCITY'
      comm(2)=' in km/s'

      v0(0)=0.0
      dv(0)=1.0
      v0(1)=doff*convert
      dv(1)=dtrans*convert
      v0(2)=-nff*dny
      dv(2)=dny

c     Resort frequencies, put to linear field for writing
      do 100 j=-nff,nff
      do 100 i=1,nbeam
      daten((j+nff)*nbeam+i)=ibeam(i,-j)
100   continue

      kerr=0
      call savefits(2,nax,daten,v0,dv,units,comm,info,nu0,lfile,kerr)
c     Correct file name
      lend=index(lfile,'.'//transname(1:mend)//'.'//useext)
      if (lend.gt.0) then
       newname=lfile(1:lend)//useext
       lfile=newname
      endif
      return
      endif
c     *************************************************************
c     ASCII file - ask for comment line structure
c     *************************************************************
      if (precomm.ge.0) then
      icomm=precomm
      else
      write(6,15)
      read (5,*) icomm
      if ((icomm.gt.2).or.(icomm.lt.0)) icomm=0
      endif
      open(14,err=130,file=lfile,status='new')
250   continue
c     *************************************************************
c     Real output
c     *************************************************************
c     The frequency scale is inverted here. All other parts of the
c     program take higher frequencies as positive, in the output
c     blue shifted frequencies will be given as negative.
c     *************************************************************
c     Header
      if (icomm.gt.0) then
      write(14,85) comml(icomm)
      write(14,83) comml(icomm),mesg
      write(14,84) comml(icomm),freqstring
      write(14,89) comml(icomm),tau
      write(14,86) comml(icomm),nbeam
      write(14,87) comml(icomm),2*nff+1
      write(14,88) comml(icomm)
      endif
c     Block structure
      if (istruc.eq.2) then
      if (icomm.gt.0) then
      write(14,81) comml(icomm),doff
      write(14,91) comml(icomm),dtrans
      write(14,82) comml(icomm),-dny*nff
      write(14,92) comml(icomm),dny
      write(14,88) comml(icomm)
      endif
      do 450 i=1,nbeam
      write(14,46) (ibeam(i,-j),j=-nff,nff)
450   continue
c     3-column structure
      else
      do 200 i=1,nbeam
      do 200 j= -nff,nff
      write(14,45) doff+(i-1)*dtrans,j*dny,ibeam(i,-j)
200   continue
      endif
      close(14)
c     Correct file name
      lend=index(lfile,'.'//transname(1:mend)//'.'//useext)
      if (lend.gt.0) then
       newname=lfile(1:lend)//useext
       lfile=newname
      endif
      return

c     *************************************************************
c     Error handling
c     *************************************************************
130   goto (250,260) ioerror(lfile,3)+1
      return
260   lfile=leer
      goto 350
      end

c     *************************************************************
c     Subroutine to write out the line profiles in a form that
c     can be read by greg or gnuplot
c     ***********************************************************
      subroutine savefits(ndim,nx,daten,v0,dv,units,comments,
     %                                          info,nu0,lfile,kerr)
      implicit none
      integer ndim,nx(3),kerr
      real v0(0:3),dv(0:3),nu0
      character*8 units(0:3)
      character*40 comments(0:3)
      real daten(*)
      character*250 lfile,leer
      character*210 info

      integer i,ioerror
      integer status,blocksize,bitpix,naxis,naxes(3)
      integer group,fpixel,nelements
      logical simple,extend
      character*6 ctype(3),crpix(3),crval(3),cdelt(3)

      data ctype /'CTYPE1','CTYPE2','CTYPE3'/
      data crpix /'CRPIX1','CRPIX2','CRPIX3'/
      data crval /'CRVAL1','CRVAL2','CRVAL3'/
      data cdelt /'CDELT1','CDELT2','CDELT3'/

c     Input formats
13      format(/,' Name of the output FITS file: ',$)
42    format(A250)

c     *************************************************************
c     Open the file for writing
c     *************************************************************
      blocksize=1
      leer=' '
      if (lfile.ne.leer) goto 200
350   write(6,13)
      read (5,42) lfile
200   status=0
      call ftinit(14,lfile,blocksize,status)
      if (status.ne.0) then
          call ftclos(14,status)
          goto (150,250,350) ioerror(lfile,3)+2
        endif
c     *************************************************************
c     Required header
c     *************************************************************
      simple=.true.
      bitpix=-32
      extend=.true.
      naxis=ndim
      do 210 i=1,ndim
        naxes(i)=nx(i)
210   continue
      call ftphpr(14,simple,bitpix,naxis,naxes,0,1,extend,status)
c     *************************************************************
c     Extended header
c     *************************************************************
      call ftpkys(14,'BUNIT ',units(0),comments(0),status)
      call ftpkye(14,'BZERO ',v0(0),4,' ',status)
      call ftpkye(14,'BSCALE',dv(0),4,' ',status)
      do 220 i=1,ndim
       call ftpkys(14,ctype(i),units(i),comments(i),status)
       call ftpkyj(14,crpix(i),1,' ',status)
       call ftpkye(14,crval(i),v0(i),4,' ',status)
       call ftpkye(14,cdelt(i),dv(i),4,' ',status)
220   continue
      call ftpkye(14,'RESTFREQ',nu0,5,' ',status)
      call ftpcom(14,info,status)
      call ftpdat(14,status)
c     ************************************************************
c     Write field
c     ************************************************************
      group=1
      fpixel=1
      nelements=1
      do 110 i=1,ndim
        nelements=nelements*nx(i)
110   continue
      call ftppre(14,group,fpixel,nelements,daten,status)
      call ftclos(14,status)
      if (status.gt.0) then
       call printerror(status)
       kerr=1
      endif
      return

c     *************************************************************
c     Error handling
c     *************************************************************
150   kerr=1
      return
250   status=0
      close(14)
      call deletefile(lfile,status)
      goto 200
      end

c     *************************************************************
c     String processing - find name without extension and blanks
c     *************************************************************
      subroutine clearstring(string,length,ignoreext)
      implicit none
      character*(*) string
      integer length
      logical ignoreext

      character*255 sbuffer,tbuffer
      integer tlen,bpos,lpos,last

c     Remove leading blanks
100   bpos=index(string,' ')
      if (bpos.eq.1) then
       sbuffer=string(2:)
       string=sbuffer
       goto 100
      endif
c     Remove trailing blanks
      tlen=len(string)
300   sbuffer=string(tlen:)
      bpos=index(sbuffer,' ')
      if (bpos.eq.1) then
         tlen=tlen-bpos
       goto 300
      endif
      last=tlen
c     All if extension not to be removed
      if (ignoreext) then
       length=last
       return
      endif
c     Find last dot
      sbuffer=string
      lpos=0
200   bpos=index(sbuffer,'.')
      lpos=lpos+bpos
      if (bpos.gt.0) then
       tbuffer=sbuffer(bpos+1:)
       sbuffer=tbuffer
       goto 200
      endif
c     Exclude dots in path
      if (lpos.eq.last) then
        length=lpos-1
      else if (string(lpos+1:lpos+1).eq.'/') then
        length=last
      else
        length=lpos-1
      endif
      return
      end

c     ***********************************************************
c     Subroutine for general file I/O error handling
c     *************************************************************
      integer function ioerror(lfile,inum)
      implicit none
      character*250 lfile
      character yn
      integer inum

      integer precomm,prestruc
      logical overwr,fileinp,numdef,anames,allstop
      common /options/ precomm,prestruc,overwr,fileinp,numdef,
     %                   anames,allstop
      save /options/

41    format(A1)
51    format(/,' I cannot open the input file ',A250)
52    format(//,' Format error in input file ',A250)
53    format(/,' This file exists already: ',A250)
54    format(' Shall I replace it ? [y/n]: ',$)
56    format(' I cannot open the output file ',A250)
57    format(' Use another file ? [y/n]: ',$)
55    format(' Continue ? [y/n]: ',$)
59    format(//,' Stop on all errors requested. Program halted.')

      ioerror=-1
      if (inum.eq.1) then
c     Error opening the input file
      write(6,51) lfile
      if (allstop) then
       write(6,59)
       stop 2
      endif
      write(6,55)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) ioerror=1
      return
      else if (inum.eq.2) then
c     The file has not the structure required
      write(6,52) lfile
      close(14)
      if (allstop) then
       write(6,59)
       stop 3
      endif
      write(6,55)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) ioerror=1
      return
      else if (inum.eq.3) then
c     Error opening the output file - Test if the file exists
      open(14,err=140,file=lfile,status='old')
      close(14)
      if (.not.overwr) then
       write(6,53) lfile
       if (allstop) then
        write(6,59)
        stop 4
       endif
       write(6,54)
       read (5,41) yn
       if ((yn.ne.'y').and.(yn.ne.'Y')) goto 150
      endif
      open(14,err=140,file=lfile,status='old')
      ioerror=0
      return
      endif

c     The file is not writable
140   write(6,56) lfile
      if (allstop) then
       write(6,59)
       stop 4
      endif
150   write(6,57)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) ioerror=1
      return
      end

c     ***********************************************************
c     Set predefined constants for the output file handling and 
c     the numerical parameters via compiler option or command
c     line argument. ! This part is not standard-FORTRAN77 and 
c     may not compile everywhere ! Taking defaults, you will be 
c     asked for the options at all file output operations.
c     ***********************************************************
      subroutine predefine(allexit,part1only,tauwrite)
      integer nmax,narg,iargc,i
      integer precomm,prestruc
      logical overwr,numdef,fileinp,anames,allstop,allexit,
     %  part1only,tauwrite
      character*20 string
      character*20 stover(2),stblock(2),stlines(2),stnocom(2),
     %  stunixcom(2),stvmscom(2),stnumdef(2),stfile(2),stfits(2),
     %  stanames(2),ststop(2),stpopulate(2),sttauwrite(2)

      common /options/ precomm,prestruc,overwr,fileinp,numdef,
     %                   anames,allstop
      save /options/

c     All allowed possibilities of command line parameters
      data stfits /'-fits','-FITS'/
      data stblock /'-block','-BLOCK'/
      data stlines /'-lines','-LINES'/
      data stnocom /'-nocomments','-NOCOMMENTS'/
      data stunixcom /'-unixcomments','-UNIXCOMMENTS'/
      data stvmscom /'-vmscomments','-VMSCOMMENTS'/
      data stover /'-overwrite','-OVERWRITE'/
      data stfile /'-fileinput','-FILEINPUT'/
      data stnumdef /'-numdefaults','-NUMDEFAULTS'/
      data stanames /'-autonames','-AUTONAMES'/
      data ststop /'-stoponerror','-STOPONERROR'/
      data stpopulate /'-populateonly','-POPULATEONLY'/
      data sttauwrite /'-writetau','-WRITETAU'/

c     ***********************************************************
c     Default values (interactive mode)
c     ***********************************************************
      precomm=-1
      prestruc=-1
      overwr=.false.
      fileinp=.false.
      numdef=.false.
      anames=.false.
      allstop=.false.
      part1only=.false.
      tauwrite=.false.
c     ***********************************************************
c     Compile time defaults
c     (may be activated on systems with preprocessor cpp)
c     ***********************************************************
c     #if defined (overwrite)
c     overwr=.true.
c     #endif
c     #if defined (block)
c     prestruc=2
c     #elif defined (lines)
c     prestruc=1
c     #elif defined (fits)
c     prestruc=0
c     #endif
c     #if defined (nocomments)
c     precomm=0
c     #elif defined (unixcomments)
c     precomm=1
c     #elif defined (vmscomments)
c     precomm=2
c     #endif
c     #if defined (writetau)
c     tauwrite=.true.
c     #endif
c     ************************************************************
c     Overwrite compile time defaults by command line parameters
c     ************************************************************
      nmax=iargc()
      do 100 narg=1,nmax
      call getarg(narg,string)
      do 100 i=1,2
      if (string.eq.stblock(i)) prestruc=2
      if (string.eq.stlines(i)) prestruc=1
      if (string.eq.stfits(i)) prestruc=0
      if (string.eq.stnocom(i)) precomm=0
      if (string.eq.stunixcom(i)) precomm=1
      if (string.eq.stvmscom(i)) precomm=2
      if (string.eq.stover(i)) overwr=.true.
      if (string.eq.stfile(i)) fileinp=.true.
      if (string.eq.stnumdef(i)) numdef=.true.
      if (string.eq.stanames(i)) anames=.true.
      if (string.eq.ststop(i)) allstop=.true.
      if (string.eq.stpopulate(i)) part1only=.true.
      if (string.eq.sttauwrite(i)) tauwrite=.true.
100   continue
        allexit=allstop
      return
      end

