c     Field sizes for radiative transfer computation	
c     radial grid, frequency grid, transitions, beam grid
      integer ns,nf,nlev,ntra,nstretch
      parameter(ns=1050,nf=600,nlev=110,ntra=320,nstretch=1800)
c     nstretch should be 1.5 to 2 times ns
c
c     Field sizes for reading the molecule files
c     	  they should be larger or equal to the field sizes above
      integer NLEVMAX,NLINMAX,NCOLMAX,NTEMPMAX
      parameter(NLEVMAX=110,NLINMAX=320,NCOLMAX=3800,NTEMPMAX=45)
