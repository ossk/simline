c     ***********************************************************
c     Program to determine additional test parameters from
c     the level densities
c     ***********************************************************
      program testcloud
      integer nn,lev

      include 'fsizes.inc'
      real r(ns),vs(ns),vt(ns),dco(ns),rc(ns)
      real*8 dd(ns,nlev),probtherm,dnn
      real tex(ns,ntra)
      real cfreq, cst1,cst2,cst3,tt
      integer lfrom(ntra),lto(ntra)
      integer ilev,i,j,itm,lin,kerr
      character*250 mol
      character yn

        common /lintable/ lfrom, lto
        save /lintable/

c     cst1=h*c^2/(8*pi*k)
c     cst2=h/k
c     cst3=(8*pi/c^2)^(1/3)
      data cst1,cst2,cst3/ 1.71622e9,4.79922e-11,3.03529e-7/

11    format(/,' Do you want to compute excitation temperatures,',/,
     %  23X,' thermalization factors [1/2] ? ',$)
12    format(/,' Kinetic temperature [K]: ',$)
15      format(/,' Error reading the molecule file !')
55    format(' Continue ? [y/n]: ',$)

1000  kerr=0
      call readdens(mol,tb,r,vs,vt,rc,dco,dd,nn,lev,kerr)
      call readspecies(mol,lin,ilev,kerr)
          if (kerr.ne.0) then
          write(6,15)
          stop 3
          endif
        if (ilev.ne.lev) call linfromlev(lin,lev)

      write(6,11)
      read(5,*) itm
c     ***************************************************************
c     Compute the different parameters
c     The last level is numerically insignificant and therefore
c     omitted
c     ***************************************************************
      if (itm.eq.1) then
      do 120 j=1,lin
      l1=lfrom(j)
      l2=lto(j)
      do 120 i=1,nn
      tex(i,j)=cst2*cfreq(j)/log(bji(j)*dd(i,l2)/(bij(j)*dd(i,l1)))
120   continue
      lev=lin-1
      else 
      write(6,12)
      read(5,*) tt
      dnn=0.0
      do 210 j=1,lev-1
210   dnn=dnn+probtherm(tt,j)
      do 220 i=1,nn
      do 220 j=1,lev-1
      tex(i,j)=dd(i,j)*dnn/probtherm(tt,j)
220   continue
      lev=lev-1
      endif
c     **************************************************************
c     Output
c     **************************************************************
      call writetemp(r,tex,nn,lev)
      write(6,55)
      read (5,'(A1)') yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) goto 1000
      end

c     ***********************************************************
c     Subroutine to write out the line profiles in a form that
c     can be read by greg or gnuplot
c     ***********************************************************
      subroutine writetemp(r,ibeam,nn,lev)
      implicit none
      integer nn,lev

      include 'fsizes.inc'
      real ibeam(ns,ntra),r(ns)
      integer icomm,i,j

      character*250 lfile
      character*2 comml(2)
      character yn
      data comml / '##','!!'/

13    format(/,' Name of the output file: ',$)
15    format(' Comment line format [none = 0/UNIX = 1/VMS = 2]: ',$) 
c     Input formats
41    format(A1)
42    format(A250)
c     Output lines
86    format(A2,' Number of radial points: ',I4)
87    format(A2,' Number of transitions: ',I4)
88    format(A2,' ')
45    format(I3,1P,2E11.3)

350   write(6,13)
      read (5,42) lfile
      open(14,err=130,file=lfile,status='new')
250   continue
      write(6,15)
      read (5,*) icomm
      if ((icomm.gt.2).or.(icomm.lt.0)) goto 180
c     *************************************************************
c     Real output
c     *************************************************************
c     Header
      if (icomm.gt.0) then
      write(14,86) comml(icomm),nn
      write(14,87) comml(icomm),lev
      write(14,88) comml(icomm)
      endif
      do 200 i=1,lev
      do 200 j=1,nn
      write(14,45) i,r(j),ibeam(j,i)
200   continue
      close(14)
c     *************************************************************
c     End - exit in case of errors
c     *************************************************************
500   return

c     *************************************************************
c     Error handling
c     *************************************************************	
53    format(/,' This file exists already: ',A250)
54    format(' Shall I replace it ? [y/n]: ',$)
56    format(' I cannot open the output file ',A250)
57    format(' Use another file ? [y/n]: ',$)
55    format(' Continue ? [y/n]: ',$)
58    format(/,' You may only choose from the given numbers !',/)

c     Error opening the output file - Test if the file exists
130   open(14,err=140,file=lfile,status='old')
      close(14)
      write(6,53) lfile
      write(6,54)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) then
      open(14,err=140,file=lfile,status='old')
      goto 250
      else
      goto 150
      endif
c     The file is not writable
140   write(6,56) lfile
150   write(6,57)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) goto 350
      goto 500
c     Wrong selection number
180   write(6,58)
      goto 250
      end

