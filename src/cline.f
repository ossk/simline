c     ***********************************************************
c     Program to extract the single line with zero beam offset
c     from the map of lines produced by LTR
c     ***********************************************************
      program cline
      implicit none
      integer nff,nf

      include 'fsizes.inc'
      real ibeam(nf),ny(nf)
      character yn

55    format(' Continue ? [y/n]: ',$)

100   call readline(ny,ibeam,nff)
      call writeclin(ny,ibeam,nff)
      write(6,55)
      read (5,'(A1)') yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) goto 100
      end


c     ***********************************************************
c     Subroutine to read the results of a previous calculation
c     giving the level populations at various radii
c     ***********************************************************
      subroutine readline(ny, ibeam, nff)
      implicit none
      integer nff,nf

      include 'fsizes.inc'
      real ibeam(nf),ny(nf),adum
      integer i

      character*250 lfile,dummy
      character yn

13    format(/,' Name of the file with the line map: ',$)
c     Input formats
41    format(A1)
42    format(A250)

350   write(6,13)
      read (5,42) lfile
      open(14,err=120,file=lfile,status='old')
c     *************************************************************
c     Read everything
c     *************************************************************
      read(14,42) dummy
      read(14,42) dummy
      read(14,'(A30,I6)') dummy,nff
      read(14,42) dummy

      do 200 i=1,nff
      read(14,*,end=110,err=110) adum,ny(i),ibeam(i)
200   continue
      close(14)
      return

c     *************************************************************
c     Error handling
c     *************************************************************	
51    format(/,' I cannot open the input file ',A250)
52    format(//,' Format error in input file ',A250)
55    format(' Continue ? [y/n]: ',$)

c     Error opening the input file
120   write(6,51) lfile
      write(6,55)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) goto 350
      stop 2
c     The file has not the structure required
110   write(6,52) lfile
      write(6,55)
      read (5,41) yn
      if ((yn.eq.'y').or.(yn.eq.'Y')) goto 350
      stop 3
      end

c     ***********************************************************
c     Subroutine to write out the line profiles in a form that
c     can be read by greg or gnuplot
c     ***********************************************************
      subroutine writeclin(ny,ibeam,nff)
      implicit none
      integer nff,nf

      include 'fsizes.inc'
      real ibeam(nf),ny(nf)
      integer i

      character*250 lfile
      character yn

13    format(' Name of the output file for the line profile: ',$)
c     Input formats
41    format(A1)
42    format(A250)
43    format(1P, 2E11.4)

350   write(6,13)
      read (5,42) lfile
      open(14,err=130,file=lfile,status='new')
250   continue
c     *************************************************************
c     Real output
c     *************************************************************
      do 200 i=1,nff
      write(14,43) ny(i),ibeam(i)
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
      end
