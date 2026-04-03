c     *************************************************************
c     Subroutine to read molecular line data from a file
c     Adapted version from the Leiden workshop
c     *************************************************************
      subroutine readspecies(filename,nlin,lev,kerr)
      implicit none
      integer nlin,lev, kerr
      character*250 filename
      include 'fsizes.inc'

c     *************************************************************
c     Original arrays for reading
c     They have to large enough to read all data provided.
c     In general only a subset of them will be used
c     *************************************************************
      real energy_lev(NLEVMAX), gdeg(NLEVMAX)
      real Aud(NLINMAX), Kud(NCOLMAX,NTEMPMAX), Temp_Kud(NTEMPMAX)
        real Kud_op(NCOLMAX,NTEMPMAX), ortho(NTEMPMAX)
      integer nlines,ncoltrans,nlevels,ncttemp,ncspecies
      integer lnum(NLEVMAX),lev_up(NLINMAX),lev_down(NLINMAX),
     %        lev_ct_up(NCOLMAX),lev_ct_down(NCOLMAX)
      character*10 names(NLEVMAX)

c     Auxiliary variables for reading
      integer i,itemp,lmax,lmin,dummyi,dummyj,dummyk
      character*250 dummyline,leer,lamdastring
      logical overwr,fileinp,numdef,anames,allstop,lamdaformat
      integer precomm,prestruc,ioerror

c     ***********************************************************
c     New arrays as they are used in the code
c     ***********************************************************
      real molmass,energy(nlev),stgew(nlev)
      real fr(ntra),eina(ntra),einb(ntra),gewratio(ntra)
      real temp(NTEMPMAX),collrates((nlev*(nlev-1))/2,NTEMPMAX)
      real clight,cst1,cst2,frmin
      integer lfrom(ntra),lto(ntra)
      integer icount,ntemp
      logical downwards
      character*10 levnam(nlev)

      common /options/ precomm,prestruc,overwr,fileinp,numdef,
     %                   anames,allstop
      common /molecule/ molmass, energy, stgew, levnam
      common /coeffein/ fr,eina,einb,gewratio
      common /crates/ downwards,ntemp,temp,collrates
      common /lintable/ lfrom, lto
      save /molecule/, /coeffein/, /crates/, /lintable/,/options/

        data lamdastring / '!MOLECULE' /

c     ***********************************************************
c     Translation coefficients into my units
c     ***********************************************************
c     clight = light velocity in cm/s
c     cst1 = clight*h/kB  (Translate energies into K)
c     cst2=(8*pi/c^2)^(1/3) (Tranlation from Einstein As to Bs)
c     ***********************************************************
      data clight,cst1,cst2 / 2.997925e10, 1.4387687, 3.03529e-7/
c     ***********************************************************
c     Minimum frequency used if neighbouring levels are too close,
c     assume a 100 MHz maser there
c     ***********************************************************
      data frmin /1.0e8/

c     Error messages
13    format(' Name of the molecular data file: ',$)
12    format(/,' Reading molecular data file: ',A250)
14    format(/,' Number of levels in the molecule file too large !')
15    format(/,' Number of lines in the molecule file too large !')
16    format(/,' Number of temperatures in the molecule file too',
     %          ' large !')
17    format(/,' Number of collision rates in the molecule file',
     %          ' too large !')
18    format(/,' Inconsistency in the level numbering detected !')
19    format(/,' Wrong order of level numbering detected !')
21    format(/,' Cannot yet handle multiple collision partners !')
22    format(/,' Cannot only handle H2 as collision partner !')
23    format(/,' First collision partner must be para-H2 !')
24    format(/,' Inconsistency between para- and ortho-H2 rates !')
c     Input formats
41    format(A1)
42    format(A250)

c     ***********************************************************
c     Start reading
c     ***********************************************************
      kerr=0
      leer=' '
100   if (filename.eq.leer) then
      write(6,13)
      read (5,42) filename
      else
      write(6,12) filename
      endif
      open(unit=14,err=120,file=filename,status='old')
c       Recognize file format from the first line
      read (14,42,err=110,end=110) dummyline
        lamdaformat=(dummyline.eq.lamdastring)
c     Global parameters
      read (14,42,err=110,end=110) dummyline
      if (lamdaformat) read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) molmass
      read (14,42,err=110,end=110) dummyline
      read (14,*,err=110,end=110) nlevels
c     Check for readability
      if (nlevels.gt.nlevmax) then
       write(6,14)
       goto 500
      endif
        if (.not.lamdaformat) then 
        read (14,42,err=110,end=110) dummyline
        read (14,*,err=110,end=110) nlines
c       Check for readability
        if (nlines.gt.nlinmax) then
         write(6,15)
         goto 500
        endif
      endif
c     Read all level parameters
      read (14,42,err=110,end=110) dummyline
      do i=1,nlevels
         if (lamdaformat) then
          read (14,*,err=110,end=110)
     %     lnum(i),energy_lev(lnum(i)),gdeg(lnum(i)),names(lnum(i))
         else
          read (14,*,err=110,end=110)
     %     lnum(i),gdeg(lnum(i)),energy_lev(lnum(i)),names(lnum(i))
         endif
      enddo
      if (lamdaformat) then 
        read (14,42,err=110,end=110) dummyline
        read (14,*,err=110,end=110) nlines
c       Check for readability
        if (nlines.gt.nlinmax) then
         write(6,15)
         goto 500
        endif
      endif
c     ***********************************************************
c     Consistency check for the level numbers
c     It is assumed that levels are numbered starting 1 and
c     sequentially up to nlevels with increasing energy
c     ***********************************************************
      lmax=0
      lmin=nlevmax
      do i=1,nlevels
      lmax=max(lmax,lnum(i))
      lmin=min(lmin,lnum(i))
      enddo
      if ((lmin.ne.1).or.(lmax.ne.nlevels)) then
       write(6,18)
       goto 500
      endif
c     ***********************************************************
c     Now check for the correct order of level numbering
c     ***********************************************************
      do i=2,nlevels
      if (energy_lev(i).lt.energy_lev(i-1)) then
        write(6,19)
        goto 500
      endif
      enddo

c     ***********************************************************
c     First reduction to used values
c     ***********************************************************
      lev=min(nlevels,nlev)
      do i=1,lev
        energy(i)=energy_lev(i)*cst1
        stgew(i)=gdeg(i)
        levnam(i)=names(i)
      enddo

c     ***********************************************************
c     Further input - Einstein A`s
c     ***********************************************************
      read (14,42,err=110,end=110) dummyline
      do i=1,nlines
         if (lamdaformat) then
        read(14,*,err=110,end=110) 
     %     dummyi,lev_up(i),lev_down(i),Aud(i)
         else
        read(14,*,err=110,end=110) 
     %     lev_up(i),lev_down(i),Aud(i)
         endif
      enddo
c     ***********************************************************
c     Consistency check for the lines
c     ***********************************************************
      lmax=-nlinmax
      lmin=nlinmax
      do i=1,nlines
      lmin=min(lmin,lev_down(i))
      lmax=max(lmax,lev_down(i)-lev_up(i))
      enddo
      if ((lmin.lt.1).or.(lmax.gt.(-1))) then
       write(6,19)
       goto 500
      endif
c     **************************************************************
c     Second reduction step for Einstein coefficients
c     **************************************************************
      i=0
      nlin=0
200   continue
        i=i+1
        if (lev_up(i).le.lev) then
        nlin=nlin+1
        lfrom(nlin)=lev_up(i)
        lto(nlin)=lev_down(i)
        fr(nlin)=max((energy_lev(lev_up(i))-energy_lev(lev_down(i)))
     %               *clight,frmin)
        gewratio(nlin)=gdeg(lev_up(i))/gdeg(lev_down(i))
        eina(nlin)=Aud(i)
        einb(nlin)=eina(nlin)/(cst2*fr(nlin))**3
        endif
      if ((nlin.lt.ntra).and.(i.lt.nlines)) goto 200

c     **************************************************************
c     Now we read the collisional rate coefficients
c     **************************************************************
      read (14,42,err=110,end=110) dummyline
        if (lamdaformat) then
        read (14,*,err=110,end=110) ncspecies
          if (ncspecies.gt.2) then
            write(6,21)
          goto 500
        endif
        endif
      read (14,42,err=110,end=110) dummyline
        if (lamdaformat) then
        read (14,42,err=110,end=110) dummyline
c         Check for H2
          if (index(dummyline,'H2').eq.0) then
          write(6,22)
          goto 500
        else 
     %     if ((ncspecies.eq.2).and.(index(dummyline,'pH2').eq.0)) then
          write(6,23)
          goto 500
          endif
          read (14,42,err=110,end=110) dummyline
        read (14,*,err=110,end=110) ncoltrans
        read (14,42,err=110,end=110) dummyline
        read (14,*,err=110,end=110) ncttemp
        else
        read (14,*,err=110,end=110) ncttemp
        endif
c     Check for readability
      if (ncttemp.gt.NTEMPMAX) then
       write(6,16)
       goto 500
      endif
c     Read temperatures
      read (14,42,err=110,end=110) dummyline
        if (lamdaformat) then
          read(14,*) (temp_kud(i),i=1,ncttemp)
        else
        do i=1,ncttemp
         read (14,*,err=110,end=110) temp_kud(i)
        enddo
        endif
        if (.not.lamdaformat) then
        read (14,42,err=110,end=110) dummyline
        read (14,*,err=110,end=110) ncoltrans
        endif
      if (ncoltrans.gt.ncolmax) then
       write(6,17)
       goto 500
      endif
c     Main loop for reading collision parameters
      read (14,42,err=110,end=110) dummyline
      do i=1,ncoltrans
         if (lamdaformat) then
         read (14,*,err=110,end=110) dummyi,lev_ct_up(i),
     %        lev_ct_down(i),(Kud(i,itemp),itemp=1,ncttemp)
         else
         read (14,*,err=110,end=110) lev_ct_up(i),lev_ct_down(i),
     %        (Kud(i,itemp),itemp=1,ncttemp)
         endif
      enddo
c       Read ortho rates from lamda files
        if (lamdaformat.and.(ncspecies.eq.2)) then
        read (14,42,err=110,end=110) dummyline
        read (14,42,err=110,end=110) dummyline
        if (index(dummyline,'oH2').eq.0) then
          write(6,23)
          goto 500
          endif
          read (14,42,err=110,end=110) dummyline
        read (14,*,err=110,end=110) dummyi
        if (ncoltrans.ne.dummyi) then
          write(6,24)
          goto 500
          endif
        read (14,42,err=110,end=110) dummyline
        read (14,*,err=110,end=110) dummyi
        if (ncttemp.ne.dummyi) then
          write(6,24)
          goto 500
          endif
c         Ignore new temperature entries
        read (14,42,err=110,end=110) dummyline
        read (14,42,err=110,end=110) dummyline
c         Create table of ortho-para rates
          do itemp=1,ncttemp
            call orthopara(temp_kud(itemp),ortho(itemp))
          enddo
c       Main loop for reading collision parameters
        read (14,42,err=110,end=110) dummyline
        do i=1,ncoltrans
         read (14,*,err=110,end=110) dummyi,dummyj,
     %        dummyk,(Kud_op(i,itemp),itemp=1,ncttemp)
         if ((lev_ct_up(i).ne.dummyj).or.
     %                    (lev_ct_down(i).ne.dummyk)) then
            write(6,24)
            goto 500
           else
             do itemp=1,ncttemp
               Kud(i,itemp)=Kud(i,itemp)*(1.0-ortho(itemp))
     %                               +ortho(itemp)*Kud_op(i,itemp)
             enddo
           endif
        enddo
        endif
      close(14)

c     **************************************************************
c     Final reduction
c     **************************************************************
      downwards=(lev_ct_up(1).gt.lev_ct_down(1))
      ntemp=ncttemp
      do itemp=1,ntemp
         temp(itemp)=Temp_Kud(itemp)
      enddo
      do i=1,lev*(lev-1)/2
      do itemp=1,ntemp
         collrates(i,itemp)=0.0
      enddo
      enddo
      if (downwards) then
       do i=1,ncoltrans
         if ((lev_ct_up(i).le.lev).and.
     %                      (lev_ct_up(i).gt.lev_ct_down(i))) then
         icount=((lev_ct_up(i)-1)*(lev_ct_up(i)-2))/2+lev_ct_down(i)
         do itemp=1,ntemp
           collrates(icount,itemp)=Kud(i,itemp)
         enddo
         endif
       enddo
      else
       do i=1,ncoltrans
         if ((lev_ct_down(i).le.lev).and.
     %                      (lev_ct_down(i).gt.lev_ct_up(i))) then
         icount=((lev_ct_down(i)-1)*(lev_ct_down(i)-2))/2+lev_ct_up(i)
         do itemp=1,ntemp
           collrates(icount,itemp)=Kud(i,itemp)
         enddo
         endif
       enddo
      endif
      return

c     *************************************************************
c     Error handling
c     *************************************************************
500   filename=leer
      close(14)
      if (allstop) then
       kerr=1
       return
      endif
      goto 100
120   filename=leer
      goto (100) ioerror(filename,1)
      kerr=1
      return
110   filename=leer
      goto (100) ioerror(filename,2)
      kerr=1
      return

      end

c     ***********************************************************
c     Calling routine for the computation of the collisional
c     rate coefficients
c     The calling program does not know about special molecules
c     ***********************************************************
      real function cij(T,J,K)
      implicit none
      real T
      integer j,k

      include 'fsizes.inc'
      real molmass, energy(nlev),stgew(nlev),coldown,boltzmann
      real temp(NTEMPMAX), coll((nlev*(nlev-1))/2,NTEMPMAX)
      logical downwards
      integer jt,ilocat,jup,jlow,icount,ntemp
      character*10 levnam(nlev)

        common /molecule/ molmass, energy, stgew, levnam
        common /crates/ downwards,ntemp,temp,coll
        save /molecule/, /crates/

c     **************************************************************
c     Throw away nontransitional collisions
c     **************************************************************
      if (j.eq.k)  then
      cij = 0.
      return
      endif
c     **************************************************************
c     Level numbers
c     **************************************************************
      jup   = max(j,k)
      jlow  = min(j,k)
      icount=((jup-1)*(jup-2))/2+jlow
c     **************************************************************
c     Square root extrapolation below the lowest rate point
c     ************************************************************** 
      if (t.lt.temp(1)) then
      coldown = coll(icount,1)*sqrt(t/temp(1))
c     **************************************************************
c     Linear interpolation for upward transition from the table
c     assumption of detailed balance for downward transitions
c     cul=clu*gl/gu*exp(h*ny/kT)
c     **************************************************************
      else
      jt = ilocat(temp, ntemp, T)
      coldown = coll(icount,jt) +(t-temp(jt))/(temp(jt+1)-temp(jt))
     %               *(coll(icount,jt+1)-coll(icount,jt))
      endif
      if (coldown.gt.0.0) then
c     transitions in the opposite direction
c     Prevent overshooting for irrelevant rates exceeding real*4
      if (downwards.eqv.(j.eq.jlow)) then
        boltzmann=min(60.0,(energy(j)-energy(k))/T)
        cij = coldown*stgew(k)/stgew(j)*exp(boltzmann)
      else
        cij = coldown
      endif
      else
        cij = 0.0
      endif
c     **************************************************************
c     To get cm^3/s all values have to be multiplied by 1e-10
c     **************************************************************
      return
      end

c     *************************************************************
c     function to select the next lower index in a field of
c     parameters for a given parameter
c     x < xy(1) < xy(n) -> j = 1
c     xy(j) < x < xy(j+1) -> j = j
c     x > xy(n) ->  j = n-1
c     It has to be guaranteed that xy is ordered positively
c     x will be increased by a small factor to find close matches
c     *************************************************************
      integer function ilocat(xy, n, x)
      implicit none
      real xy(*),x,xh
      integer n,jlow,jup,jm

      jlow = 1
      jup  = n
      xh = x*1.000002

100   jm = (jup + jlow) / 2
      if (xh.gt.xy(jm)) then
       jlow = jm
      else
       jup  = jm
      end if
      if (jup-jlow.gt.1) goto 100
       
      if (jlow.eq.n) then
      ilocat = jlow - 1
      else
      ilocat = jlow
      endif
      return
      end

c     *************************************************************
c       Function to determine the otho/para-H2 ratio in the gas
c       at a given temerature. Here, I use a close approximation
c       to the chemical model of Le Bourlot (1990, A&A 242, 135)
c       ***********************************************************
        subroutine orthopara(t,op1)
        implicit none
        real t,op1,op2
        real t1,t2

c       thresholds for the thermochemical balance
        data t1,t2 /18.73, 155.2 /

        if (t.lt.t1) then
        op1=1e-3
        else if (t.gt.t2) then
        op1=3.0
        else
        op1=9.0*exp(-170.5/t)
        endif
        op1=op1/(1.0+op1)
        op2=1.0-op1
        return
        end

c     *************************************************************
c     Small routine for string processing - remove leading blanks
c     *************************************************************
      subroutine getfirstentry(string,returnstring)
      implicit none
      character*(*) string,returnstring

      character*255 sbuffer
      integer bpos

c     Remove leading blanks
140   bpos=index(string,' ')
      if (bpos.eq.1) then
       sbuffer=string(2:)
       string=sbuffer
       goto 140
      endif
        if (bpos.gt.1) then
         returnstring=string(1:bpos-1)
       sbuffer=string(bpos:)
       string=sbuffer
        else
         returnstring=''
        endif
        return
        end

