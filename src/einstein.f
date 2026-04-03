c     ***********************************************************
c     Functions to provide the Einstein coefficients
c     for a given molecule and transition
c     J is the starting level of a transition,i.e. aij(1)=A(1->0)
c     The level counting starts from 0!
c     The fields read must be initialized by initeinstein before
c     these functions are called
c     ***********************************************************
      real function aij(J)
      implicit none
      integer j

      include 'fsizes.inc'
      real fr(ntra), eina(ntra),einb(ntra),stgew(ntra)
      common /coeffein/ fr,eina,einb,stgew
      save /coeffein/
      
      aij=eina(j)
      return
      end

c     **************************************************************
c     Coefficients for induced transitions
c     **************************************************************
      real*8 function bij(j)
      implicit none
      integer j

      include 'fsizes.inc'
      real fr(ntra), eina(ntra),einb(ntra),stgew(ntra)
      common /coeffein/ fr,eina,einb,stgew
      save /coeffein/
      
      bij=dble(einb(j))
      return
      end

c     **************************************************************
      real*8 function bji(j)
      implicit none
      integer j

      include 'fsizes.inc'
      real fr(ntra), eina(ntra),einb(ntra),stgew(ntra)
      common /coeffein/ fr,eina,einb,stgew
      save /coeffein/
      
      bji=dble(einb(j))*dble(stgew(j))
      return
      end
 
c     ***********************************************************
c     Function to provide the frequency of a transition
c     ***********************************************************
      real function cfreq(j)
      implicit none
      integer j

      include 'fsizes.inc'
      real fr(ntra), eina(ntra),einb(ntra),stgew(ntra)
      common /coeffein/ fr,eina,einb,stgew
      save /coeffein/
      
      cfreq=fr(j)
      return
      end

c     ***********************************************************
c     Function giving the thermal population density
c     ! Not yet normalized !
c     ***********************************************************
      real*8 function probtherm(t,l)
      implicit none
      real t
      integer l

      include 'fsizes.inc'
      real molmass, energy(nlev),stgew(nlev)
      character*10 levnam(nlev)

      common /molecule/ molmass, energy, stgew, levnam
      save /molecule/

      probtherm=stgew(l)*exp(-energy(l)/T)
      return
      end

c     ***********************************************************
c     Square of the thermal velocity of the molecules depending 
c     on the temperature for a given molecule (in km^2/s^2)
c     ***********************************************************
      real function vtherm(t)
      implicit none
      real t,pre

      include 'fsizes.inc'
      real molmass, energy(nlev),stgew(nlev)
      character*10 levnam(nlev)

      common /molecule/ molmass, energy, stgew, levnam
      save /molecule/

c     Prefactors 2k/m_H
c     mH = 1.672623e-27 kg
c     k = 1.38066e-23 J/K = 1.38066e-29 km^2*kg/(s^2*K)
      data pre / 1.65089e-2/

      vtherm=pre*t/molmass
      return
      end

c     ***********************************************************
c     Planck function
c     ***********************************************************
      real function planck(ny,tb)
      implicit none
      real ny,tb
      real cst1,cst2

c     Prefactors for the Planck function
c     The intensity units are divided by h/(4pi)
c     Cst1=(8*pi/c^2)^(1/3)
c     Cst2=h/k
c     c in km/s
      data  cst1, cst2 / 3.03529e-7, 4.79922e-11 /

c     For machines that cannot deal with floating overflows
c     planck = (cst1*ny)**3/(exp(min(cst2*ny/tb,85.0))-1.0)
        planck = (cst1*ny)**3/(exp(cst2*ny/tb)-1.0)
      return
      end

c     ***********************************************************
c     Subroutine to compute the energy density of the background
c     radiation for all transitions
c     ***********************************************************
      subroutine uback(uext,tbb,lin)
      implicit none
      integer lin

      include 'fsizes.inc'
      real tbb,uext(ntra)
      real cfreq,planck
      integer l

      do 100 l=1,lin
        uext(l) = planck(cfreq(l),tbb)
100   continue
      return
      end

c     ************************************************************
c     Compute the opacity and source function of an HII region 
c     from electron density and temperature for all transitions
c     ************************************************************
      subroutine hiiinit(ed,et,lin)
      implicit none
      integer lin

      include 'fsizes.inc'
      real kp(ntra),s(ntra)
      real ed,et,cst1,cst2,planck,cfreq,nu,tr
      integer l

      common /hiidata/ kp,s
      save /hiidata/

c     Constants for the absorptivity and the gaunt factor
c     cst1=8/(3*sqrt(2*pi)) * e^6/(4*pi*eps0*me)^3 * (me/k)^(3/2) / c
c     cst2=1/(d*pi) * (2k/(d*me))^(3/2) * 4*pi*eps0/e^2; d=exp(0.577) 
      data cst1,cst2 / 3.014548e16, 4.95734e7 /

      tr = sqrt(et**3)
      do 100 l=1,lin
      nu = cfreq(l)
      kp(l) = cst1*(ed/nu)**2/tr*max(log(cst2*tr/nu),0.0)
      s(l) = planck(nu,et)
100   continue
      return
      end

