!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!> \file atincl.f90
!> Module for atmospheric models

module atincl

!=============================================================================

use mesh
use ppppar
use ppincl

implicit none

!=============================================================================

! 1. Arrays specific to the atmospheric physics

! 1.1 Arrays specific to the input meteo profile
!-----------------------------------------------

!   Arrays specific to values read in the input meteo file:
!                    tmmet ---> time (in sec) of the meteo profile
!                    zdmet ---> altitudes of the dynamic profiles
!                    ztmet ---> altitudes of the temperature profile
!                    umet, vmet, wmet  --> meteo u, v, w profiles
!                    ekmet --->  meteo turbulent kinetic energy profile
!                    epmet ---> meteo turbulent dissipation profile
!                    ttmet --->  meteo temperature (Celsius) profile
!                    qvmet ---> meteo specific humidity profile
!                    ncmet ---> meteo specific humidity profile
!                    pmer  ---> Sea level pressure
!                    xmet, ymet --> cooordinates of the meteo profile

!   Arrays specific to values calculated from the meteo file
!   (cf atlecm.f90):
!                    rmet -->  density profile
!                    tpmet -->  potential temperature profile
!                    phmet -->  hydro. pressure from Laplace integration

double precision, allocatable, dimension(:) :: tmmet, zdmet, ztmet
double precision, allocatable, dimension(:,:) :: umet, vmet, wmet
double precision, allocatable, dimension(:,:) :: ekmet, epmet, ttmet, qvmet, ncmet
double precision, allocatable, dimension(:) :: pmer
double precision, allocatable, dimension(:) :: xmet, ymet
double precision, allocatable, dimension(:,:) :: rmet, tpmet, phmet

double precision, allocatable, dimension(:) :: nebdia, nn

! 1.2 Pointers for the positions of the variables (in rtp, rtpa)
!---------------------------------------------------------------

!   Variables specific to the atmospheric physics:
!   ippmod(iatmos) = 1 (Dry atmosphere):
!   ippmod(iatmos) = 2 (Humid atmosphere):
!                    itotwt---> total water content
!                    intdrp---> total number of droplets

integer, save :: itotwt, intdrp

! 1.3 Pointers for the positions of the properties for the specific phys.
!      (ipproc in propce)
!------------------------------------------------------------------------

!   Properties specific to the atmospheric physics:
!   ippmod(iatmos) = 1 or 2 (Dry or Humid atmosphere):
!                    itempc---> temperature (in celsius)
!   ippmod(iatmos) = 2 (Humid atmosphere):
!                    iliqwt---> liquid water content

integer, save :: itempc, iliqwt

!----------------------------------------------------------------------------

! 2. Data specific to the atmospheric physics

! 2.1 Data specific to the input meteo profile
!----------------------------------------------
!                   imeteo --> flag for reading the meteo input file
!                               = 0 -> no reading
!                               = 1 -> reading
!                   nbmetd --> numbers of altitudes for the dynamics
!                   nbmett --> numbers of altitudes for the temperature
!                                and specific humidity
!                   nbmetm --> numbers of time steps for the meteo profiles
!                   iprofm --> read zone boundary conditions from profile
!                initmeteo --> use meteo profile for variables initialization
!                              (0: not used; 1: used (default))

integer, save :: imeteo, nbmetd, nbmett, nbmetm, iprofm(nozppm)
integer, save :: initmeteo

! 2.1 Constant specific to the physics (defined in atini1.f90)
!-------------------------------------------------------------------------------
!                  ps     --> reference pressure (to comput potential temp: 1.0d+5
!                  rvsra  --> ratio gaz constant h2o/ dry air: 1.608d0
!                  cpvcpa --> ratio Cp h2o/ dry air: 1.866d0
!                  clatev --> latent heat of evaporation: 2.501d+6
!                  gammat --> temperature gradient for the standard atmosphere
!                              ( -6.5d-03 K/m)

double precision, save:: ps, rvsra, cpvcpa, clatev, gammat,rvap

! 2.2. Space and time reference of the run
!-------------------------------------------------------------------------------
!                     syear --> starting year
!                     squant --> starting quantile
!                     shour --> starting hour
!                     smin --> starting min
!                     ssec --> starting second
!                     xlon --> longitude of the domain origin
!                     xlat --> latitude of the domain origin
!                            defined in usati1

integer, save:: syear, squant, shour, smin
double precision, save :: ssec
double precision, save:: xlon, xlat

! 2.3 Data specific to the meteo profile above the domain
!--------------------------------------------------------

integer, save:: nbmaxt, ihpm


! 2.4 Data specific to the 1D vertical grid:
!-------------------------------------------
!                  ivert  --> flag for the definition of the vertical grid
!                              -1: no vertical grid (default)
!                              0 : automatic definition
!                              1 : definition by the user in usdefv
!                   nvert  --> number of vertical arrays
!                   kvert  --> number of levels (up to the top of the domain)
!                   kmx    --> max number of levels (up to 11000 m if ray1d used)
!                              (automatically computed)

integer, save:: ivert, nvert, kvert, kmx

! 2.5 Data specific to the 1d atmospheric radiative module:
!-------------------------------------------------------------------------------
!           iatra1 -->  flag for the use of the 1d atmo radiative model
!                              = 0 no use (default)
!                               = 1 use
!                    nfatr1 --> 1d radiative model pass frequency
!                   iqv0   --> flag for the standard atmo humidity profile
!                              = 0 : q = 0 (default)
!                              = 1 : q = decreasing exponential
!                   should be defined in usati1

integer, save:: iatra1, nfatr1, iqv0
integer, save:: idrayi, idrayst, igrid

! 2.6 Arrays specific to the 1d atmospheric radiative module
!-------------------------------------------------------------------------------

double precision, allocatable :: xyvert(:,:), zvert(:)
double precision, allocatable :: acinfe(:), dacinfe(:), aco2(:, :)
double precision, allocatable :: daco2(:,:), acsup(:), dacsup(:)
double precision, allocatable :: tauzq(:), tauz(:), zq(:)
double precision, save :: tausup
double precision, allocatable :: zray(:), rayi(:,:),rayst(:,:)

! 3.0 Data specific to the ground model
!-------------------------------------------------------------------------------
! iatsoil  --> flag to use the ground model

integer, save:: iatsoil
double precision, save:: w1ini,w2ini

!  -------------------------------------------------------------------------------
! 4.0 Microphysics parameterization options
!  -------------------------------------------------------------------------------
!   --> Option for subgrid models
!   modsub = 0 : the simplest parameterization (for numerical verifications)
!   modsub = 1 : Bechtold et al. 1995 (Luc Musson-Genon)
!   modsub = 2 : Bouzereau et al. 2004
!   modsub = 3 : Cuijpers and Duynkerke 1993, Deardorff 1976, Sommeria and
!                Deardorff 1977
!
!   --> Option for liquid water content distribution models
!   moddis = 1 : all or nothing
!   moddis = 2 : Gaussian distribution
!
!   modnuc = 0 : without nucleation
!   modnuc = 1 : Pruppacher and Klett 1997
!   modnuc = 2 : Cohard et al. 1998,1999
!   modnuc = 3 : Abdul-Razzak et al. 1998,2000 NOT IMPLEMENTED YET
!
!  logaritmic standard deviation of the log-normal law of the droplet spectrum
!  adimensional  : sigc
!
!
! sedimentation flag  : modsedi
!
! logaritmic standard deviation of the log-normal law of the droplet spectrum
! adimensional :  sigc=0.53 other referenced values are 0.28, 0.15

integer, save::  modsub,moddis,modnuc,modsedi
double precision, save:: sigc


contains

  !=============================================================================

subroutine init_meteo

use atsoil

implicit none

integer :: imode

! Allocate additional arrays for Water Microphysics

if (ippmod(iatmos).ge.2) then
  allocate(nebdia(ncelet))
  allocate(nn(ncelet))
endif

if (imeteo.gt.0) then

  imode = 0

  call atlecm &
       !==========
      ( imode )

  ! NB : only ztmet,ttmet,qvmet,ncmet are extended to 11000m if iatr1=1
  !           rmet,tpmet,phmet
  allocate(tmmet(nbmetm), zdmet(nbmetd), ztmet(nbmaxt))
  allocate(umet(nbmetd,nbmetm), vmet(nbmetd,nbmetm), wmet(nbmetd,nbmetm))
  allocate(ekmet(nbmetd,nbmetm), epmet(nbmetd,nbmetm))
  allocate(ttmet(nbmaxt,nbmetm), qvmet(nbmaxt,nbmetm), ncmet(nbmaxt,nbmetm))
  allocate(pmer(nbmetm))
  allocate(xmet(nbmetm), ymet(nbmetm))
  allocate(rmet(nbmaxt,nbmetm), tpmet(nbmaxt,nbmetm), phmet(nbmaxt,nbmetm))

  ! Allocate additional arrays for 1D radiative model

  if (iatra1.eq.1) then

    imode = 0

    call usatdv &
        !==========
        ( imode )

    allocate(xyvert(nvert,3), zvert(kmx))
    allocate(acinfe(kmx), dacinfe(kmx), aco2(kmx,kmx))
    allocate(daco2(kmx,kmx), acsup(kmx), dacsup(kmx))
    allocate(tauzq(kmx+1), tauz(kmx+1), zq(kmx+1))
    allocate(rayi(kmx,nvert), rayst(kmx,nvert), zray(kmx))

    allocate(soilvert(nvert))

    call mestcr &
         !============
         ("rayi", len("rayi"), 1, 0, idrayi)
    call mestcr &
         !============
         ("rayst", len("rayst"), 1, 0, idrayst)
    call gridcr &
         !============
         ("int_grid", len("int_grid"), igrid)

  endif

endif

end subroutine init_meteo

!=============================================================================

subroutine finalize_meteo

use atsoil

implicit none

if (ippmod(iatmos).ge.2) then
  deallocate(nebdia)
  deallocate(nn)
endif

if (imeteo.gt.0) then

  deallocate(tmmet, zdmet, ztmet)
  deallocate(umet, vmet, wmet)
  deallocate(ekmet, epmet)
  deallocate(ttmet, qvmet, ncmet)
  deallocate(pmer)
  deallocate(xmet, ymet)
  deallocate(rmet, tpmet, phmet)

  if (iatra1.eq.1) then

    deallocate(xyvert, zvert)
    deallocate(acinfe, dacinfe, aco2)
    deallocate(daco2, acsup, dacsup)
    deallocate(tauzq, tauz, zq)
    deallocate(rayi, rayst, zray)

    deallocate(soilvert)

    call mestde ()
    !============

    call grides ()
    !============

  endif

endif

end subroutine finalize_meteo

end module atincl
