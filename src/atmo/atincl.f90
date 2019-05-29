!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
!> \brief Module for atmospheric models - main variables
!>-   Nota : ippmod(iatmos) = 0 constante density, =1 --> Dry atmosphere,
!>          = 2 --> Humid atmosphere
!>-     A separate vertical grid is used for 1D radiative scheme
module atincl
!> \defgroup at_main
!=============================================================================

use mesh
use ppppar
use ppincl

implicit none
!> \addtogroup at_main
!> \{

!=============================================================================

! 1. Arrays specific to the atmospheric physics

! 1.1 Arrays specific to the input meteo profile
!-----------------------------------------------

!   Arrays specific to values read in the input meteo file:

!> time (in sec) of the meteo profile
double precision, allocatable, dimension(:) :: tmmet

!> altitudes of the dynamic profiles (read in the input meteo file)
double precision, allocatable, dimension(:) :: zdmet

!> altitudes of the temperature profile (read in the input meteo file)
double precision, allocatable, dimension(:) :: ztmet

!> meteo u  profiles (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: umet

!> meteo  v profiles (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: vmet

!> meteo w profiles - unused
double precision, allocatable, dimension(:,:) :: wmet

!> meteo turbulent kinetic energy profile (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: ekmet

!> meteo turbulent dissipation profile (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: epmet

!> meteo temperature (Celsius) profile (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: ttmet

!> meteo specific humidity profile (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: qvmet

!> meteo specific drplet number profile (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: ncmet

!> Sea level pressure (read in the input meteo file)
double precision, allocatable, dimension(:) :: pmer

!> X axis cooordinates of the meteo profile (read in the input meteo file)
double precision, allocatable, dimension(:) :: xmet

!> Y axis cooordinates of the meteo profile (read in the input meteo file)
double precision, allocatable, dimension(:) :: ymet

!   Arrays specific to values calculated from the meteo file (cf atlecm.f90):

!> density profile
double precision, allocatable, dimension(:,:) :: rmet

!> potential temperature profile
double precision, allocatable, dimension(:,:) :: tpmet

!> hydrostatic pressure from Laplace integration
double precision, allocatable, dimension(:,:) :: phmet

!> Diagnosed nebulosity
double precision, allocatable, dimension(:) :: nebdia

!> fractional nebulosity
double precision, allocatable, dimension(:) :: nn

! 1.2 Pointers for the positions of the variables
!------------------------------------------------
!   Variables specific to the atmospheric physics:
!> itotwt---> total water content (for humid atmosphere)
integer, save :: itotwt
!> intdrp---> total number of droplets (for humid atmosphere)
integer, save :: intdrp

! 1.3 Pointers for the positions of the properties for the specific phys.
!------------------------------------------------------------------------
!   Properties specific to the atmospheric physics:

!> temperature (in Celsius)
integer, save :: itempc

!> liquid water content
integer, save :: iliqwt

!> momentum source term field id (useful when iatmst > 0)
integer, save :: imomst

!----------------------------------------------------------------------------

! 2. Data specific to the atmospheric physics

! 2.1 Data specific to the input meteo profile
!----------------------------------------------

!>flag for reading the meteo input file
!> - = 0 -> no reading
!> - = 1 -> reading
integer, save :: imeteo

!> numbers of altitudes for the dynamics
integer, save :: nbmetd

!> numbers of altitudes for the temperature and specific humidity
integer, save :: nbmett

!> numbers of time steps for the meteo profiles
integer, save :: nbmetm

!> read zone boundary conditions from profile
integer, save :: iprofm(nozppm)

!> automatic inlet/outlet boundary condition flag
!> (0: not auto (default); 1,2: auto)
!> When meteo momentum source terms are activated (iatmst > 0),
!> iautom = 1 corresponds to a Dirichlet on the pressure and a
!> Neumann on the velocity, whereas iautom = 2 imposes a Dirichlet
!> on both pressure and velocity
integer, allocatable, dimension(:) :: iautom

!> use meteo profile for variables initialization
!> (0: not used; 1: used (default))
integer, save :: initmeteo

!> add a momentum source term based on the meteo profile
integer, save :: iatmst

!> flag for meteo velocity field interpolation
!> - 0: linear interpolation of the meteo profile
!> - 1: the user can directly impose the exact meteo velocity
!> by declaring the 'meteo_velocity' field
!> Useful for iatmst = 1
integer, save :: theo_interp

! 2.1 Constant specific to the physics (defined in atini1.f90)
!-------------------------------------------------------------------------------

!> reference pressure (to compute potential temp: 1.0d+5)
double precision, save:: ps

!> ratio gaz constant h2o/ dry air: 1.608d0
double precision, save:: rvsra

!> ratio Cp h2o/ dry air: 1.866d0
double precision, save:: cpvcpa

!> latent heat of evaporation: 2.501d+6
double precision, save:: clatev

!> temperature gradient for the standard atmosphere (-6.5d-03 K/m)
double precision, save:: gammat

!> rvsra*rair
double precision, save:: rvap

! 2.2. Space and time reference of the run
!-------------------------------------------------------------------------------

!> starting year
integer, save:: syear

!> starting quantile
integer, save:: squant

!> starting hour
integer, save:: shour

!> starting min
integer, save:: smin

!> starting second
double precision, save :: ssec

!> longitude of the domain origin
double precision, save:: xlon

!> latitude of the domain origin
double precision, save:: xlat

! 2.3 Data specific to the meteo profile above the domain
!--------------------------------------------------------
!> Number of vertical levels (cf. 1D radiative scheme
integer, save:: nbmaxt

!> flag to compute the hydrostastic pressure by Laplace integration
!> in the meteo profiles
!> = 0 : bottom to top Laplace integration, based on P(sea level) (default)
!> = 1 : top to bottom Laplace integration based on P computed for
!>            the standard atmosphere at z(nbmaxt)
integer, save:: ihpm


! 2.4 Data specific to the 1D vertical grid:
!-------------------------------------------

!> flag for the definition of the vertical grid
!> - -1 : no vertical grid (default)
!> -  0 : automatic definition
!> -  1 : definition by the user in cs_user_atmospheric_model.f90
integer, save:: ivert

!> number of vertical arrays
integer, save:: nvert

!> number of levels (up to the top of the domain)
integer, save:: kvert

!> Number of levels (up to 11000 m if ray1d used)
!> (automatically computed)
integer, save:: kmx

! 2.5 Data specific to the 1d atmospheric radiative module:
!-------------------------------------------------------------------------------
!> flag for the use of the 1d atmo radiative model
!> - 0 no use (default)
!> - 1 use
integer, save:: iatra1

!> 1d radiative model pass frequency
integer, save:: nfatr1

!> flag for the standard atmo humidity profile
!> - 0: q = 0 (default)
!> - 1: q = decreasing exponential
integer, save:: iqv0

!> pointer for 1D infrared profile
integer, save:: idrayi

!> pointer for 1D solar profile
integer, save:: idrayst

!> grid formed by 1D profiles
integer, save:: igrid

!> Flag for the computation of downward and upward infrared radiative fluxes
!> 0: disabled
!> 1: enabled
integer, save:: irdu

!> Flag for storage of downward and upward solar radiative fluxes
!> 0: disabled
!> 1: enabled
integer, save:: soldu

! 2.6 Arrays specific to the 1d atmospheric radiative module
!-------------------------------------------------------------------------------

!> horizontal coordinates of the vertical grid
double precision, allocatable, dimension(:,:) :: xyvert

!> vertical grid for 1D radiative scheme initialize in
!>       cs_user_atmospheric_model.f90
double precision, allocatable, dimension(:) :: zvert

!> absorption for CO2 + 03
double precision, allocatable, dimension(:) :: acinfe

!> differential absorption for CO2 + 03
double precision, allocatable, dimension(:) :: dacinfe

!> absorption for CO2 only
double precision, allocatable, dimension(:,:) :: aco2, aco2s

!> differential absorption for CO2 only
double precision, allocatable, dimension(:,:) :: daco2, daco2s

!> idem acinfe, flux descendant
double precision, allocatable, dimension(:) :: acsup, acsups

!> internal variable for 1D radiative model
double precision, allocatable, dimension(:) :: dacsup, dacsups

!> internal variable for 1D radiative model
double precision, allocatable, dimension(:) :: tauzq

!> internal variable for 1D radiative model
double precision, allocatable, dimension(:) :: tauz

!> internal variable for 1D radiative model
double precision, allocatable, dimension(:) :: zq

!> internal variable for 1D radiative model
double precision, save :: tausup

!> internal variable for 1D radiative model
double precision, allocatable, dimension(:) :: zray
double precision, allocatable, dimension(:,:) :: rayi, rayst

!> Upward and downward radiative fluxes (infrared, solar) along each vertical
double precision, allocatable, dimension(:,:) :: iru, ird, solu, sold

! 3.0 Data specific to the ground model
!-------------------------------------------------------------------------------
!> iatsoil  --> flag to use the ground model
integer, save:: iatsoil

!> Water content of the first ground reservoir
double precision, save:: w1ini

!> Water content of the second ground reservoir
double precision, save:: w2ini

!> Do we compute z ground every where?
logical(c_bool), pointer, save :: compute_z_ground

!  -------------------------------------------------------------------------------
! 4.0 Microphysics parameterization options
!  -------------------------------------------------------------------------------

!> Option for subgrid models
!>  - modsub = 0 : the simplest parameterization (for numerical verifications)
!>  - modsub = 1 : Bechtold et al. 1995 (Luc Musson-Genon)
!>  - modsub = 2 : Bouzereau et al. 2004
!>  - modsub = 3 : Cuijpers and Duynkerke 1993, Deardorff 1976, Sommeria and
!>                Deardorff 1977
integer, save::  modsub

!> Option for liquid water content distribution models
!>  - moddis = 1 : all or nothing
!>  - moddis = 2 : Gaussian distribution
integer, save::  moddis

!> Option for nucleation
!>  - modnuc = 0 : without nucleation
!>  - modnuc = 1 : Pruppacher and Klett 1997
!>  - modnuc = 2 : Cohard et al. 1998,1999
!>  - modnuc = 3 : Abdul-Razzak et al. 1998,2000
!>  logaritmic standard deviation of the log-normal law of the droplet spectrum
integer, save::  modnuc

!> sedimentation flag
integer, save::  modsedi

!> deposition flag
integer, save:: moddep

!> adimensional :  sigc=0.53 other referenced values are 0.28, 0.15
double precision, save:: sigc

!> force initilization in case of restart (this option is
!> automatically set in lecamp)
integer, save :: init_at_chem

!> key id for optimal interpolation
integer, save :: kopint

!> \}

contains

  !=============================================================================

  !> \brief Map fortran to C variables

  subroutine atmo_init

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    implicit none

    ! Local variables
    type(c_ptr) :: c_compute_z_ground

    call cs_f_atmo_get_pointers(c_compute_z_ground)

    call c_f_pointer(c_compute_z_ground, compute_z_ground)

    return

  end subroutine atmo_init

  !=============================================================================
!> \brief Initialisation of meteo data
subroutine init_meteo

use atsoil

implicit none

integer :: imode, ifac

! Allocate additional arrays for Water Microphysics

if (ippmod(iatmos).ge.2) then
  allocate(nebdia(ncelet))
  allocate(nn(ncelet))
endif

if (imeteo.gt.0) then

  imode = 0

  call atlecm ( imode )

  ! NB : only ztmet,ttmet,qvmet,ncmet are extended to 11000m if iatr1=1
  !           rmet,tpmet,phmet
  allocate(tmmet(nbmetm), zdmet(nbmetd), ztmet(nbmaxt))
  allocate(umet(nbmetd,nbmetm), vmet(nbmetd,nbmetm), wmet(nbmetd,nbmetm))
  allocate(ekmet(nbmetd,nbmetm), epmet(nbmetd,nbmetm))
  allocate(ttmet(nbmaxt,nbmetm), qvmet(nbmaxt,nbmetm), ncmet(nbmaxt,nbmetm))
  allocate(pmer(nbmetm))
  allocate(xmet(nbmetm), ymet(nbmetm))
  allocate(rmet(nbmaxt,nbmetm), tpmet(nbmaxt,nbmetm), phmet(nbmaxt,nbmetm))

  ! Allocate and initialize auto inlet/outlet flag

  allocate(iautom(nfabor))
  do ifac = 1, nfabor
    iautom(ifac) = 0
  enddo

  ! Allocate additional arrays for 1D radiative model

  if (iatra1.eq.1) then

    imode = 0

    call usatdv ( imode )

    allocate(xyvert(nvert,3), zvert(kmx))
    allocate(acinfe(kmx), dacinfe(kmx), aco2(kmx,kmx), aco2s(kmx,kmx))
    allocate(daco2(kmx,kmx), daco2s(kmx,kmx), acsup(kmx), dacsup(kmx))
    allocate(acsups(kmx), dacsups(kmx))
    allocate(tauzq(kmx+1), tauz(kmx+1), zq(kmx+1))
    allocate(rayi(kmx,nvert), rayst(kmx,nvert), zray(kmx))

    allocate(soilvert(nvert))

    call mestcr  ("rayi",  len("rayi"), 1, 0, idrayi)
    call mestcr  ("rayst", len("rayst"), 1, 0, idrayst)
    call gridcr  ("int_grid", len("int_grid"), igrid)

    if (irdu.eq.1) then
      allocate(iru(kmx,nvert), ird(kmx,nvert))
    endif

    if (soldu.eq.1) then
      allocate(sold(kmx,nvert), solu(kmx,nvert))
    endif

  endif

endif

end subroutine init_meteo

!=============================================================================
!> \brief Final step for deallocation
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

  deallocate(iautom)

  if (iatra1.eq.1) then

    deallocate(xyvert, zvert)
    deallocate(acinfe, dacinfe, aco2, aco2s)
    deallocate(daco2, daco2s, acsup, dacsup)
    deallocate(acsups, dacsups)
    deallocate(tauzq, tauz, zq)
    deallocate(rayi, rayst, zray)

    deallocate(soilvert)

    call mestde ()

    call grides ()

    if (irdu.eq.1) then
      deallocate(iru, ird)
    endif

    if (soldu.eq.1) then
      deallocate(solu, sold)
    endif

  endif

endif

end subroutine finalize_meteo

end module atincl
