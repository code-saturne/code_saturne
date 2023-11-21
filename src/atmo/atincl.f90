!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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
double precision, dimension(:), pointer :: tmmet

!> altitudes of the dynamic profiles (read in the input meteo file)
double precision, dimension(:), pointer :: zdmet

!> Pressure drop integrated over a time step (used for automatic open boundaries)
double precision, allocatable, dimension(:) :: dpdt_met

!> Momentum for each level (used for automatic open boundaries)
double precision, allocatable, dimension(:,:) :: mom_met
double precision, allocatable, dimension(:,:) :: mom

!> altitudes of the temperature profile (read in the input meteo file)
double precision, dimension(:), pointer :: ztmet

!> meteo u  profiles (read in the input meteo file)
double precision, dimension(:,:), pointer :: umet

!> meteo  v profiles (read in the input meteo file)
double precision, dimension(:,:), pointer :: vmet

!> meteo w profiles - unused
double precision, dimension(:,:), pointer :: wmet

!> meteo turbulent kinetic energy profile (read in the input meteo file)
double precision, dimension(:,:), pointer :: ekmet

!> meteo turbulent dissipation profile (read in the input meteo file)
double precision, dimension(:,:), pointer :: epmet

!> meteo temperature (Celsius) profile (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: ttmet

!> meteo specific humidity profile (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: qvmet

!> meteo specific droplet number profile (read in the input meteo file)
double precision, allocatable, dimension(:,:) :: ncmet

!> Sea level pressure (read in the input meteo file)
double precision, allocatable, dimension(:) :: pmer

!> X axis coordinates of the meteo profile (read in the input meteo file)
double precision, allocatable, dimension(:) :: xmet

!> Y axis coordinates of the meteo profile (read in the input meteo file)
double precision, allocatable, dimension(:) :: ymet

! Arrays specific to values calculated from the meteo file (cf atlecm.f90):

!> density profile
double precision, allocatable, dimension(:,:) :: rmet

!> potential temperature profile
double precision, dimension(:,:), pointer :: tpmet

!> hydrostatic pressure from Laplace integration
double precision, dimension(:,:), pointer :: phmet

! 1.2 Pointers for the positions of the variables
!------------------------------------------------
!   Variables specific to the atmospheric physics:
!> total water content (for humid atmosphere)
integer, save :: iymw
!> intdrp---> total number of droplets (for humid atmosphere)
integer, save :: intdrp = -1

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
integer(c_int), pointer, save :: imeteo

!> numbers of altitudes for the dynamics
integer(c_int), pointer, save :: nbmetd

!> numbers of altitudes for the temperature and specific humidity
integer(c_int), pointer, save :: nbmett

!> numbers of time steps for the meteo profiles
integer(c_int), pointer, save :: nbmetm

!> read zone boundary conditions from profile
integer(c_int), pointer, save :: iprofm(:)

!> automatic inlet/outlet boundary condition flag
!> (0: not auto (default); 1,2: auto)
!> When meteo momentum source terms are activated (iatmst > 0),
!> iautom = 1 corresponds to a Dirichlet on the pressure and a
!> Neumann on the velocity, whereas iautom = 2 imposes a Dirichlet
!> on both pressure and velocity
integer(c_int), pointer, save :: iautom(:)

!> use meteo profile for variables initialization
!> (0: not used; 1: used (default))
integer, save :: initmeteo

!> add a momentum source term based on the meteo profile
!> for automatic open boundaries
integer(c_int), pointer, save :: iatmst

!> flag for meteo velocity field interpolation
!> - 0: linear interpolation of the meteo profile
!> - 1: the user can directly impose the exact meteo velocity
!> by declaring the 'meteo_velocity' field
!> Useful for iatmst = 1
!> Note: deprecated, imeteo=2 can be used instead.
integer(c_int), pointer, save :: theo_interp

! 2.1 Constant specific to the physics (defined in atini1.f90)
!-------------------------------------------------------------------------------

!> reference pressure (to compute potential temp: 1.0d+5)
real(c_double), pointer, save:: ps

! 2.2. Space and time reference of the run
!-------------------------------------------------------------------------------

!> starting year
integer(c_int), pointer, save:: syear

!> starting quantile
integer(c_int), pointer, save:: squant

!> starting hour
integer(c_int), pointer, save:: shour

!> starting min
integer(c_int), pointer, save:: smin

!> starting second
real(c_double), pointer, save :: ssec

!> longitude of the domain origin
real(c_double), pointer, save:: xlon

!> latitude of the domain origin
real(c_double), pointer, save:: xlat

!> x coordinate of the domain origin in Lambert-93
real(c_double), pointer, save:: xl93

!> y coordinate of the domain origin in Lambert-93
real(c_double), pointer, save:: yl93

! 2.3 Data specific to the meteo profile above the domain
!--------------------------------------------------------
!> Number of vertical levels (cf. 1-D radiative scheme)
integer(c_int), pointer, save:: nbmaxt

!> flag to compute the hydrostatic pressure by Laplace integration
!> in the meteo profiles
!> = 0 : bottom to top Laplace integration, based on P(sea level) (default)
!> = 1 : top to bottom Laplace integration based on P computed for
!>            the standard atmosphere at z(nbmaxt)
integer, save:: ihpm

! 2.4 Data specific to the 1-D vertical grid:
!-------------------------------------------

!> number of vertical arrays
integer(c_int), pointer, save:: nvert

!> number of levels (up to the top of the domain)
integer(c_int), pointer, save:: kvert

!> Number of levels (up to 11000 m if 1-D radiative transfer used)
!> (automatically computed)
integer(c_int), pointer, save:: kmx

!> Height of the boundary layer
real(c_double), pointer, save :: meteo_zi

! 2.5 Data specific to the 1-D atmospheric radiative module:
!-------------------------------------------------------------------------------
!> flag for the use of the 1-D atmo radiative model
!> - 0 no use (default)
!> - 1 use
integer(c_int), pointer, save :: iatra1

!> 1D radiative model pass frequency
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


! 2.6 Arrays specific to the 1D atmospheric radiative module
!-------------------------------------------------------------------------------

!> horizontal coordinates of the vertical grid
double precision, dimension(:,:), pointer :: xyvert

!> vertical grid for 1D radiative scheme initialize in
!>       cs_user_atmospheric_model.f90
double precision, dimension(:), pointer :: zvert

!> absorption for CO2 + 03
double precision, dimension(:), pointer :: acinfe

!> differential absorption for CO2 + 03
double precision, dimension(:), pointer :: dacinfe

!> absorption for CO2 only
double precision, dimension(:,:), pointer :: aco2, aco2s

!> differential absorption for CO2 only
double precision, dimension(:,:), pointer :: daco2, daco2s

!> idem acinfe, flux descendant
double precision, dimension(:), pointer :: acsup, acsups

!> internal variable for 1D radiative model
double precision, dimension(:), pointer :: dacsup, dacsups

!> internal variable for 1D radiative model
double precision, dimension(:), pointer :: tauzq

!> internal variable for 1D radiative model
double precision, dimension(:), pointer :: tauz

!> internal variable for 1D radiative model
double precision, dimension(:), pointer :: zq

!> internal variable for 1D radiative model
double precision, save :: tausup

!> internal variable for 1D radiative model
double precision, dimension(:), pointer  :: zray
double precision, dimension(:,:), pointer  :: rayi, rayst

!> Upward and downward radiative fluxes (infrared, solar) along each vertical
double precision, dimension(:,:), pointer  :: iru, ird, solu, sold

! 3.0 Data specific to the soil model
!-------------------------------------------------------------------------------
!> Option for soil model
!>  - iatsoil = 0 : Do not use soil model (soil structures are not created)
!>  - iatsoil = 1 : Deardorff Force Restore model (Deardorff 1973)
!>  - iatsoil = 2 : User inputs (soil structures are created)

integer(c_int), pointer, save:: iatsoil

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
integer(c_int), pointer, save:: modsub

!> Option for liquid water content distribution models
!>  - moddis = 1 : all or nothing
!>  - moddis = 2 : Gaussian distribution
integer(c_int), pointer, save:: moddis

!> Option for nucleation
!>  - modnuc = 0 : without nucleation
!>  - modnuc = 1 : Pruppacher and Klett 1997
!>  - modnuc = 2 : Cohard et al. 1998,1999
!>  - modnuc = 3 : Abdul-Razzak et al. 1998,2000
!>  logaritmic standard deviation of the log-normal law of the droplet spectrum
integer(c_int), pointer, save:: modnuc

!> sedimentation flag
integer(c_int), pointer, save:: modsedi

!> deposition flag
integer(c_int), pointer, save:: moddep

!> adimensional :  sigc=0.53 other referenced values are 0.28, 0.15
double precision, save:: sigc

!> force initilization in case of restart (this option is
!> automatically set in lecamp)
integer, save :: init_at_chem

!> key id for optimal interpolation
integer, save :: kopint

!> Aerosol optical properties

! Aerosol optical depth
!> adimensional :  aod_o3_tot=0.2 other referenced values are  0.10, 0.16
double precision, save:: aod_o3_tot
!> adimensional :  aod_h2o_tot=0.10 other referenced values are  0.06, 0.08
double precision, save:: aod_h2o_tot

!> Asymmetry factor for O3 (non-dimensional)
!> climatic value gaero_o3=0.66
double precision, save:: gaero_o3
!> Asymmetry factor for H2O (non-dimensional)
!> climatic value gaero_h2o=0.64
double precision, save:: gaero_h2o

!> Single scattering albedo for O3 (non-dimensional)
!> climatic value piaero_o3=0.84, other referenced values are 0.963
double precision, save:: piaero_o3
!> Single scattering albedo for H2O (non-dimensional)
!> climatic value piaero_h2o=0.84, other referenced values are 0.964
double precision, save:: piaero_h2o

!> Fraction of Black carbon (non-dimensional): black_carbon_frac=1.d-8 for no BC
double precision, save:: black_carbon_frac

!> Maximal height for aerosol distribution on the vertical
!> important should be <= zqq(kmray-1);
!> in meters : referenced value: zaero=6000
double precision, save:: zaero
!> \}

!=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function returning a meteo file name

    subroutine cs_f_atmo_get_meteo_file_name(f_name_max, f_name, f_name_len)  &
      bind(C, name='cs_f_atmo_get_meteo_file_name')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value       :: f_name_max
      type(c_ptr), intent(out)    :: f_name
      integer(c_int), intent(out) :: f_name_len
    end subroutine cs_f_atmo_get_meteo_file_name

    !---------------------------------------------------------------------------

    !> \brief Return pointers to atmo includes

    subroutine cs_f_atmo_get_pointers(ps,                               &
        syear, squant, shour, smin, ssec,                               &
        longitude, latitude,                                            &
        x_l93, y_l93,                                                   &
        compute_z_ground, iatmst, theo_interp,                          &
        sedimentation_model, deposition_model, nucleation_model,        &
        subgrid_model, distribution_model,                              &
        ichemistry, isepchemistry, nespg, nrg, chem_with_photo,         &
        iaerosol, frozen_gas_chem, init_gas_with_lib,                   &
        init_aero_with_lib, n_aero, n_sizebin, imeteo,                  &
        nbmetd, nbmett, nbmetm, iatra1, nbmaxt,                         &
        meteo_zi, iatsoil,                                              &
        nvertv, kvert, kmx )                                            &
      bind(C, name='cs_f_atmo_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ps
      type(c_ptr), intent(out) :: compute_z_ground, iatmst, theo_interp
      type(c_ptr), intent(out) :: ichemistry, isepchemistry, nespg, nrg
      type(c_ptr), intent(out) :: sedimentation_model, deposition_model
      type(c_ptr), intent(out) :: nucleation_model
      type(c_ptr), intent(out) :: subgrid_model, distribution_model
      type(c_ptr), intent(out) :: syear, squant, shour, smin, ssec
      type(c_ptr), intent(out) :: longitude, latitude
      type(c_ptr), intent(out) :: x_l93, y_l93
      type(c_ptr), intent(out) :: iaerosol, frozen_gas_chem
      type(c_ptr), intent(out) :: init_gas_with_lib, init_aero_with_lib
      type(c_ptr), intent(out) :: n_aero, n_sizebin, chem_with_photo
      type(c_ptr), intent(out) :: imeteo
      type(c_ptr), intent(out) :: nbmetd, nbmett, nbmetm, iatra1, nbmaxt
      type(c_ptr), intent(out) :: meteo_zi
      type(c_ptr), intent(out) :: iatsoil
      type(c_ptr), intent(out) :: nvertv, kvert, kmx
    end subroutine cs_f_atmo_get_pointers

    !---------------------------------------------------------------------------

    !> \brief Return pointers to atmo arrays

    subroutine cs_f_atmo_arrays_get_pointers(p_zdmet, p_ztmet, p_umet, p_vmet, &
         p_wmet  , p_tmmet, p_phmet, p_tpmet, p_ekmet, p_epmet,                &
         p_xyvert, p_zvert, p_acinfe,                                          &
         p_dacinfe, p_aco2, p_aco2s,                                           &
         p_daco2, p_daco2s,                                                    &
         p_acsup, p_acsups,                                                    &
         p_dacsup, p_dacsups,                                                  &
         p_tauzq, p_tauz, p_zq,                                                &
         p_zray, p_rayi, p_rayst,                                              &
         p_iru, p_ird, p_solu, p_sold,                                         &
         dim_pumet, dim_phmet,                                                 &
         dim_tpmet, dim_ekmet, dim_epmet,                                      &
         dim_xyvert, dim_kmx2, dim_kmx_nvert )                                 &
         bind(C, name='cs_f_atmo_arrays_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), dimension(2) :: dim_phmet, dim_pumet, dim_tpmet
      integer(c_int), dimension(2) ::  dim_ekmet,  dim_epmet
      integer(c_int), dimension(2) ::  dim_xyvert, dim_kmx2, dim_kmx_nvert
      type(c_ptr), intent(out) :: p_zdmet, p_ztmet, p_umet, p_vmet, p_tmmet
      type(c_ptr), intent(out) :: p_wmet
      type(c_ptr), intent(out) :: p_phmet, p_tpmet, p_ekmet, p_epmet
      type(c_ptr), intent(out) :: p_xyvert, p_zvert, p_acinfe
      type(c_ptr), intent(out) :: p_dacinfe, p_aco2, p_aco2s
      type(c_ptr), intent(out) :: p_daco2, p_daco2s
      type(c_ptr), intent(out) :: p_acsup, p_acsups
      type(c_ptr), intent(out) :: p_dacsup, p_dacsups
      type(c_ptr), intent(out) :: p_tauzq, p_tauz, p_zq
      type(c_ptr), intent(out) :: p_zray, p_rayi, p_rayst
      type(c_ptr), intent(out) :: p_iru, p_ird, p_solu, p_sold
    end subroutine cs_f_atmo_arrays_get_pointers

    !---------------------------------------------------------------------------

    !> \brief Initialize meteo profiles if no meteo file is given

    subroutine cs_atmo_init_meteo_profiles() &
        bind(C, name='cs_atmo_init_meteo_profiles')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_atmo_init_meteo_profiles

    !---------------------------------------------------------------------------

    !> \brief Compute meteo profiles if no meteo file is given

    subroutine cs_atmo_compute_meteo_profiles() &
        bind(C, name='cs_atmo_compute_meteo_profiles')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_atmo_compute_meteo_profiles

    !---------------------------------------------------------------------------

    !> \brief Calculation of the specific enthalpy of liquid water

    !> \return specific enthalpy of liquid water

    !> \param[in]  t_l  liquid water temperature (in Celsius)

    function cs_liq_t_to_h(t_l) result(h_l) &
        bind(C, name='cs_liq_t_to_h')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: t_l
      real(c_double) :: h_l
    end function cs_liq_t_to_h

    !---------------------------------------------------------------------------

    !> \brief Calculation of the absolute humidity at saturation
    !>        for a given temperature.

    !> \param[in]  t_c  temperature (in Celsius)
    !> \param[in]  p    pressure

    function cs_air_x_sat(t_c, p) result(x_s) &
        bind(C, name='cs_air_x_sat')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: t_c, p
      real(c_double) :: x_s
    end function cs_air_x_sat

    !---------------------------------------------------------------------------

    !> \brief Calculation of the air water mass fraction at saturation
    !>        for a given temperature.

    !> \param[in]  t_c  temperature (in Celsius)
    !> \param[in]  p    pressure

    function cs_air_yw_sat(t_c, p) result(x_s) &
        bind(C, name='cs_air_yw_sat')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: t_c, p
      real(c_double) :: x_s
    end function cs_air_yw_sat

    !---------------------------------------------------------------------------

    !> \brief Computes the saturation water vapor pressure function
    !> of the temperature (C).

    !> \param[in]  t_c  temperature (in Celsius)

    function cs_air_pwv_sat(t_c) result(x_s) &
        bind(C, name='cs_air_pwv_sat')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: t_c
      real(c_double) :: x_s
    end function cs_air_pwv_sat

    !---------------------------------------------------------------------------

    !> \brief Convert the absolute humidity of humid air to the
    !>        air water mass fraction.

    !> \param[in]  x  absolute humidity of humid air

    function cs_air_x_to_yw(x) result(qw) &
        bind(C, name='cs_air_x_to_yw')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: x
      real(c_double) :: qw
    end function cs_air_x_to_yw

    !---------------------------------------------------------------------------

    !> \brief Convert the air water mass fraction to the
    !>        absolute humidity of humid air.

    !> \param[in]  qw  air water mass fraction

    function cs_air_yw_to_x(qw) result(x) &
        bind(C, name='cs_air_yw_to_x')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: qw
      real(c_double) :: x
    end function cs_air_yw_to_x

    !---------------------------------------------------------------------------

    !> \brief Calculation of the density of humid air.

    !> \param[in]  ywm           air water mass fraction
    !> \param[in]  t_liq         liquid temperature
    !> \param[in]  p             pressure
    !> \param[out] yw_liq        liquid water mass fraction
    !> \param[out] t_h           temperature of humid air in Celsius
    !> \param[out] rho_h         density of humid air

    subroutine cs_rho_humidair(ywm, t_liq, p, yw_liq, t_h, rho_h) &
        bind(C, name='cs_rho_humidair')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: ywm, t_liq, p
      real(c_double), intent(out) :: yw_liq, t_h, rho_h
    end subroutine cs_rho_humidair

    !=============================================================================

    ! Switch universal functions

    ! Derivative functions

    function cs_mo_phim(z,dlmo) result(coef) &
        bind(C, name='cs_mo_phim')
      use, intrinsic :: iso_c_binding

      implicit none

      real(c_double), value :: z,dlmo
      real(c_double) :: coef

    end function cs_mo_phim

    function cs_mo_phih(z,dlmo) result(coef) &
        bind(C, name='cs_mo_phih')
      use, intrinsic :: iso_c_binding

      implicit none

      real(c_double), value :: z,dlmo
      real(c_double) :: coef

    end function cs_mo_phih

    ! Integrated version from z0 to z

    function cs_mo_psim(z,z0,dlmo) result(coef) &
        bind(C, name='cs_mo_psim')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: z,z0,dlmo
      real(c_double) :: coef

    end function cs_mo_psim

    function cs_mo_psih (z,z0,dlmo) result(coef) &
        bind(C, name='cs_mo_psih')
      use, intrinsic :: iso_c_binding

      implicit none

      real(c_double), value :: z,z0,dlmo
      real(c_double) :: coef

    end function cs_mo_psih

    subroutine cs_f_atmo_get_soil_zone(n_faces, n_soil_cat, face_ids)  &
        bind(C, name='cs_f_atmo_get_soil_zone')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(out) :: n_faces
      integer(c_int), intent(out) :: n_soil_cat
      type(c_ptr), intent(out) :: face_ids
    end subroutine cs_f_atmo_get_soil_zone

    subroutine cs_f_boundary_conditions_get_atincl_pointers(p_iprofm, p_iautom) &
      bind(C, name='cs_f_boundary_conditions_get_atincl_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_iprofm, p_iautom
    end subroutine cs_f_boundary_conditions_get_atincl_pointers

  end interface

contains

  !=============================================================================

  subroutine at_models_bc_map(nfabor)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments
    integer :: nfabor

    ! Local variables
    type(c_ptr) :: p_iprofm, p_iautom

    call cs_f_boundary_conditions_get_atincl_pointers(p_iprofm, p_iautom)
    call c_f_pointer(p_iprofm, iprofm, [nozppm])
    call c_f_pointer(p_iautom, iautom, [nfabor])

  end subroutine at_models_bc_map

  !=============================================================================

  subroutine atmo_get_soil_zone(n_faces, n_soil_cat, face_ids_p)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments
    integer(c_int), intent(out) :: n_faces
    integer(c_int), intent(out) :: n_soil_cat
    integer, dimension(:), pointer, intent(out) ::face_ids_p

    ! Local variables
    type(c_ptr) :: c_p

    call cs_f_atmo_get_soil_zone(n_faces, n_soil_cat, c_p)
    call c_f_pointer(c_p, face_ids_p, [n_faces])

  end subroutine atmo_get_soil_zone

  !=============================================================================

  !> \brief Return meteo file name

  !> \param[out]  name   meteo file name

  subroutine atmo_get_meteo_file_name(name)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(out) :: name

    ! Local variables

    integer :: i
    integer(c_int) :: name_max, c_name_len
    type(c_ptr) :: c_name_p
    character(kind=c_char, len=1), dimension(:), pointer :: c_name

    name_max = len(name)

    call cs_f_atmo_get_meteo_file_name(name_max, c_name_p, c_name_len)
    call c_f_pointer(c_name_p, c_name, [c_name_len])

    do i = 1, c_name_len
      name(i:i) = c_name(i)
    enddo
    do i = c_name_len + 1, name_max
      name(i:i) = ' '
    enddo

    return

  end subroutine atmo_get_meteo_file_name

  !=============================================================================

  !> \brief Map Fortran to C variables

  subroutine atmo_init

    use, intrinsic :: iso_c_binding
    use atchem, only: nrg, nespg, ichemistry, isepchemistry, photolysis
    use sshaerosol
    use cs_c_bindings

    implicit none

    ! Local variables
    type(c_ptr) :: c_ps
    type(c_ptr) :: c_compute_z_ground, c_iatmst, c_model, c_nrg, c_nespg
    type(c_ptr) :: c_sedimentation_model, c_deposition_model, c_nucleation_model
    type(c_ptr) :: c_subgrid_model
    type(c_ptr) :: c_distribution_model
    type(c_ptr) :: c_syear, c_squant, c_shour, c_smin, c_ssec
    type(c_ptr) :: c_longitude, c_latitude
    type(c_ptr) :: c_xl93, c_yl93
    type(c_ptr) :: c_modelaero, c_frozen_gas_chem, c_nlayer, c_nsize
    type(c_ptr) :: c_init_gas_with_lib, c_init_aero_with_lib, c_chem_with_photo
    type(c_ptr) :: c_imeteo
    type(c_ptr) :: c_nbmetd, c_nbmett, c_nbmetm, c_iatra1, c_nbmaxt
    type(c_ptr) :: c_meteo_zi
    type(c_ptr) :: c_iatsoil, c_isepchemistry
    type(c_ptr) :: c_nvert, c_kvert, c_kmx, c_theo_interp

    call cs_f_atmo_get_pointers(c_ps,               &
      c_syear, c_squant, c_shour, c_smin, c_ssec,   &
      c_longitude, c_latitude,                      &
      c_xl93, c_yl93,                               &
      c_compute_z_ground, c_iatmst, c_theo_interp,  &
      c_sedimentation_model, c_deposition_model,    &
      c_nucleation_model, c_subgrid_model,          &
      c_distribution_model,                         &
      c_model, c_isepchemistry, c_nespg, c_nrg,     &
      c_chem_with_photo,  c_modelaero,              &
      c_frozen_gas_chem, c_init_gas_with_lib,       &
      c_init_aero_with_lib, c_nlayer,               &
      c_nsize, c_imeteo,                            &
      c_nbmetd, c_nbmett, c_nbmetm, c_iatra1, c_nbmaxt, &
      c_meteo_zi, c_iatsoil,                        &
      c_nvert, c_kvert, c_kmx )

    call c_f_pointer(c_ps, ps)
    call c_f_pointer(c_syear, syear)
    call c_f_pointer(c_squant, squant)
    call c_f_pointer(c_shour, shour)
    call c_f_pointer(c_smin, smin)
    call c_f_pointer(c_ssec, ssec)

    call c_f_pointer(c_longitude, xlon)
    call c_f_pointer(c_latitude, xlat)
    call c_f_pointer(c_xl93, xl93)
    call c_f_pointer(c_yl93, yl93)

    call c_f_pointer(c_compute_z_ground, compute_z_ground)
    call c_f_pointer(c_iatmst, iatmst)
    call c_f_pointer(c_theo_interp, theo_interp)

    call c_f_pointer(c_sedimentation_model, modsedi)
    call c_f_pointer(c_deposition_model, moddep)
    call c_f_pointer(c_nucleation_model, modnuc)
    call c_f_pointer(c_subgrid_model, modsub)
    call c_f_pointer(c_distribution_model, moddis)

    call c_f_pointer(c_model, ichemistry)
    call c_f_pointer(c_isepchemistry, isepchemistry)
    call c_f_pointer(c_nespg, nespg)
    call c_f_pointer(c_nrg, nrg)
    call c_f_pointer(c_chem_with_photo, photolysis)
    call c_f_pointer(c_modelaero, iaerosol)
    call c_f_pointer(c_frozen_gas_chem, nogaseouschemistry)
    call c_f_pointer(c_init_gas_with_lib, init_gas_with_lib)
    call c_f_pointer(c_init_aero_with_lib, init_aero_with_lib)
    call c_f_pointer(c_nlayer, nlayer_aer)
    call c_f_pointer(c_nsize, n_aer)
    call c_f_pointer(c_imeteo, imeteo)
    call c_f_pointer(c_nbmetd, nbmetd)
    call c_f_pointer(c_nbmett, nbmett)
    call c_f_pointer(c_nbmetm, nbmetm)
    call c_f_pointer(c_iatra1, iatra1)
    call c_f_pointer(c_nbmaxt, nbmaxt)
    call c_f_pointer(c_meteo_zi, meteo_zi)
    call c_f_pointer(c_iatsoil, iatsoil)

    call c_f_pointer(c_nvert, nvert)
    call c_f_pointer(c_kvert, kvert)
    call c_f_pointer(c_kmx, kmx)

    return

  end subroutine atmo_init

!===============================================================================

!> \brief Allocate and map to C meteo data
subroutine allocate_map_atmo

use cs_c_bindings
use atsoil

implicit none

procedure() :: atlecm

! Local variables
integer :: n_level, n_times, n_level_t

type(c_ptr) :: c_z_dyn_met, c_z_temp_met, c_u_met, c_v_met, c_time_met
type(c_ptr) :: c_w_met
type(c_ptr) :: c_hyd_p_met, c_pot_t_met, c_ek_met, c_ep_met
type(c_ptr) :: c_xyvert, c_zvert, c_acinfe
type(c_ptr) :: c_dacinfe, c_aco2, c_aco2s
type(c_ptr) :: c_daco2, c_daco2s
type(c_ptr) :: c_acsup, c_acsups
type(c_ptr) :: c_dacsup, c_dacsups
type(c_ptr) :: c_tauzq, c_tauz, c_zq
type(c_ptr) :: c_zray, c_rayi, c_rayst
type(c_ptr) :: c_iru, c_ird, c_solu, c_sold

integer(c_int), dimension(2) :: dim_hyd_p_met, dim_u_met, dim_pot_t_met
integer(c_int), dimension(2) :: dim_ek_met, dim_ep_met
integer(c_int), dimension(2) :: dim_xyvert, dim_kmx2, dim_kmx_nvert

if (imeteo.eq.1) then
  call atlecm(0)
endif
if (imeteo.eq.2) then
  call cs_atmo_init_meteo_profiles()
endif

! Allocate additional arrays for 1D radiative model
if (iatra1.eq.1) then

  ! Allocate additional arrays for 1D radiative model
  allocate(soilvert(nvert))

endif

call cs_f_atmo_arrays_get_pointers(c_z_dyn_met, c_z_temp_met,     &
                                   c_u_met, c_v_met, c_w_met,     &
                                   c_time_met,                    &
                                   c_hyd_p_met, c_pot_t_met,      &
                                   c_ek_met, c_ep_met,            &
                                   c_xyvert, c_zvert, c_acinfe,   &
                                   c_dacinfe, c_aco2, c_aco2s,    &
                                   c_daco2, c_daco2s,             &
                                   c_acsup, c_acsups,             &
                                   c_dacsup, c_dacsups,           &
                                   c_tauzq, c_tauz, c_zq,         &
                                   c_zray, c_rayi, c_rayst,       &
                                   c_iru, c_ird, c_solu, c_sold,  &
                                   dim_u_met, dim_hyd_p_met,      &
                                   dim_pot_t_met, dim_ek_met,     &
                                   dim_ep_met,                    &
                                   dim_xyvert, dim_kmx2, dim_kmx_nvert)

call c_f_pointer(c_z_dyn_met, zdmet, [nbmetd])
call c_f_pointer(c_z_temp_met, ztmet, [nbmaxt])
call c_f_pointer(c_u_met, umet, [dim_u_met])
call c_f_pointer(c_v_met, vmet, [dim_u_met])
call c_f_pointer(c_w_met, wmet, [dim_u_met])
call c_f_pointer(c_time_met, tmmet, [nbmetm])
call c_f_pointer(c_hyd_p_met, phmet, [dim_hyd_p_met])
call c_f_pointer(c_pot_t_met, tpmet, [dim_pot_t_met])
call c_f_pointer(c_ek_met, ekmet, [dim_ek_met])
call c_f_pointer(c_ep_met, epmet, [dim_ep_met])

call c_f_pointer(c_xyvert , xyvert , [dim_xyvert])
call c_f_pointer(c_zvert  , zvert  , [kmx])
call c_f_pointer(c_acinfe , acinfe , [kmx])
call c_f_pointer(c_dacinfe, dacinfe, [kmx])
call c_f_pointer(c_aco2   , aco2   , [dim_kmx2])
call c_f_pointer(c_aco2s  , aco2s  , [dim_kmx2])
call c_f_pointer(c_daco2  , daco2  , [dim_kmx2])
call c_f_pointer(c_daco2s , daco2s , [dim_kmx2])
call c_f_pointer(c_acsup  , acsup  , [kmx])
call c_f_pointer(c_acsups , acsups , [kmx])
call c_f_pointer(c_dacsup , dacsup , [kmx])
call c_f_pointer(c_dacsups, dacsups, [kmx])
call c_f_pointer(c_tauzq  , tauzq  , [kmx+1])!TODO check +1 works???
call c_f_pointer(c_tauz   , tauz   , [kmx+1])
call c_f_pointer(c_zq     , zq     , [kmx+1])
call c_f_pointer(c_rayi   , rayi   , [dim_kmx_nvert])
call c_f_pointer(c_zray   , zray   , [kmx])
call c_f_pointer(c_rayst  , rayst  , [dim_kmx_nvert])
call c_f_pointer(c_iru    , iru    , [dim_kmx_nvert])
call c_f_pointer(c_ird    , ird    , [dim_kmx_nvert])
call c_f_pointer(c_solu   , solu   , [dim_kmx_nvert])
call c_f_pointer(c_sold   , sold   , [dim_kmx_nvert])

! Allocate additional arrays for Water Microphysics

if (imeteo.gt.0) then

  ! NB : only ztmet,ttmet,qvmet,ncmet are extended to 11000m if iatra1=1
  !           rmet,tpmet,phmet
  n_level = max(1, nbmetd)
  n_times = max(1, nbmetm)
  n_level_t = max(1, nbmaxt)

  if (iatmst.ge.1) then
    allocate(dpdt_met(n_level))
    allocate(mom(3, n_level))
    allocate(mom_met(3, n_level))
  endif

  allocate(ttmet(n_level_t,n_times), qvmet(n_level_t,n_times),  &
           ncmet(n_level_t,n_times))
  allocate(pmer(n_times))
  allocate(xmet(n_times), ymet(n_times))
  allocate(rmet(n_level_t,n_times))

endif

end subroutine allocate_map_atmo

!==============================================================================

!> \brief Initialisation of meteo data
subroutine init_meteo

use cs_c_bindings

implicit none

procedure() :: mestcr, gridcr, mestde

! Prepare interpolation for 1D radiative model
if (imeteo.gt.0) then

  if (iatra1.eq.1) then
    call mestcr("rayi"//c_null_char, 1, 0, idrayi)
    call mestcr("rayst"//c_null_char, 1, 0, idrayst)
    call gridcr("int_grid"//c_null_char, igrid)
  endif

endif

end subroutine init_meteo

!===============================================================================

!> \brief Final step for deallocation
subroutine finalize_meteo

use cs_c_bindings
use atsoil

implicit none

procedure() :: mestcr, grides, mestde

if (imeteo.gt.0) then

  if (allocated(mom)) then
    deallocate(mom, mom_met, dpdt_met)
  endif
  deallocate(ttmet, qvmet, ncmet)
  deallocate(pmer)
  deallocate(xmet, ymet)
  deallocate(rmet)

  if (iatra1.eq.1) then

    deallocate(soilvert)

    call mestde ()

    call grides ()

  endif

endif

end subroutine finalize_meteo

!---------------------------------------------------------------------------

!> \brief Compute LMO, friction velocity ustar, friction temperature
!>        tstar from a thermal flux using Monin Obukhov

!> \param[in]  z             altitude
!> \param[in]  z0
!> \param[in]  du            velocity difference
!> \param[in]  flux          thermal flux
!> \param[in]  tm
!> \param[in]  gredu
!> \param[out] dlmo          Inverse Monin Obukhov length
!> \param[out] ustar         friction velocity

subroutine mo_compute_from_thermal_flux(z,z0,du,flux,tm,gredu,dlmo,ustar) &
  bind(C, name='cs_f_mo_compute_from_thermal_flux')

  use cstphy
  use cstnum
  use entsor

  implicit none

  ! Arguments
  real(c_double) :: z,z0,du,tm,gredu,dlmo,ustar,flux

  ! Local variables
  double precision tstar
  double precision coef_mom
  double precision coef_mom_old
  double precision prec_ustar
  double precision num, denom
  double precision :: dlmoclip = 1.0d0
  integer icompt

  ! Precision initialisation
  prec_ustar=1.d-2

  icompt=0

  ! Initial inverse LMO (neutral)
  dlmo = 0.d0

  ! Call universal functions
  coef_mom = cs_mo_psim((z+z0),z0,dlmo)

  ! Initial ustar and tstar
  ustar = xkappa * du / coef_mom
  tstar = flux / ustar

123 continue

  icompt=icompt+1

  ! Storage previous values
  coef_mom_old = coef_mom

  ! Update LMO
  num = coef_mom**3.d0 * gredu * flux
  denom = du**3.d0 * xkappa**2.d0 * tm
  if (abs(denom).gt.(epzero*num)) then
    dlmo = num / denom
  else
    dlmo = 0.d0 !FIXME other clipping ?
  endif

  ! Clipping dlmo (we want |LMO| > 1 m  ie 1/|LMO| < 1 m^-1)
  if (abs(dlmo).ge.dlmoclip) then
    if (dlmo.ge.0.d0) dlmo =   dlmoclip
    if (dlmo.le.0.d0) dlmo = - dlmoclip
  endif

  ! Evaluate universal functions
  coef_mom = cs_mo_psim((z+z0),z0,dlmo)

  ! Update ustar,tstar
  ustar = xkappa*du/coef_mom
  tstar = flux/ustar

  ! Convergence test
  if (icompt.le.1000) then
    if (abs((coef_mom-coef_mom_old)).ge.prec_ustar) go to 123
  endif

  return

end subroutine mo_compute_from_thermal_flux

!---------------------------------------------------------------------------

!> \brief Compute LMO, friction velocity ustar, friction temperature
!>        tstar from a thermal difference using Monin Obukhov

!> \param[in]  z             altitude
!> \param[in]  z0
!> \param[in]  du            velocity difference
!> \param[in]  dt            thermal difference
!> \param[in]  tm
!> \param[in]  gredu
!> \param[out] dlmo          Inverse Monin Obukhov length
!> \param[out] ustar         friction velocity

subroutine mo_compute_from_thermal_diff(z,z0,du,dt,tm,gredu,dlmo,ustar) &
  bind(C, name='cs_f_mo_compute_from_thermal_diff')

  use cstphy
  use cstnum
  use entsor

  implicit none

  ! Arguments
  real(c_double) :: z,z0,du,dt,tm,gredu,dlmo,ustar

  ! Local variables
  double precision tstar
  double precision coef_mom,coef_moh
  double precision coef_mom_old,coef_moh_old
  double precision prec_ustar,prec_tstar
  double precision num, denom
  real(c_double) :: zref
  double precision :: dlmoclip = 1.0d0
  integer icompt

  ! Precision initialisation
  prec_ustar=1.d-2
  prec_tstar=1.d-2

  icompt=0

  ! Initial LMO
  dlmo = 0.d0

  ! Call universal functions
  zref = z+z0
  coef_mom = cs_mo_psim(zref,z0,dlmo)
  coef_moh = cs_mo_psih(zref,z0,dlmo)

  ! Initial ustar and tstar
  ustar = xkappa * du / coef_mom
  if (abs(coef_moh).gt.epzero) then
    tstar = xkappa*dt/coef_moh
  else
    tstar = 0.d0
  endif


123 continue

  icompt=icompt+1

  ! Storage previous values
  coef_mom_old = coef_mom
  coef_moh_old = coef_moh

  ! Update LMO
  num = coef_mom**2.d0 * gredu * dt
  denom = (du**2.d0)*tm*coef_moh
  if (abs(denom).gt.(epzero* abs(num))) then
    dlmo = num / denom
  else
    dlmo = 0.d0 !FIXME
  endif

  ! Clipping dlmo (we want |LMO| > 1 m  ie 1/|LMO| < 1 m^-1)
  if (abs(dlmo).ge.dlmoclip) then
    if (dlmo.ge.0.d0) dlmo =   dlmoclip
    if (dlmo.le.0.d0) dlmo = - dlmoclip
  endif

  ! Evaluate universal functions
  coef_mom = cs_mo_psim((z+z0),z0,dlmo)
  coef_moh = cs_mo_psih(z+z0,z0,dlmo)

  ! Update ustar,tstar
  ustar = xkappa*du/coef_mom
  if (abs(coef_moh).gt.epzero) then
    tstar=xkappa*dt/coef_moh
  else
    tstar=0.d0
  endif

  ! Convergence test
  if (icompt.le.1000) then !FIXME compteur max 1000 a mettre en param
    if (abs((coef_mom-coef_mom_old)).ge.prec_ustar) go to 123
    if (abs((coef_moh-coef_moh_old)).ge.prec_tstar) go to 123
  endif

  return

end subroutine mo_compute_from_thermal_diff

end module atincl
