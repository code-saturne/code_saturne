!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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
double precision, dimension(:), pointer :: dpdt_met

!> Momentum for each level (used for automatic open boundaries)
double precision, dimension(:,:), pointer :: mom_met
double precision, dimension(:,:), pointer :: mom

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
double precision, dimension(:,:), pointer :: ttmet

!> meteo specific humidity profile (read in the input meteo file)
double precision, dimension(:,:), pointer :: qvmet

!> meteo specific droplet number profile (read in the input meteo file)
double precision, dimension(:,:), pointer :: ncmet

!> X, Y coordinates and sea level pressure of the meteo profile
!> (read in the input meteo file)
double precision, dimension(:,:), pointer :: xyp_met

! Arrays specific to values calculated from the meteo file (cf atlecm.f90):

!> density profile
double precision, dimension(:,:), pointer :: rmet

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
real(c_double), pointer, save :: ps

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
integer(c_int), pointer, save:: ihpm

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
integer(c_int), pointer, save :: nfatr1

!> flag for the standard atmo humidity profile
!> - 0: q = 0 (default)
!> - 1: q = decreasing exponential
integer(c_int), pointer, save :: iqv0

!> pointer for 1D infrared profile
integer(c_int), pointer, save :: idrayi

!> pointer for 1D solar profile
integer(c_int), pointer, save :: idrayst

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

!> Defines the soil constants and variables of the vertical arrays
!> used for the 1D radiative model
!> soil albedo
double precision, dimension(:), pointer :: soil_albedo
!> emissivity
double precision, dimension(:), pointer :: soil_emissi
!> soil thermo temperature
double precision, dimension(:), pointer :: soil_ttsoil
!> soil potential temperature
double precision, dimension(:), pointer :: soil_tpsoil
!> total water content
double precision, dimension(:), pointer :: soil_totwat
!> surface pressure
double precision, dimension(:), pointer :: soil_pressure
!> density
double precision, dimension(:), pointer :: soil_density

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

!> logaritmic standard deviation of the log-normal law of the droplet spectrum
!> adimensional:  sigc=0.53 other referenced values are 0.28, 0.15
real(c_double), pointer, save :: sigc

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
        ichemistry, isepchemistry, nespg, nrg,                          &
        chem_with_photo,                                                &
        iaerosol, frozen_gas_chem, init_gas_with_lib,                   &
        init_aero_with_lib, n_aero, n_sizebin, imeteo,                  &
        nbmetd, nbmett, nbmetm, iatra1, nbmaxt,                         &
        meteo_zi, iatsoil,                                              &
        nvertv, kvert, kmx, tsini, tprini, qvsini, ihpm, iqv0,          &
        nfatr1, w1ini, w2ini, sigc, idrayi, idrayst)          &
      bind(C, name='cs_f_atmo_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ps, sigc
      type(c_ptr), intent(out) :: compute_z_ground, iatmst, theo_interp
      type(c_ptr), intent(out) :: ichemistry, isepchemistry, nespg, nrg
      type(c_ptr), intent(out) :: sedimentation_model, deposition_model
      type(c_ptr), intent(out) :: nucleation_model
      type(c_ptr), intent(out) :: subgrid_model, distribution_model
      type(c_ptr), intent(out) :: syear, squant, shour, smin, ssec
      type(c_ptr), intent(out) :: longitude, latitude
      type(c_ptr), intent(out) :: x_l93, y_l93, idrayi, idrayst
      type(c_ptr), intent(out) :: iaerosol, frozen_gas_chem
      type(c_ptr), intent(out) :: init_gas_with_lib, init_aero_with_lib
      type(c_ptr), intent(out) :: n_aero, n_sizebin, chem_with_photo
      type(c_ptr), intent(out) :: imeteo
      type(c_ptr), intent(out) :: nbmetd, nbmett, nbmetm, iatra1, nbmaxt
      type(c_ptr), intent(out) :: meteo_zi, iqv0, w1ini, w2ini
      type(c_ptr), intent(out) :: iatsoil, tsini, tprini, qvsini
      type(c_ptr), intent(out) :: nvertv, kvert, kmx, ihpm, nfatr1
    end subroutine cs_f_atmo_get_pointers

    !---------------------------------------------------------------------------

    !> \brief Return pointers to atmo arrays

    subroutine cs_f_atmo_arrays_get_pointers(p_zdmet, p_ztmet, p_xyp_met,      &
         p_umet, p_vmet,                                                       &
         p_wmet  , p_tmmet, p_phmet, p_tpmet, p_ekmet, p_epmet,                &
         p_ttmet , p_rmet , p_qvmet, p_ncmet,                                  &
         p_dpdt_met,                                                           &
         p_mom_met ,                                                           &
         p_mom_cs  ,                                                           &
         p_xyvert, p_zvert, p_acinfe,                                          &
         p_dacinfe, p_aco2, p_aco2s,                                           &
         p_daco2, p_daco2s,                                                    &
         p_acsup, p_acsups,                                                    &
         p_dacsup, p_dacsups,                                                  &
         p_tauzq, p_tauz, p_zq,                                                &
         p_zray, p_rayi, p_rayst,                                              &
         p_iru, p_ird, p_solu, p_sold,                                         &
         p_soil_albedo,                                                        &
         p_soil_emissi,                                                        &
         p_soil_ttsoil,                                                        &
         p_soil_tpsoil,                                                        &
         p_soil_totwat,                                                        &
         p_soil_pressure,                                                      &
         p_soil_density,                                                       &
         dim_nd_nt, dim_ntx_nt,                                                &
         dim_nd_3, dim_nt_3,                                                   &
         dim_xyvert, dim_kmx2, dim_kmx_nvert )                                 &
         bind(C, name='cs_f_atmo_arrays_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), dimension(2) :: dim_nd_nt, dim_ntx_nt, dim_nd_3, dim_nt_3
      integer(c_int), dimension(2) ::  dim_xyvert, dim_kmx2, dim_kmx_nvert
      type(c_ptr), intent(out) :: p_zdmet, p_ztmet, p_xyp_met
      type(c_ptr), intent(out) :: p_umet, p_vmet, p_tmmet, p_wmet
      type(c_ptr), intent(out) :: p_phmet, p_tpmet, p_ekmet, p_epmet
      type(c_ptr), intent(out) :: p_ttmet, p_rmet, p_qvmet, p_ncmet
      type(c_ptr), intent(out) :: p_dpdt_met, p_mom_met, p_mom_cs
      type(c_ptr), intent(out) :: p_xyvert, p_zvert, p_acinfe
      type(c_ptr), intent(out) :: p_dacinfe, p_aco2, p_aco2s
      type(c_ptr), intent(out) :: p_daco2, p_daco2s
      type(c_ptr), intent(out) :: p_acsup, p_acsups
      type(c_ptr), intent(out) :: p_dacsup, p_dacsups
      type(c_ptr), intent(out) :: p_tauzq, p_tauz, p_zq
      type(c_ptr), intent(out) :: p_zray, p_rayi, p_rayst
      type(c_ptr), intent(out) :: p_iru, p_ird, p_solu, p_sold
      type(c_ptr), intent(out) :: p_soil_albedo
      type(c_ptr), intent(out) :: p_soil_emissi
      type(c_ptr), intent(out) :: p_soil_ttsoil
      type(c_ptr), intent(out) :: p_soil_tpsoil
      type(c_ptr), intent(out) :: p_soil_totwat
      type(c_ptr), intent(out) :: p_soil_pressure
      type(c_ptr), intent(out) :: p_soil_density
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
    !> \param[in]  theta_l       potential liquide temperature
    !> \param[in]  p             pressure
    !> \param[out] yw_liq        liquid water mass fraction
    !> \param[out] t_h           temperature of humid air in Celsius
    !> \param[out] rho_h         density of humid air
    !> \param[out] beta_h

    subroutine cs_rho_humidair(ywm, theta_l, p, yw_liq, t_h, rho_h, beta_h) &
        bind(C, name='cs_rho_humidair')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: ywm, theta_l, p
      real(c_double), intent(out) :: yw_liq, t_h, rho_h, beta_h
    end subroutine cs_rho_humidair

    !=============================================================================

    subroutine cs_f_atmo_get_soil_zone(n_faces, n_soil_cat, face_ids)  &
        bind(C, name='cs_f_atmo_get_soil_zone')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(out) :: n_faces
      integer(c_int), intent(out) :: n_soil_cat
      type(c_ptr), intent(out) :: face_ids
    end subroutine cs_f_atmo_get_soil_zone

    subroutine cs_f_boundary_conditions_get_atincl_pointers(p_iautom) &
      bind(C, name='cs_f_boundary_conditions_get_atincl_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_iautom
    end subroutine cs_f_boundary_conditions_get_atincl_pointers

  end interface

contains

  !=============================================================================

  subroutine at_models_bc_map() &
    bind(C, name='cs_f_atmo_models_boundary_conditions_map')

    use, intrinsic :: iso_c_binding
    use mesh, only: nfabor
    implicit none

    ! Arguments

    ! Local variables
    type(c_ptr) :: p_iautom

    call cs_f_boundary_conditions_get_atincl_pointers(p_iautom)
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
    use atchem
    use atsoil
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
    type(c_ptr) :: c_xl93, c_yl93, c_sigc
    type(c_ptr) :: c_modelaero, c_frozen_gas_chem, c_nlayer, c_nsize
    type(c_ptr) :: c_init_gas_with_lib, c_init_aero_with_lib, c_chem_with_photo
    type(c_ptr) :: c_imeteo, c_iqv0, c_idrayi, c_idrayst
    type(c_ptr) :: c_nbmetd, c_nbmett, c_nbmetm, c_iatra1, c_nbmaxt
    type(c_ptr) :: c_meteo_zi, c_tprini, c_qvsini, c_nfatr1
    type(c_ptr) :: c_iatsoil, c_tsini, c_isepchemistry, c_w1ini, c_w2ini
    type(c_ptr) :: c_nvert, c_kvert, c_kmx, c_theo_interp, c_ihpm

    call cs_f_atmo_get_pointers(c_ps,               &
      c_syear, c_squant, c_shour, c_smin, c_ssec,   &
      c_longitude, c_latitude,                      &
      c_xl93, c_yl93,                               &
      c_compute_z_ground, c_iatmst, c_theo_interp,  &
      c_sedimentation_model, c_deposition_model,    &
      c_nucleation_model, c_subgrid_model,          &
      c_distribution_model,                         &
      c_model, c_isepchemistry,                     &
      c_nespg, c_nrg, c_chem_with_photo,            &
      c_modelaero, c_frozen_gas_chem,               &
      c_init_gas_with_lib,                          &
      c_init_aero_with_lib, c_nlayer,               &
      c_nsize, c_imeteo,                            &
      c_nbmetd, c_nbmett, c_nbmetm, c_iatra1,       &
      c_nbmaxt, c_meteo_zi, c_iatsoil,              &
      c_nvert, c_kvert, c_kmx, c_tsini, c_tprini,   &
      c_qvsini, c_ihpm, c_iqv0, c_nfatr1, c_w1ini,  &
      c_w2ini, c_sigc, c_idrayi, c_idrayst)

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
    call c_f_pointer(c_nfatr1, nfatr1)
    call c_f_pointer(c_nbmaxt, nbmaxt)
    call c_f_pointer(c_meteo_zi, meteo_zi)
    call c_f_pointer(c_iatsoil, iatsoil)

    call c_f_pointer(c_nvert, nvert)
    call c_f_pointer(c_kvert, kvert)
    call c_f_pointer(c_kmx, kmx)
    call c_f_pointer(c_tsini, tsini)
    call c_f_pointer(c_tprini, tprini)
    call c_f_pointer(c_qvsini, qvsini)
    call c_f_pointer(c_ihpm, ihpm)
    call c_f_pointer(c_iqv0, iqv0)
    call c_f_pointer(c_w1ini, w1ini)
    call c_f_pointer(c_w2ini, w2ini)
    call c_f_pointer(c_sigc, sigc)
    call c_f_pointer(c_idrayi, idrayi)
    call c_f_pointer(c_idrayst, idrayst)

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
type(c_ptr) :: c_z_dyn_met, c_z_temp_met, c_xyp_met
type(c_ptr) :: c_u_met, c_v_met, c_time_met
type(c_ptr) :: c_w_met
type(c_ptr) :: c_hyd_p_met, c_pot_t_met, c_ek_met, c_ep_met
type(c_ptr) :: c_temp_met, c_rho_met, c_qw_met, c_ndrop_met
type(c_ptr) :: c_dpdt_met, c_mom_met, c_mom_cs
type(c_ptr) :: c_xyvert, c_zvert, c_acinfe
type(c_ptr) :: c_dacinfe, c_aco2, c_aco2s
type(c_ptr) :: c_daco2, c_daco2s
type(c_ptr) :: c_acsup, c_acsups
type(c_ptr) :: c_dacsup, c_dacsups
type(c_ptr) :: c_tauzq, c_tauz, c_zq
type(c_ptr) :: c_zray, c_rayi, c_rayst
type(c_ptr) :: c_iru, c_ird, c_solu, c_sold
type(c_ptr) :: c_soil_albedo
type(c_ptr) :: c_soil_emissi
type(c_ptr) :: c_soil_ttsoil
type(c_ptr) :: c_soil_tpsoil
type(c_ptr) :: c_soil_totwat
type(c_ptr) :: c_soil_pressure
type(c_ptr) :: c_soil_density

integer(c_int), dimension(2) :: dim_ntx_nt, dim_nd_nt
integer(c_int), dimension(2) :: dim_nd_3, dim_nt_3
integer(c_int), dimension(2) :: dim_xyvert, dim_kmx2, dim_kmx_nvert

if (imeteo.eq.1) then
  call atlecm(0)
endif
if (imeteo.eq.2) then
  call cs_atmo_init_meteo_profiles()
endif

call cs_f_atmo_arrays_get_pointers(c_z_dyn_met, c_z_temp_met,     &
                                   c_xyp_met,                     &
                                   c_u_met, c_v_met, c_w_met,     &
                                   c_time_met,                    &
                                   c_hyd_p_met, c_pot_t_met,      &
                                   c_ek_met, c_ep_met,            &
                                   c_temp_met,                    &
                                   c_rho_met,                     &
                                   c_qw_met,                      &
                                   c_ndrop_met,                   &
                                   c_dpdt_met,                    &
                                   c_mom_met ,                    &
                                   c_mom_cs  ,                    &
                                   c_xyvert, c_zvert, c_acinfe,   &
                                   c_dacinfe, c_aco2, c_aco2s,    &
                                   c_daco2, c_daco2s,             &
                                   c_acsup, c_acsups,             &
                                   c_dacsup, c_dacsups,           &
                                   c_tauzq, c_tauz, c_zq,         &
                                   c_zray, c_rayi, c_rayst,       &
                                   c_iru, c_ird, c_solu, c_sold,  &
                                   c_soil_albedo,                 &
                                   c_soil_emissi,                 &
                                   c_soil_ttsoil,                 &
                                   c_soil_tpsoil,                 &
                                   c_soil_totwat,                 &
                                   c_soil_pressure,               &
                                   c_soil_density,                &
                                   dim_nd_nt, dim_ntx_nt,         &
                                   dim_nd_3, dim_nt_3,            &
                                   dim_xyvert, dim_kmx2, dim_kmx_nvert)

call c_f_pointer(c_z_dyn_met, zdmet, [nbmetd])
call c_f_pointer(c_z_temp_met, ztmet, [nbmaxt])
call c_f_pointer(c_xyp_met, xyp_met, [dim_nt_3])
call c_f_pointer(c_u_met, umet, [dim_nd_nt])
call c_f_pointer(c_v_met, vmet, [dim_nd_nt])
call c_f_pointer(c_w_met, wmet, [dim_nd_nt])
call c_f_pointer(c_time_met, tmmet, [nbmetm])
call c_f_pointer(c_hyd_p_met, phmet, [dim_ntx_nt])
call c_f_pointer(c_pot_t_met, tpmet, [dim_ntx_nt])
call c_f_pointer(c_ek_met, ekmet, [dim_nd_nt])
call c_f_pointer(c_ep_met, epmet, [dim_nd_nt])
call c_f_pointer(c_temp_met, ttmet, [dim_ntx_nt])
call c_f_pointer(c_rho_met, rmet, [dim_ntx_nt])
call c_f_pointer(c_qw_met, qvmet, [dim_ntx_nt])
call c_f_pointer(c_ndrop_met, ncmet, [dim_ntx_nt])
call c_f_pointer(c_dpdt_met, dpdt_met, [nbmetd])
call c_f_pointer(c_mom_met, mom_met, [dim_nd_3])
call c_f_pointer(c_mom_cs, mom, [dim_nd_3])

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
call c_f_pointer(c_tauzq  , tauzq  , [kmx+1])
call c_f_pointer(c_tauz   , tauz   , [kmx+1])
call c_f_pointer(c_zq     , zq     , [kmx+1])
call c_f_pointer(c_rayi   , rayi   , [dim_kmx_nvert])
call c_f_pointer(c_zray   , zray   , [kmx])
call c_f_pointer(c_rayst  , rayst  , [dim_kmx_nvert])
call c_f_pointer(c_iru    , iru    , [dim_kmx_nvert])
call c_f_pointer(c_ird    , ird    , [dim_kmx_nvert])
call c_f_pointer(c_solu   , solu   , [dim_kmx_nvert])
call c_f_pointer(c_sold   , sold   , [dim_kmx_nvert])

call c_f_pointer(c_soil_albedo   , soil_albedo  , [nvert])
call c_f_pointer(c_soil_emissi   , soil_emissi  , [nvert])
call c_f_pointer(c_soil_ttsoil   , soil_ttsoil  , [nvert])
call c_f_pointer(c_soil_tpsoil   , soil_tpsoil  , [nvert])
call c_f_pointer(c_soil_totwat   , soil_totwat  , [nvert])
call c_f_pointer(c_soil_pressure , soil_pressure, [nvert])
call c_f_pointer(c_soil_density  , soil_density , [nvert])
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
subroutine finalize_meteo() &
  bind(C, name='cs_f_finalize_meteo')

use cs_c_bindings
use atsoil

implicit none

procedure() :: mestcr, grides, mestde

if (imeteo.gt.0) then
  if (iatra1.eq.1) then

    call mestde ()

    call grides ()

  endif

endif

end subroutine finalize_meteo

end module atincl
