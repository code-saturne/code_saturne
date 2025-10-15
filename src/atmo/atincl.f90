!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2025 EDF S.A.
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

! Arrays specific to values calculated from the meteo file

!> density profile
double precision, dimension(:,:), pointer :: rmet

!> potential temperature profile
double precision, dimension(:,:), pointer :: tpmet

!> hydrostatic pressure from Laplace integration
double precision, dimension(:,:), pointer :: phmet

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

! 2.5 Data specific to the 1-D atmospheric radiative module:
!-------------------------------------------------------------------------------
!> flag for the use of the 1-D atmo radiative model
!> - 0 no use (default)
!> - 1 use
integer(c_int), pointer, save :: iatra1

!> 1D radiative model pass frequency
integer(c_int), pointer, save :: nfatr1

!> pointer for 1D infrared profile
integer(c_int), pointer, save :: idrayi

!> pointer for 1D solar profile
integer(c_int), pointer, save :: idrayst

!> grid formed by 1D profiles
integer(c_int), pointer, save :: igrid

! 2.6 Arrays specific to the 1D atmospheric radiative module
!-------------------------------------------------------------------------------

!> horizontal coordinates of the vertical grid
double precision, dimension(:,:), pointer :: xyvert

!> vertical grid for 1D radiative scheme initialize in
!>       cs_user_atmospheric_model.f90
double precision, dimension(:), pointer  :: zray

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

!> Option for liquid water content distribution models
!>  - moddis = 1 : all or nothing
!>  - moddis = 2 : Gaussian distribution
integer(c_int), pointer, save:: moddis

!> logaritmic standard deviation of the log-normal law of the droplet spectrum
!> adimensional:  sigc=0.53 other referenced values are 0.28, 0.15
real(c_double), pointer, save :: sigc

!> Aerosol optical properties

! Aerosol optical depth
!> adimensional :  aod_o3_tot=0.2 other referenced values are  0.10, 0.16
real(c_double), pointer, save :: aod_o3_tot
!> adimensional :  aod_h2o_tot=0.10 other referenced values are  0.06, 0.08
real(c_double), pointer, save :: aod_h2o_tot

!> total aerosol optical depth in the IR domain for thermal radiation
!> deduced from aeronet data
double precision, save:: aod_ir = 0.1d0

!> CO2 concentration in cm NTP with correction to the ratio of
!> molar masses for Mco2 and Mair
!> default is 350ppm
double precision, save:: conco2 = 3.5d-2*44.d0/29.d0

!> Asymmetry factor for O3 (non-dimensional)
!> climatic value gaero_o3=0.66
double precision, save:: gaero_o3 = 0.66d0
!> Asymmetry factor for H2O (non-dimensional)
!> climatic value gaero_h2o=0.64
double precision, save:: gaero_h2o = 0.64d0

!> Single scattering albedo for O3 (non-dimensional)
!> climatic value piaero_o3=0.84, other referenced values are 0.963
double precision, save:: piaero_o3 = 0.84d0
!> Single scattering albedo for H2O (non-dimensional)
!> climatic value piaero_h2o=0.84, other referenced values are 0.964
double precision, save:: piaero_h2o = 0.84d0

!> Fraction of Black carbon (non-dimensional): black_carbon_frac=1.d-8 for no BC
double precision, save:: black_carbon_frac = 0.d0

!> Maximal height for aerosol distribution on the vertical
!> important should be <= zqq(kmray-1);
!> in meters : referenced value: zaero=6000
double precision, save:: zaero = 6000d0

!> Cp of dry air
real(c_double), pointer, save :: cp_a

!> Cp of water vapor
real(c_double), pointer, save :: cp_v

!> Atmospheric radiation model:
!> - Direct Solar the first bit (H2O band)
!> - Direct Solar O3 band for the second bit
!> - diFfuse Solar for the third bit (SIR H2O band)
!> - diFfuse Solar O3 for the fourth bit (SUV O3 band)
!> - InfraRed for the fifth bitPeriod of the radiation module.
integer(c_int), pointer, save :: rad_atmo_model

!> \}

!=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \brief Return pointers to atmo includes

    subroutine cs_f_atmo_get_pointers(ps,                               &
        syear, squant, shour, smin, ssec,                               &
        longitude, latitude,                                            &
        compute_z_ground, iatmst, theo_interp,                          &
        sedimentation_model, deposition_model, nucleation_model,        &
        subgrid_model, distribution_model,                              &
        imeteo, nbmetd, nbmett, nbmetm, iatra1, nbmaxt,                 &
        iatsoil,                                                        &
        nvertv, kvert, kmx, ihpm,                                       &
        nfatr1, sigc, idrayi, idrayst, igrid,                           &
        aod_o3_tot, aod_h2o_tot)                                        &
      bind(C, name='cs_f_atmo_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ps, sigc
      type(c_ptr), intent(out) :: compute_z_ground, iatmst, theo_interp
      type(c_ptr), intent(out) :: sedimentation_model, deposition_model
      type(c_ptr), intent(out) :: nucleation_model
      type(c_ptr), intent(out) :: subgrid_model, distribution_model
      type(c_ptr), intent(out) :: syear, squant, shour, smin, ssec
      type(c_ptr), intent(out) :: longitude, latitude
      type(c_ptr), intent(out) :: idrayi, idrayst, igrid
      type(c_ptr), intent(out) :: imeteo
      type(c_ptr), intent(out) :: nbmetd, nbmett, nbmetm, iatra1, nbmaxt
      type(c_ptr), intent(out) :: iatsoil
      type(c_ptr), intent(out) :: nvertv, kvert, kmx, ihpm, nfatr1
      type(c_ptr), intent(out) :: aod_o3_tot, aod_h2o_tot
    end subroutine cs_f_atmo_get_pointers

    !---------------------------------------------------------------------------

    !> \brief Return pointers to atmo arrays

    subroutine cs_f_atmo_arrays_get_pointers(p_ztmet, p_xyp_met,               &
         p_umet, p_vmet,                                                       &
         p_wmet  , p_tmmet, p_phmet, p_tpmet, p_ekmet, p_epmet,                &
         p_ttmet , p_rmet , p_qvmet, p_ncmet,                                  &
         p_xyvert, p_zray , p_acinfe,                                          &
         p_dacinfe, p_aco2, p_aco2s,                                           &
         p_daco2, p_daco2s,                                                    &
         p_acsup, p_acsups,                                                    &
         p_dacsup, p_dacsups,                                                  &
         p_tauzq, p_tauz, p_zq,                                                &
         p_rayi, p_rayst,                                                      &
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
      type(c_ptr), intent(out) :: p_ztmet, p_xyp_met
      type(c_ptr), intent(out) :: p_umet, p_vmet, p_tmmet, p_wmet
      type(c_ptr), intent(out) :: p_phmet, p_tpmet, p_ekmet, p_epmet
      type(c_ptr), intent(out) :: p_ttmet, p_rmet, p_qvmet, p_ncmet
      type(c_ptr), intent(out) :: p_xyvert, p_zray , p_acinfe
      type(c_ptr), intent(out) :: p_dacinfe, p_aco2, p_aco2s
      type(c_ptr), intent(out) :: p_daco2, p_daco2s
      type(c_ptr), intent(out) :: p_acsup, p_acsups
      type(c_ptr), intent(out) :: p_dacsup, p_dacsups
      type(c_ptr), intent(out) :: p_tauzq, p_tauz, p_zq
      type(c_ptr), intent(out) :: p_rayi, p_rayst
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

    subroutine cs_f_atmo_get_soil_zone(n_faces, n_soil_cat, face_ids)  &
        bind(C, name='cs_f_atmo_get_soil_zone')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(out) :: n_faces
      integer(c_int), intent(out) :: n_soil_cat
      type(c_ptr), intent(out) :: face_ids
    end subroutine cs_f_atmo_get_soil_zone

    !---------------------------------------------------------------------------

    subroutine cs_air_glob_properties_get_pointers(cp_a, cp_v) &
        bind(C, name='cs_air_glob_properties_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: cp_a, cp_v
    end subroutine cs_air_glob_properties_get_pointers

    !---------------------------------------------------------------------------

    subroutine cs_rad_transfer_get_pointers(p_rad_atmo_model)    &
      bind(C, name='cs_rad_transfer_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_rad_atmo_model
    end subroutine cs_rad_transfer_get_pointers

  end interface

contains

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

  !> \brief Map Fortran to C variables

  subroutine atmo_init

    use, intrinsic :: iso_c_binding
    use cs_c_bindings

    implicit none

    ! Local variables
    type(c_ptr) :: c_ps
    type(c_ptr) :: c_compute_z_ground, c_iatmst
    type(c_ptr) :: c_sedimentation_model, c_deposition_model, c_nucleation_model
    type(c_ptr) :: c_subgrid_model
    type(c_ptr) :: c_distribution_model
    type(c_ptr) :: c_syear, c_squant, c_shour, c_smin, c_ssec
    type(c_ptr) :: c_longitude, c_latitude
    type(c_ptr) :: c_sigc
    type(c_ptr) :: c_imeteo, c_idrayi, c_idrayst, c_igrid
    type(c_ptr) :: c_nbmetd, c_nbmett, c_nbmetm, c_iatra1, c_nbmaxt
    type(c_ptr) :: c_nfatr1
    type(c_ptr) :: c_iatsoil
    type(c_ptr) :: c_nvert, c_kvert, c_kmx, c_theo_interp, c_ihpm
    type(c_ptr) :: c_aod_o3_tot, c_aod_h2o_tot
    type(c_ptr) :: c_cp_a, c_cp_v
    type(c_ptr) :: p_rad_atmo_model

    call cs_f_atmo_get_pointers(c_ps,               &
      c_syear, c_squant, c_shour, c_smin, c_ssec,   &
      c_longitude, c_latitude,                      &
      c_compute_z_ground, c_iatmst, c_theo_interp,  &
      c_sedimentation_model, c_deposition_model,    &
      c_nucleation_model, c_subgrid_model,          &
      c_distribution_model, c_imeteo,               &
      c_nbmetd, c_nbmett, c_nbmetm, c_iatra1,       &
      c_nbmaxt, c_iatsoil,                          &
      c_nvert, c_kvert, c_kmx,                      &
      c_ihpm, c_nfatr1,                             &
      c_sigc, c_idrayi, c_idrayst, c_igrid,         &
      c_aod_o3_tot, c_aod_h2o_tot)

    call c_f_pointer(c_ps, ps)
    call c_f_pointer(c_syear, syear)
    call c_f_pointer(c_squant, squant)
    call c_f_pointer(c_shour, shour)
    call c_f_pointer(c_smin, smin)
    call c_f_pointer(c_ssec, ssec)

    call c_f_pointer(c_longitude, xlon)
    call c_f_pointer(c_latitude, xlat)

    call c_f_pointer(c_compute_z_ground, compute_z_ground)
    call c_f_pointer(c_iatmst, iatmst)
    call c_f_pointer(c_theo_interp, theo_interp)

    call c_f_pointer(c_distribution_model, moddis)

    call c_f_pointer(c_imeteo, imeteo)
    call c_f_pointer(c_nbmetd, nbmetd)
    call c_f_pointer(c_nbmett, nbmett)
    call c_f_pointer(c_nbmetm, nbmetm)
    call c_f_pointer(c_iatra1, iatra1)
    call c_f_pointer(c_nfatr1, nfatr1)
    call c_f_pointer(c_nbmaxt, nbmaxt)
    call c_f_pointer(c_iatsoil, iatsoil)

    call c_f_pointer(c_nvert, nvert)
    call c_f_pointer(c_kvert, kvert)
    call c_f_pointer(c_kmx, kmx)
    call c_f_pointer(c_ihpm, ihpm)
    call c_f_pointer(c_sigc, sigc)
    call c_f_pointer(c_idrayi, idrayi)
    call c_f_pointer(c_idrayst, idrayst)
    call c_f_pointer(c_igrid, igrid)

    call c_f_pointer(c_aod_o3_tot,  aod_o3_tot)
    call c_f_pointer(c_aod_h2o_tot, aod_h2o_tot)

    call cs_air_glob_properties_get_pointers(c_cp_a, c_cp_v)

    call c_f_pointer(c_cp_a, cp_a)
    call c_f_pointer(c_cp_v, cp_v)

    call cs_rad_transfer_get_pointers(p_rad_atmo_model)

    call c_f_pointer(p_rad_atmo_model, rad_atmo_model)

    return

  end subroutine atmo_init

!===============================================================================

!> \brief Allocate and map to C meteo data
subroutine allocate_map_atmo () &
  bind(C, name='cs_f_allocate_map_atmo')

  use cs_c_bindings

  implicit none

! Local variables
type(c_ptr) :: c_z_temp_met, c_xyp_met
type(c_ptr) :: c_u_met, c_v_met, c_time_met
type(c_ptr) :: c_w_met
type(c_ptr) :: c_hyd_p_met, c_pot_t_met, c_ek_met, c_ep_met
type(c_ptr) :: c_temp_met, c_rho_met, c_qw_met, c_ndrop_met
type(c_ptr) :: c_xyvert, c_zray , c_acinfe
type(c_ptr) :: c_dacinfe, c_aco2, c_aco2s
type(c_ptr) :: c_daco2, c_daco2s
type(c_ptr) :: c_acsup, c_acsups
type(c_ptr) :: c_dacsup, c_dacsups
type(c_ptr) :: c_tauzq, c_tauz, c_zq
type(c_ptr) :: c_rayi, c_rayst
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

call cs_f_atmo_arrays_get_pointers(c_z_temp_met,                  &
                                   c_xyp_met,                     &
                                   c_u_met, c_v_met, c_w_met,     &
                                   c_time_met,                    &
                                   c_hyd_p_met, c_pot_t_met,      &
                                   c_ek_met, c_ep_met,            &
                                   c_temp_met,                    &
                                   c_rho_met,                     &
                                   c_qw_met,                      &
                                   c_ndrop_met,                   &
                                   c_xyvert, c_zray , c_acinfe,   &
                                   c_dacinfe, c_aco2, c_aco2s,    &
                                   c_daco2, c_daco2s,             &
                                   c_acsup, c_acsups,             &
                                   c_dacsup, c_dacsups,           &
                                   c_tauzq, c_tauz, c_zq,         &
                                   c_rayi, c_rayst,               &
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

call c_f_pointer(c_xyvert , xyvert , [dim_xyvert])
call c_f_pointer(c_zray   , zray   , [kmx])
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

!===============================================================================

end module atincl
