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

!> \file coincl.f90
!> Module for gas combustion

module coincl

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use ppincl
  use radiat, only : nwsgg

  implicit none

  !=============================================================================

  !> maximal number of global species
  integer    ngazgm
  parameter( ngazgm = 25 )

  !--> MODELES DE COMBUSTION

  integer(c_int), pointer, save :: cmtype

  !--> MODELE FLAMME DE DIFFUSION: CHIMIE 3 POINTS

  ! ---- Grandeurs deduites

  !       HINOXY       --> Enthalpie massique d'entree pour l'oxydant
  !       HINFUE       --> Enthalpie massique d'entree pour le fuel

  real(c_double), pointer, save :: hinfue, hinoxy

  !--> MODELE FLAMME DE DIFFUSION: Steady laminar flamelet
  !=============================================================================
  ! Grandeurs flamelets
  ! ngazfl : nombre d'espèces dans la librarie de flammelette
  !          YFUE YOXY YCO2 YH2O YCO YH2
  ! nki    : nombre de flamelets (strain rate)
  ! nxr    : discrétisation de l'enthalpie defect
  ! nzm    : discretisation de la fraction de mélange
  ! nzvar  : Discretisation de la variance
  ! nlibvar: Nombre des variables stockés dans la librairie
  ! ikimid : Indice pour la flammelette sur la branche middle

  integer(c_int), pointer, save :: ngazfl  ! ngazfl <= ngazgm - 1
  integer(c_int), pointer, save :: nki, nxr, nzm, nzvar, nlibvar
  integer(c_int), pointer, save :: ikimid

  ! Manière de calculer la variance de fraction de mélange
  ! mode_fp2m 0: Variance transport equation(VTE)
  ! mode_fp2m 1: 2nd moment of mixture fraction transport equation (STE)
  integer(c_int), pointer, save :: mode_fp2m

  ! Coef. of SGS kinetic energy used for the variance dissipation calculation
  double precision, save :: coef_k = 7.d-2

  ! Column index for each variable in the look-up table
  integer, save :: flamelet_zm,    flamelet_zvar,  flamelet_ki,  flamelet_xr
  integer, save :: flamelet_temp,  flamelet_rho,   flamelet_vis, flamelet_dt
  integer, save :: flamelet_temp2, flamelet_hrr
  integer, save :: flamelet_species(ngazgm)
  integer, save :: flamelet_c,     flamelet_omg_C

  !========================================================================

  !> Library for thermochemical properties in SLFM
  real(c_double), pointer, save :: flamelet_library(:,:,:,:,:) => null()

  !========================================================================
  ! Rayonnement

  !> Library for radiative properties in SLFM
  double precision, allocatable, dimension(:,:,:,:,:,:) :: radiation_library

  !========================================================================

  ! --- Soot model

  real(c_double), pointer, save :: lsp_fuel

  !--> POINTEURS VARIABLES COMBUSTION GAZ

  !> \defgroup gas_combustion Gaz combustion variables pointers

  !> \addtogroup gas_combustion
  !> \{

  ! ---- Variables transportees

  !> id of mixing rate field
  integer, save :: ifm = -1

  !> id of mixing rate variance field
  integer, save :: ifp2m = -1

  !> id of field specifying the second moment of the mixing rate:
  integer, save :: ifsqm = -1

  !> id of transported progress variable field (for ippmod(islfm) >= 2)
  integer, save :: ipvm = -1

  ! ---- Variables d'etat

  !> mass fractions :
  !>  - iym(1): is fuel mass fraction
  !>  - iym(2): oxidiser mass fraction
  !>  - iym(3): product mass fraction
  integer, save :: iym(ngazgm)

  !> state variable (temperature)
  integer, save :: itemp = -1

  !> state variable: Pointer to the reconstructed variance in case of mode_fp2m = 1
  integer, save :: irecvr = -1

  !> state variable: Pointer to the total scalar dissipation rate
  integer, save :: itotki = -1

  !> state variable: Pointer to volumetric heat release rate
  integer, save :: ihrr = -1

  !> state variable: Pointer to defect enthalpy
  integer, save :: ixr = -1

  !> state variable: Pointer to Omega C
  integer, save :: iomgc = -1

  !> state variable: absorption coefficient, when the radiation modelling is activated
  integer, save :: ickabs = -1

  !> state variable:  \f$T^2\f$ term
  integer, save :: it2m = -1
  !> state variable:  \f$T^3\f$ term, when the radiation modelling is activated
  integer, save :: it3m = -1
  !> state variable:  \f$T^4\f$ term, when the radiation modelling is activated
  integer, save :: it4m = -1

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global combustion model flags

    subroutine cs_f_coincl_get_pointers(p_cmtype,              &
                                        p_ngazfl, p_nki,       &
                                        p_nxr, p_nzm,          &
                                        p_nzvar, p_nlibvar,    &
                                        p_ikimid, p_mode_fp2m, &
                                        p_lsp_fuel,            &
                                        p_hinfue, p_hinoxy)    &
      bind(C, name='cs_f_coincl_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_cmtype, p_ngazfl, p_nki, p_nxr, p_nzm
      type(c_ptr), intent(out) :: p_nzvar, p_nlibvar, p_ikimid, p_mode_fp2m
      type(c_ptr), intent(out) :: p_lsp_fuel
      type(c_ptr), intent(out) :: p_hinfue, p_hinoxy
    end subroutine cs_f_coincl_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global combustion model flags

    subroutine cs_f_init_steady_laminar_flamelet_library(p_radiation_library)  &
      bind(C, name='cs_f_init_steady_laminar_flamelet_library')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_radiation_library
    end subroutine cs_f_init_steady_laminar_flamelet_library

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran combustion models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine co_models_init() &
    bind(C, name='cs_f_co_models_init')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_cmtype, c_ngazfl, c_nki, c_nxr, c_nzm,            &
                   c_nzvar, c_nlibvar, c_ikimid, c_mode_fp2m,          &
                   c_lsp_fuel, c_hinfue, c_hinoxy

    call cs_f_coincl_get_pointers(c_cmtype, c_ngazfl, c_nki,           &
                                  c_nxr, c_nzm, c_nzvar,               &
                                  c_nlibvar, c_ikimid,                 &
                                  c_mode_fp2m,                         &
                                  c_lsp_fuel,                          &
                                  c_hinfue, c_hinoxy)

    call c_f_pointer(c_cmtype, cmtype)
    call c_f_pointer(c_ngazfl, ngazfl)
    call c_f_pointer(c_nki, nki)
    call c_f_pointer(c_nxr, nxr)
    call c_f_pointer(c_nzm, nzm)
    call c_f_pointer(c_nzvar, nzvar)
    call c_f_pointer(c_nlibvar, nlibvar)
    call c_f_pointer(c_ikimid, ikimid)
    call c_f_pointer(c_mode_fp2m, mode_fp2m)
    call c_f_pointer(c_lsp_fuel, lsp_fuel)
    call c_f_pointer(c_hinfue, hinfue)
    call c_f_pointer(c_hinoxy, hinoxy)

  end subroutine co_models_init

  !=============================================================================

  subroutine init_steady_laminar_flamelet_library

    use radiat
    implicit none

    type(c_ptr) :: p_flamelet_library

    call cs_f_init_steady_laminar_flamelet_library(p_flamelet_library)

    call c_f_pointer(p_flamelet_library, flamelet_library, &
                     [nlibvar, nxr, nki, nzvar, nzm])

    if (iirayo.eq.1) then
      if(.not.allocated(radiation_library)) then
        allocate(radiation_library(2, nwsgg, nxr, nki, nzvar, nzm))
        radiation_library = 0.d0
      endif
    endif

    return

  end subroutine init_steady_laminar_flamelet_library

  !=============================================================================

  ! Free related arrays
  subroutine finalize_steady_laminar_flamelet_library() &
    bind(C, name='cs_f_finalize_steady_laminar_flamelet_library')

    implicit none

    flamelet_library => null()
    if(allocated(radiation_library)) deallocate(radiation_library)

    return

  end subroutine finalize_steady_laminar_flamelet_library

  !=============================================================================

  !> \brief Initialize Fortran combustion models properties API.
  !> This maps Fortran pointers to global C variables.

  function cs_f_flamelet_rho_idx() result(idx)  &
    bind(C, name='cs_f_flamelet_rho_idx')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int) :: idx

    idx = flamelet_rho - 1  ! C index is zero-based, so shift by 1
  end function cs_f_flamelet_rho_idx

  !=============================================================================

  !> \brief Map Fortran combustion model field ids

  subroutine cs_f_combustion_map_variables()  &
    bind(C, name='cs_f_combustion_map_variables')
    use, intrinsic :: iso_c_binding
    use paramx
    use numvar
    use ppincl
    use field
    use cs_c_bindings
    implicit none

    call field_get_id_try('enthalpy', ihm)
    call field_get_id_try('mixture_fraction', ifm)
    call field_get_id_try('mixture_fraction_variance', ifp2m)
    call field_get_id_try('mixture_fraction_2nd_moment', ifsqm)
    call field_get_id_try('progress_variable', ipvm)

  end subroutine cs_f_combustion_map_variables

  !=============================================================================

  !> \brief Map Fortran combustion model field ids

  subroutine cs_f_combustion_map_properties(iym_c)  &
    bind(C, name='cs_f_combustion_map_properties')
    use, intrinsic :: iso_c_binding
    use paramx
    use numvar
    use ppincl
    use field
    use cs_c_bindings
    implicit none

    integer(c_int), dimension(*) :: iym_c

    integer :: ifm

    call field_get_id_try('temperature', itemp)

    do ifm = 1, ngazgm
      iym(ifm) = iym_c(ifm)
    enddo

    if (ippmod(islfm).ge.0) then

      call field_get_id_try('temperature_2', it2m)

      call field_get_id_try('heat_loss', ixr)
      call field_get_id_try('heat_release_rate', ihrr)
      call field_get_id_try('omega_c', iomgc)
      call field_get_id_try('total_dissipation', itotki)
      call field_get_id_try('reconstructed_fp2m', irecvr)

    endif

    call field_get_id_try('kabs', ickabs)
    call field_get_id_try('temperature_4', it4m)
    call field_get_id_try('temperature_3', it3m)

  end subroutine cs_f_combustion_map_properties

end module coincl
