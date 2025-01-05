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

  use ppppar
  use ppincl
  use radiat, only : nwsgg

  implicit none

  !=============================================================================

  !--> MODELES DE COMBUSTION

  integer(c_int), pointer, save :: cmtype

  ! ---- Grandeurs communes

  ! combustible name
  character(len=12) :: namgas

  ! combustible reaction enthalpy (Pouvoir Calorifique Inferieur)
  double precision, save, pointer :: pcigas

  ! mass fraction conversion coefficients
  ! from global species to elementary species
  double precision, pointer, save :: coefeg(:,:)

  ! mole fraction conversion coefficients
  ! from global species to elementary species
  double precision, pointer, save :: compog(:,:)

  !--> MODELE FLAMME DE DIFFUSION: CHIMIE 3 POINTS

  ! ---- Grandeurs fournies par l'utilisateur dans usd3pc.f90

  !       TINOXY       --> Temperature d'entree pour l'oxydant en K
  !       TINFUE       --> Temperature d'entree pour le fuel en K

  ! ---- Grandeurs deduites

  !       HINOXY       --> Enthalpie massique d'entree pour l'oxydant
  !       HINFUE       --> Enthalpie massique d'entree pour le fuel
  !       HSTOEA       --> Temperature a la stoechiometrie adiabatique
  !       NMAXF        --> Nb de points de tabulation en F
  !       NMAXFM       --> Nb maximal de points de tabulation en F
  !       NMAXH        --> Nb de points de tabulation en H
  !       NMAXHM       --> Nb maximal de points de tabulation en H
  !       HH           --> Enthalpie stoechiometrique tabulee
  !       FF           --> Richesse tabulee
  !       TFH(IF,IH)   --> Tabulation richesse - enthalpie stoechiometrique

  integer    nmaxf, nmaxfm, nmaxh, nmaxhm
  parameter( nmaxfm = 15 , nmaxhm = 15)

  double precision, save :: hstoea
  double precision, save :: hh(nmaxhm), ff(nmaxfm), tfh(nmaxfm,nmaxhm)
  real(c_double), pointer, save :: tinfue, tinoxy
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

  character(len=12) :: flamelet_species_name(ngazgm)

  !========================================================================

  !> Library for thermochemical properties in SLFM
  real(c_double), pointer, save :: flamelet_library(:,:,:,:,:) => null()

  !========================================================================
  ! Rayonnement

  !> Library for radiative properties in SLFM
  double precision, allocatable, dimension(:,:,:,:,:,:) :: radiation_library

  !========================================================================

  !--> MODELE FLAMME DE PREMELANGE (MODELE EBU)

  ! ---- Grandeurs fournies par l'utilisateur dans usebuc.f90

  !       FRMEL        --> Taux de melange constant pour modeles 0 et 1
  !       TGF          --> Temperature gaz frais en K identique
  !                        pour premelange frais et dilution
  !       CEBU         --> Constante Eddy break-Up

  real(c_double), pointer, save :: frmel, tgf, cebu

  !--> MODELE DE FLAMME DE PREMELANGE LWC

  !       NDRACM : nombre de pics de Dirac maximum
  !       NDIRAC : nombre de Dirac (en fonction du modele)

  integer ndracm
  parameter (ndracm = 5)

  integer, save :: ndirac

  ! --- Grandeurs fournies par l'utilisateur dans uslwc1.f90

  !       VREF : Vitesse de reference
  !       LREF : Longueur de reference
  !         TA : Temperature d'activation
  !      TSTAR : Temperature de cross-over

  integer, save :: irhol(ndracm), iteml(ndracm), ifmel(ndracm)
  integer, save :: ifmal(ndracm), iampl(ndracm), itscl(ndracm)
  integer, save :: imaml(ndracm), ihhhh(ndracm), imam

  real(c_double), pointer, save :: vref, lref, ta, tstar
  real(c_double), pointer, save :: fmin, fmax, hmin, hmax
  real(c_double), pointer, save :: coeff1, coeff2, coeff3

  ! --- Soot model

  !     XSOOT : soot fraction production (isoot = 0)
  !     ROSOOT: soot density

  integer(c_int), pointer, save :: isoot
  real(c_double), pointer, save :: xsoot, rosoot, lsp_fuel

  ! --- Temperature/Enthalpy conversion

  !> use JANAF or not
  logical(c_bool), pointer, save :: use_janaf

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

  !> id of fresh gas mass fraction field
  integer, save :: iygfm = -1

  !> id of mass fraction field
  integer, save :: iyfm = -1

  !> id of mass fraction variance field
  integer, save :: iyfp2m = -1

  !> id of mass fraction covariance field
  integer, save :: icoyfp = -1

  ! ---- Variables d'etat

  !> mass fractions :
  !>  - iym(1): is fuel mass fraction
  !>  - iym(2): oxidiser mass fraction
  !>  - iym(3): product mass fraction
  !> ibym() contains the matching field ids at boundary faces
  integer, save :: iym(ngazgm)
  integer, save :: ibym(ngazgm)

  !> state variable (temperature)
  integer, save :: itemp = -1
  !> state variable
  integer, save :: ifmin
  !> state variable
  integer, save :: ifmax

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

  ! pointer for source term in combustion
  ! TODO
  integer, save :: itsc = -1

  !> pointer for soot precursor number in isca (isoot = 1)
  integer, save :: inpm = -1

  !> pointer for soot mass fraction in isca (isoot = 1)
  integer, save :: ifsm = -1

  !> Burke Schumann combustion model constants
  integer, parameter :: n_z = 80, n_xr = 5, n_zvar = 10
  integer, parameter :: n_var_bsh = 7, nvar_turb = 10

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global combustion model flags

    subroutine cs_f_coincl_get_pointers(p_cmtype, p_isoot,     &
                                        p_ngazfl, p_nki,       &
                                        p_nxr, p_nzm,          &
                                        p_nzvar, p_nlibvar,    &
                                        p_ikimid, p_mode_fp2m, &
                                        p_use_janaf,           &
                                        p_coefeg, p_compog,    &
                                        p_xsoot, p_rosoot,     &
                                        p_lsp_fuel,            &
                                        p_hinfue, p_hinoxy,    &
                                        p_pcigas, p_tinfue,    &
                                        p_tinoxy,              &
                                        p_vref, p_lref,        &
                                        p_ta, p_tstar,         &
                                        p_fmin, p_fmax,        &
                                        p_hmin, p_hmax,        &
                                        p_coeff1, p_coeff2,    &
                                        p_coeff3,              &
                                        p_tgf, p_frmel,        &
                                        p_cebu)                &
      bind(C, name='cs_f_coincl_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_cmtype, p_isoot, p_ngazfl, p_nki, p_nxr, p_nzm
      type(c_ptr), intent(out) :: p_nzvar, p_nlibvar, p_ikimid, p_mode_fp2m
      type(c_ptr), intent(out) :: p_use_janaf
      type(c_ptr), intent(out) :: p_coefeg, p_compog, p_xsoot, p_rosoot, p_lsp_fuel
      type(c_ptr), intent(out) :: p_hinfue, p_hinoxy, p_pcigas, p_tinfue
      type(c_ptr), intent(out) :: p_tinoxy, p_vref, p_lref, p_ta, p_tstar
      type(c_ptr), intent(out) :: p_fmin, p_fmax, p_hmin, p_hmax
      type(c_ptr), intent(out) :: p_coeff1, p_coeff2, p_coeff3
      type(c_ptr), intent(out) :: p_tgf, p_frmel, p_cebu
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

    type(c_ptr) :: c_cmtype, c_isoot, c_ngazfl, c_nki, c_nxr, c_nzm,   &
                   c_nzvar, c_nlibvar, c_ikimid, c_mode_fp2m,          &
                   c_use_janaf, c_coefeg, c_compog, c_xsoot,           &
                   c_rosoot, c_lsp_fuel, c_hinfue, c_hinoxy,           &
                   c_pcigas, c_tinfue, c_tinoxy,                       &
                   c_vref, c_lref, c_ta, c_tstar,                      &
                   c_fmin, c_fmax, c_hmin, c_hmax,                     &
                   c_coeff1, c_coeff2, c_coeff3,                       &
                   c_tgf, c_frmel, c_cebu

    call cs_f_coincl_get_pointers(c_cmtype, c_isoot, c_ngazfl, c_nki,  &
                                  c_nxr, c_nzm, c_nzvar,               &
                                  c_nlibvar, c_ikimid,                 &
                                  c_mode_fp2m, c_use_janaf,            &
                                  c_coefeg, c_compog,                  &
                                  c_xsoot,  c_rosoot,                  &
                                  c_lsp_fuel,                          &
                                  c_hinfue, c_hinoxy,                  &
                                  c_pcigas, c_tinfue, c_tinoxy,        &
                                  c_vref, c_lref, c_ta, c_tstar,       &
                                  c_fmin, c_fmax, c_hmin, c_hmax,      &
                                  c_coeff1, c_coeff2, c_coeff3,        &
                                  c_tgf, c_frmel, c_cebu)

    call c_f_pointer(c_cmtype, cmtype)
    call c_f_pointer(c_isoot, isoot)
    call c_f_pointer(c_ngazfl, ngazfl)
    call c_f_pointer(c_nki, nki)
    call c_f_pointer(c_nxr, nxr)
    call c_f_pointer(c_nzm, nzm)
    call c_f_pointer(c_nzvar, nzvar)
    call c_f_pointer(c_nlibvar, nlibvar)
    call c_f_pointer(c_ikimid, ikimid)
    call c_f_pointer(c_mode_fp2m, mode_fp2m)
    call c_f_pointer(c_use_janaf, use_janaf)
    call c_f_pointer(c_coefeg, coefeg, [ngazem, ngazgm])
    call c_f_pointer(c_compog, compog, [ngazem, ngazgm])
    call c_f_pointer(c_xsoot, xsoot)
    call c_f_pointer(c_rosoot, rosoot)
    call c_f_pointer(c_lsp_fuel, lsp_fuel)
    call c_f_pointer(c_hinfue, hinfue)
    call c_f_pointer(c_hinoxy, hinoxy)
    call c_f_pointer(c_pcigas, pcigas)
    call c_f_pointer(c_tinfue, tinfue)
    call c_f_pointer(c_tinoxy, tinoxy)
    call c_f_pointer(c_vref, vref)
    call c_f_pointer(c_lref, lref)
    call c_f_pointer(c_ta, ta)
    call c_f_pointer(c_tstar, tstar)
    call c_f_pointer(c_fmin, fmin)
    call c_f_pointer(c_fmax, fmax)
    call c_f_pointer(c_hmin, hmin)
    call c_f_pointer(c_hmax, hmax)
    call c_f_pointer(c_coeff1, coeff1)
    call c_f_pointer(c_coeff2, coeff2)
    call c_f_pointer(c_coeff3, coeff3)
    call c_f_pointer(c_tgf, tgf);
    call c_f_pointer(c_frmel, frmel);
    call c_f_pointer(c_cebu, cebu);

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

end module coincl
