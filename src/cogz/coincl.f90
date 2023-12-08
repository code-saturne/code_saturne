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
  !       IENTOX       --> indicateur oxydant par type de facette d'entree
  !       IENTFU       --> indicateur fuel    par type de facette d'entree

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

  integer(c_int), pointer, save :: ientox(:), ientfu(:)

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

  integer, save :: ngazfl  ! ngazfl <= ngazgm - 1
  integer, save :: nki, nxr, nzm, nzvar, nlibvar
  integer, save :: ikimid = 1

  ! Manière de calculer la variance de fraction de mélange
  ! mode_fp2m 0: Variance transport equation(VTE)
  ! mode_fp2m 1: 2nd moment of mixture fraction transport equation (STE)
  integer, save :: mode_fp2m = 1

  ! Coef. of SGS kinetic energy used for the variance dissipation calculation
  double precision, save :: coef_k = 7.d-2

  ! Column index for each variable in the look-up table
  integer, save :: FLAMELET_ZM,    FLAMELET_ZVAR,  FLAMELET_KI,  FLAMELET_XR
  integer, save :: FLAMELET_TEMP,  FLAMELET_RHO,   FLAMELET_VIS, FLAMELET_DT
  integer, save :: FLAMELET_TEMP2, FLAMELET_HRR
  integer, save :: FLAMELET_SPECIES(ngazgm)
  integer, save :: FLAMELET_C,     FLAMELET_OMG_C

  character(len=12) :: FLAMELET_SPECIES_NAME(ngazgm)

  !========================================================================

  !> Library for thermochemical properties in SLFM
  double precision, allocatable, dimension(:,:,:,:,:) :: flamelet_library

  !========================================================================
  ! Rayonnement

  !> Library for radiative properties in SLFM
  double precision, allocatable, dimension(:,:,:,:,:,:) :: radiation_library

  !========================================================================

  !--> MODELE FLAMME DE PREMELANGE (MODELE EBU)

  ! ---- Grandeurs fournies par l'utilisateur dans usebuc.f90

  !       IENTGF       --> indicateur gaz frais  par type de facette d'entree
  !       IENTGB       --> indicateur gaz brules par type de facette d'entree
  !       QIMP         --> Debit impose en kg/s
  !       FMENT        --> Taux de melange par type de facette d'entree
  !       TKENT        --> Temperature en K par type de facette d'entree
  !       FRMEL        --> Taux de melange constant pour modeles 0 et 1
  !       TGF          --> Temperature gaz frais en K identique
  !                        pour premelange frais et dilution
  !       CEBU         --> Constante Eddy break-Up

  ! ---- Grandeurs deduites

  !       HGF          --> Enthalpie massique gaz frais identique
  !                        pour premelange frais et dilution
  !       TGBAD        --> Temperature adiabatique gaz brules en K

  !integer, save ::          ientgf(nozppm), ientgb(nozppm)
  integer(c_int), pointer, save :: ientgf(:), ientgb(:)
  real(c_double), pointer, save :: qimp(:), fment(:), tkent(:)
  !double precision, save :: fment(nozppm), tkent(nozppm)
  double precision, save :: frmel, tgf, cebu, hgf, tgbad

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

  double precision, save :: vref, lref, ta, tstar
  double precision, save :: fmin, fmax, hmin, hmax
  double precision, save :: coeff1, coeff2, coeff3

  ! --- Soot model

  !     XSOOT : soot fraction production (isoot = 0)
  !     ROSOOT: soot density

  double precision, pointer, save :: xsoot, rosoot

  !=============================================================================

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global combustion model flags

    subroutine cs_f_coincl_get_pointers(p_coefeg, p_compog,   &
                                        p_xsoot, p_rosoot,    &
                                        p_hinfue, p_hinoxy,   &
                                        p_pcigas, p_tinfue,   &
                                        p_tinoxy)             &
      bind(C, name='cs_f_coincl_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_coefeg, p_compog, p_xsoot, p_rosoot
      type(c_ptr), intent(out) :: p_hinfue, p_hinoxy, p_pcigas, p_tinfue
      type(c_ptr), intent(out) :: p_tinoxy
    end subroutine cs_f_coincl_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving BC zone array pointers

    subroutine cs_f_boundary_conditions_get_coincl_pointers(p_ientfu, p_ientox, &
                                                            p_ientgb, p_ientgf, &
                                                            p_tkent,  p_fment,  &
                                                            p_qimp) &
      bind(C, name='cs_f_boundary_conditions_get_coincl_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_ientfu, p_ientox, p_ientgb, p_ientgf
      type(c_ptr), intent(out) :: p_tkent,  p_fment, p_qimp
    end subroutine cs_f_boundary_conditions_get_coincl_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran combustion models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine co_models_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_coefeg, c_compog, c_xsoot,  &
                   c_rosoot, c_hinfue, c_hinoxy, &
                   c_pcigas, c_tinfue, c_tinoxy

    call cs_f_coincl_get_pointers(c_coefeg, c_compog,  &
                                  c_xsoot,  c_rosoot,  &
                                  c_hinfue, c_hinoxy,  &
                                  c_pcigas, c_tinfue,  &
                                  c_tinoxy)

    call c_f_pointer(c_coefeg, coefeg, [ngazem, ngazgm])
    call c_f_pointer(c_compog, compog, [ngazem, ngazgm])
    call c_f_pointer(c_xsoot, xsoot)
    call c_f_pointer(c_rosoot, rosoot)
    call c_f_pointer(c_hinfue, hinfue)
    call c_f_pointer(c_hinoxy, hinoxy)
    call c_f_pointer(c_pcigas, pcigas)
    call c_f_pointer(c_tinfue, tinfue)
    call c_f_pointer(c_tinoxy, tinoxy)

  end subroutine co_models_init


  !> \brief Map Fortran physical models boundary condition info.
  !> This maps Fortran pointers to global C variables.

  subroutine co_models_bc_map

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_ientfu, p_ientox, p_ientgb, p_ientgf
    type(c_ptr) :: p_tkent,  p_fment,  p_qimp

    call cs_f_boundary_conditions_get_coincl_pointers(p_ientfu, p_ientox, &
                                                      p_ientgb, p_ientgf, &
                                                      p_tkent,  p_fment,  &
                                                      p_qimp)
    call c_f_pointer(p_ientfu, ientfu, [nozppm])
    call c_f_pointer(p_ientox, ientox, [nozppm])
    call c_f_pointer(p_ientgb, ientgb, [nozppm])
    call c_f_pointer(p_ientgf, ientgf, [nozppm])
    call c_f_pointer(p_tkent,  tkent,  [nozppm])
    call c_f_pointer(p_fment,  fment,  [nozppm])
    call c_f_pointer(p_qimp,   qimp,   [nozppm])

  end subroutine co_models_bc_map

  !=============================================================================

  subroutine init_steady_laminar_flamelet_library

    use radiat
    implicit none

    if(.not.allocated(flamelet_library)) then
      allocate(flamelet_library(nlibvar, nxr, nki, nzvar, nzm))
      flamelet_library = 0.d0
    endif

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
  subroutine finalize_steady_laminar_flamelet_library

    implicit none

    if(allocated(flamelet_library)) deallocate(flamelet_library)
    if(allocated(radiation_library)) deallocate(radiation_library)

    return

  end subroutine finalize_steady_laminar_flamelet_library

end module coincl
