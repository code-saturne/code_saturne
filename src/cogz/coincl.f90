!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

  implicit none

  !=============================================================================

  !--> MODELES DE COMBUSTION

  ! ---- Grandeurs communes

  ! combustible name
  character(len=12) :: namgas

  ! combustible reaction enthalpy (Pouvoir Calorifique Inferieur)
  double precision, save :: pcigas

  ! mass fraction conversion coefficients
  ! from global species to elementary species
  double precision, pointer, save :: coefeg(:,:)

  ! mole fraction conversion coefficients
  ! from global species to elementary species
  double precision, pointer, save :: compog(:,:)

  !--> MODELE FLAMME DE DIFFUSION (CHIMIE 3 POINTS)

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

  integer, save ::         ientox(nozppm), ientfu(nozppm)

  double precision, save :: tinoxy, tinfue, hinfue, hinoxy, hstoea
  double precision, save :: hh(nmaxhm), ff(nmaxfm), tfh(nmaxfm,nmaxhm)

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

  integer, save ::          ientgf(nozppm), ientgb(nozppm)
  double precision, save :: fment(nozppm), tkent(nozppm), qimp(nozppm)
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

    subroutine cs_f_coincl_get_pointers(p_coefeg, p_compog, &
                                        p_xsoot, p_rosoot)  &
      bind(C, name='cs_f_coincl_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_coefeg, p_compog, p_xsoot, p_rosoot
    end subroutine cs_f_coincl_get_pointers

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

    type(c_ptr) :: c_coefeg, c_compog, c_xsoot, c_rosoot

    call cs_f_coincl_get_pointers(c_coefeg, c_compog, c_xsoot, c_rosoot)

    call c_f_pointer(c_coefeg, coefeg, [ngazem, ngazgm])
    call c_f_pointer(c_compog, compog, [ngazem, ngazgm])
    call c_f_pointer(c_xsoot, xsoot)
    call c_f_pointer(c_rosoot, rosoot)

  end subroutine co_models_init

  !=============================================================================

end module coincl
