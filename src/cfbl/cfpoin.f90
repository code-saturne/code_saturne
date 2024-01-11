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

!> \file cfpoin.f90
!> Module for fuel combustion

module cfpoin

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  !> \addtogroup compressible
  !> \{

  !> indicator of equation of state mapping \ref cs_cf_model_t::ieos
  integer(c_int), pointer, save :: ieos

  !> thermodynamic variables indicator for initialization
  !> mapping cs_cf_model_t::ithvar
  integer(c_int), pointer, save :: ithvar

  !> indicator for thermodynamic variables initialization
  integer(c_int), pointer, save :: icfgrp

  !> imposed thermal flux indicator at the boundary
  !> (some boundary contributions of the total energy eq. have to be cancelled)
  integer, allocatable, dimension(:), target :: ifbet

  !> boundary convection flux indicator of a Rusanov or an analytical flux
  !> (some boundary contributions of the momentum eq. have to be cancelled)
  integer, allocatable, dimension(:), target :: icvfli

  !> Stiffened gas limit pressure (Pa) for single phase model
  !> Equal to zero in perfect gas
  !> mapping cs_cf_model_t::psginf
  real(c_double), pointer, save :: psginf

  !> Stiffened gas polytropic coefficient (dimensionless) for single phase model
  !> mapping cs_cf_model_t::gammasg
  real(c_double), pointer, save :: gammasg

  !> \addtogroup comp_homogeneous
  !> \{

  !> \anchor hgn_relax_eq_st
  !> homogeneous two-phase flow model indicator for source terms
  !>    -  -1: source terms are disabled
  !>    -   0: source terms are enabled
  !> mapping cs_cf_model_t::hgn_relax_eq_st
  integer(c_int), pointer, save :: hgn_relax_eq_st

  !> \}
  !> \}

  !=============================================================================

  type(c_ptr) :: p_icvfli
  bind(C, name='cs_glob_cf_icvfli') :: p_icvfli

  type(c_ptr) :: p_ifbet
  bind(C, name='cs_glob_cf_ifbet') :: p_ifbet

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global compressible model structure

    subroutine cs_f_cf_model_get_pointers(ieos,            &
                                          ithvar,          &
                                          icfgrp,          &
                                          psginf,          &
                                          gammasg,         &
                                          hgn_relax_eq_st) &
      bind(C, name='cs_f_cf_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ieos, ithvar, icfgrp, psginf, &
                                  gammasg, hgn_relax_eq_st
    end subroutine cs_f_cf_model_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran compressible model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine cf_model_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ieos, c_ithvar, c_icfgrp
    type(c_ptr) :: c_psginf, c_gammasg, c_hgn_relax_eq_st

    call cs_f_cf_model_get_pointers(c_ieos,           &
                                    c_ithvar,         &
                                    c_icfgrp,        &
                                    c_psginf,         &
                                    c_gammasg,        &
                                    c_hgn_relax_eq_st)

    call c_f_pointer(c_ieos, ieos)
    call c_f_pointer(c_ithvar, ithvar)
    call c_f_pointer(c_icfgrp, icfgrp)
    call c_f_pointer(c_psginf, psginf)
    call c_f_pointer(c_gammasg, gammasg)
    call c_f_pointer(c_hgn_relax_eq_st, hgn_relax_eq_st)

  end subroutine cf_model_init

  !> \brief Allocate boundary flux indicators array

  subroutine init_compf (nfabor)

    implicit none

    integer nfabor

    allocate(ifbet(nfabor))
    allocate(icvfli(nfabor))

    ! Map pointers to C
    p_icvfli = c_loc(icvfli)
    p_ifbet = c_loc(ifbet)

  end subroutine init_compf

  !> \brief Deallocate boundary flux indicators array

  subroutine finalize_compf

    implicit none

    deallocate(ifbet, icvfli)

  end subroutine finalize_compf

  !=============================================================================

end module cfpoin
