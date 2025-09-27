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

!> \file ppincl.f90
!> General module for specific physics

module ppincl

  !===========================================================================

  use, intrinsic :: iso_c_binding

  !=============================================================================

  implicit none

  !===========================================================================

  !> \defgroup ppincl General module for specific physics

  !> \addtogroup ppincl
  !> \{

  !> \defgroup choice Indicator table for specific physics

  !> \addtogroup choice
  !> \{

  !----------------------------------------------------------------------------
  !--> TABLEAU INDICATEURS DU CHOIX DE LA PHYSIQUE PARTICULIERE CHOISIE

  !> number of specific physics
  integer   nmodmx
  parameter(nmodmx = 15)

  !> global indicator for speciphic physics
  !> By default, all the indicators ippmod(i.....) are initialized to -1,
  !> which means that no specific physics is activated.
  !>
  !> Only the atmospheric flow module still includes Fortran calls and bindings.

  integer(c_int), pointer, save :: ippmod(:)

  !> pointer to specify atmospheric flow module with indicator ippmod(iatmos)
  !> - ippmod(iatmos) =-1 module not activated
  !> - ippmod(iatmos) = 0 standard modelling
  !> - ippmod(iatmos) = 1 dry atmosphere
  !> - ippmod(iatmos) = 2 humid atmosphere
  integer ::  iatmos

  parameter  (iatmos = 10)

  !> \}

  !=============================================================================

  !> \}

  !=============================================================================

  interface

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_physical_model_get_pointers(p_ippmod)    &
      bind(C, name='cs_f_physical_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_ippmod
    end subroutine cs_f_physical_model_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran physical models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine pp_models_init() &
     bind(C, name='cs_f_pp_models_init')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_ippmod

    call cs_f_physical_model_get_pointers(p_ippmod)
    call c_f_pointer(p_ippmod, ippmod, [nmodmx])

  end subroutine pp_models_init

  !=============================================================================

end module ppincl
