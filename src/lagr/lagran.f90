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

!> \file lagran.f90
!> \brief Module for Lagrangian model.

module lagran

  !===========================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !> \defgroup lagran Module for Lagrangian model

  !> \addtogroup lagran
  !> \{

  !=============================================================================

  !> \defgroup base Base

  !> \addtogroup base
  !> \{

  !> \anchor iilagr
  !> activates (>0) or deactivates (=0) the Lagrangian module
  !> the different values correspond to the following modellings:
  !> - = 1 Lagrangian two-phase flow in one-way coupling (no influence of
  !> the particles on the continuous phase)
  !> - = 2 Lagrangian two-phase flow with two-way coupling (influence of
  !> the particles on the dynamics of the continuous phase).
  !> Dynamics, temperature and mass may be coupled independently.
  !> - = 3 Lagrangian two-phase flow on frozen continuous phase. This option can
  !> only be used in case of a calculation restart. All the
  !> Eulerian fields are frozen (including the scalar fields). This option
  !> automatically implies \ref iccvfg = 1
  integer(c_int), pointer, save :: iilagr

  !> - 0: no head losses calculation for influence of the deposit on the flow
  !> - 1: head losses calculation for influence of the deposit on the flow
  integer(c_int), pointer, save :: iflow

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function returning particle attribute pointers

    subroutine cs_f_lagr_params_pointers(p_iilagr, p_iflow) &
      bind(C, name='cs_f_lagr_params_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_iilagr, p_iflow
    end subroutine cs_f_lagr_params_pointers

  !=============================================================================

  end interface

contains

  !=============================================================================

  subroutine cs_f_lagran_init_map () &
    bind(C, name='cs_f_lagran_init_map')

    implicit none

    ! lagr_option_base
    type(c_ptr) :: p_iilagr

    ! lagr_params
    type(c_ptr) :: p_iflow

    call cs_f_lagr_params_pointers(p_iilagr, p_iflow)
    call c_f_pointer(p_iilagr, iilagr)
    call c_f_pointer(p_iflow , iflow)

    return

  end subroutine cs_f_lagran_init_map

  !=============================================================================

end module lagran
