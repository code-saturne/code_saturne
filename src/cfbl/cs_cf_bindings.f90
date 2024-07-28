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

!> \file cs_cf_bindings.f90
!> Definition of C functions and subroutine bindings for compressible
!> flow module.

module cs_cf_bindings

  !=============================================================================

  use, intrinsic :: iso_c_binding

  implicit none

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function allowing to call the appropriate
    ! thermodynamical functions depending on the quantities provided
    ! by the user.

    subroutine cs_cf_thermo(iccfth, face_id, bc_en, bc_pr, bc_tk, bc_vel) &
      bind(C, name='cs_cf_thermo')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: iccfth, face_id
      real(kind=c_double), dimension(*) :: bc_en, bc_pr, bc_tk
      real(kind=c_double), dimension(3, *) :: bc_vel
    end subroutine cs_cf_thermo

    !---------------------------------------------------------------------------

    ! Interface to C function allowing to set the variability of the isobaric
    ! specific heat and the isochoric specific heat according to the
    ! thermodynamic law chosen by the user.

    subroutine cs_cf_set_thermo_options() &
      bind(C, name='cs_cf_set_thermo_options')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_cf_set_thermo_options

    !---------------------------------------------------------------------------

    ! Interface to C function initializing density, total energy and isochoric
    ! specific heat.

    subroutine cs_cf_thermo_default_init() &
      bind(C, name='cs_cf_thermo_default_init')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_cf_thermo_default_init

  end interface

  !=============================================================================

  end module cs_cf_bindings
