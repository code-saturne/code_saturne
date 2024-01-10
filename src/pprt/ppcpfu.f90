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

!> \file ppcpfu.f90
!> Module for specific physics for pulverized coal

module ppcpfu

  !=============================================================================

  use ppppar
  use ppthch

  implicit none

  !===========================================================================

  ! Oxydants
  integer   mnoxyd
  parameter(mnoxyd = 3)

  real(c_double), pointer, save :: oxyo2(:), oxyn2(:), oxyh2o(:), oxyco2(:)
  !
  !=============================================================================

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical model flags

    subroutine cs_f_ppcpfu_get_pointers(p_oxyo2, p_oxyn2, p_oxyh2o, p_oxyco2)  &
      bind(C, name='cs_f_ppcpfu_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: p_oxyo2, p_oxyn2, p_oxyh2o, p_oxyco2
    end subroutine cs_f_ppcpfu_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran combustion models properties API.
  !> This maps Fortran pointers to global C variables.

  subroutine ppcpfu_models_init() &
    bind(C, name='cs_f_ppcpfu_models_init')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: p_oxyo2, p_oxyn2, p_oxyh2o, p_oxyco2

    call cs_f_ppcpfu_get_pointers(p_oxyo2, p_oxyn2, p_oxyh2o, p_oxyco2)

    call c_f_pointer(p_oxyo2, oxyo2, [mnoxyd])
    call c_f_pointer(p_oxyn2, oxyn2, [mnoxyd])
    call c_f_pointer(p_oxyh2o, oxyh2o, [mnoxyd])
    call c_f_pointer(p_oxyco2, oxyco2, [mnoxyd])

  end subroutine ppcpfu_models_init

  !=============================================================================

end module ppcpfu
