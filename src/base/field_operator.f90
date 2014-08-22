!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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


!> \file field_operator.f90
!> Module for field-based algebraic operations

module field_operator

  !=============================================================================

  implicit none

  !=============================================================================

  integer :: FIELD_INTENSIVE, FIELD_EXTENSIVE
  integer :: FIELD_VARIABLE, FIELD_PROPERTY
  integer :: FIELD_POSTPROCESS, FIELD_ACCUMULATOR, FIELD_USER

  integer :: FIELD_OK, FIELD_INVALID_KEY_NAME, FIELD_INVALID_KEY_ID,   &
             FIELD_INVALID_CATEGORY, FIELD_INVALID_TYPE

  parameter (FIELD_INTENSIVE=1)
  parameter (FIELD_EXTENSIVE=2)
  parameter (FIELD_VARIABLE=4)
  parameter (FIELD_PROPERTY=8)
  parameter (FIELD_POSTPROCESS=16)
  parameter (FIELD_ACCUMULATOR=32)
  parameter (FIELD_USER=64)

  parameter (FIELD_OK=0)
  parameter (FIELD_INVALID_KEY_NAME=1)
  parameter (FIELD_INVALID_KEY_ID=2)
  parameter (FIELD_INVALID_CATEGORY=3)
  parameter (FIELD_INVALID_TYPE=4)

  !=============================================================================

  interface

    !> \brief  Compute cell gradient of scalar field or component of vector or
    !>         tensor field.

    !> \param[in]   f_id             field id
    !> \param[in]   use_previous_t   1 if values at previous time step should
    !>                               be used, 0 otherwise
    !> \param[in]   imrgra           gradient computation mode
    !> \param[in]   inc              0: increment; 1: do not increment
    !> \param[in]   recompute_cocg   1 or 0: recompute COCG or not
    !> \param[out]  grad             gradient

    subroutine field_gradient_scalar(f_id, use_previous_t, imrgra, inc,        &
                                     recompute_cocg,                           &
                                     grad)                                     &
      bind(C, name='cs_f_field_gradient_scalar')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: f_id, use_previous_t, imrgra, inc
      integer(c_int), value             :: recompute_cocg
      real(kind=c_double), dimension(*) :: grad
    end subroutine field_gradient_scalar

    !---------------------------------------------------------------------------

    !> \brief  Compute cell gradient of potential field

    !> \param[in]   f_id             field id
    !> \param[in]   use_previous_t   1 if values at previous time step should
    !>                               be used, 0 otherwise
    !> \param[in]   imrgra           gradient computation mode
    !> \param[in]   inc              0: increment; 1: do not increment
    !> \param[in]   recompute_cocg   1 or 0: recompute COCG or not
    !> \param[in]   hyd_p_flag       flag for hydrostatic pressure
    !> \param[in]   f_ext            exterior force generating
    !>                               the hydrostatic pressure
    !> \param[out]  grad             gradient

    subroutine field_gradient_potential(f_id, use_previous_t, imrgra, inc,     &
                                        recompute_cocg,                        &
                                        hyd_p_flag,                            &
                                        f_ext, grad)                           &
      bind(C, name='cs_f_field_gradient_potential')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value                :: f_id, use_previous_t, imrgra, inc
      integer(c_int), value                :: recompute_cocg
      integer(c_int), value                :: hyd_p_flag
      real(kind=c_double), dimension(3, *) :: f_ext
      real(kind=c_double), dimension(*)    :: grad
    end subroutine field_gradient_potential

    !---------------------------------------------------------------------------

    !> \brief  Compute cell gradient of vector field.

    !> \param[in]   f_id             field id
    !> \param[in]   use_previous_t   1 if values at previous time step should
    !>                               be used, 0 otherwise
    !> \param[in]   imrgra           gradient computation mode
    !> \param[in]   inc              0: increment; 1: do not increment
    !> \param[out]  grad             gradient

    subroutine field_gradient_vector(f_id, use_previous_t, imrgra, inc,        &
                                     grad)                                     &
      bind(C, name='cs_f_field_gradient_vector')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value                 :: f_id, use_previous_t, imrgra, inc
      real(kind=c_double), dimension(3,3,*) :: grad
    end subroutine field_gradient_vector

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

end module field_operator
