!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

  use field

  !=============================================================================

  implicit none

  !=============================================================================

  interface

    !> \brief  Compute cell gradient of scalar field or component of vector or
    !>         tensor field.

    !> \param[in]   f_id             field id
    !> \param[in]   use_previous_t   1 if values at previous time step should
    !>                               be used, 0 otherwise
    !> \param[in]   imrgrO           unused
    !> \param[in]   inc              0: increment; 1: do not increment
    !> \param[in]   recompute_cocg   1 or 0: recompute COCG or not
    !> \param[out]  grad             gradient

    subroutine field_gradient_scalar(f_id, use_previous_t, imrgr0, inc,        &
                                     recompute_cocg,                           &
                                     grad)                                     &
      bind(C, name='cs_f_field_gradient_scalar')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: f_id, use_previous_t, imrgr0, inc
      integer(c_int), value             :: recompute_cocg
      real(kind=c_double), dimension(*) :: grad
    end subroutine field_gradient_scalar

    !---------------------------------------------------------------------------

    !> \brief  Compute cell gradient of potential field

    !> \param[in]   f_id             field id
    !> \param[in]   use_previous_t   1 if values at previous time step should
    !>                               be used, 0 otherwise
    !> \param[in]   imrgr0           ignored
    !> \param[in]   inc              0: increment; 1: do not increment
    !> \param[in]   recompute_cocg   1 or 0: recompute COCG or not
    !> \param[in]   hyd_p_flag       flag for hydrostatic pressure
    !> \param[in]   f_ext            exterior force generating
    !>                               the hydrostatic pressure
    !> \param[out]  grad             gradient

    subroutine field_gradient_potential(f_id, use_previous_t, imrgr0, inc,     &
                                        recompute_cocg,                        &
                                        hyd_p_flag,                            &
                                        f_ext, grad)                           &
      bind(C, name='cs_f_field_gradient_potential')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value                :: f_id, use_previous_t, imrgr0, inc
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
    !> \param[in]   imrgr0           ignored
    !> \param[in]   inc              0: increment; 1: do not increment
    !> \param[out]  grad             gradient

    subroutine field_gradient_vector(f_id, use_previous_t, imrgr0, inc,        &
                                     grad)                                     &
      bind(C, name='cs_f_field_gradient_vector')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value                 :: f_id, use_previous_t, imrgr0, inc
      real(kind=c_double), dimension(3,3,*) :: grad
    end subroutine field_gradient_vector

    !---------------------------------------------------------------------------

    !> \brief  Compute cell gradient of tensor field.

    !> \param[in]   f_id             field id
    !> \param[in]   use_previous_t   1 if values at previous time step should
    !>                               be used, 0 otherwise
    !> \param[in]   imrgr0           ignored
    !> \param[in]   inc              0: increment; 1: do not increment
    !> \param[out]  grad             gradient

    subroutine field_gradient_tensor(f_id, use_previous_t, imrgr0, inc,        &
                                     grad)                                     &
      bind(C, name='cs_f_field_gradient_tensor')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value                 :: f_id, use_previous_t, imrgr0, inc
      real(kind=c_double), dimension(6,3,*) :: grad
    end subroutine field_gradient_tensor

    !---------------------------------------------------------------------------

    !> \brief  Shift field values in order to set its spatial average to a given
    !>         value.

    !> \param[in]   f_id  field id
    !> \param[in]   va    real value of volume average to be set

    subroutine field_set_volume_average(f_id, va)                              &
      bind(C, name='cs_f_field_set_volume_average')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value      :: f_id
      real(kind=c_double), value :: va
    end subroutine field_set_volume_average

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

end module field_operator
