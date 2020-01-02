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

!> \file rotation.f90
!> \brief Module for rotating zones (rotors)

module rotation

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use paramx

  implicit none

  !=============================================================================
  !> \defgroup at_rotation Module for rotation zones (rotors)

  !> \addtogroup at_rotation
  !> \{

  !> Maximum number of rotors

  integer :: nrotmx

  parameter(nrotmx=20)

  !> Rotation axis

  double precision, save :: rotax(3, nrotmx)

  !> Rotation origin coordinates

  double precision, save :: rotcen(3, nrotmx)

  !> \}

  !=============================================================================

  interface

    !> \addtogroup at_rotation
    !> \{

    !---------------------------------------------------------------------------

    !> \brief Define a global rotation.

    !> \param[in]  omega_x      rotation vector x component
    !> \param[in]  omega_y      rotation vector y component
    !> \param[in]  omega_z      rotation vector z component
    !> \param[in]  invariant_x  invariant point x component
    !> \param[in]  invariant_y  invariant point y component
    !> \param[in]  invariant_z  invariant point z component

    subroutine rotation_define(omega_x, omega_y, omega_z,              &
                               invariant_x, invariant_y, invariant_z)  &
      bind(C, name='cs_rotation_define')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: omega_x, omega_y, omega_z
      real(c_double), value :: invariant_x, invariant_y, invariant_z
    end subroutine rotation_define

    !---------------------------------------------------------------------------

    !> \brief Update coordinates based on a global rotation and time.

    !> \param[in]       n_coords  number of coordinates
    !> \param[in]       t_rot     time since rotation start
    !> \param[in, out]  coords    coordinates array

    subroutine rotation_update_coords(n_coords, t_rot, coords)  &
      bind(C, name='cs_rotation_update_coords')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n_coords
      real(c_double), value :: t_rot
      real(c_double), dimension(3,*), intent(inout) :: coords
    end subroutine rotation_update_coords

    !---------------------------------------------------------------------------

    !> \brief Return angular velocity associated with a rotation.

    !> \param[in]    r_num   rotation number (0 for none, > 0 otherwise)
    !> \param[out]   omegea  angular velocity

    subroutine angular_velocity(r_num, omega)  &
      bind(C, name='cs_f_rotation_angular_velocity')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: r_num
      real(c_double), intent(out) :: omega
    end subroutine angular_velocity

    !---------------------------------------------------------------------------

    !> \brief Compute rotation velocity at given point coordinates
    !> \f[
    !> \vect{v}=\vect{\Omega}\wedge\vect{x},\quad
    !> \text{with }\vect{\Omega}\textnormal{ the rotation vector and }
    !> \vect{x}\text{ the point location.}
    !> \f]

    !> \param[in]    r_num   rotation number (0 for none, > 0 otherwise)
    !> \param[in]    coords  coordinates at point
    !> \param[out]   vr      rotation velocity

    subroutine rotation_velocity(r_num, coords, vr)  &
      bind(C, name='cs_f_rotation_velocity')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: r_num
      real(c_double), dimension(3), intent(in) :: coords
      real(c_double), dimension(3), intent(out) :: vr
    end subroutine rotation_velocity

    !---------------------------------------------------------------------------

    !> \brief Add a Coriolis term to a local vector
    !> \f[
    !> \vect{w}=\vect{w}+C\vect{\Omega}\wedge\vect{u},\quad
    !> \text{with }C\text{ a multiplicative coefficient, }
    !> \vect{\Omega}\text{ the rotation vector and }
    !> \vect{u}\text{ the velocity.}
    !> \f]

    !> \param[in]       r_num  rotation number (0 for none, > 0 otherwise)
    !> \param[in]       c      multiplicative coefficient
    !> \param[in]       v      velocity
    !> \param[in, out]  vr     resulting Coriolis term

    subroutine add_coriolis_v(r_num, c, v, vr) &
      bind(C, name='cs_f_rotation_add_coriolis_v')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: r_num
      real(c_double), value :: c
      real(c_double), dimension(3), intent(in) :: v
      real(c_double), dimension(3), intent(inout) :: vr
    end subroutine add_coriolis_v

    !---------------------------------------------------------------------------

    !> \brief Compute a Coriolis term for a vector
    !> \f[
    !> \vect{w}=C\vect{\Omega}\wedge\vect{u},\quad
    !> \text{with }C\text{ a multiplicative coefficient, }
    !> \vect{\Omega}\text{ the rotation vector and }
    !> \vect{u}\text{ the velocity.}
    !> \f]

    !> \param[in]   r_num  rotation number (0 for none, > 0 otherwise)
    !> \param[in]   c      multiplicative coefficient
    !> \param[in]   v      velocity
    !> \param[out]  vr     resulting Coriolis term

    subroutine coriolis_v(r_num, c, v, vr) &
      bind(C, name='cs_f_rotation_coriolis_v')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: r_num
      real(c_double), value :: c
      real(c_double), dimension(3), intent(in) :: v
      real(c_double), dimension(3), intent(out) :: vr
    end subroutine coriolis_v

    !---------------------------------------------------------------------------

    !> \brief Add the dual tensor of a rotation vector to a tensor
    !> \f[
    !> T_{ij}=T_{ij}+\epsilon_{imj}\Omega_m,\quad
    !> \text{for }\Omega_m\text{ the rotation vector components and }
    !> \epsilon_{imj}\text{ the permutation tensor components.}
    !> \f]

    !> \param[in]      r_num  rotation number (0 for none, > 0 otherwise)
    !> \param[in]      c      multiplicative coefficient
    !> \param[in,out]  tr     tensor to which dual tensor of rotation is added

    subroutine add_coriolis_t(r_num, c, tr) &
      bind(C, name='cs_f_rotation_add_coriolis_t')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: r_num
      real(c_double), value :: c
      real(c_double), dimension(3,3), intent(inout) :: tr
    end subroutine add_coriolis_t

    !---------------------------------------------------------------------------

    !> \brief Compute the dual tensor of a rotation vector
    !> \f[
    !> T_{ij}=\epsilon_{imj}\Omega_m,\quad
    !> \text{for }\Omega_m\text{ the rotation vector components and }
    !> \epsilon_{imj}\text{ the permutation tensor components.}
    !> \f]

    !> \param[in]   r_num  rotation number (0 for none, > 0 otherwise)
    !> \param[in]   c      multiplicative coefficient
    !> \param[out]  tr     dual tensor of rotation

    subroutine coriolis_t(r_num, c, tr) &
      bind(C, name='cs_f_rotation_coriolis_t')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: r_num
      real(c_double), value :: c
      real(c_double), dimension(3,3), intent(inout) :: tr
    end subroutine coriolis_t

    !> \}

    !---------------------------------------------------------------------------

    !> \brief Copy rotation structure values to an array

    !> This may be useful to avoid requiring specific type mappings for MPI or
    !> other programming languages.

    !> \param[in]   r_num  rotation number (0 for none, > 0 otherwise)
    !> \param[out]  fra    flat rotation array: axis (1-3), invariant(4-6),
    !>                     omega (7), angle(8)

    subroutine rotation_to_array(r_num, fra) &
      bind(C, name='cs_rotation_to_array')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: r_num
      real(c_double), dimension(8), intent(out) :: fra
    end subroutine rotation_to_array

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

!contains

  !=============================================================================


  !=============================================================================

end module rotation
