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

!> \file post.f90
!> Module for post-processing related operations

module post

  !=============================================================================

  implicit none

  !=============================================================================

  integer :: POST_ON_LOCATION, POST_BOUNDARY_NR, POST_MONITOR

  parameter (POST_ON_LOCATION=1)
  parameter (POST_BOUNDARY_NR=2)
  parameter (POST_MONITOR=4)

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function outputting a floating point variable defined at
    ! cells or faces of post-processing mesh using associated writers.

    ! If the field id is not valid, a fatal error is provoked.

    subroutine cs_f_post_write_var(mesh_id, var_name, var_dim, interlace,  &
                                   use_parent, nt_cur_abs, t_cur_abs,      &
                                   cel_vals, i_face_vals, b_face_vals)     &
      bind(C, name='cs_f_post_write_var')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value                                   :: mesh_id
      character(kind=c_char, len=1), dimension(*), intent(in) :: var_name
      integer(c_int), value                                   :: var_dim
      logical(c_bool), value                                  :: interlace
      logical(c_bool), value                                  :: use_parent
      integer(c_int), value                                   :: nt_cur_abs
      real(c_double), value                                   :: t_cur_abs
      real(c_double), dimension(*), intent(in)                :: cel_vals
      real(c_double), dimension(*), intent(in)                :: i_face_vals
      real(c_double), dimension(*), intent(in)                :: b_face_vals
    end subroutine cs_f_post_write_var

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

contains

  !=============================================================================

  !> \brief  Output a variable defined at cells or faces of a post-processing
  !>         mesh using associated writers.

  !> \param[in]  mesh_id        id of associated mesh
  !> \param[in]  var_name       name of variable to output
  !> \param[in]  var_dim        1 for scalar, 3 for vector, 6/9 for tensor
  !> \param[in]  interleaved    .true. if values interleaved
  !>                            (ignored if scalar)
  !> \param[in]  use_parent     .true. if values are defined on "parent" mesh,
  !>                            .false. if values are defined directly on
  !>                            post-processing mesh
  !> \param[in]  nt_cur_abs     current time step number, or -1 if
  !>                            time-independent
  !> \param[in]  t_cur_abs      current physical time
  !> \param[in]  cel_vals       cell values array
  !> \param[in]  i_face_vals    interior face values array
  !> \param[in]  b_face_vals    boundary face values array

  subroutine post_write_var(mesh_id, var_name, var_dim, interleaved,  &
                            use_parent, nt_cur_abs, t_cur_abs,        &
                            cel_vals, i_face_vals, b_face_vals)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                      :: mesh_id
    character(len=*), intent(in)             :: var_name
    integer, intent(in)                      :: var_dim
    logical, intent(in)                      :: interleaved
    logical, intent(in)                      :: use_parent
    integer, intent(in)                      :: nt_cur_abs
    double precision, intent(in)             :: t_cur_abs
    real(c_double), dimension(*), intent(in) :: cel_vals
    real(c_double), dimension(*), intent(in) :: i_face_vals
    real(c_double), dimension(*), intent(in) :: b_face_vals

    ! Local variables

    character(len=len_trim(var_name)+1, kind=c_char) :: c_var_name
    integer(c_int) :: c_mesh_id
    integer(c_int) :: c_var_dim
    logical(c_bool) :: c_interleave
    logical(c_bool) :: c_use_parent
    integer(c_int) :: c_nt_cur_abs
    real(c_double) :: c_t_cur_abs

    c_mesh_id = mesh_id
    c_var_name = trim(var_name)//c_null_char
    c_var_dim = var_dim
    c_interleave = interleaved
    c_use_parent = use_parent
    c_nt_cur_abs = nt_cur_abs
    c_t_cur_abs = t_cur_abs

    call cs_f_post_write_var(c_mesh_id, c_var_name, c_var_dim, c_interleave,  &
                             c_use_parent, c_nt_cur_abs, c_t_cur_abs,         &
                             cel_vals, i_face_vals, b_face_vals)

  end subroutine post_write_var

  !=============================================================================

end module post
