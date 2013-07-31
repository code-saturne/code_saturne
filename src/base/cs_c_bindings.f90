!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file cs_c_bindings.f90
!> Definition of structure, subroutine, and function C bindings.

module cs_c_bindings

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use field

  !=============================================================================

  integer :: MESH_LOCATION_NONE, MESH_LOCATION_CELLS
  integer :: MESH_LOCATION_INTERIOR_FACES, MESH_LOCATION_BOUNDARY_FACES
  integer :: MESH_LOCATION_VERTICES, MESH_LOCATION_PARTICLES
  integer :: MESH_LOCATION_OTHER

  parameter (MESH_LOCATION_NONE=0)
  parameter (MESH_LOCATION_CELLS=1)
  parameter (MESH_LOCATION_INTERIOR_FACES=2)
  parameter (MESH_LOCATION_BOUNDARY_FACES=3)
  parameter (MESH_LOCATION_VERTICES=4)
  parameter (MESH_LOCATION_PARTICLES=5)
  parameter (MESH_LOCATION_OTHER=6)

    !---------------------------------------------------------------------------

  type, bind(c)  :: var_cal_opt
    integer(c_int) :: iwarni
    integer(c_int) :: iconv
    integer(c_int) :: istat
    integer(c_int) :: idiff
    integer(c_int) :: idifft
    integer(c_int) :: idften
    integer(c_int) :: iswdyn
    integer(c_int) :: ischcv
    integer(c_int) :: isstpc
    integer(c_int) :: nswrgr
    integer(c_int) :: nswrsm
    integer(c_int) :: imligr
    integer(c_int) :: ircflu
    integer(c_int) :: inc
    real(c_double) :: thetav
    real(c_double) :: blencv
    real(c_double) :: epsilo
    real(c_double) :: epsrsm
    real(c_double) :: epsrgr
    real(c_double) :: climgr
    real(c_double) :: extrag
    real(c_double) :: relaxv
  end type var_cal_opt

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function logging field and other array statistics
    ! at relevant time steps.

    !> \brief Log field and other array statistics for a given time step.

    subroutine log_iteration()  &
      bind(C, name='cs_log_iteration')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine log_iteration

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function checking the presence of a control file
    ! and dealing with the interactive control.

    subroutine cs_control_check_file()  &
      bind(C, name='cs_control_check_file')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_control_check_file

    !---------------------------------------------------------------------------

    ! Interface to C function initializoing post-processing of moments

    subroutine cs_log_init_moments(cumulative_time)     &
      bind(C, name='cs_log_init_moments')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), dimension(*), intent(in) :: cumulative_time
    end subroutine cs_log_init_moments

    !---------------------------------------------------------------------------

    ! Interface to C function adding an array not saved as a permanent field
    ! to logging of fields

    subroutine cs_log_iteration_add_array(name, category, ml, is_intensive,    &
                                          dim, val)                            &
      bind(C, name='cs_log_iteration_add_array')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      character(kind=c_char, len=1), dimension(*), intent(in) :: category
      integer(c_int), value :: ml
      logical(c_bool), value :: is_intensive
      integer(c_int), value :: dim
      real(kind=c_double), dimension(*) :: val
    end subroutine cs_log_iteration_add_array

    !---------------------------------------------------------------------------

    ! Interface to C function adding an array not saved as a permanent field
    ! to logging of fields

    subroutine cs_log_iteration_clipping(name, dim, n_clip_min, n_clip_max,    &
                                         min_pre_clip, max_pre_clip)           &
      bind(C, name='cs_log_iteration_clipping')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      integer(c_int), value :: dim, n_clip_max, n_clip_min
      real(kind=c_double), dimension(*) :: min_pre_clip, max_pre_clip
    end subroutine cs_log_iteration_clipping

    !---------------------------------------------------------------------------

    ! Interface to C function adding an array not saved as a permanent field
    ! to logging of fields

    subroutine cs_log_iteration_clipping_field(f_id, n_clip_min, n_clip_max,   &
                                               min_pre_clip, max_pre_clip)     &
      bind(C, name='cs_log_iteration_clipping_field')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id, n_clip_max, n_clip_min
      real(kind=c_double), dimension(*) :: min_pre_clip, max_pre_clip
    end subroutine cs_log_iteration_clipping_field

    !---------------------------------------------------------------------------

    !> \endcond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Assign a var_cal_opt for a cs_var_cal_opt_t key to a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_value  structure associated with key

  subroutine field_set_key_struct_var_cal_opt (f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                   :: f_id
    type(var_cal_opt), intent(in), target :: k_value

    ! Local variables

    integer(c_int)                 :: c_f_id, c_k_id
    type(var_cal_opt),pointer      :: p_k_value
    type(c_ptr)                    :: c_k_value
    character(len=11+1, kind=c_char) :: c_name

    c_name = "var_cal_opt"//c_null_char
    c_k_id = cs_f_field_key_id_try(c_name)

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_set_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_set_key_struct_var_cal_opt

  !=============================================================================

  !> \brief Return a pointer to the var_cal_opt structure for cs_var_cal_opt key
  !> associated with a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[out]  k_value  integer value associated with key id for this field

  subroutine field_get_key_struct_var_cal_opt (f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                      :: f_id
    type(var_cal_opt), intent(inout), target :: k_value

    ! Local variables

    integer(c_int)                 :: c_f_id, c_k_id
    type(var_cal_opt),pointer      :: p_k_value
    type(c_ptr)                    :: c_k_value
    character(len=11+1, kind=c_char) :: c_name

    c_name = "var_cal_opt"//c_null_char
    c_k_id = cs_f_field_key_id_try(c_name)

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_get_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_get_key_struct_var_cal_opt

  !=============================================================================

  ! Interface to C function adding an array not saved as a permanent field
  ! to logging of fields

  !> \brief Add array not saved as permanent field to logging of fields.

  !> \param[in]  name         array name
  !> \param[in]  category     category name
  !> \param[in]  location     associated mesh location
  !> \param[in]  is_intensive associated mesh location
  !> \param[in]  dim          associated dimension (interleaved)
  !> \param[in]  val          associated values

  subroutine log_iteration_add_array(name, category, location, is_intensive,   &
                                     dim, val)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)      :: name, category
    integer, intent(in)               :: location, dim
    logical, intent(in)               :: is_intensive
    real(kind=c_double), dimension(*) :: val

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    character(len=len_trim(category)+1, kind=c_char) :: c_cat
    integer(c_int) :: c_ml, c_dim
    logical(c_bool) :: c_inten

    c_name = trim(name)//c_null_char
    c_cat = trim(category)//c_null_char
    c_ml = location
    c_inten = is_intensive
    c_dim = dim

    call cs_log_iteration_add_array(c_name, c_cat, c_ml, c_inten, c_dim, val)

    return

  end subroutine log_iteration_add_array

  !=============================================================================

  ! Interface to C function adding an array not saved as a permanent field
  ! to logging of fields

  !> \brief Add array not saved as permanent field to logging of fields.

  !> \param[in]  name          array name
  !> \param[in]  dim           associated dimension (interleaved)
  !> \param[in]  n_clip_min    local number of clipped to min values
  !> \param[in]  n_clip_max    local number of clipped to max values
  !> \param[in]  min_pre_clip  min local value prior to clip
  !> \param[in]  max_pre_clip  max local value prior to clip

  subroutine log_iteration_clipping(name, dim, n_clip_min, n_clip_max,        &
                                    min_pre_clip, max_pre_clip)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)      :: name
    integer, intent(in)               :: dim, n_clip_min, n_clip_max
    real(kind=c_double), dimension(*) :: min_pre_clip, max_pre_clip

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: c_dim, c_clip_min, c_clip_max

    c_name = trim(name)//c_null_char
    c_dim = dim
    c_clip_min = n_clip_min
    c_clip_max = n_clip_max

    call cs_log_iteration_clipping(c_name, c_dim, c_clip_min, c_clip_max, &
                                   min_pre_clip, max_pre_clip)

    return

  end subroutine log_iteration_clipping

  !=============================================================================

  ! Interface to C function adding an array not saved as a permanent field
  ! to logging of fields

  !> \brief Add array not saved as permanent field to logging of fields.

  !> \param[in]  f_id          associated dimension (interleaved)
  !> \param[in]  n_clip_min    local number of clipped to min values
  !> \param[in]  n_clip_max    local number of clipped to max values
  !> \param[in]  min_pre_clip  min local value prior to clip
  !> \param[in]  max_pre_clip  max local value prior to clip

  subroutine log_iteration_clipping_field(f_id, n_clip_min, n_clip_max,        &
                                          min_pre_clip, max_pre_clip)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)               :: f_id, n_clip_min, n_clip_max
    real(kind=c_double), dimension(*) :: min_pre_clip, max_pre_clip

    ! Local variables

    integer(c_int) :: c_f_id, c_clip_min, c_clip_max

    c_f_id = f_id
    c_clip_min = n_clip_min
    c_clip_max = n_clip_max

    call cs_log_iteration_clipping_field(c_f_id, c_clip_min, c_clip_max, &
                                         min_pre_clip, max_pre_clip)

    return

  end subroutine log_iteration_clipping_field

  !=============================================================================

end module cs_c_bindings
