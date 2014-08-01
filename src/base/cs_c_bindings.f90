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

!> \file cs_c_bindings.f90
!> Definition of structure, subroutine, and function C bindings.

module cs_c_bindings

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use field

  implicit none

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
    integer(c_int) :: imrgra
    integer(c_int) :: imligr
    integer(c_int) :: ircflu
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

    !> \brief  Read restart metadata.

    subroutine parameters_read_restart_info()  &
      bind(C, name='cs_parameters_read_restart_info')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine parameters_read_restart_info

    !---------------------------------------------------------------------------

    !> \brief  Update temporal moments.

    subroutine time_moment_update_all()  &
      bind(C, name='cs_time_moment_update_all')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine time_moment_update_all

    !---------------------------------------------------------------------------

    !> \brief  Log temporal moments initialization

    subroutine time_moment_log_iteration()  &
      bind(C, name='cs_time_moment_log_iteration')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine time_moment_log_iteration

    !---------------------------------------------------------------------------

    !> \brief  Get field id associated with a given moment.

    !> \param[in]   m_num   moment number (based on moment creation order,
    !>                      1 to n numbering)

    !> \return      f_id    associated field id

    function time_moment_field_id(m_num) result(f_id)  &
      bind(C, name='cs_f_time_moment_field_id')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: m_num
      integer(c_int)        :: f_id
    end function time_moment_field_id

    !---------------------------------------------------------------------------

    !> \brief  Read temporal moments checkpoint information.

    subroutine time_moment_restart_read(r_num)  &
      bind(C, name='cs_f_time_moment_restart_read')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: r_num
    end subroutine time_moment_restart_read

    !---------------------------------------------------------------------------

    !> \brief  Checkpoint temporal moments.

    subroutine time_moment_restart_write(r_num)  &
      bind(C, name='cs_f_time_moment_restart_write')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: r_num
    end subroutine time_moment_restart_write

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

    ! Interface to C function initializing gradient rotational periodicity
    ! computation API.

    subroutine cs_gradient_perio_initialize()  &
      bind(C, name='cs_gradient_perio_initialize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gradient_perio_initialize

    !---------------------------------------------------------------------------

    ! Interface to C function finalizing gradient rotational periodicity
    ! computation API.

    subroutine cs_gradient_perio_finalize()  &
      bind(C, name='cs_gradient_perio_finalize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gradient_perio_finalize

    !---------------------------------------------------------------------------

    !> \brief  Define user variables through the GUI.

    subroutine cs_gui_user_variables()  &
      bind(C, name='cs_gui_user_variables')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gui_user_variables

    !---------------------------------------------------------------------------

    ! Interface to C function adding an array not saved as a permanent field
    ! to logging of fields

    subroutine cs_log_iteration_add_array(name, category, ml, is_intensive,   &
                                          dim, val)                           &
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

    subroutine cs_log_iteration_clipping(name, dim, n_clip_min, n_clip_max,   &
                                         min_pre_clip, max_pre_clip)          &
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

    subroutine cs_log_iteration_clipping_field(f_id, n_clip_min, n_clip_max,  &
                                               min_pre_clip, max_pre_clip)    &
      bind(C, name='cs_log_iteration_clipping_field')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id, n_clip_max, n_clip_min
      real(kind=c_double), dimension(*) :: min_pre_clip, max_pre_clip
    end subroutine cs_log_iteration_clipping_field

    !---------------------------------------------------------------------------

    ! Interface to C function to compute properties with Freesteam in a
    ! defined thermal plane.

    subroutine phys_prop_freesteam(thermo_plane, property, n_vals,            &
                                   var1, var2, val)                           &
      bind(C, name='cs_phys_prop_freesteam')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: thermo_plane, property, n_vals
      real(kind=c_double), dimension(*), intent(in) :: var1, var2
      real(kind=c_double), dimension(*), intent(out) :: val
    end subroutine phys_prop_freesteam

    !---------------------------------------------------------------------------

    ! Interface to C function reading a cs_int_t vector section from a
    ! restart file, when that section may have used a different name and
    ! been non-interleaved in a previous version.

    subroutine cs_f_restart_read_int_t_compat(file_num, sec_name,            &
                                              old_name, location_id,         &
                                              n_location_vals,               &
                                              val, ierror)                   &
      bind(C, name='cs_f_restart_read_int_t_compat')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: file_num
      character(kind=c_char, len=1), dimension(*), intent(in) :: sec_name
      character(kind=c_char, len=1), dimension(*), intent(in) :: old_name
      integer(c_int), value :: location_id, n_location_vals
      integer(kind=c_int), dimension(*) :: val
      integer(c_int) :: ierror
    end subroutine cs_f_restart_read_int_t_compat

    !---------------------------------------------------------------------------

    ! Interface to C function reading a cs_real_t vector section from a
    ! restart file, when that section may have used a different name and
    ! been non-interleaved in a previous version.

    subroutine cs_f_restart_read_real_t_compat(file_num, sec_name,           &
                                               old_name, location_id,        &
                                               n_location_vals,              &
                                               val, ierror)                  &
      bind(C, name='cs_f_restart_read_real_t_compat')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: file_num
      character(kind=c_char, len=1), dimension(*), intent(in) :: sec_name
      character(kind=c_char, len=1), dimension(*), intent(in) :: old_name
      integer(c_int), value :: location_id, n_location_vals
      real(kind=c_double), dimension(*) :: val
      integer(c_int) :: ierror
    end subroutine cs_f_restart_read_real_t_compat

    !---------------------------------------------------------------------------

    ! Interface to C function reading a cs_real_3_t vector section from a
    ! restart file, when that section may have used a different name and
    ! been non-interleaved in a previous version.

    subroutine cs_f_restart_read_real_3_t_compat(file_num, sec_name,           &
                                                 old_name_x, old_name_y,       &
                                                 old_name_z, location_id,      &
                                                 val, ierror)                  &
      bind(C, name='cs_f_restart_read_real_3_t_compat')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: file_num
      character(kind=c_char, len=1), dimension(*), intent(in) :: sec_name
      character(kind=c_char, len=1), dimension(*), intent(in) :: old_name_x
      character(kind=c_char, len=1), dimension(*), intent(in) :: old_name_y
      character(kind=c_char, len=1), dimension(*), intent(in) :: old_name_z
      integer(c_int), value :: location_id
      real(kind=c_double), dimension(*) :: val
      integer(c_int) :: ierror
    end subroutine cs_f_restart_read_real_3_t_compat

    !---------------------------------------------------------------------------

    ! Interface to C function reading a cs_real_3_t vector section from a
    ! restart file, when that section may have used a different name and
    ! been non-interleaved in a previous version.

    function cs_f_restart_ptr(file_num) result(r)  &
      bind(C, name='cs_f_restart_ptr')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: file_num
      type(c_ptr)           :: r
    end function cs_f_restart_ptr

    !---------------------------------------------------------------------------

    ! Interface to C function reading BC coefficeints for all fields.

    subroutine cs_restart_read_bc_coeffs(restart)  &
      bind(C, name='cs_restart_read_bc_coeffs')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: restart
    end subroutine cs_restart_read_bc_coeffs

    !---------------------------------------------------------------------------

    ! Interface to C function writing BC coefficeints for all fields.

    subroutine cs_restart_write_bc_coeffs(restart)  &
      bind(C, name='cs_restart_write_bc_coeffs')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: restart
    end subroutine cs_restart_write_bc_coeffs

    !---------------------------------------------------------------------------

    ! Define a moment of a product of existing fields components.

    ! Moments will involve the tensor products of their component fields,
    ! and only scalar, vector, or rank-2 tensors are handled (for
    ! post-processing output reasons), so a moment may not involve more than
    ! 2 vectors or 1 tensor, unless single components are specified.

    function cs_f_time_moment_define_by_field_ids(name, n_fields,              &
                                                  field_id, component_id,      &
                                                  type, nt_start, t_start,     &
                                                  restart_mode_or_id)          &
      result(moment_id)                                                        &
      bind(C, name='cs_f_time_moment_define_by_field_ids')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      integer(c_int), value :: n_fields, restart_mode_or_id
      integer(kind=c_int), dimension(*) :: field_id, component_id
      integer(c_int), value :: type, nt_start
      real(kind=c_double), value :: t_start
      integer(c_int) :: moment_id
    end function cs_f_time_moment_define_by_field_ids

    !---------------------------------------------------------------------------

    ! Add terms from backward differentiation in time.

    subroutine cs_backward_differentiation_in_time(field_id,                  &
                                                   exp_part, imp_part)        &
      bind(C, name='cs_backward_differentiation_in_time')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: field_id
      real(kind=c_double), dimension(*), intent(inout) :: exp_part, imp_part
    end subroutine cs_backward_differentiation_in_time

    !---------------------------------------------------------------------------

    ! Interface to C user function for extra operations

    subroutine cs_user_extra_operations()  &
      bind(C, name='cs_user_extra_operations')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_user_extra_operations

    !---------------------------------------------------------------------------

    ! Interface to C user function for physical model options

    subroutine cs_user_model()  &
      bind(C, name='cs_user_model')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_user_model

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

  !---------------------------------------------------------------------------

  !> \brief Read a section of integers from a restart file,
  !> when that section may have used a different name in a previous version.

  !> \param[in]   f_num         restart file number
  !> \param[in]   sec_name      name of section
  !> \param[in]   old_name      old name of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[in]   n_loc_vals    number of valeus per location
  !> \param[out]  val           min local value prior to clip
  !> \param[out]  ierror        0: success, < 0: error code

  subroutine restart_read_int_t_compat(f_num, sec_name, old_name,            &
                                       location_id, n_loc_vals, val,         &
                                       ierror)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)               :: f_num
    character(len=*), intent(in)      :: sec_name
    character(len=*), intent(in)      :: old_name
    integer, intent(in)               :: location_id, n_loc_vals
    integer, dimension(*)             :: val
    integer, intent(out)              :: ierror

    ! Local variables

    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_n
    character(len=len_trim(old_name)+1, kind=c_char) :: c_o_n
    integer(c_int) :: c_f_num, c_loc_id, c_n_l_vals, c_ierror

    c_f_num = f_num
    c_s_n = trim(sec_name)//c_null_char
    c_o_n = trim(old_name)//c_null_char
    c_loc_id = location_id
    c_n_l_vals = n_loc_vals

    call cs_f_restart_read_int_t_compat(c_f_num, c_s_n, c_o_n,      &
                                        c_loc_id, c_n_l_vals, val,  &
                                        c_ierror)

    ierror = c_ierror

  end subroutine restart_read_int_t_compat

  !---------------------------------------------------------------------------

  !> \brief Read a section of double precision reals from a restart file,
  !> when that section may have used a different name in a previous version.

  !> \param[in]   f_num         restart file number
  !> \param[in]   sec_name      name of section
  !> \param[in]   old_name      old name of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[in]   n_loc_vals    number of valeus per location
  !> \param[out]  val           min local value prior to clip
  !> \param[out]  ierror        0: success, < 0: error code

  subroutine restart_read_real_t_compat(f_num, sec_name, old_name,           &
                                        location_id, n_loc_vals, val,        &
                                        ierror)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)               :: f_num
    character(len=*), intent(in)      :: sec_name
    character(len=*), intent(in)      :: old_name
    integer, intent(in)               :: location_id, n_loc_vals
    real(kind=c_double), dimension(*) :: val
    integer, intent(out)              :: ierror

    ! Local variables

    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_n
    character(len=len_trim(old_name)+1, kind=c_char) :: c_o_n
    integer(c_int) :: c_f_num, c_loc_id, c_n_l_vals, c_ierror

    c_f_num = f_num
    c_s_n = trim(sec_name)//c_null_char
    c_o_n = trim(old_name)//c_null_char
    c_loc_id = location_id
    c_n_l_vals = n_loc_vals

    call cs_f_restart_read_real_t_compat(c_f_num, c_s_n, c_o_n,      &
                                         c_loc_id, c_n_l_vals, val,  &
                                         c_ierror)

    ierror = c_ierror

  end subroutine restart_read_real_t_compat

  !---------------------------------------------------------------------------

  !> \brief Read a vector of double precision reals of dimension (3,*) from a
  !> restart file, when that section may have used a different name and
  !> been non-interleaved in a previous version.

  !> \param[in]   f_num         restart file number
  !> \param[in]   sec_name      name of section
  !> \param[in]   old_name_x    old name of component x of section
  !> \param[in]   old_name_y    old name of component y of section
  !> \param[in]   old_name_z    old name of component z of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[out]  val           min local value prior to clip
  !> \param[out]  ierror        0: success, < 0: error code

  subroutine restart_read_real_3_t_compat(f_num, sec_name,                     &
                                          old_name_x, old_name_y, old_name_z,  &
                                          location_id, val, ierror)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)               :: f_num
    character(len=*), intent(in)      :: sec_name
    character(len=*), intent(in)      :: old_name_x, old_name_y, old_name_z
    integer, intent(in)               :: location_id
    real(kind=c_double), dimension(*) :: val
    integer, intent(out)              :: ierror

    ! Local variables

    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_n
    character(len=len_trim(old_name_x)+1, kind=c_char) :: c_o_n_x
    character(len=len_trim(old_name_y)+1, kind=c_char) :: c_o_n_y
    character(len=len_trim(old_name_z)+1, kind=c_char) :: c_o_n_z
    integer(c_int) :: c_f_num, c_loc_id, c_ierror

    c_f_num = f_num
    c_s_n = trim(sec_name)//c_null_char
    c_o_n_x = trim(old_name_x)//c_null_char
    c_o_n_y = trim(old_name_y)//c_null_char
    c_o_n_z = trim(old_name_z)//c_null_char
    c_loc_id = location_id

    call cs_f_restart_read_real_3_t_compat(c_f_num, c_s_n, c_o_n_x, c_o_n_y,   &
                                           c_o_n_z, c_loc_id, val, c_ierror)

    ierror = c_ierror

  end subroutine restart_read_real_3_t_compat

  !---------------------------------------------------------------------------

  !> \brief Read boundary condition coefficients for all fields from
  !>        checkpoint.

  !> \param[in]   f_num     restart file number

  subroutine restart_read_bc_coeffs(f_num)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)   :: f_num

    ! Local variables

    integer(c_int) :: c_f_num
    type(c_ptr) :: r

    c_f_num = f_num
    r = cs_f_restart_ptr(c_f_num)

    call cs_restart_read_bc_coeffs(r)

  end subroutine restart_read_bc_coeffs

  !---------------------------------------------------------------------------

  !> \brief Write boundary condition coefficients for all fields to
  !>        checkpoint.

  !> \param[in]   f_num     restart file number

  subroutine restart_write_bc_coeffs(f_num)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)   :: f_num

    ! Local variables

    integer(c_int) :: c_f_num
    type(c_ptr) :: r

    c_f_num = f_num
    r = cs_f_restart_ptr(c_f_num)

    call cs_restart_write_bc_coeffs(r)

  end subroutine restart_write_bc_coeffs

  !---------------------------------------------------------------------------

  !> \brief Define a moment of a product of existing fields components.

  !> Moment ids are based on the moment creation order, starting at 0.

  !> Moments will involve the tensor products of their component fields,
  !> and only scalar, vector, or rank-2 tensors are handled (for
  !> post-processing output reasons), so a moment may not involve more than
  !> 2 vectors or 1 tensor, unless single components are specified.
  !>
  !> \param[in]   name          name of associated moment
  !> \param[in]   n_fields      number of associated fields
  !> \param[in]   field_id      ids of associated fields
  !> \param[in]   component_id  ids of matching field components (-1 for all)
  !> \param[in]   type          moment type
  !> \param[in]   nt_start      starting time step (or -1 to use t_start)
  !> \param[in]   t_start       starting time
  !> \param[in]   restart_id    -2: automatic, -1: reset, >= 0: id of
  !>                            matching moment in restart data
  !> \param[out]  moment_id     id of new moment in case of success,
  !>                            -1 in case of error.

  subroutine time_moment_define_by_field_ids(name, n_fields,                   &
                                             field_id, component_id,           &
                                             type, nt_start, t_start,          &
                                             restart_id,                       &
                                             moment_id)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)             :: name
    integer, intent(in)                      :: n_fields
    integer(c_int), intent(in), dimension(*) :: field_id, component_id
    integer, intent(in)                      :: type, nt_start, restart_id
    double precision, intent(in)             :: t_start
    integer, intent(out)                     :: moment_id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: c_n_fields, c_type, c_nt_start, c_moment_id
    integer(c_int) :: c_restart_id
    real(c_double) :: c_t_start

    c_name = trim(name)//c_null_char
    c_n_fields = n_fields
    c_type = type
    c_nt_start = nt_start
    c_t_start = t_start
    c_restart_id = restart_id

    c_moment_id = cs_f_time_moment_define_by_field_ids(c_name, c_n_fields,    &
                                                       field_id,              &
                                                       component_id, c_type,  &
                                                       c_nt_start, c_t_start, &
                                                       c_restart_id)

    moment_id = c_moment_id

  end subroutine time_moment_define_by_field_ids

  !=============================================================================

  end module cs_c_bindings
