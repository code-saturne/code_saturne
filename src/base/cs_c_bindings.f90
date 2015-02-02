!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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
!> Definition of C function and subroutine bindings.

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

  integer :: RESTART_VAL_TYPE_INT_T, RESTART_VAL_TYPE_REAL_T

  parameter (MESH_LOCATION_NONE=0)
  parameter (MESH_LOCATION_CELLS=1)
  parameter (MESH_LOCATION_INTERIOR_FACES=2)
  parameter (MESH_LOCATION_BOUNDARY_FACES=3)
  parameter (MESH_LOCATION_VERTICES=4)
  parameter (MESH_LOCATION_PARTICLES=5)
  parameter (MESH_LOCATION_OTHER=6)

  parameter (RESTART_VAL_TYPE_INT_T=1)
  parameter (RESTART_VAL_TYPE_REAL_T=3)

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
    integer(c_int) :: ibdtso
    integer(c_int) :: isstpc
    integer(c_int) :: nswrgr
    integer(c_int) :: nswrsm
    integer(c_int) :: imrgra
    integer(c_int) :: imligr
    integer(c_int) :: ircflu
    integer(c_int) :: iwgrec
    real(c_double) :: thetav
    real(c_double) :: blencv
    real(c_double) :: epsilo
    real(c_double) :: epsrsm
    real(c_double) :: epsrgr
    real(c_double) :: climgr
    real(c_double) :: extrag
    real(c_double) :: relaxv
  end type var_cal_opt

  !---------------------------------------------------------------------------

  type, bind(c)  :: solving_info
    integer(c_int) :: nbivar
    real(c_double) :: rnsmbr
    real(c_double) :: resvar
    real(c_double) :: dervar
    real(c_double) :: l2residual
  end type solving_info

  !---------------------------------------------------------------------------

  type, bind(c)  :: gas_mix_species_prop
    real(c_double) :: mol_mas
    real(c_double) :: cp
    real(c_double) :: vol_dif
    real(c_double) :: mu_a
    real(c_double) :: mu_b
    real(c_double) :: lambda_a
    real(c_double) :: lambda_b
  end type gas_mix_species_prop

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    ! Interface to C function logging field and other array statistics
    ! at relevant time steps.

    ! \brief Log field and other array statistics for a given time step.

    subroutine log_iteration()  &
      bind(C, name='cs_log_iteration')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine log_iteration

    !---------------------------------------------------------------------------

    !> \brief Compute filters for dynamic models.

    !> \param[in]   dim            stride of array to filter
    !> \param[in]   val            array of values to filter
    !> \param[out]  f_val          array of filtered values

    subroutine les_filter(stride, val, f_val)  &
      bind(C, name='cs_les_filter')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: stride
      real(kind=c_double), dimension(*) :: val
      real(kind=c_double), dimension(*), intent(out) :: f_val
    end subroutine les_filter

    !---------------------------------------------------------------------------

    !> \brief  Destroy name to id map structure.

    !> \param[in, out] m pointer to map structure

    subroutine cs_map_name_to_id_destroy(m)  &
      bind(C, name='cs_map_name_to_id_destroy')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: m
    end subroutine cs_map_name_to_id_destroy

    !---------------------------------------------------------------------------

    !> \brief  Read restart metadata.

    subroutine parameters_read_restart_info()  &
      bind(C, name='cs_parameters_read_restart_info')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine parameters_read_restart_info

    !---------------------------------------------------------------------------

    !> \brief  Destroy structure associated with a restart file
    !>         (and close the file).

    !> \param[in, out] r pointer to map structure

    subroutine restart_destroy(r)  &
      bind(C, name='cs_restart_destroy')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(inout) :: r
    end subroutine restart_destroy

    !---------------------------------------------------------------------------

    !> \brief  Check the locations associated with a restart file.

    !> For each type of entity, return .true. if the associated number
    !> of entities matches the current value (and so that we consider the
    !> mesh locations, false otherwise.

    !> \param[in]   r     restart structure pointer
    !> \param[out]  lcel  match for cells
    !> \param[out]  lfac  match for interior faces
    !> \param[out]  lfbr  match for boundary faces
    !> \param[out]  lsom  match for vertices

    subroutine restart_check_base_location(r, lcel, lfac, lfbr, lsom)  &
      bind(C, name='cs_restart_check_base_location')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      logical(kind=c_bool), intent(out) :: lcel, lfac, lfbr, lsom
    end subroutine restart_check_base_location

    !---------------------------------------------------------------------------

    !> \brief Read field metadata from checkpoint.

    !> \param[in]   r              restart structure pointer
    !> \param[in]   old_field_map  old field map pointer

    subroutine restart_read_field_info(r, old_field_map)  &
      bind(C, name='cs_restart_read_field_info')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      type(c_ptr), intent(out) :: old_field_map
    end subroutine restart_read_field_info

    !---------------------------------------------------------------------------

    !> \brief Write field metadata to checkpoint.

    !> \param[in]   r  restart structure pointer

    subroutine restart_write_field_info(r)  &
      bind(C, name='cs_restart_write_field_info')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
    end subroutine restart_write_field_info

    !---------------------------------------------------------------------------

    !> \brief Read boundary condition coefficients for all fields from
    !>        checkpoint.

    !> \param[in]   r  pointer to restart structure

    subroutine restart_read_bc_coeffs(r)  &
      bind(C, name='cs_restart_read_bc_coeffs')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
    end subroutine restart_read_bc_coeffs

    !---------------------------------------------------------------------------

    !> \brief Write boundary condition coefficients for all fields to
    !>        checkpoint.

    !> \param[in]   r  pointer to restart structure

    subroutine restart_write_bc_coeffs(r)  &
      bind(C, name='cs_restart_write_bc_coeffs')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
    end subroutine restart_write_bc_coeffs

    !---------------------------------------------------------------------------

    !> \brief  Return the number of temporal moments.

    !> \return number of defined moments

    function cs_time_moment_n_moments() result(n_moments)  &
      bind(C, name='cs_time_moment_n_moments')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: n_moments
    end function cs_time_moment_n_moments

    !---------------------------------------------------------------------------

    !> \brief  Return if moment is active (1) or not (0).

    !> \return 1 if moment is active, 0 if not

    function cs_time_moment_is_active(m_id) result(is_active)  &
      bind(C, name='cs_time_moment_is_active')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: m_id
      integer(c_int) :: is_active
    end function cs_time_moment_is_active

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

    !> For moments not defined by the user, but defined automatically so as
    !> to allow computation of higher order moments (i.e. variances), no field
    !> is associated, so the returned value is -1.

    !> \param[in]   m_num   moment number (based on moment creation order,
    !>                      1 to n numbering)

    !> \return      f_id    associated field id, or -1

    function time_moment_field_id(m_num) result(f_id)  &
      bind(C, name='cs_f_time_moment_field_id')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: m_num
      integer(c_int)        :: f_id
    end function time_moment_field_id

    !---------------------------------------------------------------------------

    !> \brief  Read temporal moments checkpoint information.

    subroutine time_moment_restart_read(r)  &
      bind(C, name='cs_time_moment_restart_read')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
    end subroutine time_moment_restart_read

    !---------------------------------------------------------------------------

    !> \brief  Checkpoint temporal moments.

    subroutine time_moment_restart_write(r)  &
      bind(C, name='cs_time_moment_restart_write')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
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

    ! Interface to C function returning the global dot product of 2 vectors

    function cs_gdot(n, x, y) result(gdot) &
      bind(C, name='cs_gdot')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n
      real(kind=c_double), dimension(*), intent(in) :: x, y
      real(kind=c_double) :: gdot
    end function cs_gdot

    !---------------------------------------------------------------------------

    ! Interface to C function returning the global residual of 2 vectors

    function cs_gres(n, vol, x, y) result(gres) &
      bind(C, name='cs_gres')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: n
      real(kind=c_double), dimension(*), intent(in) :: vol, x, y
      real(kind=c_double) :: gres
    end function cs_gres

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

    ! Interface to C function initializing ghost cell values
    ! for Reynolds stress tensor gradient.

    subroutine cs_gradient_perio_init_rij(f, idimtr, grad) &
      bind(C, name='cs_gradient_perio_init_rij')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value                                  :: f
      integer(c_int), intent(out)                         :: idimtr
      real(kind=c_double), dimension(3, *), intent(inout) :: grad
    end subroutine cs_gradient_perio_init_rij

    !---------------------------------------------------------------------------

    !> \brief  Define user variables through the GUI.

    subroutine cs_gui_user_variables()  &
      bind(C, name='cs_gui_user_variables')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gui_user_variables

    !---------------------------------------------------------------------------

    !> \brief  Define time moments through the GUI.

    subroutine cs_gui_time_moments()  &
      bind(C, name='cs_gui_time_moments')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gui_time_moments

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

    ! Interface to C function initializing codensation-related field key.

    subroutine cs_parameters_define_field_key_gas_mix()  &
      bind(C, name='cs_parameters_define_field_key_gas_mix')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_parameters_define_field_key_gas_mix

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

    ! Interface to C function initializing a restart file

    function cs_restart_create(name, path, mode) result(r) &
      bind(C, name='cs_restart_create')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name, path
      integer(c_int), value :: mode
      type(c_ptr) :: r
    end function cs_restart_create

    !---------------------------------------------------------------------------

    ! Interface to C function reading a section from a restart file.

    function cs_restart_read_section(r, sec_name,                           &
                                     location_id, n_location_vals,          &
                                     val_type, val) result(error)           &
      bind(C, name='cs_restart_read_section')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      character(kind=c_char, len=1), dimension(*), intent(in) :: sec_name
      integer(c_int), value :: location_id, n_location_vals, val_type
      type(c_ptr), value :: val
      integer(c_int) :: error
    end function cs_restart_read_section

    !---------------------------------------------------------------------------

    ! Interface to C function reading a section from a restart file, when
    ! that section may have used a different name in a previous version.

    function cs_restart_read_section_compat(r, sec_name, old_name,          &
                                            location_id, n_location_vals,   &
                                            val_type, val) result(error)    &
      bind(C, name='cs_restart_read_section_compat')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      character(kind=c_char, len=1), dimension(*), intent(in) :: sec_name
      character(kind=c_char, len=1), dimension(*), intent(in) :: old_name
      integer(c_int), value :: location_id, n_location_vals, val_type
      type(c_ptr), value :: val
      integer(c_int) :: error
    end function cs_restart_read_section_compat

    !---------------------------------------------------------------------------

    ! Interface to C function writing a section to a checkpoint file.

    subroutine cs_restart_write_section(r, sec_name,                        &
                                        location_id, n_location_vals,       &
                                        val_type, val)                      &
      bind(C, name='cs_restart_write_section')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      character(kind=c_char, len=1), dimension(*), intent(in) :: sec_name
      integer(c_int), value :: location_id, n_location_vals, val_type
      type(c_ptr), value :: val
      integer(c_int) :: error
    end subroutine cs_restart_write_section

    !---------------------------------------------------------------------------

    ! Interface to C function reading variables

    subroutine cs_restart_read_variables(r, old_field_map,                   &
                                         t_id_flag, read_flag)               &
      bind(C, name='cs_restart_read_variables')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      type(c_ptr), value :: old_field_map
      integer(kind=c_int), value :: t_id_flag
      type(c_ptr), value :: read_flag
      ! integer(kind=c_int), dimension(*) :: read_flag ! (swap below to use)
    end subroutine cs_restart_read_variables

    !---------------------------------------------------------------------------

    ! Interface to C function writing variables

    subroutine cs_restart_write_variables(r, t_id_flag, write_flag)          &
      bind(C, name='cs_restart_write_variables')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      integer(kind=c_int), value :: t_id_flag
      type(c_ptr), value :: write_flag
      ! integer(kind=c_int), dimension(*) :: write_flag ! (swap below to use)
    end subroutine cs_restart_write_variables

    !---------------------------------------------------------------------------

    ! Interface to C function reading a cs_real_3_t vector section from a
    ! restart file, when that section may have used a different name and
    ! been non-interleaved in a previous version.

    function cs_restart_read_real_3_t_compat(r, sec_name,                     &
                                             old_name_x, old_name_y,          &
                                             old_name_z, location_id,         &
                                             val) result(ierror)              &
      bind(C, name='cs_restart_read_real_3_t_compat')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      character(kind=c_char, len=1), dimension(*), intent(in) :: sec_name
      character(kind=c_char, len=1), dimension(*), intent(in) :: old_name_x
      character(kind=c_char, len=1), dimension(*), intent(in) :: old_name_y
      character(kind=c_char, len=1), dimension(*), intent(in) :: old_name_z
      integer(c_int), value :: location_id
      real(kind=c_double), dimension(*) :: val
      integer(c_int) :: ierror
    end function cs_restart_read_real_3_t_compat

    !---------------------------------------------------------------------------

    ! Interface to C function reading field values from a restart file,
    ! when that section may have used a different name and
    ! been non-interleaved in a previous version.

    function cs_restart_read_field_vals(r, f_id, t_id) result(ierr)  &
      bind(C, name='cs_restart_read_field_vals')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: r
      integer(c_int), value :: f_id, t_id
      integer(c_int)        :: ierr
    end function cs_restart_read_field_vals

    !---------------------------------------------------------------------------

    ! Interface to C function writing field values to a restart file.

    subroutine cs_restart_write_field_vals(r, f_id, t_id)  &
      bind(C, name='cs_restart_write_field_vals')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value    :: r
      integer(c_int), value :: f_id, t_id
    end subroutine cs_restart_write_field_vals

    !---------------------------------------------------------------------------

    ! Interface to C function reading fields depending on others from checkpoint

    function cs_restart_read_linked_fields(r, old_field_map, key, read_flag) &
      result(n)  &
      bind(C, name='cs_restart_read_linked_fields')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      type(c_ptr), value :: old_field_map
      character(kind=c_char, len=1), dimension(*), intent(in) :: key
      ! integer(kind=c_int), dimension(*) :: read_flag ! (swap below to use)
      type(c_ptr), value :: read_flag
      integer(c_int)     :: n
    end function cs_restart_read_linked_fields

    !---------------------------------------------------------------------------

    ! Interface to C function writing fields depending on others to a checkpoint

    function cs_restart_write_linked_fields(r, key, write_flag) result(n)  &
      bind(C, name='cs_restart_write_linked_fields')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      character(kind=c_char, len=1), dimension(*), intent(in) :: key
      ! integer(kind=c_int), dimension(*) :: write_flag ! (swap below to use)
      type(c_ptr), value :: write_flag
      integer(c_int)     :: n
    end function cs_restart_write_linked_fields

    !---------------------------------------------------------------------------

    ! Interface to C function calling sparse linear equation solver
    ! using native matrix arrays.

    function cs_sles_solve_native(f_id, name, symmetric,                      &
                                  diag_block_size, extra_diag_block_size,     &
                                  da, xa, rotation_mode, precision, r_norm,   &
                                  n_iter, residue, rhs, vx) result(state)     &
      bind(C, name='cs_sles_solve_native')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      logical(kind=c_bool), value :: symmetric
      integer(c_int), value :: rotation_mode
      integer(c_int), dimension(*) :: diag_block_size, extra_diag_block_size
      real(kind=c_double), value :: precision, r_norm
      integer(c_int), intent(out) :: n_iter
      real(kind=c_double), intent(out) :: residue
      real(kind=c_double), dimension(*), intent(in) :: da, xa, rhs
      real(kind=c_double), dimension(*), intent(inout) :: vx
      integer(c_int) :: state
    end function cs_sles_solve_native

    !---------------------------------------------------------------------------

    ! Interface to C function freeing sparse linear equation solver setup
    ! using native matrix arrays.

    subroutine cs_sles_free_native(f_id, name)                                &
      bind(C, name='cs_sles_free_native')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
    end subroutine cs_sles_free_native

    !---------------------------------------------------------------------------

    ! Temporarily replace field id with name for matching calls
    ! to cs_sles_solve_native.

    subroutine cs_sles_push(f_id, name)                                       &
      bind(C, name='cs_sles_push')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
    end subroutine cs_sles_push

    !---------------------------------------------------------------------------

    ! Revert to normal behavior of field id for matching calls
    ! to cs_sles_solve_native.

    subroutine cs_sles_pop(f_id)                                             &
      bind(C, name='cs_sles_pop')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
    end subroutine cs_sles_pop

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
    ! Interface to C function for balance computation

    subroutine cs_balance_by_zone(itypfb, selection_crit, scalar_name)  &
      bind(C, name='cs_balance_by_zone')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), dimension(*), intent(in) :: itypfb
      character(kind=c_char, len=1), dimension(*), intent(in) :: selection_crit
      character(kind=c_char, len=1), dimension(*), intent(in) :: scalar_name
    end subroutine cs_balance_by_zone

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

    !> Interface to C user function for time moments

    subroutine cs_user_time_moments()  &
      bind(C, name='cs_user_time_moments')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_user_time_moments

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

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

    integer(c_int)                 :: c_f_id
    type(var_cal_opt),pointer      :: p_k_value
    type(c_ptr)                    :: c_k_value
    character(len=11+1, kind=c_char) :: c_name

    integer(c_int), save           :: c_k_id = -1

    if (c_k_id .eq. -1) then
      c_name = "var_cal_opt"//c_null_char
      c_k_id = cs_f_field_key_id(c_name)
    endif

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_set_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_set_key_struct_var_cal_opt

  !=============================================================================

  !> \brief Assign a solving_info for a cs_solving_info_t key to a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_value  structure associated with key

  subroutine field_set_key_struct_solving_info (f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                    :: f_id
    type(solving_info), intent(in), target :: k_value

    ! Local variables

    integer(c_int)                   :: c_f_id
    type(solving_info), pointer      :: p_k_value
    type(c_ptr)                      :: c_k_value
    character(len=12+1, kind=c_char) :: c_name

    integer(c_int), save           :: c_k_id = -1

    if (c_k_id .eq. -1) then
      c_name = "solving_info"//c_null_char
      c_k_id = cs_f_field_key_id(c_name)
    endif

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_set_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_set_key_struct_solving_info

  !=============================================================================

  !> \brief Assign a gas_mix_species_prop for a cs_gas_mix_species_prop_t
  !> key to a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_value  structure associated with key

  subroutine field_set_key_struct_gas_mix_species_prop(f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                               :: f_id
    type(gas_mix_species_prop), intent(in), target :: k_value

    ! Local variables

    integer(c_int)                             :: c_f_id
    type(gas_mix_species_prop),pointer      :: p_k_value
    type(c_ptr)                                :: c_k_value
    character(len=23+1, kind=c_char)           :: c_name

    integer(c_int), save           :: c_k_id = -1

    if (c_k_id .eq. -1) then
      c_name = "gas_mix_species_prop"//c_null_char
      c_k_id = cs_f_field_key_id(c_name)
    endif

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_set_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_set_key_struct_gas_mix_species_prop

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

    integer(c_int)                 :: c_f_id
    type(var_cal_opt),pointer      :: p_k_value
    type(c_ptr)                    :: c_k_value
    character(len=11+1, kind=c_char) :: c_name

    integer(c_int), save           :: c_k_id = -1

    if (c_k_id .eq. -1) then
      c_name = "var_cal_opt"//c_null_char
      c_k_id = cs_f_field_key_id(c_name)
    endif

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_get_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_get_key_struct_var_cal_opt

  !=============================================================================

  !> \brief Return a pointer to the solving_info structure for
  !>        cs_solving_info_t key associated with a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[out]  k_value  integer value associated with key id for this field

  subroutine field_get_key_struct_solving_info (f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                       :: f_id
    type(solving_info), intent(inout), target :: k_value

    ! Local variables

    integer(c_int)                   :: c_f_id, c_k_id
    type(solving_info), pointer      :: p_k_value
    type(c_ptr)                      :: c_k_value
    character(len=12+1, kind=c_char) :: c_name

    c_name = "solving_info"//c_null_char
    c_k_id = cs_f_field_key_id(c_name)

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_get_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_get_key_struct_solving_info

  !=============================================================================

  !> \brief Return a pointer to the gas_mix_species_prop structure for
  !>        cs_gas_mix_species_prop_t key associated with a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[out]  k_value  integer value associated with key id for this field

  subroutine field_get_key_struct_gas_mix_species_prop (f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                                  :: f_id
    type(gas_mix_species_prop), intent(inout), target :: k_value

    ! Local variables

    integer(c_int)                             :: c_f_id, c_k_id
    type(gas_mix_species_prop),pointer      :: p_k_value
    type(c_ptr)                                :: c_k_value
    character(len=23+1, kind=c_char)           :: c_name

    c_name = "gas_mix_species_prop"//c_null_char
    c_k_id = cs_f_field_key_id(c_name)

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_get_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_get_key_struct_gas_mix_species_prop

  !=============================================================================

  !> \brief  Compute cell gradient

  !> \param[in]       f_id             field id, or -1
  !> \param[in]       imrgra           gradient computation mode
  !> \param[in]       inc              0: increment; 1: do not increment
  !> \param[in]       recompute_cocg   1 or 0: recompute COCG or not
  !> \param[in]       nswrgp           number of sweeps for reconstruction
  !> \param[in]       imligp           gradient limitation method:
  !>                                     < 0 no limitation
  !>                                     = 0 based on neighboring gradients
  !>                                     = 1 based on mean gradient
  !> \param[in]       iwarnp           verbosity
  !> \param[in]       epsgrp           relative precision for reconstruction
  !> \param[in]       climgp           limiter coefficient for imligp
  !> \param[in]       extrap           gradient extrapolation coefficient
  !> \param[in, out]  pvar             cell values whose gradient is computed
  !> \param[in]       coefap           boundary coefap coefficients
  !> \param[in]       coefbp           boundary coefap coefficients
  !> \param[out]      grad             resulting gradient

  subroutine gradient_s(f_id, imrgra, inc, recompute_cocg, nswrgp,             &
                        imligp, iwarnp, epsrgp, climgp, extrap,                &
                        pvar, coefap, coefbp, grad)

    use, intrinsic :: iso_c_binding
    use paramx
    use mesh
    use field
    use period

    implicit none

    ! Arguments

    integer, intent(in) :: f_id, imrgra, inc, recompute_cocg , nswrgp
    integer, intent(in) :: imligp, iwarnp
    double precision, intent(in) :: epsrgp, climgp, extrap
    real(kind=c_double), dimension(nfabor), intent(in) :: coefap, coefbp
    real(kind=c_double), dimension(ncelet), intent(inout) :: pvar
    real(kind=c_double), dimension(3, ncelet), intent(out) :: grad

    ! Local variables

    integer        :: hyd_p_flag
    integer        :: idimtr, ipond
    type(c_ptr)    :: f

    ! Preparation for periodicity of rotation

    ! By default, the gradient will be treated as a vector ...
    !   (i.e. we assume it is the gradient of a scalar field)

    ! If rotational periodicities are present,
    !   we determine if the variable is a tensor (Reynolds stresses)
    !   so as to apply the necessary treatment.
    !   We set idimtr and we retrieve the matching gradient.
    ! Note that if halo gradients have not been saved before, they cannot be
    !   retrieved here (...)
    !   So this subroutine is called by phyvar (in perinr)
    !   to compute gradients at the beginning of the time step and save them
    !   in dudxyz et drdxyz

    ! It is necessary for idimtr to always be initialized, even with no
    !   periodicity of rotation, so it's default value is set.

    idimtr = 0

    if (iperot.eq.1 .and. f_id.gt.-1) then
      f = cs_field_by_id(f_id)
      call cs_gradient_perio_init_rij(f, idimtr, grad)
    endif

    ! The gradient of a potential (pressure, ...) is a vector

    hyd_p_flag = 0
    ipond = 0

    call cgdcel(f_id, imrgra, inc, recompute_cocg, nswrgp,                     &
                idimtr, hyd_p_flag, ipond, iwarnp, imligp, epsrgp, extrap,     &
                climgp, c_null_ptr, coefap, coefbp,                            &
                pvar, c_null_ptr, grad)

  end subroutine gradient_s

  !=============================================================================

  !> \brief  Compute cell gradient of potential-type values

  !> \param[in]       f_id             field id, or -1
  !> \param[in]       imrgra           gradient computation mode
  !> \param[in]       inc              0: increment; 1: do not increment
  !> \param[in]       recompute_cocg   1 or 0: recompute COCG or not
  !> \param[in]       nswrgp           number of sweeps for reconstruction
  !> \param[in]       imligp           gradient limitation method:
  !>                                     < 0 no limitation
  !>                                     = 0 based on neighboring gradients
  !>                                     = 1 based on mean gradient
  !> \param[in]       hyd_p_flag       flag for hydrostatic pressure
  !> \param[in]       iwarnp           verbosity
  !> \param[in]       epsgrp           relative precision for reconstruction
  !> \param[in]       climgp           limiter coefficient for imligp
  !> \param[in]       extrap           gradient extrapolation coefficient
  !> \param[in]       f_ext            exterior force generating
  !>                                   the hydrostatic pressure
  !> \param[in, out]  pvar             cell values whose gradient is computed
  !> \param[in]       coefap           boundary coefap coefficients
  !> \param[in]       coefbp           boundary coefap coefficients
  !> \param[out]      grad             resulting gradient

  subroutine gradient_potential_s(f_id, imrgra, inc, recompute_cocg, nswrgp,   &
                                  imligp, hyd_p_flag, iwarnp, epsrgp, climgp,  &
                                  extrap, f_ext, pvar, coefap, coefbp, grad)

    use, intrinsic :: iso_c_binding
    use paramx
    use mesh
    use field

    implicit none

    ! Arguments

    integer, intent(in) :: f_id, imrgra, inc, recompute_cocg , nswrgp
    integer, intent(in) :: imligp, hyd_p_flag, iwarnp
    double precision, intent(in) :: epsrgp, climgp, extrap
    real(kind=c_double), dimension(nfabor), intent(in) :: coefap, coefbp
    real(kind=c_double), dimension(ncelet), intent(inout) :: pvar
    real(kind=c_double), dimension(3, *), intent(in) :: f_ext
    real(kind=c_double), dimension(3, ncelet), intent(out) :: grad

    ! Local variables

    integer          :: imrgrp
    integer          :: idimtr, ipond

    ! Use iterative gradient

    if (imrgra.lt.0) then
      imrgrp = 0
    else
      imrgrp = imrgra
    endif

    ! The gradient of a potential (pressure, ...) is a vector

    idimtr = 0
    ipond = 0

    call cgdcel(f_id, imrgrp, inc, recompute_cocg, nswrgp,                     &
                idimtr, hyd_p_flag, ipond, iwarnp, imligp, epsrgp, extrap,     &
                climgp, f_ext, coefap, coefbp,                                 &
                pvar, c_null_ptr, grad)

  end subroutine gradient_potential_s

  !=============================================================================

  !> \brief  Compute cell gradient of a scalar with weighting

  !> \param[in]       f_id             field id, or -1
  !> \param[in]       imrgra           gradient computation mode
  !> \param[in]       inc              0: increment; 1: do not increment
  !> \param[in]       recompute_cocg   1 or 0: recompute COCG or not
  !> \param[in]       nswrgp           number of sweeps for reconstruction
  !> \param[in]       imligp           gradient limitation method:
  !>                                     < 0 no limitation
  !>                                     = 0 based on neighboring gradients
  !>                                     = 1 based on mean gradient
  !> \param[in]       iwarnp           verbosity
  !> \param[in]       epsgrp           relative precision for reconstruction
  !> \param[in]       climgp           limiter coefficient for imligp
  !> \param[in]       extrap           gradient extrapolation coefficient
  !> \param[in, out]  pvar             cell values whose gradient is computed
  !> \param[in, out]  c_weight         cell weighting coefficient
  !> \param[in]       coefap           boundary coefap coefficients
  !> \param[in]       coefbp           boundary coefap coefficients
  !> \param[out]      grad             resulting gradient

  subroutine gradient_weighted_s(f_id, imrgra, inc, recompute_cocg, nswrgp,   &
                                 imligp, iwarnp, epsrgp, climgp, extrap,      &
                                 pvar, c_weight, coefap, coefbp, grad)

    use, intrinsic :: iso_c_binding
    use paramx
    use mesh
    use field

    implicit none

    ! Arguments

    integer, intent(in) :: f_id, imrgra, inc, recompute_cocg , nswrgp
    integer, intent(in) :: imligp, iwarnp
    double precision, intent(in) :: epsrgp, climgp, extrap
    real(kind=c_double), dimension(nfabor), intent(in) :: coefap, coefbp
    real(kind=c_double), dimension(ncelet), intent(inout) :: pvar
    real(kind=c_double), dimension(*), intent(in) :: c_weight
    real(kind=c_double), dimension(3, ncelet), intent(out) :: grad

    ! Local variables

    integer          :: hyd_p_flag
    integer          :: idimtr, ipond

    ! The current variable is a scalar
    idimtr = 0

    ! the gradient is computed with no extern hydrostatic force
    hyd_p_flag = 0

    ! the pressure gradient coefficient weighting is used
    ipond = 1

    call cgdcel(f_id, imrgra, inc, recompute_cocg, nswrgp,                     &
                idimtr, hyd_p_flag, ipond, iwarnp, imligp, epsrgp, extrap,     &
                climgp, c_null_ptr, coefap, coefbp,                            &
                pvar, c_weight, grad)

  end subroutine gradient_weighted_s

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

  !> \brief Initialize a restart file

  !> \param[in]   name  file name
  !> \param[in]   path  optional directory name for output
  !>                    (automatically created if necessary)
  !> \param[in]   mode  read (0) or write (1)
  !> \param[out]  r     pointer to restart structure

  subroutine restart_create(name, path, mode, r)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name, path
    integer, intent(in)          :: mode
    type(c_ptr), intent(out)     :: r

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    character(len=len_trim(path)+1, kind=c_char) :: c_path
    integer(c_int) :: c_mode

    c_name = trim(name)//c_null_char
    c_path = trim(path)//c_null_char
    c_mode = mode

    r = cs_restart_create(c_name, c_path, c_mode)

  end subroutine restart_create

  !---------------------------------------------------------------------------

  !> \brief Read variables from checkpoint.

  !> \param[in]   r              pointer to restart structure
  !> \param[in]   old_field_map  old field map pointer
  !> \param[in]   t_id_flag      -1: all time values; 0: current values;
  !>                             > 0: previous values

  subroutine restart_read_variables(r, old_field_map, t_id_flag)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in) :: r
    integer, intent(in)     :: t_id_flag
    type(c_ptr), intent(in) :: old_field_map

    ! Local variables

    integer(c_int) :: c_t_id_flag

    c_t_id_flag = t_id_flag

    call cs_restart_read_variables(r, old_field_map, c_t_id_flag, c_null_ptr)

  end subroutine restart_read_variables

  !-----------------------------------------------------------------------------

  !> \brief Write variables to checkpoint

  !> \param[in]   r          pointer to restart structure
  !> \param[in]   t_id_flag  -1: all time values; 0: current values;
  !>                         > 0: previous values

  subroutine restart_write_variables(r, t_id_flag)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in) :: r
    integer, intent(in)     :: t_id_flag

    ! Local variables

    integer(c_int) :: c_t_id_flag

    c_t_id_flag = t_id_flag

    call cs_restart_write_variables(r, c_t_id_flag, c_null_ptr)

  end subroutine restart_write_variables

  !---------------------------------------------------------------------------

  !> \brief Read a section of integers from a restart file.

  !> \param[in]   r             pointer to restart structure
  !> \param[in]   sec_name      name of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[in]   n_loc_vals    number of values per location
  !> \param[out]  val           values array
  !> \param[out]  ierror        0: success, < 0: error code

  subroutine restart_read_section_int_t(r, sec_name,                       &
                                        location_id, n_loc_vals, val,      &
                                        ierror)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in)           :: r
    character(len=*), intent(in)      :: sec_name
    integer, intent(in)               :: location_id, n_loc_vals
    integer, dimension(*), target     :: val
    integer, intent(out)              :: ierror

    ! Local variables

    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_n
    integer(c_int) :: c_loc_id, c_n_l_vals, c_val_type, c_ierror
    type(c_ptr) :: c_val

    c_s_n = trim(sec_name)//c_null_char
    c_loc_id = location_id
    c_n_l_vals = n_loc_vals
    c_val_type = RESTART_VAL_TYPE_INT_T
    c_val = c_loc(val)

    c_ierror = cs_restart_read_section(r, c_s_n, c_loc_id,         &
                                       c_n_l_vals, c_val_type,     &
                                       c_val)

    ierror = c_ierror

  end subroutine restart_read_section_int_t

  !---------------------------------------------------------------------------

  !> \brief Read a section of integers from a restart file,
  !> when that section may have used a different name in a previous version.

  !> \param[in]   r             pointer to restart structure
  !> \param[in]   sec_name      name of section
  !> \param[in]   old_name      old name of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[in]   n_loc_vals    number of values per location
  !> \param[out]  val           values array
  !> \param[out]  ierror        0: success, < 0: error code

  subroutine restart_read_int_t_compat(r, sec_name, old_name,                &
                                       location_id, n_loc_vals, val,         &
                                       ierror)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in)           :: r
    character(len=*), intent(in)      :: sec_name
    character(len=*), intent(in)      :: old_name
    integer, intent(in)               :: location_id, n_loc_vals
    integer, dimension(*), target     :: val
    integer, intent(out)              :: ierror

    ! Local variables

    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_n
    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_o
    integer(c_int) :: c_loc_id, c_n_l_vals, c_val_type, c_ierror
    type(c_ptr) :: c_val

    c_s_n = trim(sec_name)//c_null_char
    c_s_o = trim(old_name)//c_null_char
    c_loc_id = location_id
    c_n_l_vals = n_loc_vals
    c_val_type = RESTART_VAL_TYPE_INT_T
    c_val = c_loc(val)

    c_ierror = cs_restart_read_section_compat(r, c_s_n, c_s_o,        &
                                              c_loc_id, c_n_l_vals,   &
                                              c_val_type, c_val)

    ierror = c_ierror

  end subroutine restart_read_int_t_compat

  !---------------------------------------------------------------------------

  !> \brief Write a section of integers to a checkpoint file.

  !> \param[in]   r             pointer to restart structure
  !> \param[in]   sec_name      name of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[in]   n_loc_vals    number of values per location
  !> \param[in]   val           values array

  subroutine restart_write_section_int_t(r, sec_name,                      &
                                         location_id, n_loc_vals, val)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in)           :: r
    character(len=*), intent(in)      :: sec_name
    integer, intent(in)               :: location_id, n_loc_vals
    integer, dimension(*), intent(in), target :: val

    ! Local variables

    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_n
    integer(c_int) :: c_loc_id, c_n_l_vals, c_val_type
    type(c_ptr) :: c_val

    c_s_n = trim(sec_name)//c_null_char
    c_loc_id = location_id
    c_n_l_vals = n_loc_vals
    c_val_type = RESTART_VAL_TYPE_INT_T
    c_val = c_loc(val)

    call cs_restart_write_section(r, c_s_n, c_loc_id,         &
                                  c_n_l_vals, c_val_type,     &
                                  c_val)

  end subroutine restart_write_section_int_t

  !---------------------------------------------------------------------------

  !> \brief Read a section of doubles from a restart file.

  !> \param[in]   r             pointer to restart structure
  !> \param[in]   sec_name      name of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[in]   n_loc_vals    number of values per location
  !> \param[out]  val           values array
  !> \param[out]  ierror        0: success, < 0: error code

  subroutine restart_read_section_real_t(r, sec_name,                      &
                                         location_id, n_loc_vals, val,     &
                                         ierror)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in)                   :: r
    character(len=*), intent(in)              :: sec_name
    integer, intent(in)                       :: location_id, n_loc_vals
    real(kind=c_double), dimension(*), target :: val
    integer, intent(out)                      :: ierror

    ! Local variables

    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_n
    integer(c_int) :: c_loc_id, c_n_l_vals, c_val_type, c_ierror
    type(c_ptr) :: c_val

    c_s_n = trim(sec_name)//c_null_char
    c_loc_id = location_id
    c_n_l_vals = n_loc_vals
    c_val_type = RESTART_VAL_TYPE_REAL_T
    c_val = c_loc(val)

    c_ierror = cs_restart_read_section(r, c_s_n, c_loc_id,   &
                                       c_n_l_vals, c_val_type,     &
                                       c_val)

    ierror = c_ierror

  end subroutine restart_read_section_real_t

  !---------------------------------------------------------------------------

  !> \brief Read a section of double precision reals from a restart file,
  !> when that section may have used a different name in a previous version.

  !> \param[in]   r             pointer to restart structure
  !> \param[in]   sec_name      name of section
  !> \param[in]   old_name      old name of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[in]   n_loc_vals    number of values per location
  !> \param[out]  val           values array
  !> \param[out]  ierror        0: success, < 0: error code

  subroutine restart_read_real_t_compat(r, sec_name, old_name,               &
                                        location_id, n_loc_vals, val,        &
                                        ierror)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in)           :: r
    character(len=*), intent(in)      :: sec_name
    character(len=*), intent(in)      :: old_name
    integer, intent(in)               :: location_id, n_loc_vals
    real(kind=c_double), dimension(*), target :: val
    integer, intent(out)              :: ierror

    ! Local variables

    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_n
    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_o
    integer(c_int) :: c_loc_id, c_n_l_vals, c_val_type, c_ierror
    type(c_ptr) :: c_val

    c_s_n = trim(sec_name)//c_null_char
    c_s_o = trim(old_name)//c_null_char
    c_loc_id = location_id
    c_n_l_vals = n_loc_vals
    c_val_type = RESTART_VAL_TYPE_REAL_T
    c_val = c_loc(val)

    c_ierror = cs_restart_read_section_compat(r, c_s_n, c_s_o,        &
                                              c_loc_id, c_n_l_vals,   &
                                              c_val_type, c_val)

    ierror = c_ierror

  end subroutine restart_read_real_t_compat

  !---------------------------------------------------------------------------

  !> \brief Read a vector of double precision reals of dimension (3,*) from a
  !> restart file, when that section may have used a different name and
  !> been non-interleaved in a previous version.

  !> \param[in]   r             pointer to restart structure
  !> \param[in]   sec_name      name of section
  !> \param[in]   old_name_x    old name of component x of section
  !> \param[in]   old_name_y    old name of component y of section
  !> \param[in]   old_name_z    old name of component z of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[out]  val           values array
  !> \param[out]  ierror        0: success, < 0: error code

  subroutine restart_read_real_3_t_compat(r, sec_name,                         &
                                          old_name_x, old_name_y, old_name_z,  &
                                          location_id, val, ierror)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in)           :: r
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
    integer(c_int) :: c_loc_id, c_ierror

    c_s_n = trim(sec_name)//c_null_char
    c_o_n_x = trim(old_name_x)//c_null_char
    c_o_n_y = trim(old_name_y)//c_null_char
    c_o_n_z = trim(old_name_z)//c_null_char
    c_loc_id = location_id

    c_ierror = cs_restart_read_real_3_t_compat(r, c_s_n, c_o_n_x, c_o_n_y,     &
                                               c_o_n_z, c_loc_id, val)

    ierror = c_ierror

  end subroutine restart_read_real_3_t_compat

  !---------------------------------------------------------------------------

  !> \brief write a section of doubles to a checkpoint file.

  !> \param[in]   r             pointer to restart structure
  !> \param[in]   sec_name      name of section
  !> \param[in]   location_id   id of associated mesh location
  !> \param[in]   n_loc_vals    number of values per location
  !> \param[in]   val           values array
  !> \param[out]  ierror        0: success, < 0: error code

  subroutine restart_write_section_real_t(r, sec_name,                     &
                                          location_id, n_loc_vals, val)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in)                       :: r
    character(len=*), intent(in)                  :: sec_name
    integer, intent(in)                           :: location_id, n_loc_vals
    real(kind=c_double), dimension(*), target, intent(in) :: val

    ! Local variables

    character(len=len_trim(sec_name)+1, kind=c_char) :: c_s_n
    integer(c_int) :: c_loc_id, c_n_l_vals, c_val_type
    type(c_ptr) :: c_val

    c_s_n = trim(sec_name)//c_null_char
    c_loc_id = location_id
    c_n_l_vals = n_loc_vals
    c_val_type = RESTART_VAL_TYPE_REAL_T
    c_val = c_loc(val)

    call cs_restart_write_section(r, c_s_n, c_loc_id,         &
                                  c_n_l_vals, c_val_type,     &
                                  c_val)

  end subroutine restart_write_section_real_t

  !---------------------------------------------------------------------------

  !> \brief Read field values from checkpoint.

  !> If the values are not found using the default rules based on the
  !> field's name, its name itself, or a "restart_rename" keyed string value,
  !> an old name may be used for compatibility with older files.
  !> For cell-based fields, the old name base is appended automatically with
  !> "_ce_phase01", except for scalars, where the name uses a different scheme,
  !> based on "scalaire_ce_%04" % s_num;

  !> \param[in]   r       pointer to restart structure
  !> \param[in]   f_id    field id
  !> \param[in]   t_id    time id (0 for current, 1 for previous, ...)
  !> \param[out]  ierror  return code

  subroutine restart_read_field_vals(r, f_id, t_id, ierror)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in) :: r
    integer, intent(in)     :: f_id, t_id
    integer, intent(out)    :: ierror

    ! Local variables

    integer(c_int) :: c_f_id, c_t_id, c_retcode
    c_f_id = f_id
    c_t_id = t_id

    c_retcode = cs_restart_read_field_vals(r, c_f_id, c_t_id)
    ierror = c_retcode

  end subroutine restart_read_field_vals

  !---------------------------------------------------------------------------

  !> \brief Write field values to checkpoint.

  !> \param[in]   r       pointer to restart structure
  !> \param[in]   f_id    field id
  !> \param[in]   t_id    time id (0 for current, 1 for previous, ...)
  !> \param[out]  ierror  return code

  subroutine restart_write_field_vals(r, f_id, t_id)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in) :: r
    integer, intent(in)     :: f_id, t_id

    ! Local variables

    integer(c_int) :: c_f_id, c_t_id
    c_f_id = f_id
    c_t_id = t_id

    call cs_restart_write_field_vals(r, c_f_id, c_t_id)

  end subroutine restart_write_field_vals

  !---------------------------------------------------------------------------

  !> \brief Read fields depending on others from checkpoint.

  !> \param[in]   r              pointer to restart structure
  !> \param[in]   old_field_map  pointer to old field map
  !> \param[in]   key            key for field association
  !> \param[out]  n_w            number of fields read

    ! Interface to C function writing

  subroutine restart_read_linked_fields(r, old_field_map, key, n_w)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in)      :: r
    type(c_ptr), intent(in)      :: old_field_map
    character(len=*), intent(in) :: key
    integer, intent(out)         :: n_w

    ! Local variables

    integer(c_int) :: c_n_w
    character(len=len_trim(key)+1, kind=c_char) :: c_key

    c_key = trim(key)//c_null_char

    c_n_w = cs_restart_read_linked_fields(r, old_field_map, c_key, c_null_ptr)

    n_w = c_n_w

  end subroutine restart_read_linked_fields

  !---------------------------------------------------------------------------

  !> \brief Write fields depending on others to checkpoint.

  !> \param[in]   r    pointer to restart structure
  !> \param[in]   key  key for field association
  !> \param[out]  n_w  number of fields written

    ! Interface to C function writing

  subroutine restart_write_linked_fields(r, key, n_w)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr), intent(in)      :: r
    character(len=*), intent(in) :: key
    integer, intent(out)         :: n_w

    ! Local variables

    integer(c_int) :: c_n_w
    character(len=len_trim(key)+1, kind=c_char) :: c_key

    c_key = trim(key)//c_null_char

    c_n_w = cs_restart_write_linked_fields(r, c_key, c_null_ptr)

    n_w = c_n_w

  end subroutine restart_write_linked_fields

  !=============================================================================

  !> \brief Call sparse linear equation solver using native matrix arrays.

  !> param[in]       f_id     associated field id, or < 0
  !> param[in]       name     associated name if f_id < 0, or ignored
  !> param[in]       isym     symmetry indicator: 1 symmetric, 2: not symmetric
  !> param[in]       ibsize   block sizes for diagonal
  !> param[in]       iesize   block sizes for extra diagonal
  !> param[in]       dam      matrix diagonal
  !> param[in]       xam      matrix extra-diagonal terms
  !> param[in]       iinvpe   Indicator to cancel increments in rotational
  !>                          periodicty (2) or to exchange them as scalars (1)
  !> param[in]       epsilp   precision for iterative resolution
  !> param[in]       rnorm    residue normalization
  !> param[out]      niter    number of "equivalent" iterations
  !> param[out]      residue  residue
  !> param[in]       rhs      right hand side
  !> param[in, out]  vx       system solution

  subroutine sles_solve_native(f_id, name, isym, ibsize, iesize, dam, xam,     &
                               iinvpe, epsilp, rnorm, niter, residue, rhs, vx)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)      :: name
    integer, intent(in)               :: f_id, isym, ibsize, iesize, iinvpe
    double precision, intent(in)      :: rnorm, epsilp
    integer, intent(out)              :: niter
    double precision, intent(out)     :: residue
    real(kind=c_double), dimension(*), intent(in) :: dam, xam, rhs
    real(kind=c_double), dimension(*), intent(inout) :: vx

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: rotation_mode, cvg
    integer(c_int), dimension(4) :: db_size, eb_size
    logical(kind=c_bool) :: c_sym

    c_name = trim(name)//c_null_char

    if (isym.eq.1) then
      c_sym = .true.
    else
      c_sym = .false.
    endif

    if (iinvpe.eq.2) then
      rotation_mode = 1 ! CS_HALO_ROTATION_ZERO
    else if (iinvpe.eq.3) then
      rotation_mode = 2 ! CS_HALO_ROTATION_IGNORE
    else
      rotation_mode = 0 ! CS_HALO_ROTATION_COPY, might not be called
    endif

    db_size(1) = ibsize
    db_size(2) = ibsize
    db_size(3) = ibsize
    db_size(4) = ibsize*ibsize

    eb_size(1) = iesize
    eb_size(2) = iesize
    eb_size(3) = iesize
    eb_size(4) = iesize*iesize

    cvg = cs_sles_solve_native(f_id, c_name, c_sym, db_size, eb_size,         &
                               dam, xam, rotation_mode, epsilp, rnorm,        &
                               niter, residue, rhs, vx)

    return

  end subroutine sles_solve_native

  !=============================================================================

  !> \brief Free sparse linear equation solver setup using native matrix arrays.

  !> param[in]       f_id     associated field id, or < 0
  !> param[in]       name     associated name if f_id < 0, or ignored

  subroutine sles_free_native(f_id, name)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)      :: name
    integer, intent(in)               :: f_id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_name = trim(name)//c_null_char

    call cs_sles_free_native(f_id, c_name)

    return

  end subroutine sles_free_native

  !=============================================================================

  !> \brief Temporarily replace field id with name for matching calls
  !>        to \ref sles_solve_native

  !> param[in]       f_id     associated field id, or < 0
  !> param[in]       name     associated name if f_id < 0, or ignored

  subroutine sles_push(f_id, name)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)      :: name
    integer, intent(in)               :: f_id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_name = trim(name)//c_null_char

    call cs_sles_push(f_id, c_name)

    return

  end subroutine sles_push

  !=============================================================================

  !> \brief Revert to normal behavior of field id for matching calls
  !>        to \ref sles_solve_native

  !> param[in]  f_id   associated field id, or < 0

  subroutine sles_pop(f_id)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in) :: f_id

    ! Local variables

    call cs_sles_pop(f_id)

    return

  end subroutine sles_pop

  !=============================================================================

  !> \brief Compute balance on a given zone for a given scalar

  !> param[in]       itypfb     array of boundary faces type
  !> param[in]       sel_crit   selection criterium of a volumic zone
  !> param[in]       name       scalar name

  subroutine balance_by_zone(itypfb, sel_crit, name)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer(kind=c_int), dimension(*), intent(in) :: itypfb
    character(len=*), intent(in)             :: sel_crit, name

    ! Local variables

    character(len=len_trim(sel_crit)+1, kind=c_char) :: c_sel_crit
    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_sel_crit = trim(sel_crit)//c_null_char
    c_name = trim(name)//c_null_char

    call cs_balance_by_zone(itypfb, c_sel_crit, c_name)

    return

  end subroutine balance_by_zone

  !=============================================================================

  end module cs_c_bindings
