!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

  integer :: RESTART_DISABLED, RESTART_MAIN, RESTART_AUXILIARY
  integer :: RESTART_RAD_TRANSFER, RESTART_LAGR, RESTART_LAGR_STAT
  integer :: RESTART_1D_WALL_THERMAL, RESTART_LES_INFLOW

  integer :: VOLUME_ZONE_INITIALIZATION, VOLUME_ZONE_POROSITY
  integer :: VOLUME_ZONE_HEAD_LOSS
  integer :: VOLUME_ZONE_SOURCE_TERM, VOLUME_ZONE_MASS_SOURCE_TERM

  parameter (MESH_LOCATION_NONE=0)
  parameter (MESH_LOCATION_CELLS=1)
  parameter (MESH_LOCATION_INTERIOR_FACES=2)
  parameter (MESH_LOCATION_BOUNDARY_FACES=3)
  parameter (MESH_LOCATION_VERTICES=4)
  parameter (MESH_LOCATION_PARTICLES=5)
  parameter (MESH_LOCATION_OTHER=6)

  parameter (RESTART_VAL_TYPE_INT_T=1)
  parameter (RESTART_VAL_TYPE_REAL_T=3)

  parameter (RESTART_DISABLED=-1)
  parameter (RESTART_MAIN=0)
  parameter (RESTART_AUXILIARY=1)
  parameter (RESTART_RAD_TRANSFER=2)
  parameter (RESTART_LAGR=3)
  parameter (RESTART_LAGR_STAT=4)
  parameter (RESTART_1D_WALL_THERMAL=5)
  parameter (RESTART_LES_INFLOW=6)

  parameter (VOLUME_ZONE_INITIALIZATION=1)
  parameter (VOLUME_ZONE_POROSITY=2)
  parameter (VOLUME_ZONE_HEAD_LOSS=4)
  parameter (VOLUME_ZONE_SOURCE_TERM=8)
  parameter (VOLUME_ZONE_MASS_SOURCE_TERM=16)

  !-----------------------------------------------------------------------------

  type, bind(c)  :: var_cal_opt
    integer(c_int) :: iwarni
    integer(c_int) :: iconv
    integer(c_int) :: istat
    integer(c_int) :: idircl
    integer(c_int) :: ndircl
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
    integer(c_int) :: icoupl
    real(c_double) :: thetav
    real(c_double) :: blencv
    real(c_double) :: blend_st
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

  type, bind(c)  :: gwf_soilwater_partition
    integer(c_int) :: kinetic
    integer(c_int) :: ikd
    integer(c_int) :: idel
    integer(c_int) :: ikp
    integer(c_int) :: ikm
    integer(c_int) :: imxsol
    integer(c_int) :: resol_method
  end type gwf_soilwater_partition

  !---------------------------------------------------------------------------

  type, bind(c)  :: gas_mix_species_prop
    real(c_double) :: mol_mas
    real(c_double) :: cp
    real(c_double) :: vol_dif
    real(c_double) :: mu_a
    real(c_double) :: mu_b
    real(c_double) :: lambda_a
    real(c_double) :: lambda_b
    real(c_double) :: muref
    real(c_double) :: lamref
    real(c_double) :: trefmu
    real(c_double) :: treflam
    real(c_double) :: smu
    real(c_double) :: slam
  end type gas_mix_species_prop

  !=============================================================================

  interface

    subroutine max_limiter_building(f_id, inc, rovsdt) &
    bind(C, name='cs_max_limiter_building')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      integer(c_int), value :: inc
      real(c_double), dimension(*) , intent(in) :: rovsdt
    end subroutine max_limiter_building

    !---------------------------------------------------------------------------

    !> \brief Set mapped boundary conditions for a given field and mapping
    !>        locator.

    !> param[in]       field_id         id of field whose boundary conditions
    !>                                  are set
    !> param[in]       locator          associated mapping locator, as returned
    !>                                  by \ref cs_boundary_conditions_map.
    !> param[in]       location_type    matching values location
    !>                                  (CS_MESH_LOCATION_CELLS or
    !>                                  CS_MESH_LOCATION_BOUNDARY_FACES)
    !> param[in]       normalize        normalization:
    !>                                    0: values are simply mapped
    !>                                    1: values are mapped, then multiplied
    !>                                       by a constant factor so that their
    !>                                       surface integral on selected faces
    !>                                       is preserved (relative to the
    !>                                       input values)
    !>                                    2: as 1, but with a boundary-defined
    !>                                       weight, defined by balance_w
    !>                                    3: as 1, but with a cell-defined
    !>                                       weight, defined by balance_w
    !> param[in]       interpolate      interpolation option:
    !>                                    0: values are simply based on
    !>                                       matching cell or face center values
    !>                                    1: values are based on matching cell
    !>                                       or face center values, corrected
    !>                                       by gradient interpolation
    !> param[in]       n_faces          number of selected boundary faces
    !> param[in]       faces            list of selected boundary faces (1 to n)
    !> param[in]       balance_w        optional balance weight
    !> param[in]       nvar             number of variables with BC's
    !> param[in, out]  rcodcl           boundary condition values

    subroutine boundary_conditions_mapped_set(field_id, locator,               &
                                              location_type, normalize,        &
                                              interpolate, n_faces, faces,     &
                                              balance_w, nvar, rcodcl)         &
      bind(C, name='cs_f_boundary_conditions_mapped_set')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: field_id
      type(c_ptr), value :: locator
      integer(c_int), value :: location_type, normalize, interpolate
      integer(c_int), value :: n_faces, nvar
      integer(c_int), dimension(*), intent(in) :: faces
      real(kind=c_double), dimension(*), intent(in) :: balance_w, rcodcl
    end subroutine boundary_conditions_mapped_set

    !---------------------------------------------------------------------------

    !> \brief  Return the number of fans.

    !> \return number of defined fans

    function cs_fan_n_fans() result(n_fans)  &
      bind(C, name='cs_fan_n_fans')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: n_fans
    end function cs_fan_n_fans

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

    ! Interface to C function logging L2 time residual at relevant time steps.

    ! \brief Log field and other array statistics for a given time step.

    subroutine log_l2residual()  &
      bind(C, name='cs_log_l2residual')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine log_l2residual

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

    !> \brief Loop over all fields and save them in the restart file which id is
    !>        passed in argument if it matches their "restart_file" key value.

    !> \param[in, out]   r       restart structure pointer
    !> \param[in]        r_id    key value

    subroutine restart_write_fields(r, r_id)  &
      bind(C, name='cs_restart_write_fields')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      integer(c_int), value :: r_id
    end subroutine restart_write_fields

    !---------------------------------------------------------------------------

    !> \brief Loop over all fields and read them in the restart file which id is
    !>        passed in argument if it matches their "restart_file" key value.

    !> \param[in, out]   r       restart structure pointer
    !> \param[in]        r_id    key value

    subroutine restart_read_fields(r, r_id)  &
      bind(C, name='cs_restart_read_fields')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      integer(c_int), value :: r_id
    end subroutine restart_read_fields

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

    !> \brief  Increment time step for timer statistics.

    !> \param[in]   id    id of statistic

    subroutine timer_stats_increment_time_step()  &
      bind(C, name='cs_timer_stats_increment_time_step')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine timer_stats_increment_time_step

    !---------------------------------------------------------------------------

    !> \brief  Enable or disable plotting for a timer statistic.

    !> \param[in]  id    id of statistic
    !> \param[in]  plot  0 to disable, 1 to enable

    subroutine timer_stats_set_plot(id, plot)  &
      bind(C, name='cs_timer_stats_set_plot')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: id, plot
    end subroutine timer_stats_set_plot

    !---------------------------------------------------------------------------

    !> \brief  Start a timer for a given statistic.

    !> Parents of the current statistic are also started, if not active.

    !> If a timer with the same root but different parents is active, we assume
    !> the current operation is a subset of the active timer, so the timer is
    !> not started, so as to avoid having a sum of parts larger thn the total.

    !> \param[in]   id    id of statistic

    subroutine timer_stats_start(id)  &
      bind(C, name='cs_timer_stats_start')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: id
    end subroutine timer_stats_start

    !---------------------------------------------------------------------------

    !> \brief  Stop a timer for a given statistic.

    !> \param[in]   id    id of statistic

    subroutine timer_stats_stop(id)  &
      bind(C, name='cs_timer_stats_stop')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: id
    end subroutine timer_stats_stop

    !---------------------------------------------------------------------------

    !> \brief  Start a timer for a given statistic, stopping previous timers
    !>         of the same type which are not a parent, and starting inactive
    !>         parent timers if necessary.

    !> \param[in]   id    id of statistic

    !> \return  id of previously active statistic, or -1 in case of error

    function timer_stats_switch(id)  result(old_id)  &
      bind(C, name='cs_timer_stats_switch')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: id
      integer(c_int)        :: old_id
    end function timer_stats_switch

    !---------------------------------------------------------------------------

    !> \brief Calculation of \f$ u^\star \f$, \f$ k \f$ and \f$\varepsilon \f$
    !>        from a diameter \f$ D_H \f$ and the reference velocity
    !>        \f$ U_{ref} \f$
    !>        for a circular duct flow with smooth wall
    !>        (use for inlet boundary conditions).
    !>
    !> Both \f$ u^\star \f$ and\f$ (k,\varepsilon )\f$ are returned, so that
    !> the user may compute other values of \f$ k \f$ and \f$ \varepsilon \f$
    !> with \f$ u^\star \f$.
    !>
    !> We use the laws from Idel'Cik, i.e.
    !> the head loss coefficient \f$ \lambda \f$ is defined by:
    !> \f[ |\dfrac{\Delta P}{\Delta x}| =
    !>                        \dfrac{\lambda}{D_H} \frac{1}{2} \rho U_{ref}^2 \f]
    !>
    !> then  the relation reads \f$u^\star = U_{ref} \sqrt{\dfrac{\lambda}{8}}\f$.
    !> \f$\lambda \f$ depends on the hydraulic Reynolds number
    !> \f$ Re = \dfrac{U_{ref} D_H}{ \nu} \f$ and is given by:
    !>  - for \f$ Re < 2000 \f$
    !>      \f[ \lambda = \dfrac{64}{Re} \f]
    !>
    !>  - for \f$ Re > 4000 \f$
    !>      \f[ \lambda = \dfrac{1}{( 1.8 \log_{10}(Re)-1.64 )^2} \f]
    !>
    !>  - for \f$ 2000 < Re < 4000 \f$, we complete by a straight line
    !>      \f[ \lambda = 0.021377 + 5.3115. 10^{-6} Re \f]
    !>
    !>  From \f$ u^\star \f$, we can estimate \f$ k \f$ and \f$ \varepsilon\f$
    !>  from the well known formulae of developped turbulence
    !>
    !> \f[ k = \dfrac{u^{\star 2}}{\sqrt{C_\mu}} \f]
    !> \f[ \varepsilon = \dfrac{ u^{\star 3}}{(\kappa D_H /10)} \f]
    !>
    !> \param[in]     uref2         square of the reference flow velocity
    !> \param[in]     dh            hydraulic diameter \f$ D_H \f$
    !> \param[in]     rho           mass density \f$ \rho \f$
    !> \param[in]     mu            dynamic viscosity \f$ \nu \f$
    !> \param[out]    ustar2        square of friction speed
    !> \param[out]    k             calculated turbulent intensity \f$ k \f$
    !> \param[out]    eps           calculated turbulent dissipation
    !>                               \f$ \varepsilon \f$

    subroutine turbulence_bc_ke_hyd_diam(uref2, dh, rho, mu,                   &
                                         ustar2, k, eps)                       &
      bind(C, name='cs_turbulence_bc_ke_hyd_diam')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: uref2, dh, rho, mu
      real(c_double) :: ustar2, k, eps
    end subroutine turbulence_bc_ke_hyd_diam

    !---------------------------------------------------------------------------

    !> \brief Calculation of \f$ k \f$ and \f$\varepsilon\f$
    !>        from a diameter \f$ D_H \f$, a turbulent intensity \f$ I \f$
    !>        and the reference velocity \f$ U_{ref} \f$
    !>        for a circular duct flow with smooth wall
    !>        (for inlet boundary conditions).
    !>
    !> \param[in]     uref2         square of the reference flow velocity
    !> \param[in]     t_intensity   turbulent intensity \f$ I \f$
    !> \param[in]     dh            hydraulic diameter \f$ D_H \f$
    !> \param[out]    k             calculated turbulent intensity \f$ k \f$
    !> \param[out]    eps           calculated turbulent dissipation
    !>                               \f$ \varepsilon \f$

    subroutine turbulence_bc_ke_turb_intensity(uref2, t_intensity, dh,         &
                                               k, eps)                         &
      bind(C, name='cs_turbulence_bc_ke_turb_intensity')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), value :: uref2, t_intensity, dh
      real(kind=c_double), intent(inout) :: k, eps
    end subroutine turbulence_bc_ke_turb_intensity

    !---------------------------------------------------------------------------

    !> \brief Compute matrix \f$ \tens{alpha} \f$ used in the computation of the
    !>        Reynolds stress tensor boundary conditions.
    !>
    !> \param[in]      is_sym  Constant c in description above
    !>                         (1 at a symmetry face, 0 at a wall face)
    !> \param[in]      p_lg    change of basis matrix (local to global)
    !> \param[out]     alpha   transformation matrix

    subroutine turbulence_bc_rij_transform(is_sym, p_lg, alpha)                &
      bind(C, name='cs_turbulence_bc_rij_transform')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: is_sym
      real(c_double), dimension(3,3), intent(in) :: p_lg
      real(c_double), dimension(6,6), intent(out) :: alpha
    end subroutine turbulence_bc_rij_transform

    !---------------------------------------------------------------------------

    !> \brief Set inlet boundary condition values for turbulence variables based
    !>        on a diameter \f$ D_H \f$ and the reference velocity
    !>        \f$ U_{ref} \f$
    !>        for a circular duct flow with smooth wall.
    !>
    !> We use the laws from Idel'Cik, i.e.
    !> the head loss coefficient \f$ \lambda \f$ is defined by:
    !> \f[ |\dfrac{\Delta P}{\Delta x}| =
    !>                        \dfrac{\lambda}{D_H} \frac{1}{2} \rho U_{ref}^2 \f]
    !>
    !> then  the relation reads \f$u^\star = U_{ref} \sqrt{\dfrac{\lambda}{8}}\f$.
    !> \f$\lambda \f$ depends on the hydraulic Reynolds number
    !> \f$ Re = \dfrac{U_{ref} D_H}{ \nu} \f$ and is given by:
    !>  - for \f$ Re < 2000 \f$
    !>      \f[ \lambda = \dfrac{64}{Re} \f]
    !>
    !>  - for \f$ Re > 4000 \f$
    !>      \f[ \lambda = \dfrac{1}{( 1.8 \log_{10}(Re)-1.64 )^2} \f]
    !>
    !>  - for \f$ 2000 < Re < 4000 \f$, we complete by a straight line
    !>      \f[ \lambda = 0.021377 + 5.3115. 10^{-6} Re \f]
    !>
    !>  From \f$ u^\star \f$, we can estimate \f$ k \f$ and \f$ \varepsilon\f$
    !>  from the well known formulae of developped turbulence
    !>
    !> \f[ k = \dfrac{u^{\star 2}}{\sqrt{C_\mu}} \f]
    !> \f[ \varepsilon = \dfrac{ u^{\star 3}}{(\kappa D_H /10)} \f]
    !>
    !> \param[in]     face_num   boundary face number
    !> \param[in]     uref2      square of the reference flow velocity
    !> \param[in]     dh         hydraulic diameter \f$ D_H \f$
    !> \param[in]     rho        mass density \f$ \rho \f$
    !> \param[in]     mu         dynamic viscosity \f$ \nu \f$
    !> \param[out]    rcodcl     boundary condition values

    subroutine turbulence_bc_inlet_hyd_diam(face_num, uref2, dh, rho, mu,      &
                                            rcodcl)                            &
      bind(C, name='cs_f_turbulence_bc_inlet_hyd_diam')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: face_num
      real(c_double), value :: uref2, dh, rho, mu
      real(kind=c_double), dimension(*) :: rcodcl
    end subroutine turbulence_bc_inlet_hyd_diam

    !---------------------------------------------------------------------------

    !> \brief Set inlet boundary condition values for turbulence variables based
    !>        on a diameter \f$ D_H \f$, a turbulent intensity \f$ I \f$
    !>        and the reference velocity \f$ U_{ref} \f$
    !>        for a circular duct flow with smooth wall.
    !>
    !> \param[in]     face_id       boundary face id
    !> \param[in]     uref2         square of the reference flow velocity
    !> \param[in]     t_intensity   turbulent intensity \f$ I \f$
    !> \param[in]     dh            hydraulic diameter \f$ D_H \f$
    !> \param[out]    rcodcl        boundary condition values

    subroutine turbulence_bc_inlet_turb_intensity(face_num,                    &
                                                  uref2, t_intensity, dh,      &
                                                  rcodcl)                      &
      bind(C, name='cs_f_turbulence_bc_inlet_turb_intensity')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: face_num
      real(c_double), value :: uref2, t_intensity, dh
      real(kind=c_double), dimension(*) :: rcodcl
    end subroutine turbulence_bc_inlet_turb_intensity

    !---------------------------------------------------------------------------

    !> \brief Set inlet boundary condition values for turbulence variables based
    !>        on given k and epsilon values.
    !>
    !> \param[in]     face_id       boundary face id
    !> \param[in]     k             turbulent kinetic energy
    !> \param[in]     epsilon       turbulent dissipation
    !> \param[out]    rcodcl        boundary condition values

    subroutine turbulence_bc_inlet_k_eps(face_num,                             &
                                         k, eps,                               &
                                         rcodcl)                               &
      bind(C, name='cs_f_turbulence_bc_inlet_k_eps')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: face_num
      real(c_double), value :: k, eps
      real(kind=c_double), dimension(*) :: rcodcl
    end subroutine turbulence_bc_inlet_k_eps

    !---------------------------------------------------------------------------

    !> \brief Set inlet boundary condition values for turbulence variables based
    !>        on given k and epsilon values only if not initialized already.
    !>
    !> \param[in]     face_id       boundary face id
    !> \param[in]     k             turbulent kinetic energy
    !> \param[in]     epsilon       turbulent dissipation
    !> \param[out]    rcodcl        boundary condition values

    subroutine turbulence_bc_set_uninit_inlet_k_eps(face_num,                  &
                                                    k, eps,                    &
                                                    rcodcl)                    &
      bind(C, name='cs_f_turbulence_bc_set_uninit_inlet_k_eps')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: face_num
      real(c_double), value :: k, eps
      real(kind=c_double), dimension(*) :: rcodcl
    end subroutine turbulence_bc_set_uninit_inlet_k_eps

    !---------------------------------------------------------------------------

    !> \brief  General user parameters

    subroutine user_parameters()  &
      bind(C, name='cs_user_parameters_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine user_parameters

    !---------------------------------------------------------------------------

    !> \brief  General user parameters

    subroutine user_porosity()  &
      bind(C, name='cs_user_porosity_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine user_porosity

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function handling boundary condition errors and output

    subroutine cs_boundary_conditions_error(bc_flag, type_name) &
      bind(C, name='cs_boundary_conditions_error')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), dimension(*), intent(in) :: bc_flag
      type(c_ptr), value :: type_name
    end subroutine cs_boundary_conditions_error

    !---------------------------------------------------------------------------

    ! Interface to C function locating shifted bundary face coordinates on
    ! possibly filtered cells or boundary faces for later interpolation.

    function cs_boundary_conditions_map(location_type, n_location_elts,         &
                                        n_faces, location_elts, faces,          &
                                        coord_shift, coord_stride,              &
                                        tolerance) result(l)                    &
      bind(C, name='cs_boundary_conditions_map')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: location_type, n_location_elts, n_faces
      integer(c_int), dimension(*), intent(in) :: location_elts, faces
      real(kind=c_double), dimension(*) :: coord_shift
      integer(c_int), value :: coord_stride
      real(kind=c_double), value :: tolerance
      type(c_ptr) :: l
    end function cs_boundary_conditions_map

    !---------------------------------------------------------------------------

    ! Interface to C function creating the bc type and face zone arrays

    subroutine cs_f_boundary_conditions_create() &
      bind(C, name='cs_boundary_conditions_create')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_boundary_conditions_create

    !---------------------------------------------------------------------------

    ! Interface to C function freeing the bc type and face zone arrays

    subroutine cs_f_boundary_conditions_free() &
      bind(C, name='cs_boundary_conditions_free')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_boundary_conditions_free

    !---------------------------------------------------------------------------

    ! Interface to C function to get the bc type array pointer

    subroutine cs_f_boundary_conditions_get_pointers(itypfb, izfppp) &
      bind(C, name='cs_f_boundary_conditions_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: itypfb, izfppp
    end subroutine cs_f_boundary_conditions_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function checking the presence of a control file
    ! and dealing with the interactive control.

    subroutine cs_control_check_file()  &
      bind(C, name='cs_control_check_file')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_control_check_file

    !---------------------------------------------------------------------------

    ! Interface to C function mapping field pointers

    subroutine cs_field_pointer_map_base()  &
      bind(C, name='cs_field_pointer_map_base')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_field_pointer_map_base

    !---------------------------------------------------------------------------

    ! Interface to C function mapping boundary field pointers

    subroutine cs_field_pointer_map_boundary()  &
      bind(C, name='cs_field_pointer_map_boundary')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_field_pointer_map_boundary

    !---------------------------------------------------------------------------

    ! Interface to C function returning the product of a matrix (native format)
    ! by a vector

    subroutine cs_matrix_vector_native_multiply(symmetric, db_size, eb_size,   &
                                                rotation_mode, f_id, dam, xam, &
                                                vx, vy)  &
      bind(C, name='cs_matrix_vector_native_multiply')
      use, intrinsic :: iso_c_binding
      implicit none
      logical(c_bool), value :: symmetric
      integer(c_int), dimension(4), intent(in) :: db_size, eb_size
      integer(c_int), value :: rotation_mode, f_id
      real(kind=c_double), dimension(*), intent(in) :: dam, xam, vx
      real(kind=c_double), dimension(*), intent(out) :: vy
    end subroutine cs_matrix_vector_native_multiply

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

    ! Interface to C function

    subroutine cs_bad_cells_regularisation_scalar(var) &
      bind(C, name='cs_bad_cells_regularisation_scalar')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: var
    end subroutine cs_bad_cells_regularisation_scalar

    !---------------------------------------------------------------------------

    ! Interface to C function

    subroutine cs_bad_cells_regularisation_vector(var, &
      boundary_projection) &
      bind(C, name='cs_bad_cells_regularisation_vector')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: var
      integer(c_int), value :: boundary_projection
    end subroutine cs_bad_cells_regularisation_vector

    !---------------------------------------------------------------------------

    ! Interface to C function

    subroutine cs_bad_cells_regularisation_tensor(var, &
      boundary_projection) &
      bind(C, name='cs_bad_cells_regularisation_tensor')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: var
      integer(c_int), value :: boundary_projection
    end subroutine cs_bad_cells_regularisation_tensor

    !---------------------------------------------------------------------------

    ! Interface to C function

    subroutine cs_bad_cells_regularisation_sym_tensor(var, &
      boundary_projection) &
      bind(C, name='cs_bad_cells_regularisation_sym_tensor')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: var
      integer(c_int), value :: boundary_projection
    end subroutine cs_bad_cells_regularisation_sym_tensor

    !---------------------------------------------------------------------------

    !> Interface to C function defining turbulence model through the GUI.

    subroutine cs_gui_turb_model()  &
      bind(C, name='cs_gui_turb_model')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gui_turb_model

    !---------------------------------------------------------------------------

    !> Interface to C function defining reference length and reference velocity
    !> for the initialization of the turbulence variables through the GUI.

    subroutine cs_gui_turb_ref_values()  &
      bind(C, name='cs_gui_turb_ref_values')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gui_turb_ref_values

    !---------------------------------------------------------------------------

    !> Interface to C function defining user variables through the GUI.

    subroutine cs_gui_user_variables()  &
      bind(C, name='cs_gui_user_variables')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gui_user_variables

    !---------------------------------------------------------------------------

    ! Interface to C function selecting specific physical models.

    subroutine cs_gui_physical_model_select()  &
      bind(C, name='cs_gui_physical_model_select')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gui_physical_model_select

    !---------------------------------------------------------------------------

    !> Interface to C function defining time moments through the GUI.

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
                                               min_pre_clip, max_pre_clip,    &
                                               n_clip_min_comp,               &
                                               n_clip_max_comp)               &
      bind(C, name='cs_log_iteration_clipping_field')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id, n_clip_max, n_clip_min
      real(kind=c_double), dimension(*) :: min_pre_clip, max_pre_clip
      integer(c_int), dimension(*), intent(in) :: n_clip_min_comp, n_clip_max_comp
    end subroutine cs_log_iteration_clipping_field

    !---------------------------------------------------------------------------

    ! Interface to C function initializing condensation-related field key.

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

    ! Interface to C function for uniform distribution random number

    subroutine cs_random_uniform(n, a) &
      bind(C, name='cs_random_uniform')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n
      real(kind=c_double), dimension(*), intent(out) :: a
    end subroutine cs_random_uniform

    !---------------------------------------------------------------------------

    ! Interface to C function for normal distribution random number

    subroutine cs_random_normal(n, x) &
      bind(C, name='cs_random_normal')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n
      real(kind=c_double), dimension(*), intent(out) :: x
    end subroutine cs_random_normal

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

    ! Interface to C function which checks if the restart is from NEPTUNE_CFD
    function cs_restart_check_if_restart_from_ncfd(r) result(flag) &
      bind(C, name='cs_restart_check_if_restart_from_ncfd')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
      integer(c_int) :: flag
    end function cs_restart_check_if_restart_from_ncfd

    !---------------------------------------------------------------------------

    ! Interface to C function which returns if the restart is from NEPTUNE_CFD
    function cs_restart_is_from_ncfd() result(flag) &
      bind(C, name='cs_restart_is_from_ncfd')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: flag
    end function cs_restart_is_from_ncfd

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

    ! Interface to C function defining statistic based on its name.

    function cs_timer_stats_create(parent_name, name, label) result(id) &
      bind(C, name='cs_timer_stats_create')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: parent_name
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name, label
      integer(c_int)        :: id
    end function cs_timer_stats_create

    !---------------------------------------------------------------------------

    ! Interface to C function obtaining a defined statistic based on its name.

    function cs_timer_stats_id_by_name(name) result(id) &
      bind(C, name='cs_timer_stats_id_by_name')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name
      integer(c_int)                                           :: id
    end function cs_timer_stats_id_by_name

    !---------------------------------------------------------------------------

    ! Interface to C function creating a variable field

    function cs_variable_field_create(name, label,                   &
                                      location_id, dim) result(id)   &
      bind(C, name='cs_variable_field_create')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name, label
      integer(c_int), value                                    :: location_id
      integer(c_int), value                                    :: dim
      integer(c_int)                                           :: id
    end function cs_variable_field_create

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

    subroutine cs_balance_by_zone(selection_crit, scalar_name)  &
      bind(C, name='cs_balance_by_zone')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: selection_crit
      character(kind=c_char, len=1), dimension(*), intent(in) :: scalar_name
    end subroutine cs_balance_by_zone

    !---------------------------------------------------------------------------
    ! Interface to C function for balance computation

    subroutine cs_pressure_drop_by_zone(selection_crit)  &
      bind(C, name='cs_pressure_drop_by_zone')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: selection_crit
    end subroutine cs_pressure_drop_by_zone

    !---------------------------------------------------------------------------
    ! Interface to C function for balance computation

    subroutine cs_surface_balance(selection_crit, scalar_name, normal)  &
      bind(C, name='cs_surface_balance')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: selection_crit
      character(kind=c_char, len=1), dimension(*), intent(in) :: scalar_name
      real(kind=c_double), dimension(3), intent(in) :: normal
    end subroutine cs_surface_balance

    !---------------------------------------------------------------------------

    ! Interface to C function building volume zones.

    subroutine cs_volume_zone_build_all(mesh_modified)  &
      bind(C, name='cs_volume_zone_build_all')
      use, intrinsic :: iso_c_binding
      implicit none
      logical(kind=c_bool), value :: mesh_modified
    end subroutine cs_volume_zone_build_all

    !---------------------------------------------------------------------------

    ! Interface to C function returning the number of volume zones
    ! associated with a given zone flag

    function cs_volume_zone_n_type_zones(type_flag) result(n)   &
      bind(C, name='cs_volume_zone_n_type_zones')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value                            :: type_flag
      integer(c_int)                                   :: n
    end function cs_volume_zone_n_type_zones

    !---------------------------------------------------------------------------

    ! Interface to C function returning the number of volume zone cells
    ! associated with a given zone flag

    function cs_volume_zone_n_type_cells(type_flag) result(n)   &
      bind(C, name='cs_volume_zone_n_type_cells')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value                            :: type_flag
      integer(c_int)                                   :: n
    end function cs_volume_zone_n_type_cells

    !---------------------------------------------------------------------------

    ! Interface to C function selecting cells in volume zones.

    subroutine cs_volume_zone_select_type_cells(type_flag, cell_ids)  &
      bind(C, name='cs_volume_zone_select_type_cells')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: type_flag
      type(c_ptr), value :: cell_ids
    end subroutine cs_volume_zone_select_type_cells

    !---------------------------------------------------------------------------

    ! Interface to C function building boundary zones.

    subroutine cs_boundary_zone_build_all(mesh_modified)  &
      bind(C, name='cs_boundary_zone_build_all')
      use, intrinsic :: iso_c_binding
      implicit none
      logical(kind=c_bool), value :: mesh_modified
    end subroutine cs_boundary_zone_build_all

    !---------------------------------------------------------------------------
    ! Interface to C user function for boundary conditions

    subroutine user_boundary_conditions(nvar, bc_type, icodcl, rcodcl)  &
      bind(C, name='cs_user_boundary_conditions')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: nvar
      integer(kind=c_int), dimension(*), intent(inout) :: bc_type
      integer(kind=c_int), dimension(*), intent(inout) :: icodcl
      real(kind=c_double), dimension(*), intent(inout) :: rcodcl
    end subroutine user_boundary_conditions

    !---------------------------------------------------------------------------

    ! Interface to C user function for extra operations

    subroutine user_extra_operations_initialize()  &
      bind(C, name='cs_user_extra_operations_initialize_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine user_extra_operations_initialize

    subroutine user_extra_operations()  &
      bind(C, name='cs_user_extra_operations_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine user_extra_operations

    !---------------------------------------------------------------------------

    ! Interface to C user function for initialization

    subroutine user_initialization()  &
      bind(C, name='cs_user_initialization_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine user_initialization

    !---------------------------------------------------------------------------

    ! Interface to C user function for user arrays

    subroutine cs_gui_user_arrays()  &
      bind(C, name='cs_gui_user_arrays')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_gui_user_arrays

    !---------------------------------------------------------------------------

    ! Interface to C user function for physical model options

    subroutine cs_user_model()  &
      bind(C, name='cs_user_model')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_user_model

    !---------------------------------------------------------------------------

    ! Interface to C user function for time moments

    subroutine cs_user_time_moments()  &
      bind(C, name='cs_user_time_moments')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_user_time_moments

    !---------------------------------------------------------------------------

    ! Interface to C function for the destruction of a locator structure.

    function ple_locator_destroy(this_locator) result (l) &
      bind(C, name='ple_locator_destroy')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: this_locator
      type(c_ptr) :: l
    end function ple_locator_destroy

    !---------------------------------------------------------------------------

    ! Interface to C function cs_equation_iterative_solve_scalar

    subroutine cs_equation_iterative_solve_scalar(idtvar, iterns,             &
                                                  f_id, name,                 &
                                                  iescap, imucpp, normp,      &
                                                  vcopt, pvara, pvark,        &
                                                  coefap, coefbp, cofafp,     &
                                                  cofbfp, i_massflux,         &
                                                  b_massflux, i_viscm,        &
                                                  b_viscm, i_visc, b_visc,    &
                                                  viscel, weighf, weighb,     &
                                                  icvflb, icvfli,             &
                                                  rovsdt, smbrp, pvar, dpvar, &
                                                  xcpp, eswork)               &
      bind(C, name='cs_equation_iterative_solve_scalar')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: idtvar, iterns, f_id, iescap, imucpp
      real(kind=c_double), value :: normp
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      type(c_ptr), value :: vcopt
      real(kind=c_double), dimension(*), intent(in) :: pvara, pvark, coefap
      real(kind=c_double), dimension(*), intent(in) :: coefbp, cofafp, cofbfp
      real(kind=c_double), dimension(*), intent(in) :: i_massflux, b_massflux
      real(kind=c_double), dimension(*), intent(in) :: i_viscm, b_viscm
      real(kind=c_double), dimension(*), intent(in) :: i_visc, b_visc, viscel
      real(kind=c_double), dimension(*), intent(in) :: weighf, weighb
      integer(c_int), value :: icvflb
      integer(c_int), dimension(*), intent(in) :: icvfli
      real(kind=c_double), dimension(*), intent(in) :: rovsdt, xcpp
      real(kind=c_double), dimension(*), intent(inout) :: smbrp, pvar, dpvar
      real(kind=c_double), dimension(*), intent(inout) :: eswork
    end subroutine cs_equation_iterative_solve_scalar

    !---------------------------------------------------------------------------

    ! Interface to C function cs_equation_iterative_solve_vector

    subroutine cs_equation_iterative_solve_vector(idtvar, iterns,             &
                                                  f_id, name,                 &
                                                  ivisep, iescap,             &
                                                  vcopt, pvara, pvark,        &
                                                  coefav, coefbv, cofafv,     &
                                                  cofbfv, i_massflux,         &
                                                  b_massflux, i_viscm,        &
                                                  b_viscm, i_visc, b_visc,    &
                                                  secvif, secvib,             &
                                                  viscce, weighf, weighb,     &
                                                  icvflb, icvfli, fimp,       &
                                                  smbrp, pvar, eswork)        &
      bind(C, name='cs_equation_iterative_solve_vector')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: idtvar, iterns, f_id, iescap, ivisep
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      type(c_ptr), value :: vcopt
      real(kind=c_double), dimension(*), intent(in) :: pvara, pvark, coefav
      real(kind=c_double), dimension(*), intent(in) :: coefbv, cofafv, cofbfv
      real(kind=c_double), dimension(*), intent(in) :: i_massflux, b_massflux
      real(kind=c_double), dimension(*), intent(in) :: i_visc, b_visc
      real(kind=c_double), dimension(*), intent(in) :: i_viscm, b_viscm
      real(kind=c_double), dimension(*), intent(in) :: secvif, secvib
      real(kind=c_double), dimension(*), intent(in) :: viscce
      real(kind=c_double), dimension(*), intent(in) :: weighf, weighb
      integer(c_int), value :: icvflb
      integer(c_int), dimension(*), intent(in) :: icvfli
      real(kind=c_double), dimension(*), intent(in) :: fimp
      real(kind=c_double), dimension(*), intent(inout) :: smbrp, pvar, eswork
    end subroutine cs_equation_iterative_solve_vector

    !---------------------------------------------------------------------------

    ! Interface to C function cs_equation_iterative_solve_tensor

    subroutine cs_equation_iterative_solve_tensor(idtvar, f_id, name,         &
                                                  vcopt, pvara, pvark,        &
                                                  coefats, coefbts, cofafts,  &
                                                  cofbfts, i_massflux,        &
                                                  b_massflux, i_viscm,        &
                                                  b_viscm, i_visc, b_visc,    &
                                                  viscce, weighf, weighb,     &
                                                  icvflb, icvfli,             &
                                                  fimp, smbrp, pvar)          &
      bind(C, name='cs_equation_iterative_solve_tensor')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: idtvar, f_id
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      type(c_ptr), value :: vcopt
      real(kind=c_double), dimension(*), intent(in) :: pvara, pvark, coefats
      real(kind=c_double), dimension(*), intent(in) :: coefbts, cofafts, cofbfts
      real(kind=c_double), dimension(*), intent(in) :: i_massflux, b_massflux
      real(kind=c_double), dimension(*), intent(in) :: i_visc, b_visc
      real(kind=c_double), dimension(*), intent(in) :: i_viscm, b_viscm
      real(kind=c_double), dimension(*), intent(in) :: viscce
      real(kind=c_double), dimension(*), intent(in) :: weighf, weighb
      integer(c_int), value :: icvflb
      integer(c_int), dimension(*), intent(in) :: icvfli
      real(kind=c_double), dimension(*), intent(in) :: fimp
      real(kind=c_double), dimension(*), intent(inout) :: smbrp, pvar
    end subroutine cs_equation_iterative_solve_tensor

    !---------------------------------------------------------------------------

    ! Interface to C function cs_balance_scalar

    subroutine cs_balance_scalar(idtvar, f_id , imucpp, imasac, inc,          &
                                 iccocg, vcopt , pvar , pvara,                &
                                 coefap, coefbp, cofafp, cofbfp, i_massflux,  &
                                 b_massflux, i_visc, b_visc, viscel, xcpp,    &
                                 weighf, weighb, icvflb, icvfli,              &
                                 smbrp)                                       &
      bind(C, name='cs_balance_scalar')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: idtvar, f_id, imasac, imucpp, inc
      integer(c_int), value :: iccocg
      type(c_ptr), value :: vcopt
      real(kind=c_double), dimension(*), intent(in) :: pvar, pvara, coefap
      real(kind=c_double), dimension(*), intent(in) :: coefbp, cofafp, cofbfp
      real(kind=c_double), dimension(*), intent(in) :: i_massflux, b_massflux
      real(kind=c_double), dimension(*), intent(in) :: i_visc, b_visc, viscel
      real(kind=c_double), dimension(*), intent(in) :: weighf, weighb, xcpp
      integer(c_int), value :: icvflb
      integer(c_int), dimension(*), intent(in) :: icvfli
      real(kind=c_double), dimension(*), intent(inout) :: smbrp
    end subroutine cs_balance_scalar

    !---------------------------------------------------------------------------

    ! Interface to C function cs_balance_vector

    subroutine cs_balance_vector(idtvar, f_id, imasac, inc, ivisep,          &
                                 vcopt, pvar, pvara, coefav, coefbv, cofafv, &
                                 cofbfv, i_massflux, b_massflux, i_visc,     &
                                 b_visc, secvif, secvib, viscel,             &
                                 weighf, weighb, icvflb, icvfli,             &
                                 smbrp)                                      &
      bind(C, name='cs_balance_vector')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: idtvar, f_id, imasac, inc
      integer(c_int), value :: ivisep
      type(c_ptr), value :: vcopt
      real(kind=c_double), dimension(*), intent(in) :: pvar, pvara, coefav
      real(kind=c_double), dimension(*), intent(in) :: coefbv, cofafv, cofbfv
      real(kind=c_double), dimension(*), intent(in) :: i_massflux, b_massflux
      real(kind=c_double), dimension(*), intent(in) :: i_visc, b_visc, viscel
      real(kind=c_double), dimension(*), intent(in) :: secvif, secvib
      real(kind=c_double), dimension(*), intent(in) :: weighf, weighb
      integer(c_int), value :: icvflb
      integer(c_int), dimension(*), intent(in) :: icvfli
      real(kind=c_double), dimension(*), intent(inout) :: smbrp
    end subroutine cs_balance_vector

    !---------------------------------------------------------------------------

    !> \brief Read lagrangian moments checkpoint information.

    subroutine lagr_moment_restart_read(r)  &
      bind(C, name='cs_lagr_moment_restart_read')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
    end subroutine lagr_moment_restart_read

    !---------------------------------------------------------------------------

    !> \brief  Add a species field to the gas mix (set of fields).

    subroutine gas_mix_add_species(f_id)  &
      bind(C, name='cs_gas_mix_add_species')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
    end subroutine gas_mix_add_species

    !---------------------------------------------------------------------------

    !> \brief  Free array mapping gas mix species ids to field ids.

    subroutine finalize_gas_mix()  &
      bind(C, name='cs_gas_mix_finalize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine finalize_gas_mix

    !---------------------------------------------------------------------------

    !> \brief Create global 1d wall thermal model structure.

    subroutine cs_1d_wall_thermal_create()  &
      bind(C, name='cs_1d_wall_thermal_create')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_1d_wall_thermal_create

    !---------------------------------------------------------------------------

    !> \brief Allocate the array of structures local_models.

    subroutine init_1d_wall_thermal_local_models()  &
      bind(C, name='cs_1d_wall_thermal_local_models_create')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine init_1d_wall_thermal_local_models

    !---------------------------------------------------------------------------

    !> \brief Allocate the discretization points coordinates array and
    !>        the temperature at each point of discretization.

    subroutine init_1d_wall_thermal_local_models_arrays()  &
      bind(C, name='cs_1d_wall_thermal_local_models_init')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine init_1d_wall_thermal_local_models_arrays

    !---------------------------------------------------------------------------

    !> \brief Create the 1D mesh for each face and initialize the temperature.

    subroutine cs_1d_wall_thermal_mesh_create()  &
      bind(C, name='cs_1d_wall_thermal_mesh_create')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_1d_wall_thermal_mesh_create

    !---------------------------------------------------------------------------

    !> \brief Solve the 1D equation for a given face.

    !> \param[in]   ii   face number
    !> \param[in]   tf   fluid temperature at the boundarys
    !> \param[in]   hf   exchange coefficient for the fluid

    subroutine cs_1d_wall_thermal_solve(ii, tf, hf)  &
      bind(C, name='cs_1d_wall_thermal_solve')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: ii
      real(kind=c_double), value :: tf, hf
    end subroutine cs_1d_wall_thermal_solve

    !---------------------------------------------------------------------------

    !> \brief Read the restart file of the 1D-wall thermal module.

    subroutine cs_1d_wall_thermal_read()  &
      bind(C, name='cs_1d_wall_thermal_read')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_1d_wall_thermal_read

    !---------------------------------------------------------------------------

    !> \brief Write the restart file of the 1D-wall thermal module.

    subroutine cs_1d_wall_thermal_write()  &
      bind(C, name='cs_1d_wall_thermal_write')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_1d_wall_thermal_write

    !---------------------------------------------------------------------------

    !> \brief Free members of the global 1d wall thermal structure.

    subroutine cs_1d_wall_thermal_free()  &
      bind(C, name='cs_1d_wall_thermal_free')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_1d_wall_thermal_free

    !---------------------------------------------------------------------------

    !> \brief Destroy the global 1d wall thermal structure.

    subroutine cs_1d_wall_thermal_finalize()  &
      bind(C, name='cs_1d_wall_thermal_finalize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_1d_wall_thermal_finalize

    !---------------------------------------------------------------------------

    !> \brief Data Entry of the 1D wall thermal module.

    !> \param[in]   iappel   Call number
    !> \param[in]   isuit1   Restart caculation or not

    subroutine cs_user_1d_wall_thermal(iappel, isuit1)  &
      bind(C, name='cs_user_1d_wall_thermal')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: iappel, isuit1
    end subroutine cs_user_1d_wall_thermal

    !---------------------------------------------------------------------------

    !> \brief Return pointers to nfpt1d and nfpt1t.

    !> \param[out]   nfpt1d   Pointer to nfpt1d
    !> \param[out]   nfpt1t   Pointer to nfpt1t

    subroutine cs_f_1d_wall_thermal_get_pointers(nfpt1d, nfpt1t) &
      bind(C, name='cs_f_1d_wall_thermal_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: nfpt1d, nfpt1t
    end subroutine cs_f_1d_wall_thermal_get_pointers

    !---------------------------------------------------------------------------

    !> \brief Data checking for the 1D thermal wall module.

    !> \param[in]   iappel   Call number
    !> \param[in]   isuit1   Restart caculation or not

    subroutine cs_1d_wall_thermal_check(iappel, isuit1) &
      bind(C, name='cs_1d_wall_thermal_check')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: iappel, isuit1
    end subroutine cs_1d_wall_thermal_check

    !---------------------------------------------------------------------------

    !> \brief Return a pointer to the ifpt1d array.

    !> \param[out]   ifpt1d   Pointer to ifpt1d

    subroutine cs_f_1d_wall_thermal_get_faces(ifpt1d)  &
      bind(C, name='cs_f_1d_wall_thermal_get_faces')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ifpt1d
    end subroutine cs_f_1d_wall_thermal_get_faces

    !---------------------------------------------------------------------------

    !> \brief Return a pointer to the tppt1d array.

    !> \param[out]   tppt1d   Pointer to tppt1d

    subroutine cs_f_1d_wall_thermal_get_temp(tppt1d)  &
      bind(C, name='cs_f_1d_wall_thermal_get_temp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: tppt1d
    end subroutine cs_f_1d_wall_thermal_get_temp

    !---------------------------------------------------------------------------

    !> \brief  Binding to cs_user_internal_coupling

    subroutine cs_user_internal_coupling()  &
      bind(C, name='cs_user_internal_coupling')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_user_internal_coupling

    !---------------------------------------------------------------------------

    !> \brief  Binding to cs_user_internal_coupling

    subroutine cs_internal_coupling_setup()  &
      bind(C, name='cs_internal_coupling_setup')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_internal_coupling_setup

    !---------------------------------------------------------------------------

    !> \brief Binding to cs_ic_field_set_exchcoeff

    !> \param[in]   field_id   field id
    !> \param[in]   hbnd       boundary exchange coefficients passed by face id

    subroutine cs_ic_field_set_exchcoeff(field_id, &
                                         hbnd)     &
      bind(C, name='cs_ic_field_set_exchcoeff')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: field_id
      real(kind=c_double), dimension(*), intent(in) :: hbnd
    end subroutine cs_ic_field_set_exchcoeff

    !---------------------------------------------------------------------------

    !  Binding to cs_f_ic_field_coupled_faces

    subroutine cs_f_ic_field_coupled_faces(f_id, c_p)  &
      bind(C, name='cs_f_ic_field_coupled_faces')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: f_id
      type(c_ptr), intent(out)   :: c_p
    end subroutine cs_f_ic_field_coupled_faces

    !---------------------------------------------------------------------------

    !> \brief Binding to cs_ic_field_dist_data_by_face_id

    !> \param[in]  field_id    field id
    !> \param[in]  stride      number of values (interlaced) by entity
    !> \param[in]  tab_distant exchanged data by face id
    !> \param[out] tab_local   local data by face id

    subroutine cs_ic_field_dist_data_by_face_id(field_id,    &
                                                stride,      &
                                                tab_distant, &
                                                tab_local)   &
      bind(C, name='cs_ic_field_dist_data_by_face_id')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: field_id, stride
      real(kind=c_double), dimension(*), intent(in) :: tab_distant
      real(kind=c_double), dimension(*), intent(out) :: tab_local
    end subroutine cs_ic_field_dist_data_by_face_id

    !---------------------------------------------------------------------------

    !> \brief  Binding to cs_internal_coupling_dump

    subroutine cs_internal_coupling_dump()  &
      bind(C, name='cs_internal_coupling_dump')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_internal_coupling_dump

    !---------------------------------------------------------------------------

    ! Interface to C function handling BCs for internal coupling

    subroutine cs_internal_coupling_bcs(bc_type)                &
        bind(C, name='cs_internal_coupling_bcs')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), dimension(*) :: bc_type
    end subroutine cs_internal_coupling_bcs

    !---------------------------------------------------------------------------

    !> \brief  Check calculation parameters.

    subroutine parameters_check() &
      bind(C, name='cs_parameters_check')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine parameters_check

    !---------------------------------------------------------------------------

    !> \brief Return pointers to atmo includes

    !> \param[out]   compute_z_ground     Pointer to compute_z_ground

    subroutine cs_f_atmo_get_pointers(compute_z_ground) &
      bind(C, name='cs_f_atmo_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: compute_z_ground
    end subroutine cs_f_atmo_get_pointers

    !---------------------------------------------------------------------------

    !> \brief Compute the relative ground elevation (mainly for the atmospheric
    !>  module).

    subroutine cs_atmo_z_ground_compute() &
      bind(C, name='cs_atmo_z_ground_compute')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_atmo_z_ground_compute

    !---------------------------------------------------------------------------

    ! Interface to C user function for cooling tower

    subroutine cs_ctwr_field_pointer_map()  &
      bind(C, name='cs_ctwr_field_pointer_map')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_ctwr_field_pointer_map

    !---------------------------------------------------------------------------

    ! Interface to C function for cooling towers

    subroutine cs_ctwr_init_field_vars(rho0, t0, p0, molmassrat)           &
      bind(C, name='cs_ctwr_init_field_vars')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), value :: rho0, t0, p0, molmassrat
    end subroutine cs_ctwr_init_field_vars

    !---------------------------------------------------------------------------

    ! Interface to C function for cooling towers

    subroutine cs_ctwr_restart_field_vars(rho0, t0, p0, humidity0, molmassrat)  &
      bind(C, name='cs_ctwr_restart_field_vars')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), value :: rho0, t0, p0, humidity0, molmassrat
    end subroutine cs_ctwr_restart_field_vars

    !---------------------------------------------------------------------------

    ! Interface to C function for cooling towers

    subroutine cs_ctwr_phyvar_update(rho0, t0, p0, molmassrat)             &
      bind(C, name='cs_ctwr_phyvar_update')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), value :: rho0, t0, p0, molmassrat
    end subroutine cs_ctwr_phyvar_update

    !---------------------------------------------------------------------------

    ! Interface to C function for cooling towers

    subroutine cs_ctwr_init_flow_vars(liq_mass_flow)                         &
      bind(C, name='cs_ctwr_init_flow_vars')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(inout) :: liq_mass_flow
    end subroutine cs_ctwr_init_flow_vars

    !---------------------------------------------------------------------------

    ! Interface to C function for cooling towers

    subroutine cs_ctwr_bulk_mass_source_term(p0, molmassrat,                 &
                                             mass_source)                    &
      bind(C, name='cs_ctwr_bulk_mass_source_term')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), value :: p0, molmassrat
      real(kind=c_double), dimension(*), intent(inout) :: mass_source
    end subroutine cs_ctwr_bulk_mass_source_term

    !---------------------------------------------------------------------------

    ! Interface to C function for cooling towers

    subroutine cs_ctwr_source_term(f_id, p0, molmassrat,           &
                                   exp_st, imp_st)                 &
      bind(C, name='cs_ctwr_source_term')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      real(kind=c_double), value :: p0, molmassrat
      real(kind=c_double), dimension(*), intent(inout) :: exp_st
      real(kind=c_double), dimension(*), intent(inout) :: imp_st
    end subroutine cs_ctwr_source_term

    !---------------------------------------------------------------------------

    ! Interface to C function for head losses

    subroutine cs_head_losses_compute(ckupdc)  &
      bind(C, name='cs_head_losses_compute')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: ckupdc
    end subroutine cs_head_losses_compute

    !---------------------------------------------------------------------------

    ! Interface to C function cs_f_math_sym_33_inv_cramer

    subroutine symmetric_matrix_inverse(s, sout)                    &
      bind(C, name='cs_f_math_sym_33_inv_cramer')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(in) :: s
      real(kind=c_double), dimension(*), intent(out) :: sout
    end subroutine symmetric_matrix_inverse

    !---------------------------------------------------------------------------

    ! Interface to C function cs_f_math_sym_33_product

    subroutine symmetric_matrix_product(s1, s2, sout)               &
      bind(C, name='cs_f_math_sym_33_product')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(in) :: s1, s2
      real(kind=c_double), dimension(*), intent(out) :: sout
    end subroutine symmetric_matrix_product

    !---------------------------------------------------------------------------

    ! Interface to C function cs_f_math_reduce_symprod33_to_66

    subroutine reduce_symprod33_to_66(s, sout)                      &
      bind(C, name='cs_f_math_reduce_sym_prod_33_to_66')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(in) :: s
      real(kind=c_double), dimension(*), intent(out) :: sout
    end subroutine reduce_symprod33_to_66

    !---------------------------------------------------------------------------

    ! Interface to C function cs_math_sym_33_eigen

    subroutine calc_symtens_eigvals(m, eig_vals)                   &
      bind(C, name='cs_math_sym_33_eigen')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(in) :: m
      real(kind=c_double), dimension(*), intent(out) :: eig_vals
    end subroutine calc_symtens_eigvals

    !---------------------------------------------------------------------------

    ! Interface to C function for data assimilation (atmospheric module)

    subroutine cs_at_data_assim_initialize()                        &
      bind(C, name='cs_at_data_assim_initialize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_at_data_assim_initialize

    !---------------------------------------------------------------------------

    ! Interface to C function for data assimilation (atmospheric module)

    subroutine cs_at_data_assim_build_ops()                        &
      bind(C, name='cs_at_data_assim_build_ops')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_at_data_assim_build_ops

    !---------------------------------------------------------------------------

    ! Interface to C function for data assimilation (atmospheric module)

    function cs_at_opt_interp_is_p1_proj_needed() result (ineeded)   &
      bind(C, name='cs_at_opt_interp_is_p1_proj_needed')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: ineeded
    end function cs_at_opt_interp_is_p1_proj_needed

    !---------------------------------------------------------------------------

    ! Interface to C function for data assimilation (atmospheric module)

    subroutine cs_at_data_assim_finalize()                        &
      bind(C, name='cs_at_data_assim_finalize')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_at_data_assim_finalize

    !---------------------------------------------------------------------------

    ! Interface to C function for data assimilation (atmospheric module).

    subroutine cs_at_data_assim_source_term(f_id, exp_st, imp_st)   &
      bind(C, name='cs_at_data_assim_source_term')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      real(kind=c_double), dimension(*), intent(inout) :: exp_st
      real(kind=c_double), dimension(*), intent(inout) :: imp_st
    end subroutine cs_at_data_assim_source_term

    !---------------------------------------------------------------------------

    ! Binding to cs_ic_set_temp

    subroutine cs_ic_set_temp(field_id, theipb,       &
                             temp_neig)              &
      bind(C, name='cs_ic_set_temp')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: field_id
      real(kind=c_double), dimension(*), intent(in) :: theipb, temp_neig
    end subroutine cs_ic_set_temp

    !---------------------------------------------------------------------------

    ! Set porosity model.

    subroutine cs_mesh_quantities_set_porous_model(iporos)   &
      bind(C, name='cs_mesh_quantities_set_porous_model')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: iporos
    end subroutine cs_mesh_quantities_set_porous_model

    !---------------------------------------------------------------------------

    ! Read turbomachinery metadata from restart file.

    subroutine turbomachinery_restart_read(r)  &
      bind(C, name='cs_turbomachinery_restart_read')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
    end subroutine turbomachinery_restart_read

    !---------------------------------------------------------------------------

    ! Write turbomachinery metadata from restart file.

    subroutine turbomachinery_restart_write(r)  &
      bind(C, name='cs_turbomachinery_restart_write')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: r
    end subroutine turbomachinery_restart_write

    !---------------------------------------------------------------------------

    ! Interface to C function for sorbed concentration update.

    subroutine cs_gwf_sorbed_concentration_update(f_id)  &
      bind(C, name='cs_gwf_sorbed_concentration_update')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
    end subroutine cs_gwf_sorbed_concentration_update

    !---------------------------------------------------------------------------

    ! Interface to C function for precipitation treatment.

    subroutine cs_gwf_precipitation(f_id)   &
      bind(C, name='cs_gwf_precipitation')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
    end subroutine cs_gwf_precipitation

    !---------------------------------------------------------------------------

    ! Interface to C function for decay treatment.

    subroutine cs_gwf_decay_rate(f_id, ts_imp)   &
      bind(C, name='cs_gwf_decay_rate')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      real(kind=c_double), dimension(*), intent(inout) :: ts_imp
    end subroutine cs_gwf_decay_rate

    !---------------------------------------------------------------------------

    ! Interface to C function for kinetic reaction.

    subroutine cs_gwf_kinetic_reaction(f_id, ts_imp, ts_exp)   &
      bind(C, name='cs_gwf_kinetic_reaction')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      real(kind=c_double), dimension(*), intent(inout) :: ts_imp
      real(kind=c_double), dimension(*), intent(inout) :: ts_exp
    end subroutine cs_gwf_kinetic_reaction

    !---------------------------------------------------------------------------

    ! Interface to C function to count number of buoyant scalars.

    subroutine cs_parameters_set_n_buoyant_scalars()   &
      bind(C, name='cs_parameters_set_n_buoyant_scalars')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_parameters_set_n_buoyant_scalars

    !---------------------------------------------------------------------------

    ! Interface to C function updating mesh quantities in the ALE framework.

    subroutine cs_ale_update_mesh_quantities(min_vol, max_vol, tot_vol)   &
      bind(C, name='cs_ale_update_mesh_quantities')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), intent(inout) :: min_vol, max_vol, tot_vol
    end subroutine cs_ale_update_mesh_quantities

    !---------------------------------------------------------------------------

    ! Interface to C function updating the mesh in the ALE framework.

    subroutine cs_ale_update_mesh(itrale, xyzno0)   &
      bind(C, name='cs_ale_update_mesh')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: itrale
      real(kind=c_double), dimension(*), intent(in) :: xyzno0
    end subroutine cs_ale_update_mesh

    !---------------------------------------------------------------------------

    ! Interface to C function updating BCs in ALE framework.

    subroutine cs_ale_update_bcs(ialtyb, b_fluid_vel)   &
      bind(C, name='cs_ale_update_bcs')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), dimension(*), intent(in) :: ialtyb
      real(kind=c_double), dimension(*), intent(in) :: b_fluid_vel
    end subroutine cs_ale_update_bcs

    !---------------------------------------------------------------------------

    ! Interface to C function solving mesh velocity in ALE framework.

    subroutine cs_ale_solve_mesh_velocity(iterns, impale, ialtyb)   &
      bind(C, name='cs_ale_solve_mesh_velocity')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: iterns
      integer(c_int), dimension(*), intent(in) :: impale, ialtyb
    end subroutine cs_ale_solve_mesh_velocity

    !---------------------------------------------------------------------------

    !> \brief  Binding to cs_ale_activate

    subroutine cs_ale_activate()  &
      bind(C, name='cs_ale_activate')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_ale_activate

    !---------------------------------------------------------------------------

    ! Interface to C function to get the cs_glob_ale option

    subroutine cs_f_ale_get_pointers(iale) &
      bind(C, name='cs_f_ale_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iale
    end subroutine cs_f_ale_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function computing total, min, and max cell fluid volumes

    subroutine cs_f_mesh_quantities_fluid_vol_reductions()  &
      bind(C, name='cs_f_mesh_quantities_fluid_vol_reductions')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_mesh_quantities_fluid_vol_reductions

    !---------------------------------------------------------------------------

    ! Interface to C function to get notebook parameter value

    function cs_f_notebook_parameter_value_by_name(name) result(val) &
        bind(C, name='cs_notebook_parameter_value_by_name')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name
      real(kind=c_double) :: val
    end function cs_f_notebook_parameter_value_by_name

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Compute balance on a given zone for a given scalar

  !> param[in]       sel_crit   selection criteria of a volume zone
  !> param[in]       name       scalar name

  subroutine balance_by_zone(sel_crit, name)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)             :: sel_crit, name

    ! Local variables

    character(len=len_trim(sel_crit)+1, kind=c_char) :: c_sel_crit
    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_sel_crit = trim(sel_crit)//c_null_char
    c_name = trim(name)//c_null_char

    call cs_balance_by_zone(c_sel_crit, c_name)

    return

  end subroutine balance_by_zone

  !=============================================================================

  !> \brief Compute pressure drop for a given zone

  !> param[in]       sel_crit   selection criteria of a volume zone

  subroutine pressure_drop_by_zone(sel_crit)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)             :: sel_crit

    ! Local variables

    character(len=len_trim(sel_crit)+1, kind=c_char) :: c_sel_crit

    c_sel_crit = trim(sel_crit)//c_null_char

    call cs_pressure_drop_by_zone(c_sel_crit)

    return

  end subroutine pressure_drop_by_zone

  !=============================================================================

  ! Interface to C function returning the product of a matrix (native format)
  ! by a vector

  subroutine promav(isym, ibsize, iesize, f_id, dam, xam, vx, vy)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, value :: isym, ibsize, iesize, f_id
    real(kind=c_double), dimension(*), intent(in) :: dam, xam, vx
    real(kind=c_double), dimension(*), intent(out) :: vy

    ! Local variables

    integer(c_int), dimension(4) :: c_db_size, c_eb_size
    integer(c_int) :: c_rotation_mode
    logical(c_bool) :: c_symmetric

    if (isym.eq.1) then
      c_symmetric = .true.
    else
      c_symmetric = .false.
    endif

    c_rotation_mode = 0 ! CS_HALO_ROTATION_COPY

    c_db_size(0+1) = ibsize;
    c_db_size(1+1) = ibsize;
    c_db_size(2+1) = ibsize;
    c_db_size(3+1) = ibsize*ibsize;

    c_eb_size(0+1) = iesize;
    c_eb_size(1+1) = iesize;
    c_eb_size(2+1) = iesize;
    c_eb_size(3+1) = iesize*iesize;

    call cs_matrix_vector_native_multiply(c_symmetric, c_db_size, c_eb_size, &
                                          c_rotation_mode, f_id, dam, xam, vx, vy)

    return

  end subroutine promav

  !=============================================================================

  !> \brief Compute surface scalar balance for a given surface area

  !> param[in]       sel_crit   selection criteria of a volume zone

  subroutine surface_balance(sel_crit, name, normal)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)             :: sel_crit, name
    real(kind=c_double), dimension(3), intent(in) :: normal

    ! Local variables

    character(len=len_trim(sel_crit)+1, kind=c_char) :: c_sel_crit
    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_sel_crit = trim(sel_crit)//c_null_char
    c_name = trim(name)//c_null_char

    call cs_surface_balance(c_sel_crit, c_name, normal)

    return

  end subroutine surface_balance

  !=============================================================================

  !> \brief Handle boundary condition definition errors and associated output.

  !> For each boundary face, bc_type defines the boundary condition type.
  !> As a convention here, zero values correspond to undefined types,
  !> positive values to defined types (with no error), and negative values
  !> to defined types with inconsistent or incompatible values, the
  !> absolute value indicating the original boundary condition type.

  !> param[in]  bc_type    array og BC type ids

  subroutine boundary_conditions_error(bc_type)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer(c_int), dimension(*), intent(in) :: bc_type

    ! Call C function with default name

    call cs_boundary_conditions_error(bc_type, c_null_ptr)

  end subroutine boundary_conditions_error

  !=============================================================================

  !> \brief Locate shifted boundary face coordinates on possibly filtered
  !>        cells or boundary faces for later interpolation.

  !> param[in]  location_type    matching values location (CS_MESH_LOCATION_CELLS
  !>                             or CS_MESH_LOCATION_BOUNDARY_FACES)
  !> param[in]  n_location_elts  number of selected location elements
  !> param[in]  n_faces          number of selected boundary faces
  !> param[in]  location_elts    list of selected location elements (1 to n),
  !>                             or NULL if no indirection is needed
  !> param[in]  faces            list of selected boundary faces (1 to n),
  !>                             or NULL if no indirection is needed
  !> param[in]  coord_shift      array of coordinates shift relative to selected
  !>                             boundary faces
  !> param[in]  coord_stride     access stride in coord_shift: 0 for uniform
  !>                             shift, 1 for "per face" shift.
  !> param[in]  tolerance        relative tolerance for point location.

  !> return  associated locator structure

  function boundary_conditions_map(location_type, n_location_elts,           &
                                   n_faces, location_elts, faces,            &
                                   coord_shift, coord_stride,                &
                                   tolerance) result(l)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in) :: location_type, n_location_elts, n_faces
    integer, dimension(*), intent(in) :: location_elts, faces
    real(kind=c_double), dimension(*) :: coord_shift
    integer, intent(in) :: coord_stride
    double precision, intent(in) :: tolerance
    type(c_ptr) :: l

    ! Local variables

    integer iel, ifac
    integer(c_int) :: c_loc_type, c_n_elts, c_n_faces, c_coord_stride
    integer(c_int), dimension(:), allocatable :: c_loc_elts, c_faces
    real(kind=c_double) :: c_tolerance

    c_loc_type = location_type
    c_n_elts = n_location_elts
    c_n_faces = n_faces
    c_coord_stride = coord_stride
    c_tolerance = tolerance

    allocate(c_loc_elts(n_location_elts))
    allocate(c_faces(n_faces))

    do iel = 1, n_location_elts
      c_loc_elts(iel) = location_elts(iel) - 1
    enddo
    do ifac = 1, n_faces
      c_faces(ifac) = faces(ifac) - 1
    enddo

    l = cs_boundary_conditions_map(c_loc_type, c_n_elts, c_n_faces,          &
                                   c_loc_elts, c_faces,                      &
                                   coord_shift, c_coord_stride, c_tolerance)

    deallocate(c_faces)
    deallocate(c_loc_elts)

  end function boundary_conditions_map

  !=============================================================================

  !> \brief Assign a var_cal_opt for a cs_var_cal_opt_t key to a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_value  structure associated with key

  subroutine field_set_key_struct_var_cal_opt(f_id, k_value)

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

  !> \brief Assign a gwf_soilwater_partition for a cs_gwf_soilwater_partition_t
  !> key to a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[in]   k_value  structure associated with key

  subroutine field_set_key_struct_gwf_soilwater_partition(f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                          :: f_id
    type(gwf_soilwater_partition), intent(in), target :: k_value

    ! Local variables

    integer(c_int)                        :: c_f_id
    type(gwf_soilwater_partition),pointer      :: p_k_value
    type(c_ptr)                           :: c_k_value
    character(len=34+1, kind=c_char)      :: c_name

    integer(c_int), save           :: c_k_id = -1

    if (c_k_id .eq. -1) then
      c_name = "gwf_soilwater_partition"//c_null_char
      c_k_id = cs_f_field_key_id(c_name)
    endif

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_set_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_set_key_struct_gwf_soilwater_partition

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

  subroutine field_get_key_struct_var_cal_opt(f_id, k_value)

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

  subroutine field_get_key_struct_solving_info(f_id, k_value)

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

  !> \brief Return a pointer to the gwf_soilwater_partition structure for
  !>        cs_gwf_soilwater_partition_t key associated with a field.

  !> If the field category is not compatible, a fatal error is provoked.

  !> \param[in]   f_id     field id
  !> \param[out]  k_value  integer value associated with key id for this field

  subroutine field_get_key_struct_gwf_soilwater_partition(f_id, k_value)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)                             :: f_id
    type(gwf_soilwater_partition), intent(inout), target :: k_value

    ! Local variables

    integer(c_int)                        :: c_f_id, c_k_id
    type(gwf_soilwater_partition),pointer      :: p_k_value
    type(c_ptr)                           :: c_k_value
    character(len=34+1, kind=c_char)      :: c_name

    c_name = "gwf_soilwater_partition"//c_null_char

    c_k_id = cs_f_field_key_id(c_name)
    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_get_key_struct(c_f_id, c_k_id, c_k_value)

    return

  end subroutine field_get_key_struct_gwf_soilwater_partition

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
  !> \param[in]       epsrgp           relative precision for reconstruction
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
  !> \param[in]       epsrgp           relative precision for reconstruction
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
  !> \param[in]       hyd_p_flag       flag for hydrostatic pressure
  !> \param[in]       iwarnp           verbosity
  !> \param[in]       epsrgp           relative precision for reconstruction
  !> \param[in]       climgp           limiter coefficient for imligp
  !> \param[in]       extrap           gradient extrapolation coefficient
  !> \param[in]       f_ext            exterior force generating
  !>                                   the hydrostatic pressure
  !> \param[in, out]  pvar             cell values whose gradient is computed
  !> \param[in, out]  c_weight         cell weighting coefficient
  !> \param[in]       coefap           boundary coefap coefficients
  !> \param[in]       coefbp           boundary coefap coefficients
  !> \param[out]      grad             resulting gradient

  subroutine gradient_weighted_s(f_id, imrgra, inc, recompute_cocg, nswrgp,   &
                                 imligp, hyd_p_flag, iwarnp, epsrgp, climgp,  &
                                 extrap, f_ext, pvar, c_weight, coefap,       &
                                 coefbp, grad)

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
    real(kind=c_double), dimension(*), intent(in) :: c_weight
    real(kind=c_double), dimension(3, *), intent(in) :: f_ext
    real(kind=c_double), dimension(3, ncelet), intent(out) :: grad

    ! Local variables

    integer          :: idimtr, ipond

    ! The current variable is a scalar
    idimtr = 0

    ! the pressure gradient coefficient weighting is used
    ipond = 1

    call cgdcel(f_id, imrgra, inc, recompute_cocg, nswrgp,                     &
                idimtr, hyd_p_flag, ipond, iwarnp, imligp, epsrgp, extrap,     &
                climgp, f_ext, coefap, coefbp,                                 &
                pvar, c_weight, grad)

  end subroutine gradient_weighted_s

  !=============================================================================

  !> \brief Destruction of a locator structure.

  !> \param[in, out]   this_locator

  subroutine locator_destroy(this_locator)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    type(c_ptr) :: this_locator

    ! Local variables

    this_locator = ple_locator_destroy(this_locator)

  end subroutine locator_destroy

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

  !> \param[in]  f_id            associated dimension (interleaved)
  !> \param[in]  n_clip_min      local number of clipped to min values
  !> \param[in]  n_clip_max      local number of clipped to max values
  !> \param[in]  min_pre_clip    min local value prior to clip
  !> \param[in]  max_pre_clip    max local value prior to clip
  !> \param[in]  n_clip_min_comp number of clip min by component
  !> \param[in]  n_clip_max_comp number of clip max by component

  subroutine log_iteration_clipping_field(f_id, n_clip_min, n_clip_max,        &
                                          min_pre_clip, max_pre_clip,  &
                                          n_clip_min_comp, n_clip_max_comp)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)               :: f_id, n_clip_min, n_clip_max
    real(kind=c_double), dimension(*) :: min_pre_clip, max_pre_clip
    integer(c_int), dimension(*), intent(in) :: n_clip_min_comp, n_clip_max_comp
    ! Local variables

    integer(c_int) :: c_f_id, c_clip_min, c_clip_max

    c_f_id = f_id
    c_clip_min = n_clip_min
    c_clip_max = n_clip_max

    call cs_log_iteration_clipping_field(c_f_id, c_clip_min, c_clip_max, &
                                         min_pre_clip, max_pre_clip, &
                                         n_clip_min_comp, n_clip_max_comp)

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
  !> param[in]       epsilp   precision for iterative resolution
  !> param[in]       rnorm    residue normalization
  !> param[out]      niter    number of "equivalent" iterations
  !> param[out]      residue  residue
  !> param[in]       rhs      right hand side
  !> param[in, out]  vx       system solution

  subroutine sles_solve_native(f_id, name, isym, ibsize, iesize, dam, xam,     &
                               epsilp, rnorm, niter, residue, rhs, vx)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in)      :: name
    integer, intent(in)               :: f_id, isym, ibsize, iesize
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

    rotation_mode = 0 ! CS_HALO_ROTATION_COPY, might not be called

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

  !> \brief Create a timer statistics structure.

  !> If no timer with the given name exists, -1 is returned.

  !> \param[in]  parent_name  name of parent statistic (may be empty)
  !> \param[in]  name         associated canonical name
  !> \param[in]  label        associated label (may be empty)

  function timer_stats_create (parent_name, name, label) result(id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: parent_name, name, label
    integer :: id

    ! Local variables

    character(len=len_trim(parent_name)+1, kind=c_char) :: c_p_name
    character(len=len_trim(name)+1, kind=c_char) :: c_name
    character(len=len_trim(label)+1, kind=c_char) :: c_label
    integer(c_int) :: c_id

    c_p_name = trim(parent_name)//c_null_char
    c_name = trim(name)//c_null_char
    c_label = trim(label)//c_null_char

    c_id = cs_timer_stats_create(c_p_name, c_name, c_label)
    id = c_id

  end function timer_stats_create

  !=============================================================================

  !> \brief Return the id of a defined statistic based on its name.

  !> If no timer with the given name exists, -1 is returned.

  !> \param[in]   name   statistic name

  function timer_stats_id_by_name(name) result(id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name
    integer :: id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    integer(c_int) :: c_id

    c_name = trim(name)//c_null_char

    c_id = cs_timer_stats_id_by_name(c_name)
    id = c_id

  end function timer_stats_id_by_name

  !=============================================================================

  !> \brief  Add field defining a general solved variable, with default options.

  !> \param[in]  name           field name
  !> \param[in]  label          field default label, or empty
  !> \param[in]  location_id    field location type:
  !>                              0: none
  !>                              1: cells
  !>                              2: interior faces
  !>                              3: interior faces
  !>                              4: vertices
  !> \param[in]  dim            field dimension
  !> \param[out] id             id of defined field

  subroutine variable_field_create(name, label, location_id, dim, id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name, label
    integer, intent(in)          :: location_id, dim
    integer, intent(out)         :: id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    character(len=len_trim(label)+1, kind=c_char) :: c_label
    integer(c_int) :: c_location_id, c_dim, c_id

    c_name = trim(name)//c_null_char
    c_label = trim(label)//c_null_char
    c_location_id = location_id
    c_dim = dim

    c_id = cs_variable_field_create(c_name, c_label, c_location_id, c_dim)

    id = c_id

    return

  end subroutine variable_field_create

  !=============================================================================

  !> \brief Return the number of volume zones associated with a given type flag.

  !> \param[in]   type_flag   type flag queried

  function volume_zone_n_type_zones(type_flag) result(n)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer :: type_flag, n

    ! Local variables

    integer(c_int) :: c_type_flag, c_count

    c_type_flag = type_flag
    c_count = cs_volume_zone_n_type_zones(c_type_flag)
    n = c_count

  end function volume_zone_n_type_zones

  !=============================================================================

  !> \brief Return the number of volume zone cells associated with a given
  !>        type flag.

  !> \param[in]   type_flag   type flag queried

  function volume_zone_n_type_cells(type_flag) result(n)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer :: type_flag, n

    ! Local variables

    integer(c_int) :: c_type_flag, c_count

    c_type_flag = type_flag
    c_count = cs_volume_zone_n_type_cells(c_type_flag)
    n = c_count

  end function volume_zone_n_type_cells

  !=============================================================================

  !> \brief Return the list of volume zone cells associated with a given
  !>        type flag.

  !> \param[in]   type_flag   type flag queried
  !> \param[out]  cell_list   list of cells

  subroutine volume_zone_select_type_cells(type_flag, cell_list)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer :: type_flag
    integer, dimension(*), intent(out), target :: cell_list

    ! Local variables

    integer(c_int) :: c_type_flag, c_count, i
    type(c_ptr) :: c_cell_list

    c_type_flag = type_flag
    c_cell_list = c_loc(cell_list)
    c_count = volume_zone_n_type_cells(c_type_flag)
    call cs_volume_zone_select_type_cells(c_type_flag, c_cell_list)
    do i = 1, c_count
      cell_list(i) = cell_list(i) + 1
    enddo

  end subroutine volume_zone_select_type_cells

  !=============================================================================

  !> \brief This function solves an advection diffusion equation with source
  !>        terms for one time step for the variable \f$ a \f$.

  !> The equation reads:
  !>
  !> \f[
  !> f_s^{imp}(a^{n+1}-a^n)
  !> + \divs \left( a^{n+1} \rho \vect{u} - \mu \grad a^{n+1} \right)
  !> = Rhs
  !> \f]
  !>
  !> This equation is rewritten as:
  !>
  !> \f[
  !> f_s^{imp} \delta a
  !> + \divs \left( \delta a \rho \vect{u} - \mu \grad \delta a \right)
  !> = Rhs^1
  !> \f]
  !>
  !> where \f$ \delta a = a^{n+1} - a^n\f$ and
  !> \f$ Rhs^1 = Rhs - \divs( a^n \rho \vect{u} - \mu \grad a^n)\f$
  !>
  !> It is in fact solved with the following iterative process:
  !>
  !> \f[
  !> f_s^{imp} \delta a^k
  !> + \divs \left(\delta a^k \rho \vect{u}-\mu\grad\delta a^k \right)
  !> = Rhs^k
  !> \f]
  !>
  !> where \f$Rhs^k=Rhs-f_s^{imp}(a^k-a^n)
  !> - \divs \left( a^k\rho\vect{u}-\mu\grad a^k \right)\f$

  !> Be careful, it is forbidden to modify \f$ f_s^{imp} \f$ here!

  !> \param[in]     idtvar        indicator of the temporal scheme
  !> \param[in]     iterns        external sub-iteration number
  !> \param[in]     f_id          field id (or -1)
  !> \param[in]     iconvp        indicator
  !>                               - 1 convection,
  !>                               - 0 otherwise
  !> \param[in]     idiffp        indicator
  !>                               - 1 diffusion,
  !>                               - 0 otherwise
  !> \param[in]     ndircp        indicator (0 if the diagonal is stepped aside)
  !> \param[in]     imrgra        indicator
  !>                               - 0 iterative gradient
  !>                               - 1 least squares gradient
  !> \param[in]     nswrsp        number of reconstruction sweeps for the
  !>                               Right Hand Side
  !> \param[in]     nswrgp        number of reconstruction sweeps for the
  !>                               gradients
  !> \param[in]     imligp        clipping gradient method
  !>                               - < 0 no clipping
  !>                               - = 0 thank to neighboring gradients
  !>                               - = 1 thank to the mean gradient
  !> \param[in]     ircflp        indicator
  !>                               - 1 flux reconstruction,
  !>                               - 0 otherwise
  !> \param[in]     ischcp        indicator
  !>                               - 1 centered
  !>                               - 0 2nd order
  !> \param[in]     isstpp        indicator
  !>                               - 1 without slope test
  !>                               - 0 with slope test
  !> \param[in]     iescap        compute the predictor indicator if 1
  !> \param[in]     imucpp        indicator
  !>                               - 0 do not multiply the convectiv term by Cp
  !>                               - 1 do multiply the convectiv term by Cp
  !> \param[in]     idftnp        diffusivity type indicator
  !> \param[in]     iswdyp        indicator
  !>                               - 0 no dynamic relaxation
  !>                               - 1 dynamic relaxation depending on
  !>                                 \f$ \delta \varia^k \f$
  !>                               - 2 dynamic relaxation depending on
  !>                                 \f$ \delta \varia^k \f$  and
  !>                                 \f$ \delta \varia^{k-1} \f$
  !> \param[in]     iwarnp        verbosity
  !> \param[in]     normp         (optional) norm residual
  !> \param[in]     blencp        fraction of upwinding
  !> \param[in]     epsilp        precision pour resol iter
  !> \param[in]     epsrsp        relative precision for the iterative process
  !> \param[in]     epsrgp        relative precision for the gradient
  !>                               reconstruction
  !> \param[in]     climgp        clipping coefficient for the computation of
  !>                               the gradient
  !> \param[in]     extrap        coefficient for extrapolation of the gradient
  !> \param[in]     relaxp        coefficient of relaxation
  !> \param[in]     thetap        weighting coefficient for the theta-schema,
  !>                               - thetap = 0: explicit scheme
  !>                               - thetap = 0.5: time-centered
  !>                               scheme (mix between Crank-Nicolson and
  !>                               Adams-Bashforth)
  !>                               - thetap = 1: implicit scheme
  !> \param[in]     pvara         variable at the previous time step
  !>                               \f$ a^n \f$
  !> \param[in]     pvark         variable at the previous sub-iteration
  !>                               \f$ a^k \f$.
  !>                               If you sub-iter on Navier-Stokes, then
  !>                               it allows to initialize by something else
  !>                               than pvara (usually pvar=pvara)
  !> \param[in]     coefap        boundary condition array for the variable
  !>                               (explicit part)
  !> \param[in]     coefbp        boundary condition array for the variable
  !>                               (implicit part)
  !> \param[in]     cofafp        boundary condition array for the diffusion
  !>                               of the variable (explicit part)
  !> \param[in]     cofbfp        boundary condition array for the diffusion
  !>                               of the variable (implicit part)
  !> \param[in]     i_massflux    mass flux at interior faces
  !> \param[in]     b_massflux    mass flux at boundary faces
  !> \param[in]     i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
  !>                               at interior faces for the matrix
  !> \param[in]     b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
  !>                               at boundary faces for the matrix
  !> \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
  !>                               at interior faces for the r.h.s.
  !> \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
  !>                               at boundary faces for the r.h.s.
  !> \param[in]     viscel        symmetric cell tensor
  !>                              \f$ \tens{\mu}_\celli \f$
  !> \param[in]     weighf        internal face weight between cells i j in case
  !>                               of tensor diffusion
  !> \param[in]     weighb        boundary face weight for cells i in case
  !>                               of tensor diffusion
  !> \param[in]     icvflb        global indicator of boundary convection flux
  !>                               - 0 upwind scheme at all boundary faces
  !>                               - 1 imposed flux at some boundary faces
  !> \param[in]     icvfli        boundary face indicator array of convection
  !>                              flux
  !>                               - 0 upwind scheme
  !>                               - 1 imposed flux
  !> \param[in]     rovsdt        \f$ f_s^{imp} \f$
  !> \param[in]     smbrp         Right hand side \f$ Rhs^k \f$
  !> \param[in,out] pvar          current variable
  !> \param[in,out] dpvar         last variable increment
  !> \param[in]     xcpp          array of specific heat (Cp)
  !> \param[out]    eswork        prediction-stage error estimator
  !>                              (if iescap > 0)

  subroutine codits (idtvar, iterns,                                           &
                     f_id  , iconvp, idiffp, ndircp, imrgra, nswrsp,           &
                     nswrgp, imligp, ircflp, ischcp, isstpp, iescap, imucpp,   &
                     idftnp, iswdyp, iwarnp, normp,                            &
                     blencp, epsilp, epsrsp, epsrgp,                           &
                     climgp, extrap, relaxp, thetap, pvara, pvark, coefap,     &
                     coefbp, cofafp, cofbfp, i_massflux, b_massflux, i_viscm,  &
                     b_viscm, i_visc, b_visc, viscel, weighf, weighb, icvflb,  &
                     icvfli, rovsdt, smbrp, pvar, dpvar, xcpp, eswork)

    use, intrinsic :: iso_c_binding
    use entsor, only:nomva0
    use numvar
    use cplsat

    implicit none

    ! Arguments

    integer, intent(in) :: idtvar, iterns, f_id, iconvp, idiffp, ndircp, imrgra
    integer, intent(in) :: nswrsp, nswrgp, imligp, ircflp, ischcp, isstpp
    integer, intent(in) :: iescap, imucpp, idftnp, iswdyp, iwarnp
    double precision, intent(in) :: normp
    double precision, intent(in) :: blencp, epsilp, epsrsp, epsrgp, climgp
    double precision, intent(in) :: extrap, relaxp, thetap
    real(kind=c_double), dimension(*), intent(in) :: pvara, pvark, coefap
    real(kind=c_double), dimension(*), intent(in) :: coefbp, cofafp, cofbfp
    real(kind=c_double), dimension(*), intent(in) :: i_massflux, b_massflux
    real(kind=c_double), dimension(*), intent(in) :: i_viscm, b_viscm
    real(kind=c_double), dimension(*), intent(in) :: i_visc, b_visc, viscel
    real(kind=c_double), dimension(*), intent(in) :: weighf, weighb
    integer, intent(in) :: icvflb
    integer(c_int), dimension(*), intent(in) :: icvfli
    real(kind=c_double), dimension(*), intent(in) :: rovsdt, xcpp
    real(kind=c_double), dimension(*), intent(inout) :: smbrp, pvar, dpvar
    real(kind=c_double), dimension(*), intent(inout) :: eswork

    ! Local variables
    character(len=len_trim(nomva0)+1, kind=c_char) :: c_name
    type(var_cal_opt), target   :: vcopt
    type(var_cal_opt), pointer  :: p_k_value
    type(c_ptr)                 :: c_k_value

    c_name = trim(nomva0)//c_null_char

    p_k_value => vcopt
    c_k_value = c_loc(p_k_value)

    vcopt%iwarni = iwarnp
    vcopt%iconv  = iconvp
    vcopt%istat  = -1
    vcopt%ndircl = ndircp
    vcopt%idiff  = idiffp
    vcopt%idifft = -1
    vcopt%idften = idftnp
    vcopt%iswdyn = iswdyp
    vcopt%ischcv = ischcp
    vcopt%isstpc = isstpp
    vcopt%nswrgr = nswrgp
    vcopt%nswrsm = nswrsp
    vcopt%imrgra = imrgra
    vcopt%imligr = imligp
    vcopt%ircflu = ircflp
    vcopt%iwgrec = 0 ! Warning, may be overwritten if a field
    vcopt%thetav = thetap
    vcopt%blencv = blencp
    vcopt%blend_st = 0 ! Warning, may be overwritten if a field
    vcopt%epsilo = epsilp
    vcopt%epsrsm = epsrsp
    vcopt%epsrgr = epsrgp
    vcopt%climgr = climgp
    vcopt%extrag = extrap
    vcopt%relaxv = relaxp

    call cs_equation_iterative_solve_scalar(idtvar, iterns,                    &
                                            f_id, c_name,                      &
                                            iescap, imucpp, normp, c_k_value,  &
                                            pvara, pvark,                      &
                                            coefap, coefbp, cofafp, cofbfp,    &
                                            i_massflux, b_massflux,            &
                                            i_viscm, b_viscm, i_visc, b_visc,  &
                                            viscel, weighf, weighb,            &
                                            icvflb, icvfli,                    &
                                            rovsdt, smbrp, pvar, dpvar,        &
                                            xcpp, eswork)

    return

  end subroutine codits

  !=============================================================================
  !> \brief This function solves an advection diffusion equation with source
  !>        terms for one time step for the vector variable \f$ \vect{a} \f$.

  !> The equation reads:
  !>
  !> \f[
  !> \tens{f_s}^{imp}(\vect{a}^{n+1}-\vect{a}^n)
  !> + \divv \left( \vect{a}^{n+1} \otimes \rho \vect {u}
  !>              - \mu \gradt \vect{a}^{n+1}\right)
  !> = \vect{Rhs}
  !> \f]
  !>
  !> This equation is rewritten as:
  !>
  !> \f[
  !> \tens{f_s}^{imp} \delta \vect{a}
  !> + \divv \left( \delta \vect{a} \otimes \rho \vect{u}
  !>              - \mu \gradt \delta \vect{a} \right)
  !> = \vect{Rhs}^1
  !> \f]
  !>
  !> where \f$ \delta \vect{a} = \vect{a}^{n+1} - \vect{a}^n\f$ and
  !> \f$ \vect{Rhs}^1 = \vect{Rhs}
  !> - \divv \left( \vect{a}^n \otimes \rho \vect{u}
  !>              - \mu \gradt \vect{a}^n \right)\f$
  !>
  !> It is in fact solved with the following iterative process:
  !>
  !> \f[
  !> \tens{f_s}^{imp} \delta \vect{a}^k
  !> + \divv \left( \delta \vect{a}^k \otimes \rho \vect{u}
  !>              - \mu \gradt \delta \vect{a}^k \right)
  !> = \vect{Rhs}^k
  !> \f]
  !>
  !> where \f$ \vect{Rhs}^k = \vect{Rhs}
  !> - \tens{f_s}^{imp} \left(\vect{a}^k-\vect{a}^n \right)
  !> - \divv \left( \vect{a}^k \otimes \rho \vect{u}
  !>              - \mu \gradt \vect{a}^k \right)\f$
  !>
  !> Be careful, it is forbidden to modify \f$ \tens{f_s}^{imp} \f$ here!

  !> \param[in]     idtvar        indicator of the temporal scheme
  !> \param[in]     iterns        external sub-iteration number
  !> \param[in]     f_id          field id (or -1)
  !> \param[in]     iconvp        indicator
  !>                               - 1 convection,
  !>                               - 0 otherwise
  !> \param[in]     idiffp        indicator
  !>                               - 1 diffusion,
  !>                               - 0 otherwise
  !> \param[in]     ndircp        indicator (0 if the diagonal is stepped aside)
  !> \param[in]     imrgra        indicateur
  !>                               - 0 iterative gradient
  !>                               - 1 least squares gradient
  !> \param[in]     nswrsp        number of reconstruction sweeps for the
  !>                               Right Hand Side
  !> \param[in]     nswrgp        number of reconstruction sweeps for the
  !>                               gradients
  !> \param[in]     imligp        clipping gradient method
  !>                               - < 0 no clipping
  !>                               - = 0 thank to neighboring gradients
  !>                               - = 1 thank to the mean gradient
  !> \param[in]     ircflp        indicator
  !>                               - 1 flux reconstruction,
  !>                               - 0 otherwise
  !> \param[in]     ivisep        indicator to take \f$ \divv
  !>                               \left(\mu \gradt \transpose{\vect{a}} \right)
  !>                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
  !>                               - 1 take into account,
  !>                               - 0 otherwise
  !> \param[in]     ischcp        indicator
  !>                               - 1 centered
  !>                               - 0 2nd order
  !> \param[in]     isstpp        indicator
  !>                               - 1 without slope test
  !>                               - 0 with slope test
  !> \param[in]     iescap        compute the predictor indicator if 1
  !> \param[in]     idftnp        indicator
  !>                               - 1 the diffusivity is scalar
  !>                               - 6 the diffusivity is a symmetric tensor
  !> \param[in]     iswdyp        indicator
  !>                               - 0 no dynamic relaxation
  !>                               - 1 dynamic relaxation depending on
  !>                                 \f$ \delta \vect{\varia}^k \f$
  !>                               - 2 dynamic relaxation depending on
  !>                                 \f$ \delta \vect{\varia}^k \f$  and
  !>                                 \f$ \delta \vect{\varia}^{k-1} \f$
  !> \param[in]     iwarnp        verbosity
  !> \param[in]     blencp        fraction of upwinding
  !> \param[in]     epsilp        precision pour resol iter
  !> \param[in]     epsrsp        relative precision for the iterative process
  !> \param[in]     epsrgp        relative precision for the gradient
  !>                               reconstruction
  !> \param[in]     climgp        clipping coefficient for the computation of
  !>                               the gradient
  !> \param[in]     relaxp        coefficient of relaxation
  !> \param[in]     thetap        weighting coefficient for the theta-schema,
  !>                               - thetap = 0: explicit scheme
  !>                               - thetap = 0.5: time-centered
  !>                               scheme (mix between Crank-Nicolson and
  !>                               Adams-Bashforth)
  !>                               - thetap = 1: implicit scheme
  !> \param[in]     pvara         variable at the previous time step
  !>                               \f$ \vect{a}^n \f$
  !> \param[in]     pvark         variable at the previous sub-iteration
  !>                               \f$ \vect{a}^k \f$.
  !>                               If you sub-iter on Navier-Stokes, then
  !>                               it allows to initialize by something else
  !>                               than pvara (usually pvar=pvara)
  !> \param[in]     coefav        boundary condition array for the variable
  !>                               (explicit part)
  !> \param[in]     coefbv        boundary condition array for the variable
  !>                               (implicit part)
  !> \param[in]     cofafv        boundary condition array for the diffusion
  !>                               of the variable (Explicit part)
  !> \param[in]     cofbfv        boundary condition array for the diffusion
  !>                               of the variable (Implicit part)
  !> \param[in]     i_massflux    mass flux at interior faces
  !> \param[in]     b_massflux    mass flux at boundary faces
  !> \param[in]     i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
  !>                               at interior faces for the matrix
  !> \param[in]     b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
  !>                               at boundary faces for the matrix
  !> \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
  !>                               at interior faces for the r.h.s.
  !> \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
  !>                               at boundary faces for the r.h.s.
  !> \param[in]     secvif        secondary viscosity at interior faces
  !> \param[in]     secvib        secondary viscosity at boundary faces
  !> \param[in]     viscce        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
  !> \param[in]     weighf        internal face weight between cells i j in case
  !>                               of tensor diffusion
  !> \param[in]     weighb        boundary face weight for cells i in case
  !>                               of tensor diffusion
  !> \param[in]     icvflb        global indicator of boundary convection flux
  !>                               - 0 upwind scheme at all boundary faces
  !>                               - 1 imposed flux at some boundary faces
  !> \param[in]     icvfli        boundary face indicator array of convection flux
  !>                               - 0 upwind scheme
  !>                               - 1 imposed flux
  !> \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
  !> \param[in]     smbrp         Right hand side \f$ \vect{Rhs}^k \f$
  !> \param[in,out] pvar          current variable
  !> \param[out]    eswork        prediction-stage error estimator
  !>                              (if iescap > 0)

  subroutine coditv (idtvar, iterns,                                           &
                     f_id  , iconvp, idiffp, ndircp, imrgra, nswrsp,           &
                     nswrgp, imligp, ircflp, ivisep, ischcp, isstpp, iescap,   &
                     idftnp, iswdyp, iwarnp, blencp, epsilp, epsrsp, epsrgp,   &
                     climgp, relaxp, thetap, pvara , pvark , coefav, coefbv,   &
                     cofafv, cofbfv, i_massflux, b_massflux, i_viscm,          &
                     b_viscm, i_visc, b_visc, secvif, secvib,                  &
                     viscce, weighf, weighb, icvflb, icvfli,                   &
                     fimp, smbrp, pvar, eswork)

    use, intrinsic :: iso_c_binding
    use entsor, only:nomva0
    use numvar
    use cplsat

    implicit none

    ! Arguments

    integer, intent(in) :: idtvar, iterns, f_id, iconvp, idiffp, ndircp, imrgra
    integer, intent(in) :: nswrsp, nswrgp, imligp, ircflp, ischcp, isstpp
    integer, intent(in) :: iescap, ivisep, idftnp, iswdyp, iwarnp
    double precision, intent(in) :: blencp, epsilp, epsrsp, epsrgp, climgp
    double precision, intent(in) :: relaxp, thetap
    real(kind=c_double), dimension(*), intent(in) :: pvara, pvark, coefav
    real(kind=c_double), dimension(*), intent(in) :: coefbv, cofafv, cofbfv
    real(kind=c_double), dimension(*), intent(in) :: i_massflux, b_massflux
    real(kind=c_double), dimension(*), intent(in) :: i_visc, b_visc
    real(kind=c_double), dimension(*), intent(in) :: i_viscm, b_viscm
    real(kind=c_double), dimension(*), intent(in) :: secvif, secvib
    real(kind=c_double), dimension(*), intent(in) :: viscce
    real(kind=c_double), dimension(*), intent(in) :: weighf, weighb
    integer, intent(in) :: icvflb
    integer(c_int), dimension(*), intent(in) :: icvfli
    real(kind=c_double), dimension(*), intent(in) :: fimp
    real(kind=c_double), dimension(*), intent(inout) :: smbrp, pvar, eswork

    ! Local variables
    character(len=len_trim(nomva0)+1, kind=c_char) :: c_name
    type(var_cal_opt), target   :: vcopt
    type(var_cal_opt), pointer  :: p_k_value
    type(c_ptr)                 :: c_k_value

    c_name = trim(nomva0)//c_null_char

    p_k_value => vcopt
    c_k_value = c_loc(p_k_value)

    vcopt%iwarni = iwarnp
    vcopt%iconv  = iconvp
    vcopt%istat  = -1
    vcopt%ndircl = ndircp
    vcopt%idiff  = idiffp
    vcopt%idifft = -1
    vcopt%idften = idftnp
    vcopt%iswdyn = iswdyp
    vcopt%ischcv = ischcp
    vcopt%isstpc = isstpp
    vcopt%nswrgr = nswrgp
    vcopt%nswrsm = nswrsp
    vcopt%imrgra = imrgra
    vcopt%imligr = imligp
    vcopt%ircflu = ircflp
    vcopt%iwgrec = 0
    vcopt%thetav = thetap
    vcopt%blencv = blencp
    vcopt%blend_st = 0 ! Warning, may be overwritten if a field
    vcopt%epsilo = epsilp
    vcopt%epsrsm = epsrsp
    vcopt%epsrgr = epsrgp
    vcopt%climgr = climgp
    vcopt%extrag = 0
    vcopt%relaxv = relaxp

    call cs_equation_iterative_solve_vector(idtvar, iterns,                    &
                                            f_id, c_name,                      &
                                            ivisep, iescap, c_k_value,         &
                                            pvara, pvark,                      &
                                            coefav, coefbv, cofafv, cofbfv,    &
                                            i_massflux, b_massflux, i_viscm,   &
                                            b_viscm, i_visc, b_visc, secvif,   &
                                            secvib, viscce, weighf, weighb,    &
                                            icvflb, icvfli,                    &
                                            fimp, smbrp, pvar, eswork)
    return

  end subroutine coditv

  !=============================================================================
  !> \brief This function solves an advection diffusion equation with source
  !>        terms for one time step for the symmetric tensor variable
  !>        \f$ \tens{\variat} \f$.

  !> The equation reads:
  !>
  !> \f[
  !> \tens{f_s}^{imp}(\tens{\variat}^{n+1}-\tens{\variat}^n)
  !> + \divt \left( \tens{\variat}^{n+1} \otimes \rho \vect {u}
  !>              - \mu \gradtt \tens{\variat}^{n+1}\right)
  !> = \tens{Rhs}
  !> \f]
  !>
  !> This equation is rewritten as:
  !>
  !> \f[
  !> \tens{f_s}^{imp} \delta \tens{\variat}
  !> + \divt \left( \delta \tens{\variat} \otimes \rho \vect{u}
  !>              - \mu \gradtt \delta \tens{\variat} \right)
  !> = \tens{Rhs}^1
  !> \f]
  !>
  !> where \f$ \delta \tens{\variat} = \tens{\variat}^{n+1} -\tens{\variat}^n\f$
  !> and \f$ \tens{Rhs}^1 = \tens{Rhs}
  !> - \divt \left( \tens{\variat}^n \otimes \rho \vect{u}
  !>              - \mu \gradtt \tens{\variat}^n \right)\f$
  !>
  !> It is in fact solved with the following iterative process:
  !>
  !> \f[
  !> \tens{f_s}^{imp} \delta \tens{\variat}^k
  !> + \divt \left( \delta \tens{\variat}^k \otimes \rho \vect{u}
  !>              - \mu \gradtt \delta \tens{\variat}^k \right)
  !> = \tens{Rhs}^k
  !> \f]
  !>
  !> where \f$ \tens{Rhs}^k = \tens{Rhs}
  !> - \tens{f_s}^{imp} \left(\tens{\variat}^k-\tens{\variat}^n \right)
  !> - \divt \left( \tens{\variat}^k \otimes \rho \vect{u}
  !>              - \mu \gradtt \tens{\variat}^k \right)\f$
  !>
  !> Be careful, it is forbidden to modify \f$ \tens{f_s}^{imp} \f$ here!

  !> \param[in]     idtvar        indicator of the temporal scheme
  !> \param[in]     f_id          field id (or -1)
  !> \param[in]     iconvp        indicator
  !>                               - 1 convection,
  !>                               - 0 otherwise
  !> \param[in]     idiffp        indicator
  !>                               - 1 diffusion,
  !>                               - 0 otherwise
  !> \param[in]     ndircp        indicator (0 if the diagonal is stepped aside)
  !> \param[in]     imrgra        indicateur
  !>                               - 0 iterative gradient
  !>                               - 1 least squares gradient
  !> \param[in]     nswrsp        number of reconstruction sweeps for the
  !>                               Right Hand Side
  !> \param[in]     nswrgp        number of reconstruction sweeps for the
  !>                               gradients
  !> \param[in]     imligp        clipping gradient method
  !>                               - < 0 no clipping
  !>                               - = 0 thank to neighboring gradients
  !>                               - = 1 thank to the mean gradient
  !> \param[in]     ircflp        indicator
  !>                               - 1 flux reconstruction,
  !>                               - 0 otherwise
  !> \param[in]     ischcp        indicator
  !>                               - 1 centered
  !>                               - 0 2nd order
  !> \param[in]     isstpp        indicator
  !>                               - 1 without slope test
  !>                               - 0 with slope test
  !> \param[in]     idftnp        indicator
  !>                               - 1 the diffusivity is scalar
  !>                               - 6 the diffusivity is a symmetric tensor
  !> \param[in]     iswdyp        indicator
  !>                               - 0 no dynamic relaxation
  !>                               - 1 dynamic relaxation depending on
  !>                                 \f$ \delta \vect{\varia}^k \f$
  !>                               - 2 dynamic relaxation depending on
  !>                                 \f$ \delta \vect{\varia}^k \f$  and
  !>                                 \f$ \delta \vect{\varia}^{k-1} \f$
  !> \param[in]     iwarnp        verbosity
  !> \param[in]     blencp        fraction of upwinding
  !> \param[in]     epsilp        precision pour resol iter
  !> \param[in]     epsrsp        relative precision for the iterative process
  !> \param[in]     epsrgp        relative precision for the gradient
  !>                               reconstruction
  !> \param[in]     climgp        clipping coefficient for the computation of
  !>                               the gradient
  !> \param[in]     relaxp        coefficient of relaxation
  !> \param[in]     thetap        weighting coefficient for the theta-schema,
  !>                               - thetap = 0: explicit scheme
  !>                               - thetap = 0.5: time-centered
  !>                               scheme (mix between Crank-Nicolson and
  !>                               Adams-Bashforth)
  !>                               - thetap = 1: implicit scheme
  !> \param[in]     pvara         variable at the previous time step
  !>                               \f$ \vect{a}^n \f$
  !> \param[in]     pvark         variable at the previous sub-iteration
  !>                               \f$ \vect{a}^k \f$.
  !>                               If you sub-iter on Navier-Stokes, then
  !>                               it allows to initialize by something else
  !>                               than pvara (usually pvar=pvara)
  !> \param[in]     coefats        boundary condition array for the variable
  !>                               (Explicit part)
  !> \param[in]     coefbts        boundary condition array for the variable
  !>                               (Impplicit part)
  !> \param[in]     cofafts        boundary condition array for the diffusion
  !>                               of the variable (Explicit part)
  !> \param[in]     cofbfts        boundary condition array for the diffusion
  !>                               of the variable (Implicit part)
  !> \param[in]     i_massflux    mass flux at interior faces
  !> \param[in]     b_massflux    mass flux at boundary faces
  !> \param[in]     i_viscm       \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
  !>                               at interior faces for the matrix
  !> \param[in]     b_viscm       \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
  !>                               at boundary faces for the matrix
  !> \param[in]     i_visc        \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
  !>                               at interior faces for the r.h.s.
  !> \param[in]     b_visc        \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
  !>                               at boundary faces for the r.h.s.
  !> \param[in]     viscce        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
  !> \param[in]     weighf        internal face weight between cells i j in case
  !>                               of tensor diffusion
  !> \param[in]     weighb        boundary face weight for cells i in case
  !>                               of tensor diffusion
  !> \param[in]     icvflb        global indicator of boundary convection flux
  !>                               - 0 upwind scheme at all boundary faces
  !>                               - 1 imposed flux at some boundary faces
  !> \param[in]     icvfli        boundary face indicator array of convection flux
  !>                               - 0 upwind scheme
  !>                               - 1 imposed flux
  !> \param[in]     fimp          \f$ \tens{f_s}^{imp} \f$
  !> \param[in]     smbrp         Right hand side \f$ \vect{Rhs}^k \f$
  !> \param[in,out] pvar          current variable

  subroutine coditts (idtvar, f_id  , iconvp, idiffp, ndircp, imrgra, nswrsp,  &
                      nswrgp, imligp, ircflp, ischcp, isstpp, idftnp, iswdyp,  &
                      iwarnp, blencp, epsilp, epsrsp, epsrgp, climgp, relaxp,  &
                      thetap, pvara , pvark , coefats , coefbts , cofafts ,    &
                      cofbfts , i_massflux, b_massflux, i_viscm, b_viscm,      &
                      i_visc, b_visc, viscce, weighf , weighb , icvflb,        &
                      icvfli , fimp, smbrp, pvar)

    use, intrinsic :: iso_c_binding
    use entsor, only:nomva0
    use numvar
    use cplsat

    implicit none

    ! Arguments

    integer, intent(in) :: idtvar, f_id, iconvp, idiffp, ndircp, imrgra
    integer, intent(in) :: nswrsp, nswrgp, imligp, ircflp, ischcp, isstpp
    integer, intent(in) :: idftnp, iswdyp, iwarnp
    double precision, intent(in) :: blencp, epsilp, epsrsp, epsrgp, climgp
    double precision, intent(in) :: relaxp, thetap
    real(kind=c_double), dimension(*), intent(in) :: pvara, pvark, coefats
    real(kind=c_double), dimension(*), intent(in) :: coefbts, cofafts, cofbfts
    real(kind=c_double), dimension(*), intent(in) :: i_massflux, b_massflux
    real(kind=c_double), dimension(*), intent(in) :: i_visc, b_visc
    real(kind=c_double), dimension(*), intent(in) :: i_viscm, b_viscm
    real(kind=c_double), dimension(*), intent(in) :: viscce
    real(kind=c_double), dimension(*), intent(in) :: weighf, weighb
    integer, intent(in) :: icvflb
    integer(c_int), dimension(*), intent(in) :: icvfli
    real(kind=c_double), dimension(*), intent(in) :: fimp
    real(kind=c_double), dimension(*), intent(inout) :: smbrp, pvar

    ! Local variables
    character(len=len_trim(nomva0)+1, kind=c_char) :: c_name
    type(var_cal_opt), target   :: vcopt
    type(var_cal_opt), pointer  :: p_k_value
    type(c_ptr)                 :: c_k_value

    c_name = trim(nomva0)//c_null_char

    p_k_value => vcopt
    c_k_value = c_loc(p_k_value)

    vcopt%iwarni = iwarnp
    vcopt%iconv  = iconvp
    vcopt%istat  = -1
    vcopt%ndircl = ndircp
    vcopt%idiff  = idiffp
    vcopt%idifft = -1
    vcopt%idften = idftnp
    vcopt%iswdyn = iswdyp
    vcopt%ischcv = ischcp
    vcopt%isstpc = isstpp
    vcopt%nswrgr = nswrgp
    vcopt%nswrsm = nswrsp
    vcopt%imrgra = imrgra
    vcopt%imligr = imligp
    vcopt%ircflu = ircflp
    vcopt%iwgrec = 0
    vcopt%thetav = thetap
    vcopt%blencv = blencp
    vcopt%blend_st = 0 ! Warning, may be overwritten if a field
    vcopt%epsilo = epsilp
    vcopt%epsrsm = epsrsp
    vcopt%epsrgr = epsrgp
    vcopt%climgr = climgp
    vcopt%extrag = 0
    vcopt%relaxv = relaxp

    call cs_equation_iterative_solve_tensor(idtvar, f_id, c_name,              &
                                            c_k_value,                         &
                                            pvara, pvark,                      &
                                            coefats, coefbts, cofafts, cofbfts,&
                                            i_massflux, b_massflux, i_viscm,   &
                                            b_viscm, i_visc, b_visc, viscce,   &
                                            weighf, weighb , icvflb, icvfli,   &
                                            fimp, smbrp, pvar)
    return

  end subroutine coditts

  !=============================================================================
  !> \brief Wrapper to the function which adds the explicit part of the
  !>        convection/diffusion terms of a transport equation of
  !>        a scalar field \f$ \varia \f$.

  !> More precisely, the right hand side \f$ Rhs \f$ is updated as
  !> follows:
  !> \f[
  !> Rhs = Rhs - \sum_{\fij \in \Facei{\celli}}      \left(
  !>        \dot{m}_\ij \left( \varia_\fij - \varia_\celli \right)
  !>      - \mu_\fij \gradv_\fij \varia \cdot \vect{S}_\ij  \right)
  !> \f]
  !>
  !> Warning:
  !> - \f$ Rhs \f$ has already been initialized before calling bilsca!
  !> - mind the minus sign
  !>
  !> Options for the diffusive scheme:
  !> - idftnp = 1: scalar diffusivity
  !> - idftnp = 6: symmetric tensor diffusivity
  !>
  !> Options for the convective scheme:
  !> - blencp = 0: upwind scheme for the advection
  !> - blencp = 1: no upwind scheme except in the slope test
  !> - ischcp = 0: second order
  !> - ischcp = 1: centered
  !> - imucpp = 0: do not multiply the convective part by \f$ C_p \f$
  !> - imucpp = 1: multiply the convective part by \f$ C_p \f$

  !> \param[in]     idtvar        indicator of the temporal scheme
  !> \param[in]     f_id          field id (or -1)
  !> \param[in]     iconvp        indicator
  !>                               - 1 convection,
  !>                               - 0 otherwise
  !> \param[in]     idiffp        indicator
  !>                               - 1 diffusion,
  !>                               - 0 otherwise
  !> \param[in]     nswrgp        number of reconstruction sweeps for the
  !>                               gradients
  !> \param[in]     imligp        clipping gradient method
  !>                               - < 0 no clipping
  !>                               - = 0 by neighboring gradients
  !>                               - = 1 by the mean gradient
  !> \param[in]     ircflp        indicator
  !>                               - 1 flux reconstruction,
  !>                               - 0 otherwise
  !> \param[in]     ischcp        indicator
  !>                               - 1 centered
  !>                               - 0 2nd order
  !> \param[in]     isstpp        indicator
  !>                               - 1 without slope test
  !>                               - 0 with slope test
  !> \param[in]     inc           indicator
  !>                               - 0 when solving an increment
  !>                               - 1 otherwise
  !> \param[in]     imrgra        indicator
  !>                               - 0 iterative gradient
  !>                               - 1 least squares gradient
  !> \param[in]     iccocg        indicator
  !>                               - 1 re-compute cocg matrix
  !>                                 (for iterative gradients)
  !>                               - 0 otherwise
  !> \param[in]     iwarnp        verbosity
  !> \param[in]     imucpp        indicator
  !>                               - 0 do not multiply the convective term by Cp
  !>                               - 1 do multiply the convective term by Cp
  !> \param[in]     idftnp        indicator
  !>                               - 1 scalar diffusivity
  !>                               - 6 symmetric tensor diffusivity
  !> \param[in]     imasac        take mass accumulation into account?
  !> \param[in]     blencp        fraction of upwinding
  !> \param[in]     epsrgp        relative precision for the gradient
  !>                               reconstruction
  !> \param[in]     climgp        clipping coefficient for the computation of
  !>                               the gradient
  !> \param[in]     extrap        coefficient for extrapolation of the gradient
  !> \param[in]     relaxp        coefficient of relaxation
  !> \param[in]     thetap        weighting coefficient for the theta-schema,
  !>                               - thetap = 0: explicit scheme
  !>                               - thetap = 0.5: time-centered
  !>                               scheme (mix between Crank-Nicolson and
  !>                               Adams-Bashforth)
  !>                               - thetap = 1: implicit scheme
  !> \param[in]     pvar          solved variable (current time step)
  !> \param[in]     pvara         solved variable (previous time step)
  !> \param[in]     coefap        boundary condition array for the variable
  !>                               (explicit part)
  !> \param[in]     coefbp        boundary condition array for the variable
  !>                               (implicit part)
  !> \param[in]     cofafp        boundary condition array for the diffusion
  !>                               of the variable (explicit part)
  !> \param[in]     cofbfp        boundary condition array for the diffusion
  !>                               of the variable (implicit part)
  !> \param[in]     flumas        mass flux at interior faces
  !> \param[in]     flumab        mass flux at boundary faces
  !> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
  !>                               at interior faces for the r.h.s.
  !> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
  !>                               at boundary faces for the r.h.s.
  !> \param[in]     viscce        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
  !> \param[in]     xcpp          array of specific heat (Cp)
  !> \param[in]     weighf        internal face weight between cells i j in case
  !>                               of tensor diffusion
  !> \param[in]     weighb        boundary face weight for cells i in case
  !>                               of tensor diffusion
  !> \param[in]     icvflb        global indicator of boundary convection flux
  !>                               - 0 upwind scheme at all boundary faces
  !>                               - 1 imposed flux at some boundary faces
  !> \param[in]     icvfli        boundary face indicator array of convection flux
  !>                               - 0 upwind scheme
  !>                               - 1 imposed flux
  !> \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$

  subroutine bilsca (idtvar, f_id, iconvp, idiffp, nswrgp, imligp, ircflp,    &
                     ischcp, isstpp, inc, imrgra, iccocg, iwarnp, imucpp,     &
                     idftnp, imasac, blencp, epsrgp, climgp, extrap, relaxp,  &
                     thetap, pvar, pvara, coefap, coefbp, cofafp, cofbfp,     &
                     flumas, flumab, viscf, viscb, viscce, xcpp, weighf,      &
                     weighb, icvflb, icvfli, smbrp)

    use, intrinsic :: iso_c_binding
    use numvar
    use cplsat

    implicit none

    ! Arguments

    integer, intent(in) :: idtvar, f_id, iconvp, idiffp, imrgra, imucpp
    integer, intent(in) :: imligp, ircflp, ischcp, isstpp, inc, iccocg
    integer, intent(in) :: idftnp, iwarnp, imasac, nswrgp
    double precision, intent(in) :: blencp, epsrgp, climgp
    double precision, intent(in) :: relaxp, thetap, extrap
    real(kind=c_double), dimension(*), intent(in) :: pvar, pvara, coefap
    real(kind=c_double), dimension(*), intent(in) :: coefbp, cofafp, cofbfp
    real(kind=c_double), dimension(*), intent(in) :: flumas, flumab
    real(kind=c_double), dimension(*), intent(in) :: viscf, viscb
    real(kind=c_double), dimension(*), intent(in) :: viscce, xcpp
    real(kind=c_double), dimension(*), intent(in) :: weighf, weighb
    integer, intent(in) :: icvflb
    integer(c_int), dimension(*), intent(in) :: icvfli
    real(kind=c_double), dimension(*), intent(inout) :: smbrp

    ! Local variables
    type(var_cal_opt), target   :: vcopt
    type(var_cal_opt), pointer  :: p_k_value
    type(c_ptr)                 :: c_k_value

    p_k_value => vcopt
    c_k_value = c_loc(p_k_value)

    vcopt%iwarni = iwarnp
    vcopt%iconv  = iconvp
    vcopt%istat  = -1
    vcopt%idiff  = idiffp
    vcopt%idifft = -1
    vcopt%idften = idftnp
    vcopt%iswdyn = -1
    vcopt%ischcv = ischcp
    vcopt%isstpc = isstpp
    vcopt%nswrgr = nswrgp
    vcopt%nswrsm = -1
    vcopt%imrgra = imrgra
    vcopt%imligr = imligp
    vcopt%ircflu = ircflp
    vcopt%iwgrec = 0
    vcopt%thetav = thetap
    vcopt%blencv = blencp
    vcopt%blend_st = 0 ! Warning, may be overwritten if a field
    vcopt%epsilo = -1
    vcopt%epsrsm = -1
    vcopt%epsrgr = epsrgp
    vcopt%climgr = climgp
    vcopt%extrag = extrap
    vcopt%relaxv = relaxp

    call cs_balance_scalar(idtvar, f_id , imucpp, imasac, inc, iccocg,        &
                           c_k_value, pvar , pvara , coefap, coefbp,          &
                           cofafp, cofbfp, flumas, flumab, viscf, viscb,      &
                           viscce, xcpp , weighf, weighb, icvflb, icvfli,     &
                           smbrp)

    return

  end subroutine bilsca

  !=============================================================================
  !> \brief Wrapper to the function which adds the explicit part of the
  !>        convection/diffusion terms of a transport equation of
  !>        a vector field \f$ \vect{\varia} \f$.

  !> More precisely, the right hand side \f$ \vect{Rhs} \f$ is updated as
  !> follows:
  !> \f[
  !> \vect{Rhs} = \vect{Rhs} - \sum_{\fij \in \Facei{\celli}}      \left(
  !>        \dot{m}_\ij \left( \vect{\varia}_\fij - \vect{\varia}_\celli \right)
  !>      - \mu_\fij \gradt_\fij \vect{\varia} \cdot \vect{S}_\ij  \right)
  !> \f]
  !>
  !> Remark:
  !> if ivisep = 1, then we also take \f$ \mu \transpose{\gradt\vect{\varia}}
  !> + \lambda \trace{\gradt\vect{\varia}} \f$, where \f$ \lambda \f$ is
  !> the secondary viscosity, i.e. usually \f$ -\frac{2}{3} \mu \f$.
  !>
  !> Warning:
  !> - \f$ \vect{Rhs} \f$ has already been initialized before calling bilscv!
  !> - mind the sign minus
  !>
  !> Options for the diffusive scheme:
  !> - idftnp = 1: scalar diffusivity
  !> - idftnp = 6: symmetric tensor diffusivity
  !>
  !> Options for the convective scheme:
  !> - blencp = 0: upwind scheme for the advection
  !> - blencp = 1: no upwind scheme except in the slope test
  !> - ischcp = 0: second order
  !> - ischcp = 1: centered

  !> \param[in]     idtvar        indicator of the temporal scheme
  !> \param[in]     f_id          field id (or -1)
  !> \param[in]     iconvp        indicator
  !>                               - 1 convection,
  !>                               - 0 otherwise
  !> \param[in]     idiffp        indicator
  !>                               - 1 diffusion,
  !>                               - 0 otherwise
  !> \param[in]     nswrgp        number of reconstruction sweeps for the
  !>                               gradients
  !> \param[in]     imligp        clipping gradient method
  !>                               - < 0 no clipping
  !>                               - = 0 by neighboring gradients
  !>                               - = 1 by the mean gradient
  !> \param[in]     ircflp        indicator
  !>                               - 1 flux reconstruction,
  !>                               - 0 otherwise
  !> \param[in]     ischcp        indicator
  !>                               - 1 centered
  !>                               - 0 2nd order
  !> \param[in]     isstpp        indicator
  !>                               - 1 without slope test
  !>                               - 0 with slope test
  !> \param[in]     inc           indicator
  !>                               - 0 when solving an increment
  !>                               - 1 otherwise
  !> \param[in]     imrgra        indicator
  !>                               - 0 iterative gradient
  !>                               - 1 least squares gradient
  !> \param[in]     ivisep        indicator to take \f$ \divv
  !>                               \left(\mu \gradt \transpose{\vect{a}} \right)
  !>                               -2/3 \grad\left( \mu \dive \vect{a} \right)\f$
  !>                               - 1 take into account,
  !>                               - 0 otherwise
  !> \param[in]     iwarnp        verbosity
  !> \param[in]     idftnp        indicator
  !>                               - 1 scalar diffusivity
  !>                               - 6 symmetric tensor diffusivity
  !> \param[in]     imasac        take mass accumulation into account?
  !> \param[in]     blencp        fraction of upwinding
  !> \param[in]     epsrgp        relative precision for the gradient
  !>                               reconstruction
  !> \param[in]     climgp        clipping coefficient for the computation of
  !>                               the gradient
  !> \param[in]     relaxp        coefficient of relaxation
  !> \param[in]     thetap        weighting coefficient for the theta-schema,
  !>                               - thetap = 0: explicit scheme
  !>                               - thetap = 0.5: time-centered
  !>                               scheme (mix between Crank-Nicolson and
  !>                               Adams-Bashforth)
  !>                               - thetap = 1: implicit scheme
  !> \param[in]     pvar          solved velocity (current time step)
  !> \param[in]     pvara         solved velocity (previous time step)
  !> \param[in]     coefav        boundary condition array for the variable
  !>                               (explicit part)
  !> \param[in]     coefbv        boundary condition array for the variable
  !>                               (implicit part)
  !> \param[in]     cofafv        boundary condition array for the diffusion
  !>                               of the variable (explicit part)
  !> \param[in]     cofbfv        boundary condition array for the diffusion
  !>                               of the variable (implicit part)
  !> \param[in]     flumas        mass flux at interior faces
  !> \param[in]     flumab        mass flux at boundary faces
  !> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
  !>                               at interior faces for the r.h.s.
  !> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
  !>                               at boundary faces for the r.h.s.
  !> \param[in]     secvif        secondary viscosity at interior faces
  !> \param[in]     secvib        secondary viscosity at boundary faces
  !> \param[in]     viscce        symmetric cell tensor \f$ \tens{\mu}_\celli \f$
  !> \param[in]     weighf        internal face weight between cells i j in case
  !>                               of tensor diffusion
  !> \param[in]     weighb        boundary face weight for cells i in case
  !>                               of tensor diffusion
  !> \param[in]     icvflb        global indicator of boundary convection flux
  !>                               - 0 upwind scheme at all boundary faces
  !>                               - 1 imposed flux at some boundary faces
  !> \param[in]     icvfli        boundary face indicator of convection flux
  !>                               - 0 upwind scheme
  !>                               - 1 imposed flux
  !> \param[in,out] smbrp         right hand side \f$ \vect{Rhs} \f$

  subroutine bilscv (idtvar, f_id, iconvp, idiffp, nswrgp, imligp, ircflp,    &
                     ischcp, isstpp, inc, imrgra, ivisep, iwarnp, idftnp,     &
                     imasac, blencp, epsrgp, climgp, relaxp, thetap, pvar,    &
                     pvara, coefav, coefbv, cofafv, cofbfv, flumas, flumab,   &
                     viscf, viscb, secvif, secvib, viscce, weighf, weighb,    &
                     icvflb, icvfli, smbrp)

    use, intrinsic :: iso_c_binding
    use numvar
    use cplsat

    implicit none

    ! Arguments

    integer, intent(in) :: idtvar, f_id, iconvp, idiffp, imrgra
    integer, intent(in) :: imligp, ircflp, ischcp, isstpp, inc, ivisep
    integer, intent(in) :: idftnp, iwarnp, imasac, nswrgp
    double precision, intent(in) :: blencp, epsrgp, climgp
    double precision, intent(in) :: relaxp, thetap
    real(kind=c_double), dimension(*), intent(in) :: pvar, pvara, coefav
    real(kind=c_double), dimension(*), intent(in) :: coefbv, cofafv, cofbfv
    real(kind=c_double), dimension(*), intent(in) :: flumas, flumab
    real(kind=c_double), dimension(*), intent(in) :: viscf, viscb
    real(kind=c_double), dimension(*), intent(in) :: viscce, weighf, weighb
    real(kind=c_double), dimension(*), intent(in) :: secvif, secvib
    integer, intent(in) :: icvflb
    integer(c_int), dimension(*), intent(in) :: icvfli
    real(kind=c_double), dimension(*), intent(inout) :: smbrp

    ! Local variables
    type(var_cal_opt), target   :: vcopt
    type(var_cal_opt), pointer  :: p_k_value
    type(c_ptr)                 :: c_k_value

    p_k_value => vcopt
    c_k_value = c_loc(p_k_value)

    vcopt%iwarni = iwarnp
    vcopt%iconv  = iconvp
    vcopt%istat  = -1
    vcopt%idiff  = idiffp
    vcopt%idifft = -1
    vcopt%idften = idftnp
    vcopt%iswdyn = -1
    vcopt%ischcv = ischcp
    vcopt%isstpc = isstpp
    vcopt%nswrgr = nswrgp
    vcopt%nswrsm = -1
    vcopt%imrgra = imrgra
    vcopt%imligr = imligp
    vcopt%ircflu = ircflp
    vcopt%iwgrec = 0
    vcopt%thetav = thetap
    vcopt%blencv = blencp
    vcopt%blend_st = 0 ! Warning, may be overwritten if a field
    vcopt%epsilo = -1
    vcopt%epsrsm = -1
    vcopt%epsrgr = epsrgp
    vcopt%climgr = climgp
    vcopt%extrag = -1
    vcopt%relaxv = relaxp

    call cs_balance_vector(idtvar, f_id, imasac, inc, ivisep,                &
                           c_k_value, pvar, pvara , coefav, coefbv, cofafv,  &
                           cofbfv, flumas, flumab, viscf, viscb, secvif,     &
                           secvib, viscce, weighf, weighb,                   &
                           icvflb, icvfli, smbrp)

    return

  end subroutine bilscv

  !=============================================================================

  !> \brief Return pointer to coupling face indicator for a field

  !> \param[in]     f_id      id of given field
  !> \param[out]    cpl_faces pointer to coupling face indicator

  subroutine field_get_coupled_faces(f_id, cpl_faces)

    use, intrinsic :: iso_c_binding
    use mesh, only:nfabor

    implicit none

    ! Arguments

    integer, intent(in) :: f_id
    logical(kind=c_bool), dimension(:), pointer, intent(inout) :: cpl_faces

    ! Local variables
    type(c_ptr) :: c_p

    call cs_f_ic_field_coupled_faces(f_id, c_p)
    call c_f_pointer(c_p, cpl_faces, [nfabor])

    return

  end subroutine field_get_coupled_faces

  !=============================================================================

  !> \brief Return notebook parameter value

  !> \param[in]     name      name of the notebook parameter
  !> \result        val

  function notebook_parameter_value_by_name(name) result(val)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name
    double precision :: val

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    real(kind=c_double) :: c_val

    c_name = trim(name)//c_null_char

    c_val = cs_f_notebook_parameter_value_by_name(c_name)
    val = c_val

  end function notebook_parameter_value_by_name

  !=============================================================================

  end module cs_c_bindings
