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

  integer :: POST_ON_LOCATION, POST_BOUNDARY_NR, POST_MONITOR

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

  parameter (POST_ON_LOCATION=1)
  parameter (POST_BOUNDARY_NR=2)
  parameter (POST_MONITOR=4)

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

  procedure() :: csexit
  procedure() :: cslogname

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
    integer(c_int) :: imvisf
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
    real(c_double) :: relaxv
  end type var_cal_opt

  !=============================================================================

  interface

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

    ! Interface to C function copying a var_cal_opt structure associated
    ! with a field.

    subroutine cs_f_field_set_key_struct_var_cal_opt(f_id, k_value) &
      bind(C, name='cs_f_field_set_key_struct_var_cal_opt')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: f_id
      type(c_ptr), value                :: k_value
    end subroutine cs_f_field_set_key_struct_var_cal_opt

    !---------------------------------------------------------------------------

    ! Interface to C function setting a var_cal_opt structure associated
    ! with a field.

    subroutine cs_f_field_get_key_struct_var_cal_opt(f_id, k_value) &
      bind(C, name='cs_f_field_get_key_struct_var_cal_opt')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value             :: f_id
      type(c_ptr), value                :: k_value
    end subroutine cs_f_field_get_key_struct_var_cal_opt

    !---------------------------------------------------------------------------

    ! Interface to C function returninng a pointer to a cs_equation_param_t
    ! structure based on a given var_cal_opt structure.

    function equation_param_from_vcopt(k_value) result(eqp) &
      bind(C, name='cs_f_equation_param_from_var_cal_opt')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), value :: k_value
      type(c_ptr)        :: eqp
    end function equation_param_from_vcopt

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

    !> \brief Convert enthalpy to temperature for gas combustion.

    function cs_gas_combustion_h_to_t(xespec, enthal) result(temper)  &
      bind(C, name='cs_gas_combustion_h_to_t')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: xespec
      real(c_double), value :: enthal
      real(c_double) :: temper
    end function cs_gas_combustion_h_to_t

    !---------------------------------------------------------------------------

    !> \brief Convert temperature to enthalpy for gas combustion.

    function cs_gas_combustion_t_to_h(xespec, temper) result(enthal)  &
      bind(C, name='cs_gas_combustion_t_to_h')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: xespec
      real(c_double), value :: temper
      real(c_double) :: enthal
    end function cs_gas_combustion_t_to_h

    !---------------------------------------------------------------------------

    !> \brief Define automatic turbulence values for specific physical modules.

    !> The definitions are similar to those of the standard case, though wall
    !> shear direction is not computed for second-order models, and determination
    !> of face BC types is done using the legacy physical model zone info
    !> (izfpp, ...).

    !> \param[in]  bc_type  type of boundary for each face

    subroutine cs_boundary_conditions_legacy_turbulence(bc_type) &
      bind(C, name='cs_boundary_conditions_legacy_turbulence')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), dimension(*) :: bc_type
    end subroutine cs_boundary_conditions_legacy_turbulence

    !---------------------------------------------------------------------------

    ! Interface to C function activating default log.

    subroutine cs_log_default_activate(activate)  &
      bind(C, name='cs_log_default_activate')
      use, intrinsic :: iso_c_binding
      implicit none
      logical(kind=c_bool), value :: activate
    end subroutine cs_log_default_activate

    !---------------------------------------------------------------------------

    ! Interface to C function activating default log.

    function cs_log_default_is_active() result(is_active) &
      bind(C, name='cs_log_default_is_active')
      use, intrinsic :: iso_c_binding
      implicit none
      logical(kind=c_bool) :: is_active
    end function cs_log_default_is_active

    !---------------------------------------------------------------------------

    ! Initialize turbulence model structures

    subroutine cs_turb_model_init()  &
      bind(C, name='cs_turb_model_init')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_turb_model_init

    !---------------------------------------------------------------------------

    ! Set type and order of the turbulence model

    subroutine cs_set_type_order_turbulence_model()  &
      bind(C, name='cs_set_type_order_turbulence_model')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_set_type_order_turbulence_model

    !---------------------------------------------------------------------------

    ! Scalar clipping

    subroutine cs_f_scalar_clipping(f_id)  &
      bind(C, name='cs_f_scalar_clipping')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: f_id
    end subroutine cs_f_scalar_clipping

    !---------------------------------------------------------------------------

    ! Temporal and z-axis interpolation for meteorological profiles

    function cs_intprf(nprofz, nproft, profz, proft,              &
                       profv, xz, t) result (var)                 &
         bind(C, name='cs_intprf')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: nprofz, nproft
      real(kind=c_double), dimension(nprofz), intent(in) :: profz
      real(kind=c_double), dimension(nproft), intent(in) :: proft
      real(kind=c_double), dimension(nprofz,nproft), intent(in) :: profv
      real(kind=c_double), intent(in), value :: xz, t
      real(kind=c_double) :: var
    end function cs_intprf

    !---------------------------------------------------------------------------

    ! Z-axis interpolation for meteorological profiles

    subroutine cs_intprz(nprofz, profz, profv, xz, z_lv, var)  &
      bind(C, name='cs_intprz')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: nprofz
      real(kind=c_double), dimension(nprofz), intent(in) :: profz, profv
      real(kind=c_double), intent(in), value :: xz
      integer(c_int), dimension(2), intent(out) :: z_lv
      real(kind=c_double), intent(out) :: var
    end subroutine cs_intprz

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

    ! Interface to C function returning number of SYRTHES couplingsg.

    function cs_syr_coupling_n_couplings() result(n_couplings) &
      bind(C, name='cs_syr_coupling_n_couplings')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int) :: n_couplings
    end function cs_syr_coupling_n_couplings

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
                                         vel_dir, shear_dir,                   &
                                         rcodcl)                               &
      bind(C, name='cs_f_turbulence_bc_inlet_k_eps')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: face_num
      real(c_double), value :: k, eps
      real(kind=c_double), dimension(3) :: vel_dir
      real(kind=c_double), dimension(3) :: shear_dir
      real(kind=c_double), dimension(*) :: rcodcl
    end subroutine turbulence_bc_inlet_k_eps

    !---------------------------------------------------------------------------

    !> \brief Set inlet boundary condition values for turbulence variables based
    !>        on given k and epsilon values only if not initialized already.
    !>
    !> \param[in]     face_id       boundary face id
    !> \param[in]     k             turbulent kinetic energy
    !> \param[in]     epsilon       turbulent dissipation
    !> \param[in]     vel_dir       velocity direction
    !> \param[in]     shear_dir     shear direction
    !> \param[out]    rcodcl        boundary condition values

    subroutine turbulence_bc_set_uninit_inlet_k_eps(face_num,                  &
                                                    k, eps,                    &
                                                    vel_dir, shear_dir,        &
                                                    rcodcl)                    &
      bind(C, name='cs_f_turbulence_bc_set_uninit_inlet_k_eps')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: face_num
      real(c_double), value :: k, eps
      real(kind=c_double), dimension(3) :: vel_dir
      real(kind=c_double), dimension(3) :: shear_dir
      real(kind=c_double), dimension(*) :: rcodcl
    end subroutine turbulence_bc_set_uninit_inlet_k_eps

    !---------------------------------------------------------------------------

    !> \brief Compute molar and mass fractions of elementary species Ye, Xe
    !>  (fuel, O2, CO2, H2O, N2) from global species Yg (fuel, oxidant, products)

    !> \param[in]     yg            global mass fraction
    !> \param[out]    ye            elementary mass fraction
    !> \param[out]    xe            elementary molar fraction

    subroutine yg2xye(yg, ye, xe)  &
      bind(C, name='cs_combustion_gas_yg2xye')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: yg
      real(kind=c_double), dimension(*), intent(out) :: ye, xe
    end subroutine yg2xye

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

    ! Interface to C function to get the bc type array pointer

    subroutine cs_f_boundary_conditions_get_pointers(itypfb, izfppp) &
      bind(C, name='cs_f_boundary_conditions_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: itypfb, izfppp
    end subroutine cs_f_boundary_conditions_get_pointers

    !---------------------------------------------------------------------------

    ! Interface to C function creating a directory

    subroutine cs_file_mkdir_default(path)  &
      bind(C, name='cs_file_mkdir_default')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: path
    end subroutine cs_file_mkdir_default

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

    ! Interface to C function initializing condensation-related field key.

    function cs_gas_mix_get_field_key()  &
      result(k_id) &
      bind(C, name='cs_gas_mix_get_field_key')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: k_id
    end function cs_gas_mix_get_field_key

    !---------------------------------------------------------------------------

    ! Interface to C function
    !> \brief Get reference value of a physical property

    !> \param[in] name  property name
    !> \return reference value (c_double)

    function cs_physical_property_get_ref_value(name) result(val) &
      bind(C, name='cs_physical_property_get_ref_value')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      real(c_double)                                          :: val
    end function cs_physical_property_get_ref_value

    !---------------------------------------------------------------------------

    ! Interface to C function
    !> \brief Set reference value for a physical property

    !> \param[in] name  property name
    !> \param[in] val   new value to set

    subroutine cs_physical_property_set_ref_value(name, val) &
      bind(C, name='cs_physical_property_set_ref_value')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      real(kind=c_double), value, intent(in)                  :: val
    end subroutine cs_physical_property_set_ref_value

    !---------------------------------------------------------------------------

    ! Interface to C function
    !> \brief Create physical property

    !> \param[in] name    property name
    !> \param[in] dim     property dimension
    !> \param[in] refval  reference value

    subroutine cs_physical_property_create(name, dim, refval) &
      bind(C, name='cs_physical_property_create')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      integer(c_int), value, intent(in)                       :: dim
      real(kind=c_double), value, intent(in)                  :: refval
    end subroutine cs_physical_property_create

    !---------------------------------------------------------------------------

    ! Interface to C function
    !> \brief Get property reference values for a given zone

    !> \param[in] name    property name
    !> \param[in] zname   zone name
    !> \param[in] retval  array of values to return

    subroutine cs_physical_property_get_zone_values(name, zname, retval) &
      bind(C, name='cs_physical_property_get_zone_values')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      character(kind=c_char, len=1), dimension(*), intent(in) :: zname
      real(kind=c_double), dimension(*), intent(out)          :: retval
    end subroutine cs_physical_property_get_zone_values

    !---------------------------------------------------------------------------

    ! Interface to C function
    !> \brief Update reference values for a property on a given zone

    !> \param[in] name   property name
    !> \param[in] zname  zone name
    !> \param[in] vals   array of values to set

    subroutine cs_physical_property_update_zone_values(name, zname, vals) &
      bind(C, name='cs_physical_property_update_zone_values')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      character(kind=c_char, len=1), dimension(*), intent(in) :: zname
      real(kind=c_double), dimension(*), intent(in)           :: vals
    end subroutine cs_physical_property_update_zone_values

    !---------------------------------------------------------------------------

    ! Interface to C function
    !> \brief Add a property definition on a given zone using a single value

    !> \param[in] name   property name
    !> \param[in] zname  zone name
    !> \param[in] dim    property dimension
    !> \param[in] val    reference value for the zone

    subroutine cs_physical_property_define_from_value(name, zname, dim, val) &
      bind(C, name='cs_physical_property_define_from_value')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      character(kind=c_char, len=1), dimension(*), intent(in) :: zname
      integer(c_int), value                                   :: dim
      real(kind=c_double), value, intent(in)                  :: val
    end subroutine cs_physical_property_define_from_value

    !---------------------------------------------------------------------------

    ! Interface to C function
    !> \brief Add a property multi-diemnsional definition on a given zone

    !> \param[in] name   property name
    !> \param[in] zname  zone name
    !> \param[in] dim    property dimension (>1)
    !> \param[in] vals   array of values to set

    subroutine cs_physical_property_define_from_values(name, zname, dim, vals) &
      bind(C, name='cs_physical_property_define_from_values')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      character(kind=c_char, len=1), dimension(*), intent(in) :: zname
      integer(c_int), value                                   :: dim
      real(kind=c_double), dimension(*), intent(in)           :: vals
    end subroutine cs_physical_property_define_from_values

    !---------------------------------------------------------------------------

    ! Interface to C function
    !> \brief Add a property definition based on a cs_field_t. Field is created if needed

    !> \param[in] name          property name
    !> \param[in] type_flag     field type flag
    !> \param[in] location_id   location id flag
    !> \param[in] dim           field dimension
    !> \param[in] has_previous  does the field has val_pre

    subroutine cs_physical_property_define_from_field(name, type_flag, &
      location_id, dim, has_previous) &
      bind(C, name='cs_physical_property_define_from_field')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      integer(c_int), value                                   :: type_flag
      integer(c_int), value                                   :: location_id
      integer(c_int), value                                   :: dim
      logical(c_bool), value                                  :: has_previous
    end subroutine cs_physical_property_define_from_field

    !---------------------------------------------------------------------------

    ! Interface to C function
    !> \brief Return id of field associated to property

    !> \param[in] name  property name
    !> \return field id (int)

    function cs_physical_property_field_id_by_name(name) &
      result(f_id) &
      bind(C, name='cs_physical_property_field_id_by_name')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
      integer(c_int)                                          :: f_id
    end function cs_physical_property_field_id_by_name

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

    ! Interface to C function which returns if the restart is from NEPTUNE_CFD
    function cs_restart_present() result(flag) &
      bind(C, name='cs_restart_present')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int) :: flag
    end function cs_restart_present

    !---------------------------------------------------------------------------

    ! Interface to C function to check field read status from checkpoint file

    function cs_restart_get_field_read_status(f_id) result(retval) &
      bind(C, name='cs_restart_get_field_read_status')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      integer(c_int)        :: retval
    end function cs_restart_get_field_read_status

    !---------------------------------------------------------------------------

    ! Interface to C function incrementing time step

    subroutine cs_time_step_update_dt(dt) &
      bind(C, name='cs_time_step_update_dt')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), value :: dt
    end subroutine cs_time_step_update_dt

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

    ! Interface to C function creating a CDO variable field

    function cs_variable_cdo_field_create(name, label, location_id,       &
                                          dim, has_previous) result(id)   &
      bind(C, name='cs_variable_cdo_field_create')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in)  :: name, label
      integer(c_int), value                                    :: location_id
      integer(c_int), value                                    :: has_previous
      integer(c_int), value                                    :: dim
      integer(c_int)                                           :: id
    end function cs_variable_cdo_field_create

    !---------------------------------------------------------------------------

    ! Add terms from backward differentiation in time.

    subroutine cs_backward_differentiation_in_time(field_id,                  &
                                                   exp_part, imp_part)        &
      bind(C, name='cs_f_backward_differentiation_in_time')
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

    ! Interface to C user function for boundary conditions

    subroutine user_boundary_conditions(bc_type)  &
      bind(C, name='cs_user_boundary_conditions_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), dimension(*), intent(inout) :: bc_type
    end subroutine user_boundary_conditions

    !---------------------------------------------------------------------------

    ! Interface to C user function for extra operations

    subroutine user_extra_operations_initialize()  &
      bind(C, name='cs_user_extra_operations_initialize_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine user_extra_operations_initialize

    !---------------------------------------------------------------------------

    ! Interface to C user function for initialization

    subroutine user_initialization()  &
      bind(C, name='cs_user_initialization_wrapper')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine user_initialization

    !---------------------------------------------------------------------------

    ! Interface to C function cs_clip_turbulent_fluxes

    subroutine cs_clip_turbulent_fluxes(flux_id, variance_id) &
      bind(C, name='cs_clip_turbulent_fluxes')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: flux_id
      integer(c_int), value :: variance_id
    end subroutine cs_clip_turbulent_fluxes

    !---------------------------------------------------------------------------

    ! Interface to C function cs_balance_scalar

    subroutine cs_balance_scalar(idtvar, f_id , imucpp, imasac, inc,          &
                                 vcopt , pvar , pvara,                        &
                                 coefap, coefbp, cofafp, cofbfp, i_massflux,  &
                                 b_massflux, i_visc, b_visc, viscel, xcpp,    &
                                 weighf, weighb, icvflb, icvfli,              &
                                 smbrp)                                       &
      bind(C, name='cs_balance_scalar')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: idtvar, f_id, imasac, imucpp, inc
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

    ! Interface to C soot production function
    subroutine cs_soot_production(f_id, smbrs, rovsdt) &
      bind(C, name='cs_soot_production')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
      real(kind=c_double), dimension(*) :: smbrs, rovsdt
    end subroutine cs_soot_production

    !---------------------------------------------------------------------------

    !> \brief  Add a species field to the gas mix (set of fields).

    subroutine gas_mix_add_species(f_id)  &
      bind(C, name='cs_gas_mix_add_species')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id
    end subroutine gas_mix_add_species

    !---------------------------------------------------------------------------

    !> \brief Get the aerosols concentrations and numbers from aerosol code

    subroutine cs_atmo_aerosol_get_aero(array) &
      bind(C, name='cs_atmo_aerosol_get_aero')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(out) :: array
    end subroutine cs_atmo_aerosol_get_aero

    !---------------------------------------------------------------------------

    !> \brief Get the gas concentrations from aerosol code

    subroutine cs_atmo_aerosol_get_gas(array) &
      bind(C, name='cs_atmo_aerosol_get_gas')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(out) :: array
    end subroutine cs_atmo_aerosol_get_gas

    !---------------------------------------------------------------------------

    !> \brief Return pointers to atmo chemistry arrays

    subroutine cs_f_atmo_chem_arrays_get_pointers(isca_chem, dmmk, &
        chempoint) &
      bind(C, name='cs_f_atmo_chem_arrays_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: isca_chem, dmmk, chempoint
    end subroutine cs_f_atmo_chem_arrays_get_pointers

    !---------------------------------------------------------------------------

    !> \brief Sets the chemistry concentration file name

    subroutine cs_atmo_set_chem_conc_file_name(name) &
      bind(C, name='cs_atmo_set_chem_conc_file_name')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
    end subroutine cs_atmo_set_chem_conc_file_name


    !> \brief Sets the file name used to initialize the aerosol shared library

    subroutine cs_atmo_chemistry_set_aerosol_file_name(name) &
      bind(C, name='cs_atmo_chemistry_set_aerosol_file_name')
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char, len=1), dimension(*), intent(in) :: name
    end subroutine cs_atmo_chemistry_set_aerosol_file_name

    !---------------------------------------------------------------------------

    ! Interface to C function for atmo

    subroutine raysze(xlat, xlong, jour, heurtu, imer, albe, za, muzero, omega, fo) &
      bind(C, name='cs_atmo_compute_solar_angles')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), value :: xlat, xlong, jour, heurtu
      integer(kind=c_int), value :: imer
      real(kind=c_double), intent(inout) :: albe, za, muzero, omega, fo
    end subroutine raysze

    !---------------------------------------------------------------------------

    !> \brief Initialize C chemistry structure from Fortran

    subroutine cs_f_atmo_chem_initialize_species_to_fid(species_fid) &
      bind(C, name='cs_f_atmo_chem_initialize_species_to_fid')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), dimension(*), intent(in) :: species_fid
    end subroutine cs_f_atmo_chem_initialize_species_to_fid

    !---------------------------------------------------------------------------

    ! Computes the explicit chemical source term for atmospheric chemistry
    ! in case of a semi-coupled resolution

     subroutine chem_source_terms(iscal, st_exp, st_imp) &
       bind(C, name='cs_atmo_chem_source_terms')
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int), value :: iscal
       real(kind=c_double), dimension(*), intent(inout) :: st_exp, st_imp
     end subroutine chem_source_terms

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

    subroutine cs_ctwr_phyvar_update(rho0, t0, p0)             &
      bind(C, name='cs_ctwr_phyvar_update')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), value :: rho0, t0, p0
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

    ! Interface to C function for atmospheric module

    subroutine cs_atmo_phyvar_update()             &
      bind(C, name='cs_atmo_phyvar_update')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_atmo_phyvar_update

    !---------------------------------------------------------------------------

    ! Interface to C function for head losses

    subroutine cs_head_losses_compute(ckupdc)  &
      bind(C, name='cs_head_losses_compute')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*) :: ckupdc
    end subroutine cs_head_losses_compute

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


    ! Interface to C function cs_math_3_normalize

    subroutine vector_normalize(vin, vout)                   &
      bind(C, name='cs_f_math_3_normalize')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(in) :: vin
      real(kind=c_double), dimension(*), intent(out) :: vout
    end subroutine vector_normalize

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

    ! Compute solid mesh quantities

    subroutine cs_f_mesh_quantities_solid_compute()   &
      bind(C, name='cs_f_mesh_quantities_solid_compute')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_f_mesh_quantities_solid_compute


    !---------------------------------------------------------------------------

    ! Init fluid mesh quantities

    subroutine cs_porous_model_init_fluid_quantities()   &
      bind(C, name='cs_porous_model_init_fluid_quantities')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_porous_model_init_fluid_quantities

    !---------------------------------------------------------------------------

    ! Initialize has_disable_flag

    subroutine cs_porous_model_init_disable_flag()   &
      bind(C, name='cs_porous_model_init_disable_flag')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine cs_porous_model_init_disable_flag

    !---------------------------------------------------------------------------

    ! Set porosity model.

    subroutine cs_porous_model_set_model(iporos)   &
      bind(C, name='cs_porous_model_set_model')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=c_int), value :: iporos
    end subroutine cs_porous_model_set_model

    !---------------------------------------------------------------------------

    !> \brief Return pointers

    !> \param[out]   ibm_porosity_mode  Pointer to ibm_porosity_mode

    subroutine cs_f_porosity_ibm_get_pointer(ibm_porosity_mode) &
      bind(C, name='cs_f_porosity_ibm_get_pointer')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ibm_porosity_mode
    end subroutine cs_f_porosity_ibm_get_pointer

   !---------------------------------------------------------------------------

    !> \brief Return pointers

    !> \param[out]   compute_porosity_from_scan  Pointer to
    !>                                           compute_porosity_from_scan

    subroutine cs_f_porosity_from_scan_get_pointer(compute_porosity_from_scan) &
      bind(C, name='cs_f_porosity_from_scan_get_pointer')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: compute_porosity_from_scan
    end subroutine cs_f_porosity_from_scan_get_pointer

    !---------------------------------------------------------------------------

    ! Interface to C function for scalar gradient with homogeneous Neumann BCs

    subroutine cs_f_gradient_hn_s(f_id, imrgra, inc, n_r_sweeps,               &
                                  iwarnp, imligp, epsrgp, climgp,              &
                                  pvar, grad)                                  &
      bind(C, name='cs_f_gradient_hn_s')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: f_id, imrgra, inc, n_r_sweeps
      integer(c_int), value :: iwarnp, imligp
      real(kind=c_double), value :: epsrgp, climgp
      real(kind=c_double), dimension(*), intent(inout) :: pvar
      real(kind=c_double), dimension(*), intent(inout) :: grad
    end subroutine cs_f_gradient_hn_s

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

    ! Interface to C function returning 1 for active cells

    function cs_f_porous_model_cell_is_active(cell_id) result(is_active) &
      bind(C, name='cs_f_porous_model_cell_is_active')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: cell_id
      integer(kind=c_int) :: is_active
    end function cs_f_porous_model_cell_is_active

    !---------------------------------------------------------------------------

    ! Interface to C function for temperature-enthalpy conversion at
    ! selected faces

    subroutine cs_ht_convert_t_to_h_faces_l(n_faces, face_ids, t, h)  &
      bind(C, name='cs_ht_convert_t_to_h_faces_l')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: n_faces
      integer(c_int), dimension(*), intent(in) :: face_ids
      real(kind=c_double), dimension(*), intent(in) :: t
      real(kind=c_double), dimension(*), intent(inout) :: h
    end subroutine cs_ht_convert_t_to_h_faces_l

    !---------------------------------------------------------------------------

    ! Interface to C function computing standard atmospheric profile

    subroutine atmstd(z, p, t, r) &
      bind(C, name='cs_atmo_profile_std')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), intent(in), value :: z
      real(kind=c_double), intent(out) :: p, t, r
    end subroutine atmstd

    !---------------------------------------------------------------------------

    ! Interface to C function computing etheta and eq variable
    ! knowing the saturation.

    subroutine atprke(tinstk, smbrk, smbre)  &
      bind(C, name='cs_atprke')
      use, intrinsic :: iso_c_binding
      implicit none
      real(kind=c_double), dimension(*), intent(inout) :: tinstk, smbrk, smbre
    end subroutine atprke

    !---------------------------------------------------------------------------

    ! Interface to C function solving the quadratic k-epsilon model.

    subroutine cs_turbulence_ke_q(phase_id, rij) &
      bind(C, name='cs_turbulence_ke_q')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: phase_id
      real(kind=c_double), dimension(6,*), intent(out) :: rij
    end subroutine cs_turbulence_ke_q

    !---------------------------------------------------------------------------
    !> \brief Clipping of the turbulent Reynods stress tensor and the turbulent
    !> dissipation (segregated version).
    !>
    !> \param[in]     ncel          number of cells
    !> \param[in]     iclip         indicator = 0 if viscl0 is used
    !>                              otherwise viscl is used.

    subroutine cs_turbulence_rij_clip_sg(ncel, iclip) &
      bind(C, name='cs_turbulence_rij_clip_sg')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: ncel, iclip
    end subroutine cs_turbulence_rij_clip_sg

    !---------------------------------------------------------------------------
    ! Interface to C function to compute the number of aerosols

    subroutine cs_atmo_aerosol_ssh_compute_numbers(dlconc0) &
       bind(C, name='cs_atmo_aerosol_ssh_compute_numbers')
       use, intrinsic :: iso_c_binding
       implicit none
       real(kind=c_double), dimension(*), intent(inout) :: dlconc0
    end subroutine cs_atmo_aerosol_ssh_compute_numbers

    !---------------------------------------------------------------------------
    ! Interface to C function to set the humidity in SSH

    subroutine cs_atmo_aerosol_ssh_set_t_p_h(t, p, h) &
       bind(C, name='cs_atmo_aerosol_ssh_set_t_p_h')
       use, intrinsic :: iso_c_binding
       implicit none
       real(kind=c_double), intent(inout) :: t, p, h
    end subroutine cs_atmo_aerosol_ssh_set_t_p_h

    !---------------------------------------------------------------------------
    ! Interface to C function for physical properties variable

    subroutine cs_physical_properties_update(iterns) &
      bind(C, name='cs_physical_properties_update')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: iterns
    end subroutine cs_physical_properties_update

    !---------------------------------------------------------------------------

    ! Interface to C function updating scalar array ghost values.

    subroutine synsca(var)  &
      bind(C, name='cs_mesh_sync_var_scal')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), dimension(*) :: var
    end subroutine synsca

    !---------------------------------------------------------------------------

    ! Interface to C function updating scalar array extended ghost values.

    subroutine synsce(var)  &
      bind(C, name='cs_mesh_sync_var_scal_ext')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), dimension(*) :: var
    end subroutine synsce

    !---------------------------------------------------------------------------

    ! Interface to C function updating vector array ghost values.

    subroutine synvin(var)  &
      bind(C, name='cs_mesh_sync_var_vect')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), dimension(3, *) :: var
    end subroutine synvin

    !---------------------------------------------------------------------------

    ! Interface to C function updating vector array extended ghost values.

    subroutine synvie(var)  &
      bind(C, name='cs_mesh_sync_var_vect_ext')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), dimension(3, *) :: var
    end subroutine synvie

    !---------------------------------------------------------------------------

    ! Interface to C function updating tensor array ghost values.

    subroutine syntin(var)  &
      bind(C, name='cs_mesh_sync_var_tens')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), dimension(3, 3, *) :: var
    end subroutine syntin

    !---------------------------------------------------------------------------

    ! Interface to C function updating symmetric tensor array ghost values.

    subroutine syntis(var)  &
      bind(C, name='cs_mesh_sync_var_sym_tens')
      use, intrinsic :: iso_c_binding
      implicit none
      real(c_double), dimension(6, *) :: var
    end subroutine syntis

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

  !> \brief Temporal and z-axis interpolation for meteorological profiles

  !> An optimized linear interpolation is used.

  subroutine intprf(nprofz, nproft, profz, proft,              &
                    profv, xz, temps, var)
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), intent(in), value :: nprofz, nproft
    real(kind=c_double), dimension(nprofz), intent(in) :: profz
    real(kind=c_double), dimension(nproft), intent(in) :: proft
    real(kind=c_double), dimension(nprofz, nproft), intent(in) :: profv
    real(kind=c_double), intent(in), value :: xz, temps
    real(kind=c_double), intent(out) :: var

    var = cs_intprf(nprofz, nproft, profz, proft, profv, xz, temps)

  end subroutine intprf

  !=============================================================================

  !> \brief z-axis interpolation for meteorological profiles

  !> An optimized linear interpolation is used.

  subroutine intprz(nprofz, profz, profv, xz, iz1, iz2, var)
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), intent(in), value :: nprofz
    real(kind=c_double), dimension(nprofz), intent(in) :: profz, profv
    real(kind=c_double), intent(in), value :: xz
    integer(c_int), intent(out) :: iz1, iz2
    real(kind=c_double), intent(out) :: var

    integer(c_int), dimension(2) :: z_lv

    call cs_intprz(nprofz, profz, profv, xz, z_lv, var)
    iz1 = z_lv(1) + 1
    iz2 = z_lv(2) + 1

  end subroutine intprz

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

    call cs_f_field_set_key_struct_var_cal_opt(c_f_id, c_k_value)

    return

  end subroutine field_set_key_struct_var_cal_opt

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
    type(var_cal_opt), intent(out), target :: k_value

    ! Local variables

    integer(c_int)                 :: c_f_id
    type(var_cal_opt),pointer      :: p_k_value
    type(c_ptr)                    :: c_k_value

    c_f_id = f_id

    p_k_value => k_value
    c_k_value = c_loc(p_k_value)

    call cs_f_field_get_key_struct_var_cal_opt(c_f_id, c_k_value)

    return

  end subroutine field_get_key_struct_var_cal_opt

  !=============================================================================

  !> \brief  Compute cell gradient for scalar with homegeneous gradient.

  !> \param[in]       f_id             field id, or -1
  !> \param[in]       imrgra           gradient computation mode
  !> \param[in]       inc              0: increment; 1: do not increment
  !> \param[in]       nswrgp           number of sweeps for reconstruction
  !> \param[in]       imligp           gradient limitation method:
  !>                                     < 0 no limitation
  !>                                     = 0 based on neighboring gradients
  !>                                     = 1 based on mean gradient
  !> \param[in]       iwarnp           verbosity
  !> \param[in]       epsrgp           relative precision for reconstruction
  !> \param[in]       climgp           limiter coefficient for imligp
  !> \param[in, out]  pvar             cell values whose gradient is computed
  !> \param[out]      grad             resulting gradient

  subroutine gradient_hn_s(f_id, imrgra, inc, nswrgp,                          &
                           imligp, iwarnp, epsrgp, climgp,                     &
                           pvar, grad)

    use, intrinsic :: iso_c_binding
    use paramx
    use mesh
    use field
    use period

    implicit none

    ! Arguments

    integer, intent(in) :: f_id, imrgra, inc, nswrgp
    integer, intent(in) :: imligp, iwarnp
    double precision, intent(in) :: epsrgp, climgp
    real(kind=c_double), dimension(ncelet), intent(inout) :: pvar
    real(kind=c_double), dimension(3, ncelet), intent(out) :: grad

    ! The gradient of a potential (pressure, ...) is a vector

    call cs_f_gradient_hn_s(f_id, imrgra, inc, nswrgp,                         &
                            iwarnp, imligp,                                    &
                            epsrgp, climgp, pvar, grad)

  end subroutine gradient_hn_s

  !=============================================================================

  !> \brief Get field checkpoint read status. Returns 1 if field was read, 0
  !         otherwise.

  !> \param[in]  f_id     field id
  !> \param[out] retval   return value. 1 f field was read, 0 otherwise

  subroutine restart_get_field_read_status(f_id, retval)
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, intent(in)  :: f_id
    integer, intent(out) :: retval

    ! Local variables
    integer(c_int) :: c_f_id, c_retval
    c_f_id = f_id

    c_retval = cs_restart_get_field_read_status(c_f_id)
    retval = c_retval

  end subroutine restart_get_field_read_status

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

  !> \brief  Wrapper to Fortran user boundary condition definitions.

  !> \param[in, out]  bc_type  boundary face types

  subroutine user_f_boundary_conditions(itypfb, izfppp, dt)  &
    bind(C, name='cs_f_user_boundary_conditions_wrapper')

    use dimens
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer(kind=c_int), dimension(*), intent(inout) :: itypfb, izfppp
    real(c_double), dimension(*), intent(in) :: dt

    ! Externals

    procedure() :: cs_f_user_boundary_conditions

    ! Local variables

    integer, pointer, dimension(:,:) :: icodcl
    double precision, pointer, dimension(:,:,:) :: rcodcl

    call field_build_bc_codes_all(icodcl, rcodcl) ! Get map

    call cs_f_user_boundary_conditions &
          (nvar, nscal, icodcl, c_null_ptr, itypfb, izfppp, dt, rcodcl)

  end subroutine user_f_boundary_conditions

  !=============================================================================

  !> \brief  Wrapper to Fortran user parameters, model choice

  subroutine cs_f_usppmo()  &
    bind(C, name='cs_f_usppmo')

    use dimens
    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    procedure() :: usppmo

    ! Local variables

    call usppmo(1)

  end subroutine cs_f_usppmo

  !=============================================================================

  !> \brief  Wrapper to Fortran user parameters

  subroutine cs_f_usipsu(nmodpp)  &
    bind(C, name='cs_f_usipsu')

    use dimens
    use, intrinsic :: iso_c_binding
    implicit none

    procedure() :: usipsu

    ! Arguments

    integer(c_int) :: nmodpp

    ! Local variables

    call usipsu(nmodpp)

  end subroutine cs_f_usipsu

  !=============================================================================

  !> \brief  Wrapper to Fortran user parameters, additional parameters

  subroutine cs_f_usipes(nmodpp)  &
    bind(C, name='cs_f_usipes')

    use dimens
    use, intrinsic :: iso_c_binding
    implicit none

    procedure() :: usipes

    ! Arguments

    integer(c_int) :: nmodpp

    ! Local variables

    call usipes(nmodpp)

  end subroutine cs_f_usipes

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

  !> \brief  Add a CDO field defining a general solved variable, with default
  !>         options.

  !> \param[in]  name           field name
  !> \param[in]  label          field default label, or empty
  !> \param[in]  location_id    field location type:
  !>                              0: none
  !>                              1: cells
  !>                              2: interior faces
  !>                              3: interior faces
  !>                              4: vertices
  !> \param[in]  dim            field dimension
  !> \param[in]  has_previous   if greater than 1 then store previous state
  !> \param[out] id             id of defined field

  subroutine variable_cdo_field_create(name, label, location_id, dim, &
                                       has_previous, id)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name, label
    integer, intent(in)          :: location_id, dim, has_previous
    integer, intent(out)         :: id

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name
    character(len=len_trim(label)+1, kind=c_char) :: c_label
    integer(c_int) :: c_location_id, c_dim, c_has_previous, c_id

    c_name = trim(name)//c_null_char
    c_label = trim(label)//c_null_char
    c_location_id = location_id
    c_dim = dim
    c_has_previous = has_previous;

    c_id = cs_variable_cdo_field_create(c_name, c_label, c_location_id, &
                                        c_dim, c_has_previous)

    id = c_id

    return

  end subroutine variable_cdo_field_create

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

  !> \brief Indicate of a cell is active fo the current variable

  !> \param[in]     iel       cell number (cell_id + 1)
  !> \result        is_active

  function cell_is_active(iel) result(is_active)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer :: iel
    integer :: is_active

    ! Local variables

    integer(kind=c_int) :: c_cell_id
    integer(kind=c_int) :: c_is_active


    c_cell_id = iel - 1
    c_is_active = cs_f_porous_model_cell_is_active(c_cell_id)
    is_active = c_is_active

  end function cell_is_active

  !=============================================================================

  !> \brief Sets the chemistry concentration file name

  !> \param[in]     name      name of the file

  subroutine atmo_set_chem_conc_file_name(name)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_name = trim(name)//c_null_char
    call cs_atmo_set_chem_conc_file_name(c_name)

  end subroutine atmo_set_chem_conc_file_name

  !=============================================================================

  !> \brief Sets the file name used to initialize the aerosol shared library

  !> \param[in]     name      name of the file

  subroutine atmo_chemistry_set_aerosol_file_name(name)

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    character(len=*), intent(in) :: name

    ! Local variables

    character(len=len_trim(name)+1, kind=c_char) :: c_name

    c_name = trim(name)//c_null_char
    call cs_atmo_chemistry_set_aerosol_file_name(c_name)

  end subroutine atmo_chemistry_set_aerosol_file_name

  !=============================================================================
  !> brief Clipping scalar field.
  ! \param[in]   iscal

  subroutine clpsca(iscal)

    use, intrinsic :: iso_c_binding
    use numvar
    implicit none

    ! Arguments
    integer, intent(in) :: iscal

    ! Local variables

    integer(c_int) :: f_id

    f_id = ivarfl(isca(iscal))

    call cs_f_scalar_clipping(f_id)

  end subroutine clpsca

  !=============================================================================

end module cs_c_bindings
