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

!> \file optcal.f90
!> Module for calculation options

module optcal

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup optcal Module for calculation options

  !> \addtogroup optcal
  !> \{

  !----------------------------------------------------------------------------
  ! Space discretisation
  !----------------------------------------------------------------------------

  !> Indicator of a calculation restart (=1) or not (=0).
  !> This value is set automatically by the code; depending on
  !> whether a restart directory is present, and should not be modified by
  !> the user (no need for C mapping).
  integer, save :: isuite = 0

  !----------------------------------------------------------------------------
  ! Time stepping options
  !----------------------------------------------------------------------------

  !> \defgroup time_step_options Time step options and variables

  !> \addtogroup time_step_options
  !> \{

  !> Absolute time step number for previous calculation.
  !>
  !> In the case of a restart calculation, \ref ntpabs
  !> is read from the restart file. Otherwise, it is
  !> initialised to 0 \ref ntpabs is initialised
  !> automatically by the code, its value is not to be
  !> modified by the user.
  integer(c_int), pointer, save :: ntpabs

  !> Current absolute time step number.
  !> In case of restart, this is equal to ntpabs + number of new iterations.
  integer(c_int), pointer, save :: ntcabs

  !> Maximum absolute time step number.
  !>
  !> For the restart calculations, \ref ntmabs takes into
  !> account the number of time steps of the previous calculations.
  !> For instance, after a first calculation of 3 time steps, a
  !> restart file of 2 time steps is realised by setting
  !> \ref ntmabs = 3+2 = 5
  integer(c_int), pointer, save :: ntmabs

  !> Number of time steps for initalization (for all steps between
  !> 0 and \ref ntinit, pressure is re-set to 0 before prediction
  !> correction).
  integer(c_int), pointer, save :: ntinit

  !> Absolute time value for previous calculation.
  !>
  !> In the case of a restart calculation, \ref ttpabs is read from
  !> the restart file. Otherwise it is initialised to 0.\n
  !> \ref ttpabs is initialised automatically by the code,
  !> its value is not to be modified by the user.
  real(c_double), pointer, save :: ttpabs

  !> Current absolute time.
  !>
  !> For the restart calculations, \ref ttcabs takes
  !> into account the physical time of the previous calculations.\n
  !> If the time step is uniform (\ref idtvar = 0 or 1), \ref ttcabs
  !> increases of \ref dt (value of the time step) at each iteration.
  !> If the time step is non-uniform (\ref idtvar=2), \ref ttcabs
  !> increases of \ref dtref at each time step.\n
  !> \ref ttcabs} is initialised and updated automatically by the code,
  !> its value is not to be modified by the user.
  real(c_double), pointer, save :: ttcabs

  !> Maximum absolute time.
  real(c_double), pointer, save :: ttmabs

  !> option for a variable time step
  !>    - -1: steady algorithm
  !>    -  0: constant time step
  !>    -  1: time step constant in space but variable in time
  !>    -  2: variable time step in space and in time
  !> If the numerical scheme is a second-order in time, only the
  !> option 0 is allowed.
  integer(c_int), pointer, save :: idtvar

  !> Reference time step
  !>
  !> This is the time step value used in the case of a calculation run with a
  !> uniform and constant time step, i.e. \ref idtvar =0 (restart calculation
  !> or not). It is the value used to initialize the time step in the case of
  !> an initial calculation run with a non-constant time step(\ref idtvar=1 or
  !> 2). It is also the value used to initialise the time step in the case of
  !> a restart calculation in which the type of time step has been changed
  !> (for instance, \ref idtvar=1 in the new calculation and \ref idtvar = 0 or
  !> 2 in the previous calculation).\n
  !> See \subpage user_initialization_time_step for examples.
  real(c_double), pointer, save :: dtref

  !> \}

  !----------------------------------------------------------------------------
  ! thermal model
  !----------------------------------------------------------------------------

  !> \defgroup thermal model

  !> \addtogroup thermal
  !> \{

  !> thermal model
  !>    - 0: no thermal model
  !>    - 1: temperature
  !>    - 2: enthalpy
  !>    - 3: total energy (only for compressible module)\n
  !> When a particular physics module is activated (gas combustion,
  !> pulverised coal, electricity or compressible), the user must not
  !> modify \ref itherm (the choice is made automatically: the solved
  !> variable is either the enthalpy or the total energy). The user is
  !> also reminded that, in the case of a coupling with SYRTHES, the
  !> solved thermal variable should be the temperature (\ref itherm = 1).
  !> More precisely, everything is designed in the code to allow for the
  !> running of a calculation coupled with SYRTHES with the enthalpy as
  !> thermal variable. With the compressible model, it is possible to
  !> carry out calculations coupled with SYRTHES, although the thermal
  !> scalar represents the total energy and not the temperature.
  integer(c_int), pointer, save :: itherm

  !> Temperature scale
  !> - 0: none
  !> - 1: Kelvin
  !> - 2: Celsius
  !> The distinction between \ref itpscl = 1 or 2 is useful only in case of
  !> radiation modelling. For calculations without radiation modelling,
  !> use \ref itpscl = 1 for the temperature.\n
  !> Useful if and only if \ref dimens::nscal "nscal" \f$\geqslant\f$ 1.
  integer(c_int), pointer, save :: itpscl

  !> Index of the thermal scalar (temperature, energy or enthalpy)
  !>
  !> The index of the corresponding variable is isca(iscalt)
  !> If \ref iscalt = -1, neither the temperature nor the enthalpy is
  !> represented by a scalar. When a specific physics module is activated
  !> (gas combustion, pulverised coal, electricity or compressible), the user
  !> must not modify \ref iscalt (the choice is made automatically). In the
  !> case of the compressible module, \ref iscalt does not correspond to
  !> the temperature nor enthalpy but to the total energy}.\n Useful if
  !> and only if \ref dimens::nscal "nscal" \f$\geqslant\f$ 1.
  integer, save :: iscalt = -1

  !> \}

  !----------------------------------------------------------------------------
  ! turbulence
  !----------------------------------------------------------------------------

  !> \defgroup turbulence turbulence options

  !> \addtogroup turbulence
  !> \{

  !> \anchor iturb
  !> turbulence model
  !>    - 0: no turbulence model (laminar flow)
  !>    - 10: mixing length model
  !>    - 20: standard \f$ k-\varepsilon \f$ model
  !>    - 21: \f$ k-\varepsilon \f$ model with Linear Production (LP) correction
  !>    - 30: \f$ R_{ij}-\epsilon \f$ (LRR)
  !>    - 31: \f$ R_{ij}-\epsilon \f$ (SSG)
  !>    - 32: \f$ R_{ij}-\epsilon \f$ (EBRSM)
  !>    - 40: LES (constant Smagorinsky model)
  !>    - 41: LES ("classical" dynamic Smagorisky model)
  !>    - 42: LES (WALE)
  !>    - 50: v2f phi-model
  !>    - 51: v2f \f$ BL-v^2-k \f$
  !>    - 60: \f$ k-\omega \f$ SST
  !>    - 70: Spalart-Allmaras model
  integer(c_int), pointer, save :: iturb

  !> Class of turbulence model (integer value iturb/10)
  integer(c_int), pointer, save :: itytur

  !> \}

  !----------------------------------------------------------------------------
  ! Stokes
  !----------------------------------------------------------------------------

  !> \defgroup stokes Stokes options

  !> \addtogroup stokes
  !> \{

  !> Algorithm to take into account the density variation in time
  !>    - 0: boussinesq algorithm with constant density
  !>    - 1: dilatable steady algorithm (default)
  !>    - 2: dilatable unsteady algorithm
  !>    - 3: low-Mach algorithm
  !>    - 4: algorithm for fire
  integer(c_int), pointer, save :: idilat

  !> \}

  !----------------------------------------------------------------------------
  ! Transported scalars parameters
  !----------------------------------------------------------------------------

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global time step structure

    subroutine cs_f_time_step_get_pointers(nt_prev, nt_cur, nt_max, nt_ini,  &
                                           dt_ref, t_prev, t_cur, t_max)     &
      bind(C, name='cs_f_time_step_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: nt_prev, nt_cur, nt_max, nt_ini
      type(c_ptr), intent(out) :: dt_ref, t_prev, t_cur, t_max
    end subroutine cs_f_time_step_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global time step options structure

    subroutine cs_f_time_step_options_get_pointers(idtvar)         &
      bind(C, name='cs_f_time_step_options_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: idtvar
    end subroutine cs_f_time_step_options_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global thermal model structure

    subroutine cs_f_thermal_model_get_pointers(itherm, itpscl) &
      bind(C, name='cs_f_thermal_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: itherm, itpscl
    end subroutine cs_f_thermal_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global turbulence model structure

    subroutine cs_f_turb_model_get_pointers(iturb, itytur) &
      bind(C, name='cs_f_turb_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iturb, itytur
    end subroutine cs_f_turb_model_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief If scalar iscal represents the mean of the square of a scalar
  !> k, return k; otherwise, return 0.

  function iscavr(iscal) result(iscvr)

    use field
    use numvar

    implicit none

    ! Parameters

    integer, intent(in) :: iscal
    integer             :: iscvr

    ! Local arguments

    integer :: f_id
    integer :: kscavr = -1
    integer :: keysca = -1

    ! Function body

    iscvr = 0

    if (kscavr .lt. 0) then
      call field_get_key_id("first_moment_id", kscavr)
      call field_get_key_id("scalar_id", keysca)
    endif

    if (kscavr.ge.0) then
      call field_get_key_int(ivarfl(isca(iscal)), kscavr, f_id)
      if (f_id.ge.0) call field_get_key_int(f_id, keysca, iscvr)
    endif

  end function iscavr

  !> \brief If scalar iscal represents the mean of the square of a scalar
  !> k, return k; otherwise, return 0.

  function visls0(iscal) result(visls_0)

    use field
    use numvar

    implicit none

    ! Parameters

    integer, intent(in) :: iscal
    double precision    :: visls_0

    call field_get_key_double(ivarfl(isca(iscal)), kvisl0, visls_0)

  end function visls0

  !> \brief Initialize isuite

  subroutine indsui () &
    bind(C, name='cs_f_indsui')

    use, intrinsic :: iso_c_binding
    use cs_c_bindings
    implicit none

    isuite = cs_restart_present()

  end subroutine indsui

  !> \brief Initialize Fortran time step API.
  !> This maps Fortran pointers to global C structure members.

  subroutine time_step_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ntpabs, c_ntcabs, c_ntmabs, c_ntinit
    type(c_ptr) :: c_dtref, c_ttpabs, c_ttcabs, c_ttmabs

    call cs_f_time_step_get_pointers(c_ntpabs, c_ntcabs, c_ntmabs, c_ntinit, &
                                     c_dtref, c_ttpabs, c_ttcabs, c_ttmabs)

    call c_f_pointer(c_ntpabs, ntpabs)
    call c_f_pointer(c_ntcabs, ntcabs)
    call c_f_pointer(c_ntmabs, ntmabs)
    call c_f_pointer(c_ntinit, ntinit)

    call c_f_pointer(c_dtref,  dtref)
    call c_f_pointer(c_ttpabs, ttpabs)
    call c_f_pointer(c_ttcabs, ttcabs)
    call c_f_pointer(c_ttmabs, ttmabs)

  end subroutine time_step_init

  !> \brief Initialize Fortran time step options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine time_step_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_idtvar

    call cs_f_time_step_options_get_pointers(c_idtvar)

    call c_f_pointer(c_idtvar, idtvar)

  end subroutine time_step_options_init

  !> \brief Initialize Fortran thermal model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine thermal_model_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_itherm, c_itpscl

    call cs_f_thermal_model_get_pointers(c_itherm, c_itpscl)

    call c_f_pointer(c_itherm, itherm)
    call c_f_pointer(c_itpscl, itpscl)

  end subroutine thermal_model_init

  !> \brief Initialize Fortran turbulence model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine turb_model_init

    use, intrinsic :: iso_c_binding
    use cs_c_bindings
    implicit none

    ! Local variables

    type(c_ptr) :: c_iturb, c_itytur

    call cs_f_turb_model_get_pointers(c_iturb, c_itytur)

    call c_f_pointer(c_iturb, iturb)
    call c_f_pointer(c_itytur, itytur)

  end subroutine turb_model_init

  !=============================================================================

end module optcal
