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

!> \file cstphy.f90
!> \brief Module for physical constants

module cstphy

  !=============================================================================

  use, intrinsic :: iso_c_binding

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup cstphy Module for physical constants

  !> \addtogroup cstphy
  !> \{

  !> Temperature in Kelvin correponding to 0 degrees Celsius (= +273,15)
  double precision :: tkelvi
  parameter(tkelvi = 273.15d0)

  !> Calories (1 cvar_al = xcal2j J)
  double precision :: xcal2j
  parameter(xcal2j = 4.1855d0)

  !> Stephan constant for the radiative module \f$\sigma\f$
  !> in \f$W.m^{-2}.K^{-4}\f$
  double precision :: stephn
  parameter(stephn = 5.6703d-8)

  !> Perfect gas constant for air (mixture)
  real(c_double), pointer, save :: rair

  !> Moist air gas constant (mixture)
  real(c_double), pointer, save :: rvapor

  !> ratio gas constant h2o/ dry air
  real(c_double), pointer, save :: rvsra

  !> Boltzmann constant (\f$J.K^{-1}\f$)
  double precision kboltz
  parameter(kboltz = 1.380649d-23)

  !> Ideal gas constant (\f$J.mol^{-1}.K^{-1}\f$)
  double precision cs_physical_constants_r
  parameter(cs_physical_constants_r = 8.31446261815324d0)

  !> Gravity
  real(c_double), pointer, save :: gy, gz

  !> indicates if the isobaric specific heat \f$C_p\f$ is variable:
  !>  - 0: constant, no property field is declared
  !>  - 1: variable, \f$C_p\f$ is declared as a property field\n
  !> When gas or coal combustion is activated, \ref icp is automatically set to 0
  !> (constant \f$C_p\f$). With the electric module, it is automatically set to 1.
  !> The user is not allowed to modify these default choices.\n
  !> When \ref icp = 1 is specified, the code automatically modifies this value to
  !> make \ref icp designate the effective index-number of the property "specific heat".
  !> For each cell iel, the value of \f$C_p\f$ is then specified by the user in the
  !> appropriate subroutine (\ref cs_user_physical_properties for the standard physics).\n
  !> Useful if there is at least 1 temperature scalar, or with the compressible module
  !> for non perfect gases.
  integer(c_int), pointer, save :: icp

  !> variable density field \f$ \rho \f$:
  !>    - 1: true, its variation law be given either
  !> in the GUI, or in the user subroutine
  !> \ref cs_user_physical_properties .\n
  !> See \subpage physical_properties for more informations.
  !>    - 0: false, its value is the reference density
  !> \ref ro0.
  integer(c_int), pointer, save :: irovar

  !> reference density.\n
  real(c_double), pointer, save :: ro0

  !> reference molecular dynamic viscosity.\n
  !>
  !> Negative value: not initialized.
  !> Always useful, it is the used value unless the user specifies the
  !> viscosity in the subroutine \ref cs_user_physical_properties
  real(c_double), pointer, save :: viscl0

  !> reference pressure for the total pressure.\n
  real(c_double), pointer, save :: p0

  !> reference temperature.
  !>
  !> Useful for the specific physics gas or coal combustion (initialization
  !> of the density), for the electricity modules to initialize the domain
  !> temperature and for the compressible module (initializations).
  !> It must be given in Kelvin.
  real(c_double), pointer, save :: t0

  !> Reference internal energy for the barotropic compressible module
  double precision, save :: eint0

  !> reference specific heat.
  !>
  !> Useful if there is 1 <= n <= nscaus,
  !> so that \ref optcal::iscalt "iscalt" = n and \ref optcal::itherm "itherm" = 1
  !> (there is a "temperature" scalar),
  !> unless the user specifies the specific heat in the user subroutine
  !> \ref cs_user_physical_properties (\ref cstphy::icp "icp" > 0) with the
  !> compressible module or
  !>  coal combustion, \ref cp0 is also needed even when there is no user scalar.
  !> \note None of the scalars from the specific physics is a temperature.
  !> \note When using the Graphical Interface, \ref cp0 is also used to
  !> calculate the diffusivity of the thermal scalars,
  !> based on their conductivity; it is therefore needed, unless the
  !> diffusivity is also specified in \ref cs_user_physical_properties.
  real(c_double), pointer, save :: cp0

  !> Thermodynamic pressure for the current time step
  real(c_double), pointer, save :: pther

  !> Initial reference density
  real(c_double), pointer, save :: roref

  !> \defgroup csttur Module for turbulence constants

  !> \addtogroup csttur
  !> \{

  !> \f$ \kappa \f$ Karman constant. (= 0.42)
  !> Useful if and only if \ref iturb >= 10.
  !> (mixing length, \f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$,
  !> LES, v2f or \f$k-\omega\f$)
  double precision, save :: xkappa = 0.42d0

  !> constant \f$C_\mu\f$ for all the RANS turbulence models
  !> Warning, different values for the v2f model
  !> Useful if and only if \ref iturb = 20, 21, 30, 31, 50, 51 or 60
  !> (\f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$ or \f$k-\omega\f$)
  real(c_double), pointer, save :: cmu

  !> Coefficient of interfacial coefficient in k-eps,
  !> used in Lagrange treatment
  !>

  !> constant of the Rij-epsilon EBRSM
  double precision, save :: xcl = 0.122d0

  !> constant in the expression of Ce1' for the Rij-epsilon EBRSM
  double precision, save :: xa1 = 0.1d0

  !> constant of the Rij-epsilon EBRSM
  double precision, save :: xct = 6.d0

  !> is a characteristic macroscopic
  !> length of the domain, used for the initialization of the turbulence and
  !> the potential clipping (with \ref optcal::iclkep "iclkep"=1)
  !>  - Negative value: not initialized (the code then uses the cubic root of
  !> the domain volume).
  !>
  !> Useful if and only if \ref iturb = 20, 21, 30, 31, 50 or 60 (RANS models)
  real(c_double), pointer, save :: almax

  !> the characteristic flow velocity,
  !> used for the initialization of the turbulence.
  !> Negative value: not initialized.
  !>
  !> Useful if and only if \ref iturb = 20, 21, 30, 31, 50 or 60 (RANS model)
  !> and the turbulence is not initialized somewhere
  !> else (restart file or subroutine \ref cs\_user\_initialization)
  real(c_double), pointer, save :: uref

  !> constant used in the definition of LES filtering diameter:
  !> \f$ \delta = \text{xlesfl} . (\text{ales} . volume)^{\text{bles}}\f$
  !> \ref xlesfl is a constant used to define, for
  !> each cell \f$\omega_i\f$, the width of the (implicit) filter:
  !> \f$\overline{\Delta}=xlesfl(ales*|\Omega_i|)^{bles}\f$\n
  !> Useful if and only if \ref iturb = 40 or 41
  real(c_double), pointer, save :: xlesfl

  !> constant used to define, for each cell \f$Omega_i\f$,
  !> the width of the (implicit) filter:
  !>  - \f$\overline{\Delta}=xlesfl(ales*|Omega_i|)^{bles}\f$
  !>
  !> Useful if and only if \ref iturb = 40 or 41.
  real(c_double), pointer, save :: ales

  !> constant used to define, for each cell \f$Omega_i\f$,
  !>
  !> the width of the (implicit) filter:
  !>  - \f$\overline{\Delta}=xlesfl(ales*|Omega_i|)^{bles}\f$
  !>
  !> Useful if and only if \ref iturb = 40 or 41
  real(c_double), pointer, save :: bles

  !> Smagorinsky constant used in the Smagorinsky model for LES.
  !> The sub-grid scale viscosity is calculated by
  !> \f$\displaystyle\mu_{sg}=
  !> \rho C_{smago}^2\bar{\Delta}^2\sqrt{2\bar{S}_{ij}\bar{S}_{ij}}\f$
  !> where \f$\bar{\Delta}\f$ is the width of the filter
  !>  and \f$\bar{S}_{ij}\f$ the filtered strain rate.
  !>
  !> Useful if and only if \ref iturb = 40
  !> \note In theory Smagorinsky constant is 0.18.
  !> For a planar canal plan, 0.065 value is rather taken.
  real(c_double), pointer, save :: csmago

  !> ratio between
  !> explicit and explicit filter width for a dynamic model
  !> constant used to define, for each cell \f$\Omega_i\f$,
  !> the width of the explicit filter used in the framework of
  !> the LES dynamic model:
  !> \f$\widetilde{\overline{\Delta}}=xlesfd\overline{\Delta}\f$.
  !>
  !> Useful if and only if \ref iturb = 41
  real(c_double), pointer, save :: xlesfd

  !> minimal control volume
  double precision, save :: volmin

  !> \}

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical constants structure

    subroutine cs_f_physical_constants_get_pointers(gz)     &
      bind(C, name='cs_f_physical_constants_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: gz
    end subroutine cs_f_physical_constants_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global fluid properties structure

    subroutine cs_f_fluid_properties_get_pointers(icp,     &
                                                  irovar,  &
                                                  ro0,     &
                                                  viscl0,  &
                                                  p0,      &
                                                  t0,      &
                                                  cp0,     &
                                                  rair,    &
                                                  rvapor,  &
                                                  rvsra,   &
                                                  pther,   &
                                                  roref)   &
      bind(C, name='cs_f_fluid_properties_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: icp, irovar
      type(c_ptr), intent(out) :: ro0, viscl0
      type(c_ptr), intent(out) :: p0, t0, cp0
      type(c_ptr), intent(out) :: rair, rvapor, rvsra
      type(c_ptr), intent(out) :: pther
      type(c_ptr), intent(out) :: roref
    end subroutine cs_f_fluid_properties_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! RANS turbulence model structure

    subroutine cs_f_turb_reference_values(almax, uref) &
      bind(C, name='cs_f_turb_reference_values')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: almax , uref
    end subroutine cs_f_turb_reference_values

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to constants of the
    ! turbulence model

    subroutine cs_f_turb_model_constants_get_pointers(cmu,      &
         xlesfd, xlesfl, ales, bles)                            &

      bind(C, name='cs_f_turb_model_constants_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: cmu
      type(c_ptr), intent(out) :: xlesfd, xlesfl
      type(c_ptr), intent(out) :: ales, bles
    end subroutine cs_f_turb_model_constants_get_pointers

    !---------------------------------------------------------------------------

    !> (DOXYGEN_SHOULD_SKIP_THIS) \endcond

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran physical constants API.
  !> This maps Fortran pointers to global C structure members.

  subroutine physical_constants_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_gz

    call cs_f_physical_constants_get_pointers(c_gz)

    call c_f_pointer(c_gz, gz)

  end subroutine physical_constants_init

  !> \brief Initialize Fortran fluid properties API.
  !> This maps Fortran pointers to global C structure members.

  subroutine fluid_properties_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_icp, c_irovar
    type(c_ptr) :: c_ro0, c_viscl0, c_p0
    type(c_ptr) :: c_t0, c_cp0
    type(c_ptr) :: c_rair,c_rvapor, c_rvsra
    type(c_ptr) :: c_pther
    type(c_ptr) :: c_roref

    call cs_f_fluid_properties_get_pointers(c_icp, c_irovar,                &
                                            c_ro0, c_viscl0,                &
                                            c_p0, c_t0, c_cp0,              &
                                            c_rair, c_rvapor, c_rvsra,      &
                                            c_pther, c_roref)

    call c_f_pointer(c_icp, icp)
    call c_f_pointer(c_irovar, irovar)
    call c_f_pointer(c_ro0, ro0)
    call c_f_pointer(c_viscl0, viscl0)
    call c_f_pointer(c_p0, p0)
    call c_f_pointer(c_t0, t0)
    call c_f_pointer(c_cp0, cp0)
    call c_f_pointer(c_rair, rair)
    call c_f_pointer(c_rvapor, rvapor)
    call c_f_pointer(c_rvsra, rvsra)
    call c_f_pointer(c_pther, pther)
    call c_f_pointer(c_roref, roref)

  end subroutine fluid_properties_init

  !> \brief Initialize Fortran RANS turbulence model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine turb_reference_values_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_almax , c_uref

    call cs_f_turb_reference_values(c_almax, c_uref)

    call c_f_pointer(c_almax, almax)
    call c_f_pointer(c_uref, uref)

  end subroutine turb_reference_values_init

  !> \brief Initialize Fortran turbulence model constants.
  !> This maps Fortran pointers to global C real numbers.

  subroutine turb_model_constants_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_cmu
    type(c_ptr) :: c_xlesfd, c_xlesfl, c_ales, c_bles

    call cs_f_turb_model_constants_get_pointers(c_cmu, c_xlesfd, c_xlesfl,   &
                                                c_ales, c_bles)

    call c_f_pointer(c_cmu   , cmu)
    call c_f_pointer(c_xlesfd, xlesfd)
    call c_f_pointer(c_xlesfl, xlesfl)
    call c_f_pointer(c_ales  , ales  )
    call c_f_pointer(c_bles  , bles  )

  end subroutine turb_model_constants_init

  !=============================================================================

end module cstphy
