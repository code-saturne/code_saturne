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

  !> latent heat of evaporation
  real(c_double), pointer, save :: clatev

  !> Boltzmann constant (\f$J.K^{-1}\f$)
  double precision kboltz
  parameter(kboltz = 1.380649d-23)

  !> Ideal gas constant (\f$J.mol^{-1}.K^{-1}\f$)
  double precision cs_physical_constants_r
  parameter(cs_physical_constants_r = 8.31446261815324d0)

  !> Gravity
  real(c_double), pointer, save :: gx, gy, gz

  !> Coriolis effects
  integer(c_int), pointer, save :: icorio

  !> Physical constants of the fluid
  !> filling \ref xyzp0 indicator
  integer(c_int), pointer, save :: ixyzp0

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

  !> isochoric specific heat \f$ C_v \f$
  integer(c_int), pointer, save :: icv

  !> variable density field \f$ \rho \f$:
  !>    - 1: true, its variation law be given either
  !> in the GUI, or in the user subroutine
  !> \ref cs_user_physical_properties .\n
  !> See \subpage physical_properties for more informations.
  !>    - 0: false, its value is the reference density
  !> \ref ro0.
  integer(c_int), pointer, save :: irovar

  !> variable viscosity field \f$ \mu \f$:
  !>    - 1: true, its variation law be given either
  !> in the GUI, or in the user subroutine
  !> \ref cs_user_physical_properties .\n
  !> See \subpage physical_properties for more informations.
  !>    - 0: false, its value is the reference molecular
  !> dynamic viscosity \ref viscl0
  integer(c_int), pointer, save :: ivivar

  !> Sutherland law for laminar viscosity and thermal conductivity
  !> Only useful in gas mix (igmix) specific physics
  !> - 1: Sutherland law
  !> - 0: low temperature law (linear except for helium)
  integer(c_int), pointer, save :: ivsuth

  !> reference density.\n
  !>
  !> Negative value: not initialized.
  !> Its value is not used in gas or coal combustion modelling (it will be
  !> calculated following the perfect gas law, with \f$P0\f$ and \f$T0\f$).
  !> With the compressible module, it is also not used by the code,
  !> but it may be (and often is) referenced by the user in user subroutines;
  !> it is therefore better to specify its value.
  !>
  !> Always useful otherwise, even if a law defining the density is given by
  !> the user subroutines \ref cs_user_physical_properties.
  !> indeed, except with the
  !> compressible module, CS  does not use
  !> the total pressure \f$P\f$ when solving the Navier-Stokes equation, but a
  !> reduced pressure .
  !> \f$P^*=P-\rho_0\vect{g}.(\vect{x}-\vect{x}_0)+P^*_0-P_0\f$.
  !> where
  !> \f$\vect{x_0}\f$ is a reference point (see \ref xyzp0) and \f$P^*_0\f$ and \f$P_0\f$ are
  !> reference values (see \ref pred0 and \ref p0). Hence, the term
  !> \f$-\grad{P}+\rho\vect{g}\f$ in the equation is treated as
  !> \f$-\grad{P^*}+(\rho-\rho_0)\vect{g}\f$. The closer \ref ro0 is to the value of \f$\rho\f$,
  !> the more \f$P^*\f$ will tend to represent only the dynamic part of the pressure and
  !> the faster and more precise its solution will be. Whatever the value of \ref ro0,
  !> both \f$P\f$ and \f$P^*\f$ appear in the log and the post-processing outputs..
  !> with the compressible module, the calculation is made directly on the total
  !> pressure
  real(c_double), pointer, save :: ro0

  !> reference molecular dynamic viscosity.\n
  !>
  !> Negative value: not initialized.
  !> Always useful, it is the used value unless the user specifies the
  !> viscosity in the subroutine \ref cs_user_physical_properties
  real(c_double), pointer, save :: viscl0

  !> reference pressure for the total pressure.\n
  !>
  !> except with the compressible module, the total pressure \f$P\f$ is evaluated
  !> from the reduced pressure \f$P^*\f$ so that \f$P\f$
  !> is equal to \ref p0 at the reference position \f$\vect{x}_0\f$
  !> (given by \ref xyzp0).
  !> with the compressible module, the total pressure is solved directly.
  !> always Useful
  real(c_double), pointer, save :: p0

  !> reference value for the reduced pressure \f$P^*\f$ (see \ref ro0).\n
  !>
  !> It is especially used to initialise the reduced pressure and as a reference
  !> value for the outlet boundary conditions.
  !> For an optimised precision in the resolution of \f$P^*\f$,
  !>  it is wiser to keep \ref pred0 to 0.
  !> With the compressible module, the "pressure" variable appearing in the
  !> equations directly represents the total pressure.
  !> It is therefore initialized to \ref p0 and not \ref pred0 (see \ref ro0).
  !> Always useful, except with the compressible module
  real(c_double), pointer, save :: pred0

  !> coordinates of the reference point \f$\vect{x}_0\f$ for the total pressure.
  !>  - When there are no Dirichlet conditions for the pressure (closed domain),
  !> \ref xyzp0
  !> does not need to be specified (unless the total pressure has a clear
  !> physical meaning in the configuration treated).
  !>  - When Dirichlet conditions on the pressure are specified but only through
  !> stantard outlet conditions (as it is in most configurations),
  !> \ref xyzp0 does not need to be specified by the user, since it will be
  !> set to the coordinates of the reference outlet face
  !> (i.e. the code will automatically select a
  !> reference outlet boundary face and set \ref xyzp0 so that \f$P\f$ equals
  !> \ref p0 at this face). Nonetheless, if \ref xyzp0 is specified by
  !> the user, the calculation will remain correct.
  !>  - When direct Dirichlet conditions are specified by the user (specific
  !> value set on specific boundary faces), it is better to specify the
  !> corresponding reference point (\em i.e. specify where the total pressure
  !> is \ref p0). This way, the boundary conditions for the reduced pressure
  !> will be close to \ref pred0, ensuring an optimal precision in the
  !> resolution. If \ref xyzp0 is not specified, the reduced
  !> pressure will be shifted, but the calculations will remain correct..
  !>  - With the compressible module, the "pressure" variable appearing in the
  !> equations directly represents the total pressure.
  !> \ref xyzp0 is therefore not used..
  !>
  !> Always useful, except with the compressible module.
  real(c_double), pointer, save :: xyzp0(:)

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

  !> Reference isochoric specific heat.
  !>
  !> Useful for the compressible module (J/kg/K)
  real(c_double), pointer, save :: cv0

  !> Reference isobaric specific heat for water vapor.
  !>
  !> Useful for the moist air module (J/kg/K)
  real(c_double), pointer, save :: cpv0

  !> Latent heat.
  !>
  !> Useful for the moist air module (J/kg)
  real(c_double), pointer, save :: l00

  !> Reference specific heat for liquid water.
  !>
  !> Useful for the moist air module (J/kg/K)
  real(c_double), pointer, save :: cvl

  !> Reference thermal conductivity.
  !>
  !> Useful with thermal model (W/m/K)
  real(c_double), pointer, save :: lambda0

  !> Molar mass of the perfect gas in \f$ kg/mol \f$
  !> (if \ref cstphy::ieos "ieos"=1)
  !>
  !> Always useful
  real(c_double), pointer, save :: xmasmr

  !> Uniform variable thermodynamic pressure for the low-Mach algorithm
  !>    - 1: true
  !>    - 0: false
  integer(c_int), pointer, save :: ipthrm

  !> Thermodynamic pressure for the current time step
  real(c_double), pointer, save :: pther

  !> Thermodynamic pressure for the previous time step
  real(c_double), pointer, save :: pthera

  !>   pthermax: Thermodynamic maximum pressure for user clipping,
  !>             used to model a venting effect
  real(c_double), pointer, save :: pthermax

  !> Leak surface
  real(c_double), pointer, save :: sleak

  !> Leak head loss (2.9 by default, from Idelcick)
  real(c_double), pointer, save :: kleak

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

  !> limit value of \f$y^+\f$ for the viscous sublayer.
  !> \ref ypluli depends on the chosen wall function: it is
  !> initialized to 10.88 for the scalable wall function
  !> (\ref optcal::iwallf "iwallf"=4),
  !> otherwise it is initialized to \f$1/\kappa\approx 2,38\f$.
  !> In LES, \ref ypluli is taken by default to be 10.88.
  !>
  !> Always useful
  real(c_double), pointer, save :: ypluli

  !> constant \f$C_\mu\f$ for all the RANS turbulence models
  !> Warning, different values for the v2f model
  !> Useful if and only if \ref iturb = 20, 21, 30, 31, 50, 51 or 60
  !> (\f$k-\varepsilon\f$, \f$R_{ij}-\varepsilon\f$ or \f$k-\omega\f$)
  real(c_double), pointer, save :: cmu

  !> Coefficient of interfacial coefficient in k-eps,
  !> used in Lagrange treatment
  !>

  !> constant \f$C_s\f$ for the \f$R_{ij}-\varepsilon\f$ models.
  real(c_double), pointer, save :: csrij

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

  !> van Driest constant appearing in the van Driest damping function
  !> applied to the Smagorinsky constant:
  !>  - \f$ (1-\exp^{(-y^+/cdries}) \f$.
  !>
  !> Useful if and only if \ref iturb = 40 or 41
  real(c_double), pointer, save :: cdries

  !> minimal control volume
  double precision, save :: volmin
  !> maximal control volume
  double precision, save :: volmax
  !> total domain volume
  double precision, save :: voltot

  !> constant of EB-AFM and EB-DFM
  real(c_double), pointer, save :: xclt

  !> \}

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global physical constants structure

    subroutine cs_f_physical_constants_get_pointers(gx, gy, gz,              &
                                                    icorio)                  &
      bind(C, name='cs_f_physical_constants_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: gx, gy, gz, icorio
    end subroutine cs_f_physical_constants_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global fluid properties structure

    subroutine cs_f_fluid_properties_get_pointers(ixyzp0,  &
                                                  icp,     &
                                                  icv,     &
                                                  irovar,  &
                                                  ivivar,  &
                                                  ivsuth,  &
                                                  ro0,     &
                                                  viscl0,  &
                                                  p0,      &
                                                  pred0,   &
                                                  xyzp0,   &
                                                  t0,      &
                                                  cp0,     &
                                                  cv0,     &
                                                  cpv0,    &
                                                  cvl,     &
                                                  l00,     &
                                                  lambda0, &
                                                  rair,    &
                                                  rvapor,  &
                                                  rvsra,   &
                                                  clatev,  &
                                                  xmasmr,  &
                                                  ipthrm,  &
                                                  pther,   &
                                                  pthera,  &
                                                  pthermax,&
                                                  sleak,   &
                                                  kleak,   &
                                                  roref)   &
      bind(C, name='cs_f_fluid_properties_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ixyzp0, icp, icv, irovar, ivivar, ivsuth
      type(c_ptr), intent(out) :: ro0, viscl0, p0, pred0
      type(c_ptr), intent(out) :: xyzp0, t0, cp0, cv0, cpv0, cvl, l00, lambda0
      type(c_ptr), intent(out) :: rair, rvapor, rvsra, clatev, xmasmr
      type(c_ptr), intent(out) :: ipthrm
      type(c_ptr), intent(out) :: pther, pthera, pthermax
      type(c_ptr), intent(out) :: sleak, kleak, roref
    end subroutine cs_f_fluid_properties_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! RANS turbulence model structure

    subroutine cs_f_turb_reference_values(almax, uref) &
      bind(C, name='cs_f_turb_reference_values')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: almax , uref
    end subroutine cs_f_turb_reference_values

    ! Interface to C function retrieving pointers to members of the
    ! global wall functions structure

    subroutine cs_f_wall_reference_values(ypluli) &
      bind(C, name='cs_f_wall_reference_values')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ypluli
    end subroutine cs_f_wall_reference_values

    !---------------------------------------------------------------------------

    ! Interface to C function completing the constants of the
    ! turbulence model

    subroutine cs_f_turb_complete_constants(phase_id) &
      bind(C, name='cs_turb_compute_constants')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: phase_id
    end subroutine cs_f_turb_complete_constants

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to constants of the
    ! turbulence model

    subroutine cs_f_turb_model_constants_get_pointers(cmu, csmago,             &
         xlesfd, xlesfl, ales, bles, cdries, csrij, xclt) &

      bind(C, name='cs_f_turb_model_constants_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: cmu, csmago
      type(c_ptr), intent(out) :: xlesfd, xlesfl, cdries
      type(c_ptr), intent(out) :: ales, bles, csrij, xclt
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

    type(c_ptr) :: c_gx, c_gy, c_gz, c_icorio

    call cs_f_physical_constants_get_pointers(c_gx, c_gy, c_gz, c_icorio)

    call c_f_pointer(c_gx, gx)
    call c_f_pointer(c_gy, gy)
    call c_f_pointer(c_gz, gz)
    call c_f_pointer(c_icorio, icorio)

  end subroutine physical_constants_init

  !> \brief Initialize Fortran fluid properties API.
  !> This maps Fortran pointers to global C structure members.

  subroutine fluid_properties_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ixyzp0, c_icp, c_icv, c_irovar, c_ivivar
    type(c_ptr) :: c_ivsuth, c_ro0, c_viscl0, c_p0
    type(c_ptr) :: c_pred0, c_xyzp0, c_t0, c_cp0, c_cv0
    type(c_ptr) :: c_cpv0, c_cvl, c_l00, c_lambda0
    type(c_ptr) :: c_rair,c_rvapor, c_rvsra, c_clatev, c_xmasmr
    type(c_ptr) :: c_ipthrm
    type(c_ptr) :: c_pther, c_pthera, c_pthermax
    type(c_ptr) :: c_sleak, c_kleak, c_roref

    call cs_f_fluid_properties_get_pointers(c_ixyzp0, c_icp, c_icv,         &
                                            c_irovar, c_ivivar, c_ivsuth,   &
                                            c_ro0, c_viscl0, c_p0, c_pred0, &
                                            c_xyzp0, c_t0, c_cp0, c_cv0,    &
                                            c_cpv0, c_cvl, c_l00,           &
                                            c_lambda0,                      &
                                            c_rair,c_rvapor, c_rvsra,       &
                                            c_clatev,                       &
                                            c_xmasmr,                       &
                                            c_ipthrm, c_pther, c_pthera,    &
                                            c_pthermax, c_sleak, c_kleak,   &
                                            c_roref)

    call c_f_pointer(c_ixyzp0, ixyzp0)
    call c_f_pointer(c_icp, icp)
    call c_f_pointer(c_icv, icv)
    call c_f_pointer(c_irovar, irovar)
    call c_f_pointer(c_ivivar, ivivar)
    call c_f_pointer(c_ivsuth, ivsuth)
    call c_f_pointer(c_ro0, ro0)
    call c_f_pointer(c_viscl0, viscl0)
    call c_f_pointer(c_p0, p0)
    call c_f_pointer(c_pred0, pred0)
    call c_f_pointer(c_xyzp0, xyzp0, [3])
    call c_f_pointer(c_t0, t0)
    call c_f_pointer(c_cp0, cp0)
    call c_f_pointer(c_cv0, cv0)
    call c_f_pointer(c_cpv0, cpv0)
    call c_f_pointer(c_cvl, cvl)
    call c_f_pointer(c_l00, l00)
    call c_f_pointer(c_lambda0, lambda0)
    call c_f_pointer(c_rair, rair)
    call c_f_pointer(c_rvapor, rvapor)
    call c_f_pointer(c_rvsra, rvsra)
    call c_f_pointer(c_clatev, clatev)
    call c_f_pointer(c_xmasmr, xmasmr)
    call c_f_pointer(c_ipthrm, ipthrm)
    call c_f_pointer(c_pther, pther)
    call c_f_pointer(c_pthera, pthera)
    call c_f_pointer(c_pthermax, pthermax)
    call c_f_pointer(c_sleak, sleak)
    call c_f_pointer(c_kleak, kleak)
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

    type(c_ptr) :: c_cmu, c_csmago
    type(c_ptr) :: c_xlesfd, c_xlesfl, c_ales, c_bles, c_cdries, c_csrij, c_xclt

    call cs_f_turb_model_constants_get_pointers(c_cmu, c_csmago,               &
                                                c_xlesfd, c_xlesfl,            &
                                                c_ales, c_bles,                &
                                                c_cdries, c_csrij, c_xclt)

    call c_f_pointer(c_cmu   , cmu)
    call c_f_pointer(c_csmago, csmago)
    call c_f_pointer(c_xlesfd, xlesfd)
    call c_f_pointer(c_xlesfl, xlesfl)
    call c_f_pointer(c_ales  , ales  )
    call c_f_pointer(c_bles  , bles  )
    call c_f_pointer(c_cdries, cdries)
    call c_f_pointer(c_csrij , csrij )
    call c_f_pointer(c_xclt,   xclt)

  end subroutine turb_model_constants_init

  !=============================================================================

end module cstphy
