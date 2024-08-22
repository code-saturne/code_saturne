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
  ! Time stepping
  !----------------------------------------------------------------------------

  !> \defgroup time_stepping Time stepping

  !> \addtogroup time_stepping
  !> \{

  !> Time order of time stepping
  !> (see \ref cs_time_scheme_t::time_order).
  integer(c_int), pointer, save :: ischtp

  !> Time order of the mass flux scheme.
  !> (see \ref cs_time_scheme_t::istmpf).
  integer(c_int), pointer, save :: istmpf

  !> Time scheme for source terms of momentum equations
  !> (see \ref cs_time_scheme_t::isno2t).
  integer(c_int), pointer, save :: isno2t

  !> Time scheme for source terms of turbulence equations
  !> (see \ref cs_time_scheme_t::isto2t).
  integer(c_int), pointer, save :: isto2t

  !> initro : =1 if density read from checkpoint file
  integer(c_int), pointer, save :: initro

  !> \}

  !----------------------------------------------------------------------------
  ! Space discretisation
  !----------------------------------------------------------------------------

  !> \defgroup space_discretisation Space discretisation

  !> \addtogroup space_discretisation
  !> \{

  !> \defgroup conv_scheme Convective scheme
  !> \addtogroup conv_scheme
  !> \{

  !> \anchor iflxmw
  !> method to compute interior mass flux due to ALE mesh velocity
  !>    - 1: based on cell center mesh velocity
  !>    - 0: based on nodes displacement
  integer(c_int), pointer, save :: iflxmw

  !> \}

  !> \defgroup gradient_calculation Gradient calculation
  !> \addtogroup gradient_calculation
  !> \{

  !> type of gradient reconstruction
  !>    - 0: iterative process
  !>    - 1: standard least squares method
  !>    - 2: least square method with extended neighborhood
  !>    - 3: least square method with reduced extended neighborhood
  !>    - 4: iterative process initialized by the least squares method
  integer(c_int), pointer, save :: imrgra

  !> \}

  !> \defgroup diffusive_scheme Diffusive scheme
  !> \addtogroup diffusive_scheme
  !> \{

  !> face viscosity field interpolation
  !>    - 1: harmonic
  !>    - 0: arithmetic (default)
  integer(c_int), pointer, save :: imvisf

  !> \}

 !> \}

  !> Indicator of a calculation restart (=1) or not (=0).
  !> This value is set automatically by the code; depending on
  !> whether a restart directory is present, and should not be modified by
  !> the user (no need for C mapping).
  integer, save :: isuite = 0

  !> Indicates the reading (=1) or not (=0) of the auxiliary
  !> calculation restart file\n
  !> Useful only in the case of a calculation restart
  integer(c_int), pointer, save :: ileaux

  !> \anchor isuit1
  !> For the 1D wall thermal module, activation (1) or not(0)
  !> of the reading of the mesh and of the wall temperature
  !> from the restart file
  !> Useful if nfpt1d > 0
  integer, save :: isuit1

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

  !> Activation of Hybrid RANS/LES model (only valid for iturb equal to 60 or 51)
  integer(c_int), pointer, save :: hybrid_turb

  !> Activation of rotation/curvature correction for eddy viscosity turbulence
  !> models
  !>    - 0: false
  !>    - 1: true
  integer(c_int), pointer, save :: irccor

  !> Type of rotation/curvature correction for eddy viscosity turbulence models
  !>    - 1 Cazalbou correction (default when irccor=1 and itytur=2 or 5)
  !>    - 2 Spalart-Shur correction (default when irccor=1 and iturb=60 or 70)
  integer(c_int), pointer, save :: itycor

  !> Turbulent diffusion model for second moment closure
  !>    - 0: scalar diffusivity (Shir model)
  !>    - 1: tensorial diffusivity (Daly and Harlow model, default model)
  integer(c_int), pointer, save :: idirsm

  !>  Wall functions
  !>  Indicates the type of wall function used for the velocity
  !>  boundary conditions on a frictional wall.
  !>  - 0: no wall functions
  !>  - 1: one scale of friction velocities (power law)
  !>  - 2: one scale of friction velocities (log law)
  !>  - 3: two scales of friction velocities (log law)
  !>  - 4: two scales of friction velocities (log law) (scalable wall functions)
  !>  - 5: two scales of friction velocities (mixing length based on V. Driest
  !>       analysis)
  !>  - 6: wall function unifying rough and smooth friction regimes
  !>  - 7: All \f$ y^+ \f$  for low Reynolds models\n
  !>  \ref iwallf is initialised to 2 for \ref iturb = 10, 40, 41 or 70
  !>  (mixing length, LES and Spalart Allmaras).\n
  !>  \ref iwallf is initialised to 0 for \ref iturb = 0, 32, 50 or 51\n
  !>  \ref iwallf is initialised to 3 for \ref iturb = 20, 21, 30, 31 or 60
  !>  (\f$k-\epsilon\f$, \f$R_{ij}-\epsilon\f$ LRR, \f$R_{ij}-\epsilon\f$ SSG and
  !> \f$ k-\omega\f$ SST models).\n
  !>  The v2f model (\ref iturb=50) is not designed to use wall functions
  !>  (the mesh must be low Reynolds).\n
  !>  The value \ref iwallf = 3 is not compatible with \ref iturb=0, 10, 40
  !>  or 41 (laminar, mixing length and LES).\n
  !>  Concerning the \f$k-\epsilon\f$ and \f$R_{ij}-\epsilon\f$ models, the
  !>  two-scales model is usually at least as satisfactory as the one-scale
  !>  model.\n
  !>  The scalable wall function allows to virtually shift the wall when
  !>  necessary in order to be always in a logarithmic layer. It is used to make
  !>  up for the problems related to the use of High-Reynolds models on very
  !>  refined meshes.\n
  !>  Useful if \ref iturb is different from 50.
  integer(c_int), pointer, save :: iwallf

  !>  Wall functions for scalar
  !>    - 0: three layer wall function of Arpaci and Larsen
  !>    - 1: Van Driest wall function
  integer(c_int), pointer, save :: iwalfs

  !> Indicates the clipping method used for \f$k\f$ and
  !> \f$\varepsilon\f$, for the \f$k-\epsilon\f$ and v2f models
  !> - 0: clipping in absolute value
  !> - 1: coupled clipping based on physical relationships\n
  !> Useful if and only if \ref iturb = 20, 21 or 50 (\f$k-\epsilon\f$ and
  !> v2f models). The results obtained with the method corresponding to
  !> \ref iclkep =1 showed in some cases a substantial sensitivity to the
  !> values of the length scale \ref cs_turb_ref_values_t::almax "almax".\n
  !> The option \ref iclkep = 1 is therefore not recommended, and,
  !> if chosen, must be used cautiously.
  integer(c_int), pointer, save :: iclkep

  !> Indicates if the term \f$\frac{2}{3}\grad \rho k\f$
  !> is taken into account in the velocity equation.
  !> - 1: true
  !> - 0: false in the velocity\n
  !> Useful if and only if \ref iturb = 20, 21, 50 or 60.\n
  !> This term may generate non-physical velocities at the wall.
  !> When it is not explicitly taken into account, it is
  !> implicitly included into the pressure.
  integer(c_int), pointer, save :: igrhok

  !> Indicates if the coupling of the source terms of
  !> \f$k\f$ and \f$\epsilon\f$ or \f$k\f$ and \f$\omega\f$
  !> is taken into account or not.
  !> - 1: true,
  !> - 0: false\n
  !> If \ref ikecou = 0 in \f$k-\epsilon\f$ model, the term
  !> in \f$\epsilon\f$ in the equation of \f$k\f$ is made implicit.\n
  !> \ref ikecou is initialised to 0 if \ref iturb = 21 or 60, and
  !> to 1 if \ref iturb = 20.\n
  !> \ref ikecou = 1 is forbidden when using the v2f model (\ref iturb = 50).\n
  !> Useful if and only if \ref iturb = 20, 21 or 60 (\f$k-\epsilon\f$ and
  !> \f$k-\omega\f$ models)
  integer(c_int), pointer, save :: ikecou

  !> Advanced re-init for EBRSM and k-omega models
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: reinit_turb


  !> Coupled solving of \f$ \tens{R} \f$
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: irijco

  !> pseudo eddy viscosity in the matrix of momentum equation to partially
  !> implicit \f$ \divv \left( \rho \tens{R} \right) \f$
  !>    - 1: true
  !>    - 0: false (default)
  !> The goal is to improve the stability of the calculation.
  !> The usefulness of \ref irijnu = 1 has however not been
  !> clearly demonstrated.\n Since the system is solved in
  !> incremental form, this extra turbulent viscosity does
  !> not change the final solution for steady flows. However,
  !> for unsteady flows, the parameter \ref nswrsm should be
  !> increased.\n Useful if and only if \ref iturb = 30 or 31
  !> (\f$R_{ij}-\epsilon\f$ model).
  !>    - 2: Rusanov scheme on the system momentum + Rij
  !>    - 3: Godunov scheme on the system momentum + Rij
  integer(c_int), pointer, save :: irijnu

  !> accurate treatment of \f$ \tens{R} \f$ at the boundary
  !> (see \ref cs_boundary_condition_set_coeffs)
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: irijrb

  !> Indicates if the wall echo terms in
  !> \f$R_{ij}-\epsilon\f$ LRR model are taken into account:
  !> - 1: true,
  !> - 0: false (default)\n
  !> Useful if and only if \ref iturb = 30 (\f$R_{ij}-\epsilon\f$
  !> LRR).\n It is not recommended to take these terms into account:
  !> they have an influence only near the walls, their expression is hardly
  !> justifiable according to some authors and, in the configurations
  !> studied with code_saturne, they did not bring any improvement in
  !> the results.\n
  !> In addition, their use induces an increase in the calculation time.\n
  !> The wall echo terms imply the calculation of the distance to the wall
  !> for every cell in the domain.
  integer(c_int), pointer, save :: irijec

  !> whole treatment of the diagonal part of the diffusion tensor of
  !> \f$ \tens{R} \f$ and \f$ \varepsilon \f$
  !>    - 1: true (default)
  !>    - 0: simplified treatment
  integer(c_int), pointer, save :: idifre

  !> partial implicitation of symmetry BCs of \f$ \tens{R} \f$
  !>    - 1: true (default)
  !>    - 0: false
  integer(c_int), pointer, save :: iclsyr

  !> partial implicitation of wall BCs of \f$ \tens{R} \f$
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: iclptr

  !> Activates or the van Driest wall-damping for the
  !> Smagorinsky constant (the Smagorinsky constant
  !> is multiplied by the damping function
  !> \f$1-e^{-y^+/ cdries}\f$, where \f$y^+\f$
  !> designates the non-dimensional distance to the
  !> nearest wall).
  !>    - 1: true
  !>    - 0: false
  !> The default value is 1 for the Smagorinsky model
  !> and 0 for the dynamic model.\n The van Driest
  !> wall-damping requires the knowledge of the
  !> distance to the nearest wall for each cell
  !> in the domain.
  !> Useful if and only if \ref iturb = 40 or 41
  integer(c_int), pointer, save :: idries

  !> Applied or not the Internal Consistency
  !> Constraint (ICC) for the HTLES model,
  !> in order to recover the correct RANS
  !> behavior when the energy ratio is forced
  !> to one in the RANS region:
  !>   - 1: True (default)
  !>   - 0: False
  !> Useful if and only if \ref hybrid_turb=4
  integer(c_int), pointer, save :: iicc

  !> Applied or not the two-fold shielding
  !> function (\f$f_s(\xi_K,\xi_D)\f$ of HTLES,
  !> to properly control the RANS-to-LES
  !> transition in the vicinity of the wall:
  !>    - 1: True (default)
  !>    - 0: False
  !> Useful if and only if \ref hybrid_turb=4
  integer(c_int), pointer, save :: ishield

  !> Wall boundary condition on omega in k-omega SST
  !> 0: Deprecated Neumann boundary condition
  !> 1: Dirichlet boundary condition consistent with Menter's
  !>    original model: w_wall = 60*nu/(beta*d**2)
  integer, save :: ikwcln = 1

  !> Activates or not the LES balance module
  !> - 0: false (default)
  !> - 1: true
  !> Useful if \ref iturb =40, 41 or 42\n
  integer(c_int), pointer, save :: i_les_balance

  !> number of variable (deprecated, used only for compatibility)
  integer, save :: nvarcl

  !> \}

  !----------------------------------------------------------------------------
  ! Stokes
  !----------------------------------------------------------------------------

  !> \defgroup stokes Stokes options

  !> \addtogroup stokes
  !> \{

  !> Time scheme option:
  !>    - 0: staggered time scheme. On the time grids, the velocity is
  !>         half a time step behind the density and the buoyant scalar.
  !>         (See the thesis of \cite Pierce:2004)
  !>    - 1: collocated time scheme. On the time grids, the velocity is
  !>         at the same location as the density and the buoyant scalar.
  !>         (See \cite Ma:2019)
  integer(c_int), pointer, save :: itpcol

  !> \anchor iccvfg
  !> indicates whether the dynamic field should be frozen or not:
  !>    - 1: true
  !>    - 0: false (default)\n
  !> In such a case, the values of velocity, pressure and the
  !> variables related to the potential turbulence model
  !> (\f$k\f$, \f$R_{ij}\f$, \f$\varepsilon\f$, \f$\varphi\f$,
  !> \f$\bar{f}\f$, \f$\omega\f$, turbulent viscosity) are kept
  !> constant over time and only the equations for the scalars
  !> are solved.\n Also, if \ref iccvfg = 1, the physical properties
  !> modified in \ref cs_user_physical_properties will keep being
  !> updated. Beware of non-consistencies if these properties would
  !> normally affect the dynamic field (modification of density for
  !> instance).\n Useful if and only if \ref dimens::nscal "nscal"
  !> \f$>\f$ 0 and the calculation is a restart.
  integer(c_int), pointer, save :: iccvfg

  !> Algorithm to take into account the density variation in time
  !>    - 0: boussinesq algorithm with constant density
  !>    - 1: dilatable steady algorithm (default)
  !>    - 2: dilatable unsteady algorithm
  !>    - 3: low-Mach algorithm
  !>    - 4: algorithm for fire
  integer(c_int), pointer, save :: idilat

  !> accurate treatment of the wall temperature
  !>    - 1: true
  !>    - 0: false (default)
  !> (see \ref cs_boundary_condition_set_coeffs,
  !>  useful in case of coupling with syrthes)
  integer(c_int), pointer, save :: itbrrb

  !> Improved pressure interpolation scheme.
  !> See \ref cs_velocity_pressure_param_t::iphydr.

  integer(c_int), pointer, save :: iphydr

  !> compute the hydrostatic pressure in order to compute the Dirichlet
  !> conditions on the pressure at outlets
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: icalhy

  !> \}

  !----------------------------------------------------------------------------
  ! Homogeneous mixture modelling
  !----------------------------------------------------------------------------

  !> \addtogroup homogeneous_mixture
  !> \{
  !> \addtogroup vof
  !> \{

  !> VoF model (sum of masks defining VoF model and submodels).
  !> See defined masks in \ref vof_masks.
  integer(c_int), pointer, save :: ivofmt

  !> \}
  !> \}

  !----------------------------------------------------------------------------
  ! Additional source terms
  !----------------------------------------------------------------------------

  !> \defgroup additional_source_terms Additional source terms

  !> \addtogroup additional_source_terms
  !> \{

  !> Global head losses indicator (ie number of head loss zones)
  integer, save :: ncpdct = 0

  !> Indicateur termes sources de masse global (ie somme sur les processeurs
  !>   de ncetsm)
  integer, save :: nctsmt = 0

  !> Global indicator of condensation source terms (ie. sum on the processors
  !> of nfbpcd) cells associated to the face with condensation phenomenon
  integer, save :: nftcdt = 0

  !> \anchor iporos
  !> take the porosity fomulation into account
  !>    - 1: Taking porosity into account
  !>    - 0: Standard algorithm (Without porosity)
  integer(c_int), pointer, save :: iporos

  !> \}

  !----------------------------------------------------------------------------
  ! Transported scalars parameters
  !----------------------------------------------------------------------------

  !> \defgroup scalar_params Transported scalars parameters

  !> \addtogroup scalar_params
  !> \{

  !> flag for computing the drift mass flux:
  !> (for coal classes for instance, only the first
  !>  scalar of a class compute the drift flux of the class
  !>  and the other scalars use it without recomputing it)
  integer :: DRIFT_SCALAR_ADD_DRIFT_FLUX

  !> flag for activating thermophoresis for drift scalars
  integer :: DRIFT_SCALAR_THERMOPHORESIS

  !> flag for activating turbophoresis for drift scalars
  integer :: DRIFT_SCALAR_TURBOPHORESIS

  ! flag for activating electrophoresis for drift scalars
  integer :: DRIFT_SCALAR_ELECTROPHORESIS

  !> flag for activating the centrifugal force for drift scalars
  integer :: DRIFT_SCALAR_CENTRIFUGALFORCE

  !> flag for activating imposed mass flux
  integer :: DRIFT_SCALAR_IMPOSED_MASS_FLUX

  !> flag for seting the mass flux to zero at all boundaries
  integer :: DRIFT_SCALAR_ZERO_BNDY_FLUX

  !> flag for seting the mass flux to zero at walls only
  integer :: DRIFT_SCALAR_ZERO_BNDY_FLUX_AT_WALLS

  parameter (DRIFT_SCALAR_ADD_DRIFT_FLUX=1)
  parameter (DRIFT_SCALAR_THERMOPHORESIS=2)
  parameter (DRIFT_SCALAR_TURBOPHORESIS=3)
  parameter (DRIFT_SCALAR_ELECTROPHORESIS=4)
  parameter (DRIFT_SCALAR_CENTRIFUGALFORCE=5)
  parameter (DRIFT_SCALAR_IMPOSED_MASS_FLUX=6)
  parameter (DRIFT_SCALAR_ZERO_BNDY_FLUX=7)
  parameter (DRIFT_SCALAR_ZERO_BNDY_FLUX_AT_WALLS=8)

  !> flag for isotropic diffusion
  integer :: ISOTROPIC_DIFFUSION

  !> flag for orthotropic diffusion
  integer :: ORTHOTROPIC_DIFFUSION

  !> flag for diffusion by a left-multiplied symmetric 3x3 tensor
  integer :: ANISOTROPIC_LEFT_DIFFUSION

  ! flag for diffusion by a right-multiplied symmetric 3x3 tensor
  integer :: ANISOTROPIC_RIGHT_DIFFUSION

  !> flag for diffusion by a symmetric 3x3 tensor
  integer :: ANISOTROPIC_DIFFUSION

  parameter (ISOTROPIC_DIFFUSION=1)
  parameter (ORTHOTROPIC_DIFFUSION=2)
  parameter (ANISOTROPIC_LEFT_DIFFUSION=4)
  parameter (ANISOTROPIC_RIGHT_DIFFUSION=8)
  parameter (ANISOTROPIC_DIFFUSION=12)

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function defining whether time step is variable in
    ! time or not

    subroutine time_step_define_variable(is_variable)  &
      bind(C, name='cs_time_step_define_variable')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: is_variable
    end subroutine time_step_define_variable

    ! Interface to C function defining whether time step is local in
    ! space or not

    subroutine time_step_define_local(is_local)  &
      bind(C, name='cs_time_step_define_local')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), value :: is_local
    end subroutine time_step_define_local

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

    subroutine cs_f_turb_model_get_pointers(iturb, itytur, hybrid_turb) &
      bind(C, name='cs_f_turb_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iturb, itytur, hybrid_turb
    end subroutine cs_f_turb_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global wall functions structure

    subroutine cs_f_wall_functions_get_pointers(iwallf, iwalfs, &
                                                ypluli)         &
      bind(C, name='cs_f_wall_functions_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iwallf, iwalfs, ypluli
    end subroutine cs_f_wall_functions_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! RANS turbulence model structure

    subroutine cs_f_turb_rans_model_get_pointers(irccor, itycor, idirsm, &
                                                 iclkep, igrhok,  &
                                                 ikecou, reinit_turb, &
                                                 irijco, irijnu,  &
                                                 irijrb, irijec, idifre, &
                                                 iclsyr, iclptr)         &
      bind(C, name='cs_f_turb_rans_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: irccor, itycor, idirsm, iclkep, igrhok
      type(c_ptr), intent(out) :: ikecou, reinit_turb, irijco
      type(c_ptr), intent(out) :: irijnu, irijrb
      type(c_ptr), intent(out) :: irijec, idifre, iclsyr, iclptr
    end subroutine cs_f_turb_rans_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! LES turbulence model structure

    subroutine cs_f_turb_les_model_get_pointers(idries) &
      bind(C, name='cs_f_turb_les_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: idries
    end subroutine cs_f_turb_les_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! hybrid turbulence model structure

    subroutine cs_f_turb_hybrid_model_get_pointers(iicc, ishield) &
      bind(C, name='cs_f_turb_hybrid_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iicc, ishield
    end subroutine cs_f_turb_hybrid_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! LES balance structure
    subroutine cs_f_les_balance_get_pointer(i_les_balance) &
      bind(C, name='cs_f_les_balance_get_pointer')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: i_les_balance
    end subroutine cs_f_les_balance_get_pointer

    ! Interface to C function retrieving pointers to mesh quantity options

    subroutine cs_f_porous_model_get_pointers(iporos)  &
      bind(C, name='cs_f_porous_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iporos
    end subroutine cs_f_porous_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! velocity pressure model options structure

    subroutine cs_f_velocity_pressure_model_get_pointers  &
      (idilat)  &
      bind(C, name='cs_f_velocity_pressure_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: idilat
    end subroutine cs_f_velocity_pressure_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! velocity pressure parameters structure

    subroutine cs_f_velocity_pressure_param_get_pointers  &
      (iphydr, icalhy, itpcol)  &
      bind(C, name='cs_f_velocity_pressure_param_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iphydr, icalhy, itpcol
    end subroutine cs_f_velocity_pressure_param_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global spatial discretisation options structure

    subroutine cs_f_space_disc_get_pointers(imvisf, imrgra, iflxmw, itbrrb)   &
      bind(C, name='cs_f_space_disc_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: imvisf, imrgra, iflxmw, itbrrb
    end subroutine cs_f_space_disc_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global time schemeoptions structure

    subroutine cs_f_time_scheme_get_pointers(ischtp, istmpf, isno2t, isto2t, &
                                             iccvfg, initro) &
      bind(C, name='cs_f_time_scheme_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ischtp, istmpf, isno2t, isto2t
      type(c_ptr), intent(out) :: iccvfg, initro
    end subroutine cs_f_time_scheme_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global restart_auxiliary options structure

    subroutine cs_f_restart_auxiliary_get_pointers(ileaux)          &
      bind(C, name='cs_f_restart_auxiliary_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ileaux
    end subroutine cs_f_restart_auxiliary_get_pointers

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

    type(c_ptr) :: c_iturb, c_itytur, c_hybrid_turb

    call cs_f_turb_model_get_pointers(c_iturb, c_itytur, c_hybrid_turb)

    call c_f_pointer(c_iturb, iturb)
    call c_f_pointer(c_itytur, itytur)
    call c_f_pointer(c_hybrid_turb, hybrid_turb)

  end subroutine turb_model_init

  !> \brief Initialize Fortran wall functions API.
  !> This maps Fortran pointers to global C structure members.

  subroutine wall_functions_init

    use, intrinsic :: iso_c_binding
    use cstphy, only: ypluli
    implicit none

    ! Local variables

    type(c_ptr) :: c_iwallf, c_iwalfs, c_ypluli

    call cs_f_wall_functions_get_pointers(c_iwallf, c_iwalfs, &
                                          c_ypluli)

    call c_f_pointer(c_iwallf, iwallf)
    call c_f_pointer(c_iwalfs, iwalfs)
    call c_f_pointer(c_ypluli, ypluli)

  end subroutine wall_functions_init

  !> \brief Initialize Fortran RANS turbulence model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine turb_rans_model_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_irccor, c_itycor, c_idirsm, c_iclkep, c_igrhok
    type(c_ptr) :: c_ikecou, c_reinit_turb, c_irijco, c_irijnu
    type(c_ptr) :: c_irijrb, c_irijec, c_idifre
    type(c_ptr) :: c_iclsyr, c_iclptr

    call cs_f_turb_rans_model_get_pointers( c_irccor, c_itycor, c_idirsm, &
                                            c_iclkep, c_igrhok, &
                                            c_ikecou, c_reinit_turb, &
                                            c_irijco, c_irijnu, &
                                            c_irijrb, c_irijec, c_idifre, &
                                            c_iclsyr, c_iclptr)

    call c_f_pointer(c_irccor, irccor)
    call c_f_pointer(c_itycor, itycor)
    call c_f_pointer(c_idirsm, idirsm)
    call c_f_pointer(c_iclkep, iclkep)
    call c_f_pointer(c_igrhok, igrhok)
    call c_f_pointer(c_ikecou, ikecou)
    call c_f_pointer(c_reinit_turb, reinit_turb)
    call c_f_pointer(c_irijco, irijco)
    call c_f_pointer(c_irijnu, irijnu)
    call c_f_pointer(c_irijrb, irijrb)
    call c_f_pointer(c_irijec, irijec)
    call c_f_pointer(c_idifre, idifre)
    call c_f_pointer(c_iclsyr, iclsyr)
    call c_f_pointer(c_iclptr, iclptr)

  end subroutine turb_rans_model_init

  !> \brief Initialize Fortran LES turbulence model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine turb_les_model_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_idries, c_i_les_balance

    call cs_f_turb_les_model_get_pointers(c_idries)
    call cs_f_les_balance_get_pointer(c_i_les_balance)

    call c_f_pointer(c_idries, idries)
    call c_f_pointer(c_i_les_balance, i_les_balance)

  end subroutine turb_les_model_init

  !> \brief Initialize Fortran hybrid turbulence model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine turb_hybrid_model_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_iicc, c_ishield

    call cs_f_turb_hybrid_model_get_pointers(c_iicc, c_ishield)

    call c_f_pointer(c_iicc, iicc)
    call c_f_pointer(c_ishield, ishield)

  end subroutine turb_hybrid_model_init

  !> \brief Initialize Fortran Stokes options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine velocity_pressure_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_iporos
    type(c_ptr) :: c_itpcol, c_idilat, c_iphydr, c_icalhy

    call cs_f_porous_model_get_pointers(c_iporos)

    call c_f_pointer(c_iporos, iporos)

    call cs_f_velocity_pressure_model_get_pointers(c_idilat)

    call c_f_pointer(c_idilat, idilat)

    call cs_f_velocity_pressure_param_get_pointers  &
      (c_iphydr, c_icalhy, c_itpcol)

    call c_f_pointer(c_iphydr, iphydr)
    call c_f_pointer(c_icalhy, icalhy)
    call c_f_pointer(c_itpcol, itpcol)

  end subroutine velocity_pressure_options_init

  !> \brief Initialize Fortran space discretisation options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine space_disc_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_imvisf, c_imrgra, c_iflxmw, c_itbrrb

    call cs_f_space_disc_get_pointers(c_imvisf, c_imrgra, c_iflxmw, c_itbrrb)

    call c_f_pointer(c_imvisf, imvisf)
    call c_f_pointer(c_imrgra, imrgra)
    call c_f_pointer(c_iflxmw, iflxmw)
    call c_f_pointer(c_itbrrb, itbrrb)

  end subroutine space_disc_options_init

  !> \brief Initialize Fortran time scheme options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine time_scheme_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ischtp, c_istmpf, c_isno2t, c_isto2t
    type(c_ptr) :: c_iccvfg, c_initro

    call cs_f_time_scheme_get_pointers(c_ischtp, c_istmpf, c_isno2t, c_isto2t, &
                                       c_iccvfg, c_initro)

    call c_f_pointer(c_ischtp, ischtp)
    call c_f_pointer(c_istmpf, istmpf)
    call c_f_pointer(c_isno2t, isno2t)
    call c_f_pointer(c_isto2t, isto2t)
    call c_f_pointer(c_iccvfg, iccvfg)
    call c_f_pointer(c_initro, initro)

  end subroutine time_scheme_options_init

  !> \brief Initialize Fortran auxiliary options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine restart_auxiliary_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ileaux

    call cs_f_restart_auxiliary_get_pointers(c_ileaux)

    call c_f_pointer(c_ileaux, ileaux)

  end subroutine restart_auxiliary_options_init

  !=============================================================================

end module optcal
