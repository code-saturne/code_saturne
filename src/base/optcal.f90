!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

  !> time order of time stepping
  !>    - 2: 2nd order
  !>    - 1: 1st order (default)
  integer, save ::          ischtp

  !> time order of the mass flux scheme
  !> The chosen value for \ref istmpf will automatically
  !> determine the value given to the variable \ref thetfl.
  !> - 2: theta scheme with theta > 0 (theta=0.5 means 2nd order)
  !> the mass flow used in the momentum equations is extrapolated at
  !> n+ \ref thetfl (= n+1/2) from the values at the two former time
  !> steps (Adams Bashforth); the mass flow used in the equations for
  !> turbulence and scalars is interpolated at time n+ \ref thetfl
  !> (= n+1/2) from the values at the former time step and at the
  !> newly calculated \f$n+1\f$ time step.
  !> - 0: theta scheme with theta = 0 (explicit): the mass flow
  !> calculated at the previous time step is used in the convective
  !> terms of all the equations (momentum, turbulence and scalars)
  !> - 1: implicit scheme (default) : the mass flow calculated
  !> at the previous time step is used in the convective terms of the
  !> momentum equation, and the updated mass flow is used in the
  !> equations of turbulence and scalars. By default, \ref istmpf=2
  !> is used in the case of a second-order time scheme (if \ref ischtp=2)
  !> and \ref istmpf = 1 otherwise.
  integer, save ::          istmpf

  !> number of interations on the pressure-velocity coupling on Navier-Stokes
  !> (for the PISO algorithm)
  integer(c_int), pointer, save ::          nterup

  !> \ref isno2t specifies the time scheme activated for the source
  !> terms of the momentum equation, apart from convection and
  !> diffusion (for instance: head loss, transposed gradient, ...).
  !> - 0: "standard" first-order: the terms which are linear
  !> functions of the solved variable are implicit and the others
  !> are explicit
  !> - 1: second-order: the terms of the form \f$S_i\phi\f$ which are
  !> linear functions of the solved variable \f$\phi\f$ are expressed
  !> as second-order terms by interpolation (according to the formula
  !> \f$(S_i\phi)^{n+\theta}=S_i^n[(1-\theta)\phi^n+\theta\phi^{n+1}]\f$,
  !> \f$\theta\f$ being given by the value of \ref thetav associated
  !> with the variable \f$\phi\f$); the other terms \f$S_e\f$ are
  !> expressed as second-order terms by extrapolation (according to the
  !> formula \f$(S_e)^{n+\theta}=[(1+\theta)S_e^n-\theta S_e^{n-1}]\f$,
  !> \f$\theta\f$ being given by the value of \ref thetsn = 0.5).\n
  !> - 2: the linear terms \f$S_i\phi\f$ are treated in the same
  !> way as when \ref isno2t = 1; the other terms \f$S_e\f$ are
  !> extrapolated according to the same formula as when \ref isno2t = 1,
  !> but with \f$\theta\f$= \ref thetsn = 1. By default, \ref isno2t
  !> is initialised to 1 (second-order) when the selected time scheme
  !> is second-order (\ref ischtp = 2), otherwise to 0.
  integer, save ::          isno2t

  !> \ref isto2t specifies the time scheme activated for
  !> the source terms of the turbulence equations i.e. related
  !> to \f$k\f$, \f$R_{ij}\f$, \f$\varepsilon\f$, \f$\omega\f$, \f$\varphi\f$,
  !> \f$\overline{f}\f$), apart from convection and diffusion.
  !> - 0: standard first-order: the terms which are linear
  !> functions of the solved variable are implicit and the others are explicit
  !> - 1: second-order: the terms of the form \f$S_i\phi\f$ which are linear functions
  !> of the solved variable \f$\phi\f$ are expressed as second-order terms
  !> by interpolation (according to the formula
  !> \f$(S_i\phi)^{n+\theta}=S_i^n[(1-\theta)\phi^n+\theta\phi^{n+1}]\f$,
  !> \f$\theta\f$ being given by the value of \ref thetav associated with the
  !> variable \f$\phi\f$); the other terms \f$S_e\f$ are expressed as second-order
  !> terms by extrapolation (according to the formula
  !> \f$(S_e)^{n+\theta}=[(1+\theta)S_e^n-\theta S_e^{n-1}]\f$, \f$\theta\f$ being
  !> given by the value of \ref thetst = 0.5)
  !> - 2: the linear terms \f$S_i\phi\f$ are treated in the same
  !> way as when \ref isto2t = 1; the other terms \f$S_e\f$ are
  !> extrapolated according to the same formula as when \ref isto2t = 1,
  !> but with \f$\theta\f$= \ref thetst = 1.\n
  !> Due to certain specific couplings between the turbulence equations,
  !> \ref isto2t is allowed the value 1 or 2 only for the \f$R_{ij}\f$ models
  !> (\ref iturb = 30 or 31); hence, it is always initialised to 0.
  integer, save ::          isto2t

  !> for each scalar, \ref isso2t specifies the time scheme activated
  !> for the source terms of the equation for the scalar, apart from convection and
  !> diffusion (for instance: variance production, user-specified terms, ...).
  !> - 0: "standard" first-order: the terms which are linear
  !> functions of the solved variable are implicit and the others are explicit
  !> - 1: second-order: the terms of the form \f$S_i\phi\f$ which are
  !> linear functions of the solved variable \f$\phi\f$ are expressed
  !> as second-order terms by interpolation (according to the formula
  !> \f$(S_i\phi)^{n+\theta}=S_i^n[(1-\theta)\phi^n+\theta\phi^{n+1}]\f$, \f$\theta\f$
  !> being given by the value of \ref thetav associated with the variable \f$\phi\f$);
  !> the other terms \f$S_e\f$ are expressed as second-order terms by
  !> extrapolation (according to the formula
  !> \f$(S_e)^{n+\theta}=[(1+\theta)S_e^n-\theta S_e^{n-1}]\f$, \f$\theta\f$ being
  !> given by the value of \ref thetss (iscal) = 0.5)
  !> - 2: the linear terms \f$S_i\phi\f$ are treated in the same way as
  !> when \ref isso2t = 1; the other terms \f$S_e\f$ are extrapolated
  !> according to the same formula as when \ref isso2t = 1, but with
  !> \f$\theta\f$ = \ref thetss (iscal) = 1.\n
  !> By default, \ref isso2t (iscal) is initialised to 1 (second-order)
  !> when the selected time scheme is second-order (\ref ischtp = 2),
  !> otherwise to 0.
  integer, save ::          isso2t(nscamx)

  !> \ref iroext specifies the time scheme activated
  !> for the physical property \f$\phi\f$ density.
  !> - 0: "standard" first-order: the value calculated at
  !> the beginning of the current time step (from the
  !> variables known at the end of the previous time step) is used
  !> - 1: second-order: the physical property \f$\phi\f$ is
  !> extrapolated according to the formula
  !> \f$\phi^{n+\theta}=[(1+\theta)\phi^n-\theta \phi^{n-1}]\f$, \f$\theta\f$ being
  !> given by the value of \ref thetro = 0.5
  !> - 2: first-order: the physical property \f$\phi\f$ is
  !> extrapolated at $n+1$ according to the same formula
  !> as when \ref iroext = 1 but with \f$\theta\f$ = \ref thetro = 1
  integer, save ::          iroext

  !> \ref iviext specifies the time scheme activated
  !> for the physical property \f$\phi\f$ "total viscosity"
  !> (molecular+turbulent or sub-grid viscosities).
  !> - 0: "standard" first-order: the value calculated at
  !> the beginning of the current time step (from the
  !> variables known at the end of the previous time step) is used
  !> - 1: second-order: the physical property \f$\phi\f$ is
  !> extrapolated according to the formula
  !> \f$\phi^{n+\theta}=[(1+\theta)\phi^n-\theta \phi^{n-1}]\f$, \f$\theta\f$
  !> being given by the value of \ref thetvi = 0.5
  !> - 2: first-order: the physical property \f$\phi\f$ is
  !> extrapolated at \f$n+1\f$ according to the
  !> same formula as when \ref iviext = 1, but with \f$\theta\f$= \ref thetvi = 1
  integer, save ::          iviext

  !> \ref icpext specifies the time scheme activated
  !> for the physical property \f$\phi\f$ "specific heat".
  !> - 0: "standard" first-order: the value calculated at
  !> the beginning of the current time step (from the
  !> variables known at the end of the previous time step) is used
  !> - 1: second-order: the physical property \f$\phi\f$ is
  !> extrapolated according to the formula
  !> \f$\phi^{n+\theta}=[(1+\theta)\phi^n-\theta \phi^{n-1}]\f$, \f$\theta\f$
  !> being given by the value of \ref thetcp = 0.5
  !> - 2: first-order: the physical property \f$\phi\f$ is
  !> extrapolated at \f$n+1\f$ according to the
  !> same formula as when \ref icpext = 1, but with \f$\theta\f$ = \ref thetcp = 1
  integer, save ::          icpext

  !> for each scalar iscal, \ref ivsext (iscal) specifies the time scheme
  !> activated for the physical property \f$\phi\f$ "diffusivity".
  !> - 0: "standard" first-order: the value calculated at
  !> the beginning of the current time step (from the variables known
  !> at the end of the previous time step) is used
  !> - 1: second-order: the physical property \f$\phi\f$ is
  !> extrapolated according to the formula
  !> \f$\phi^{n+\theta}=[(1+\theta)\phi^n-\theta \phi^{n-1}]\f$, \f$\theta\f$
  !> being given by the value of \ref thetvs (iscal) = 0.5
  !> - 2: first-order: the physical property \f$\phi\f$ is
  !> extrapolated at $n+1$ according to the same formula as
  !> when \ref ivsext = 1, but with \f$\theta\f$ = \ref thetvs (iscal) = 1
  integer, save ::          ivsext(nscamx)

  !> initvi : =1 if total viscosity read from checkpoint file
  integer, save ::          initvi

  !> initro : =1 if density read from checkpoint file
  integer, save ::          initro

  !> initcp : =1 if specific heat read from checkpoint file
  integer, save ::          initcp

  !> initvs : =1 if scalar diffusivity read from checkpoint file
  integer, save ::          initvs(nscamx)

  !> \f$ \theta_S \f$-scheme for the source terms \f$S_e\f$ in the
  !> Navier-Stokes equations when the source term extrapolation has
  !> been activated (see \ref isno2t), following the formula
  !> \f$(S_e)^{n+\theta}=(1+\theta)S_e^n-\theta S_e^{n-1}\f$.\n The value
  !> of \f$theta\f$ = \ref thetsn is deduced from the value chosen for
  !> \ref isno2t. Generally only the value 0.5 is used.
  !>    -  0 : second viscosity explicit
  !>    - 1/2: second viscosity extrapolated in n+1/2
  !>    -  1 : second viscosity extrapolated in n+1
  double precision, save :: thetsn

  !> \f$ \theta \f$-scheme for the extrapolation of the nonlinear
  !> explicit source terms $S_e$ of the turbulence equations when the
  !> source term extrapolation has been activated (see \ref isto2t),
  !> following the formula \f$(S_e)^{n+\theta}=(1+\theta)S_e^n-\theta S_e^{n-1}\f$.\n
  !> The value of \f$theta\f$ is deduced from the value chosen for
  !> \ref isto2t. Generally, only the value 0.5 is used.
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetst

  !> \f$ \theta \f$-scheme for the extrapolation of the nonlinear
  !> explicit source term \f$S_e\f$ of the scalar transport equation
  !> when the source term extrapolation has been activated (see
  !> \ref isso2t), following the formula
  !> \f$(S_e)^{n+\theta}=(1+\theta)S_e^n-\theta S_e^{n-1}\f$.\n
  !> The value of \f$\theta\f$ = \ref thetss is deduced from the value
  !> chosen for \ref isso2t. Generally, only the value 0.5 is used.
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetss(nscamx)

  !> \f$ \theta \f$-scheme for the mass flux when a second-order
  !> time scheme has been activated for the mass flow (see \ref istmpf).
  !>    -  0 : explicit first-order (corresponds to \ref istmpf = 0 or 1)
  !>    - 1/2: extrapolated in n+1/2 (corresponds to \ref istmpf = 2). The mass
  !> flux will be interpolated according to the formula
  !> \f$Q^{n+\theta}=\frac{1}{2-\theta}Q^{n+1}+\frac{1-\theta}{2-\theta}Q^{n+1-\theta}\f$)
  !>    -  1 : extrapolated in n+1\n
  !> Generally, only the value 0.5 is used.
  double precision, save :: thetfl

  !> \f$ \theta \f$-scheme for the extrapolation of the physical
  !> property \f$\phi\f$ "total viscosity" when the extrapolation
  !> has been activated (see \ref iviext), according to the formula
  !> \f$\phi^{n+\theta}=(1+\theta)\phi^n-\theta \phi^{n-1}\f$.\n
  !> The value of \f$\theta\f$ = \ref thetvi is deduced from the value
  !> chosen for \ref iviext. Generally, only the value 0.5 is used.
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetvi

  !> \f$ \theta \f$-scheme for the extrapolation of the physical
  !> property \f$\phi\f$ "density" when the extrapolation has been
  !> activated (see \ref iroext), according to the formula
  !> \f$\phi^{n+\theta}=(1+\theta)\phi^n-\theta \phi^{n-1}\f$.\n
  !> The value of \f$\theta\f$ = \ref thetro is deduced from the value chosen
  !> for \ref iroext. Generally, only the value 0.5 is used.
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetro

  !> \f$ \theta \f$-scheme for the extrapolation of the physical
  !> property \f$\phi\f$ "specific heat" when the extrapolation
  !> has been activated (see \ref icpext), according to the
  !> formula \f$\phi^{n+\theta}=(1+\theta)\phi^n-\theta \phi^{n-1}\f$.\n
  !> The value of \f$\theta\f$ = \ref thetcp is deduced from the value chosen for
  !> \ref icpext. Generally, only the value 0.5 is used.
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetcp

  !> \f$ \theta \f$-scheme for the extrapolation of the physical
  !> property \f$\phi\f$ "diffusivity" when the extrapolation has
  !> been activated (see \ref ivsext), according to the formula
  !> \f$\phi^{n+\theta}=(1+\theta)\phi^n-\theta \phi^{n-1}\f$.\n
  !> The value of\f$\theta\f$ = \ref thetvs is deduced from the value
  !> chosen for \ref ivsext. Generally, only the value 0.5 is used.
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetvs(nscamx)

  !> relative precision for the convergence test of the iterative process on
  !> pressure-velocity coupling (PISO)
  real(c_double), pointer, save :: epsup

  !> norm  of the increment \f$ \vect{u}^{k+1} - \vect{u}^k \f$
  !> of the iterative process on pressure-velocity coupling (PISO)
  real(c_double), pointer, save :: xnrmu

  !> norm of \f$ \vect{u}^0 \f$ (used by PISO algorithm)
  real(c_double), pointer, save :: xnrmu0

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
  !>    - 4: iterative precess initialized by the least squares method
  integer(c_int), pointer, save :: imrgra

  !> non orthogonality angle of the faces, in radians.
  !> For larger angle values, cells with one node on the wall
  !> are kept in the extended support of the neighboring cells.
  real(c_double), pointer, save :: anomax

  !> \}

  !> \defgroup diffusive_scheme Diffusive scheme
  !> \addtogroup diffusive_scheme
  !> \{

  !> face viscosity field interpolation
  !>    - 1: harmonic
  !>    - 0: arithmetic (default)
  integer(c_int), pointer, save :: imvisf

  !> \}

  !> \defgroup linear_solver Linear solver
  !> \addtogroup linear_solver
  !> \{

  !> \anchor idircl
  !> indicates whether the diagonal of the matrix should be slightly
  !> shifted or not if there is no Dirichlet boundary condition and
  !> if \ref cs_var_cal_opt_t::istat "istat" = 0.
  !>    - 0: false
  !>    - 1: true
  !> Indeed, in such a case, the matrix for the general
  !> advection/diffusion equation is singular. A slight shift in the
  !> diagonal will make it invertible again.\n By default, \ref idircl
  !> is set to 1 for all the unknowns, except \f$\overline{f}\f$ in v2f
  !> modelling, since its equation contains another diagonal term
  !> that ensures the regularity of the matrix.
  !> \remark
  !> the code computes automatically for each variable the number of Dirichlet
  !> BCs
  integer, save ::          idircl(nvarmx)

  !> number of Dirichlet BCs
  integer, save ::          ndircl(nvarmx)

  !> \}

  !> \}

  !> Indicator of a calculation restart (=1) or not (=0).
  !> This value is set automatically by the code; depending on
  !> whether a restart directory is present, and should not be modified by
  !> the user
  integer, save :: isuite

  !> Indicates the reading (=1) or not (=0) of the auxiliary
  !> calculation restart file\n
  !> Useful only in the case of a calculation restart
  integer, save :: ileaux

  !> Indicates the writing (=1) or not (=0) of the auxiliary calculation
  !> restart file.
  integer, save :: iecaux

  !> \anchor isuit1
  !> For the 1D wall thermal module, activation (1) or not(0)
  !> of the reading of the mesh and of the wall temperature
  !> from the restart file
  !> Useful if nfpt1d > 0
  integer, save :: isuit1

  !> For the vortex method, indicates whether the synthetic
  !> vortices at the inlet should be initialised or read
  !> from the restart file.
  !> Useful if \ref iturb = 40, 41, 42 and \ref ivrtex = 1
  !> - 0: initialized
  !> - 1: read
  integer, save :: isuivo

  !> Reading of the LES inflow module restart file.
  !> -0: not activated
  !> -1: activated\n
  !> If \ref isuisy = 1, synthetic fluctuations are
  !> not re-initialized in case of restart calculation.
  !> Useful if \ref iturb = 40, 41 or 42
  integer, save :: isuisy

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

  !> indicator "zero time step"
  !>    - 0: standard calculation
  !>    - 1: to simulate no time step
  !>         - for non-restarted computations:
  !>           only resolution (Navier-Stokes, turbulence, scalars) is skipped
  !>         - for restarted computations:
  !>           resolution, computation of physical properties, and definition
  !>           of boundary conditions is skipped (values are read from
  !>           checkpoint file)
  integer(c_int), pointer, save :: inpdt0

  !> Clip the time step with respect to the buoyant effects
  !>
  !> When density gradients and gravity are present, a local thermal time
  !> step can be calculated, based on the Brunt-Vaisala frequency. In
  !> numerical simulations, it is usually wise for the time step to be
  !> lower than this limit, otherwise numerical instabilities may appear.\n
  !> \ref iptlro indicates whether the time step should be limited to the
  !> local thermal time step (=1) or not (=0).\n
  !> When \ref iptlro=1, the listing shows the number of cells where the
  !> time step has been clipped due to the thermal criterion, as well as
  !> the maximum ratio between the time step and the maximum thermal time
  !> step. If \ref idtvar=0, since the time step is fixed and cannot be
  !> clipped, this ratio can be greater than 1. When \ref idtvar > 0, this
  !> ratio will be less than 1, except if the constraint \ref dtmin has
  !> prevented the code from reaching a sufficiently low value for \ref dt.
  !> Useful when density gradients and gravity are present.
  integer(c_int), pointer, save :: iptlro

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

  !> maximum Courant number (when \ref idtvar is different from 0)
  real(c_double), pointer, save :: coumax

  !> maximum Courant number for the continuity equation in compressible model
  real(c_double), pointer, save :: cflmmx

  !> maximum Fourier number (when \ref idtvar is different from 0)
  real(c_double), pointer, save :: foumax

  !> maximum allowed relative increase in the calculated time step value
  !> between two successive time steps (to ensure stability, any decrease
  !> in the time step is immediate and without limit).\n
  !> Useful when \ref idtvar is different from 0.
  real(c_double), pointer, save :: varrdt

  !> lower limit for the calculated time step when idtvar is different from 0.\n
  !> Take \ref dtmin = min (ld/ud, sqrt(lt/(gdelta rho/rho)), ...)
  real(c_double), pointer, save :: dtmin

  !> upper limit for the calculated time step when idtvar is different from 0.\n
  !> Take \ref dtmax = max (ld/ud, sqrt(lt/(gdelta rho/rho)), ...)
  real(c_double), pointer, save :: dtmax

  !> multiplicator coefficient for the time step of each variable
  !>    - useless for u,v,w,p
  !>    - for k,e     the same value is taken (value of k)
  !>    - for Rij, e  the same value is taken (value of r11)\n
  !> Hence, the time step used when solving the evolution equation for
  !> the variable is the time step used for the dynamic equations (velocity/pressure)
  !> multiplied by \ref cdtvar.
  !> The size of the array \ref cdtvar is \ref dimens::nvar "nvar". For instance, the
  !> multiplicative coefficient applied to the scalar 2 is cdtvar(isca(2))). Yet, the
  !> value of cdtvar for the velocity components and the pressure is not used. Also,
  !> although it is possible to change the value of \ref cdtvar for the turbulent
  !> variables, it is highly not recommended.
  double precision, save :: cdtvar(nvarmx)

  !> relaxation coefficient for the steady algorithm
  !> \ref relxst = 1 : no relaxation.
  real(c_double), pointer, save :: relxst

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
  !> thermal variable (the correspondence and conversion is then specified
  !> by the user in the subroutine \ref usthht). However this case has never
  !> been used in practice and has therefore not been tested. With the
  !> compressible model, it is possible to carry out calculations coupled with
  !> SYRTHES, although the thermal scalar represents the total energy and not
  !> the temperature.
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
  integer(c_int), pointer, save :: iscalt

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

  !> Activation of rotation/curvature correction for eddy viscosity turbulence models
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
  !>  Indicates the type of wall function is used for the velocity
  !>  boundary conditions on a frictional wall.
  !>  - 0: no wall functions
  !>  - 1: one scale of friction velocities (power law)
  !>  - 2: one scale of friction velocities (log law)
  !>  - 3: two scales of friction velocities (log law)
  !>  - 4: two scales of friction velocities (log law) (scalable wall functions)
  !>  - 5: two scales of friction velocities (mixing length based on V. Driest analysis)\n
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
  !>  necessary in order to be always in a logarithmic layer. It is used to make up for
  !>  the problems related to the use of High-Reynolds models on very refined
  !>  meshes.\n
  !>  Useful if \ref iturb is different from 50.
  integer(c_int), pointer, save :: iwallf

  !>  Wall functions for scalar
  !>    - 0: three layer wall function of Arpaci and Larsen
  !>    - 1: Van Driest wall function
  integer(c_int), pointer, save :: iwalfs

  !> exchange coefficient correlation
  !>    - 0: not use by default
  !>    - 1: exchange coefficient computed with a correlation
  integer(c_int), pointer, save :: iwallt

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

  !> Indicates if the terms related to gravity are taken
  !> into account in the equations of \f$k-\epsilon\f$.
  !> - 1: true (default if \f$ \rho \f$ is variable)
  !> - 0: false
  !> Useful if and only if \ref iturb = 20, 21, 50 or 60 and
  !> (\ref cstphy::gx "gx", \ref cstphy::gy "gy", \ref cstphy::gz "gz"})
  !> \f$\ne\f$ (0,0,0) and the density is not uniform.
  integer(c_int), pointer, save :: igrake

  !> Indicates if the terms related to gravity are taken
  !> into account in the equations of \f$R_{ij}-\varepsilon\f$.
  !> - 1: true (default if \f$ \rho \f$ is variable)
  !> - 0: false
  !> Useful if and only if \ref iturb = 30 or 31 and (\ref cstphy::gx "gx",
  !> \ref cstphy::gy "gy", \ref cstphy::gz "gz"}) \f$\ne\f$
  !> (0,0,0) (\f$R_{ij}-\epsilon\f$ model with gravity) and the
  !> density is not uniform.
  integer(c_int), pointer, save :: igrari

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

  !> Activation of Hybrid DDES model (only valid for iturb equal to 60)
  integer(c_int), pointer, save :: iddes

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
  integer(c_int), pointer, save :: irijnu

  !> accurate treatment of \f$ \tens{R} \f$ at the boundary (see \ref condli)
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
  !> studied with Code_Saturne, they did not bring any improvement in the results.\n
  !> In addition, their use induces an increase in the calculation time.\n
  !> The wall echo terms imply the calculation of the distance to the wall
  !> for every cell in the domain. See \ref icdpar for potential restrictions
  !> due to this.
  integer(c_int), pointer, save :: irijec

  !> whole treatment of the diagonal part of the dissusion tensor of
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
  !> in the domain. Refer to keyword \ref icdpar
  !> for potential limitations.\n
  !> Useful if and only if \ref iturb = 40 or 41
  integer(c_int), pointer, save :: idries

  !> Activates or not the generation of synthetic turbulence at the
  !> different inlet boundaries with the LES model (generation of
  !> unsteady synthetic eddies).\n
  !> - 1: true
  !> - 0: false (default)
  !> Useful if \ref iturb =40, 41 or 42\n
  !> This keyword requires the completion of the routine  \ref usvort
  integer(c_int), pointer, save :: ivrtex

  !> Wall boundary condition on omega in k-omega SST
  !> 0: Deprecated Neumann boundary condition
  !> 1: Dirichlet boundary condition consistent with Menter's
  !>    original model: w_wall = 60*nu/(beta*d**2)
  integer, save :: ikwcln

  !> turbulent flux model for \f$ \overline{\varia^\prime \vect{u}^\prime} \f$
  !> for any scalar \f$ \varia \f$, iturt(isca)
  !>    - 0: SGDH
  !>    - 10: GGDH
  !>    - 11: EB-GGDH (Elliptic Blending)
  !>    - 20: AFM
  !>    - 21: EB-AFM (Elliptic Blending)
  !>    - 30: DFM (Transport equation modelized)
  !>    - 31: EB-DFM (Elliptic Blending)
  !> GGDH, AFM and DFM are only available when a second order closure is used.
  integer, save :: iturt(nscamx)

  !> class turbulent flux model (=iturt/10)
  integer, save :: ityturt(nscamx)

  !> number of variable (deprecated, used only for compatibility)
  integer, save :: nvarcl

  !> \}

  !----------------------------------------------------------------------------
  ! Stokes
  !----------------------------------------------------------------------------

  !> \defgroup stokes Stokes options

  !> \addtogroup stokes
  !> \{

  !> Indicates whether the source terms in transposed gradient
  !> and velocity divergence should be taken into account in the
  !> momentum equation. In the compressible module, these terms
  !> also account for the volume viscosity (cf. \ref ppincl::viscv0 "viscv0"
  !> and \ref ppincl::iviscv "iviscv")
  !> \f$\partial_i \left[(\kappa -2/3\,(\mu+\mu_t))\partial_k U_k  \right]
  !> +     \partial_j \left[ (\mu+\mu_t)\partial_i U_j \right]\f$:
  !> - 0: not taken into account,
  !> - 1: taken into account.
  integer(c_int), pointer, save :: ivisse

  !> Reconstruction of the velocity field with the updated pressure option
  !>    - 0: default
  integer(c_int), pointer, save ::          irevmc

  !> Compute the pressure step thanks to the continuity equation
  !>    - 1: true (default)
  !>    - 0: false
  integer(c_int), pointer, save ::          iprco

  !> Arakawa multiplicator for the Rhie and Chow filter (1 by default)
  real(c_double), pointer, save :: arak

  !> indicates the algorithm for velocity/pressure coupling:
  !> - 0: standard algorithm,
  !> - 1: reinforced coupling in case calculation with long time steps\n
  !> Always useful (it is seldom advised, but it can prove very useful,
  !> for instance, in case of flows with weak convection effects and
  !> highly variable viscosity).
  integer(c_int), pointer, save :: ipucou

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
  !>    - 1: dilatable steady algorithm (default)
  !>    - 2: dilatable unsteady algorithm
  !>    - 3: low-Mach algorithm
  !>    - 4: algorithm for fire
  !    - 0: boussinesq algorithm with constant density
  integer(c_int), pointer, save :: idilat

  !> Option to switch on massflux predcition befor momentum solving
  !> to be fully conservative in momentum over time for variable density flows.
  !> This option is to be removed.
  integer, save :: ipredfl

  !> parameter of diagonal pressure strengthening
  real(c_double), pointer, save :: epsdp

  !TODO doxygen
  ! Type des conditions limites et index min et max
  !                 des sous listes defaces de bord
  integer, save :: idebty(ntypmx), ifinty(ntypmx)

  !> accurate treatment of the wall temperature
  !>    - 1: true
  !>    - 0: false (default)
  !> (see \ref condli, useful in case of coupling with syrthes)
  integer(c_int), pointer, save :: itbrrb

  !> improve static pressure algorithm
  !>    - 1: impose the equilibrium of the static part of the pressure
  !>         with any external force, even head losses
  !>    - 2: compute an hydrostatic pressure due to buoyancy forces before
  !>         the prediction step
  !>    - 0: no treatment (default)
  !>        When the density effects are important, the choice of \ref iphydr = 1
  !>        allows to improve the interpolation of the pressure and correct the
  !>        non-physical velocities which may appear in highly stratified areas
  !>        or near horizontal walls (thus avoiding the use of
  !>        \ref cs_var_cal_opt_t::extrag "extrag"
  !>        if the non-physical velocities are due only to gravity effects).\n
  !>        The improved algorithm also allows eradicating the velocity oscillations
  !>        which tend to appear at the frontiers of areas with high head losses.\n
  !>        In the case of a stratified flow, the calculation cost is higher when
  !>        the improved algorithm is used (about 30\% depending on the case)
  !>        because the hydrostatic pressure must be recalculated at the outlet
  !>        boundary conditions: see \ref icalhy.\n
  !>        On meshes of insufficient quality, in order to
  !>        improve the convergence, it may be useful to increase the number of
  !>        iterations for the reconstruction of the pressure right-hand side,
  !>        i.e. \ref cs_var_cal_opt_t::nswrsm "nswrsm".\n If head losses are
  !>        present just along an outlet boundary, it is necessary to specify
  !>        \ref icalhy = 0 in order to deactivate the recalculation of the
  !>        hydrostatic pressure at the boundary, which may otherwise cause
  !>        instabilities. Please refer to the
  !>    <a href="../../theory.pdf#iphydr"><b>handling of the hydrostatic pressure</b></a>
  !>      section of the theory guide for more informations.
  integer(c_int), pointer, save :: iphydr

  !> improve static pressure algorithm
  !>    - 1: take -div(rho R) in the static pressure
  !>      treatment IF iphydr=1
  !>    - 0: no treatment (default)
  integer(c_int), pointer, save :: igprij

  !> improve static pressure algorithm
  !>    - 1: take user source term in the static pressure
  !>      treatment IF iphydr=1 (default)
  !>    - 0: no treatment
  integer(c_int), pointer, save :: igpust

  !> indicates the presence of a Bernoulli boundary face (automatically computed)
  !>    - 0: no face
  !>    - 1: at least one face
  integer(c_int), pointer, save :: iifren

  !> number of the closest free standard outlet (or free inlet) face to xyzp0
  integer, save :: ifrslb

  !> max of ifrslb on all ranks, standard outlet face presence indicator
  integer, save :: itbslb

  !> compute the hydrostatic pressure in order to compute the Dirichlet
  !> conditions on the pressure at outlets
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: icalhy

  !> use interpolated face diffusion coefficient instead of cell diffusion coefficient
  !> for the mass flux reconstruction for the non-orthogonalities
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: irecmf

  !> choice the way to compute the exchange coefficient of the
  !> condensation source term used by the copain model
  !>    - 1: the turbulent exchange coefficient of the flow
  !>    - 2: the exchange coefficient of the copain correlation
  !>    - 3: the maximal value between the two previous exchange coefficients
  integer, save :: icophc

  !> choice the way to compute the thermal exchange coefficient associated
  !> to the heat transfer to wall due to the condensation phenomenon
  !>    - 2: the thermal exchange coefficient of the copain correlation
  !>    - 3: the maximal value between the current and previous thermal
  !>         exchange coefficient evaluated by the copain correlation
  integer, save :: icophg

  !> choice the way to compute the wall temperature at the solid/fluid interface
  !> coupled with condensation to the wall
  !>    - 1: the wall temperature is computed with a 1-D thermal model
  !>         with implicit numerical scheme
  !>    - 0: the wall temperature is imposed as constant by the user (default)
  !>         exchange coefficient evaluated by the copain correlation
  integer, save :: itag1d

  !> choice the way to compute the wall temperature at the solid/fluid interface
  !> coupled with condensation to the metal mass structures wall
  !>    - 1: the wall temperature is computed with a 0-D thermal model
  !>         with explicit numerical scheme
  !>    - 0: the wall temperature is imposed as constant by the user (default)
  !>         and past to the copain correlation to evaluate the exchange coefficient
  integer, save :: itagms


  !> \ref iescal indicates the calculation mode for the error estimator
  !> \ref paramx::iespre "iespre", \ref paramx::iesder "iesder",
  !> \ref paramx::iescor "iescor" or \ref paramx::iestot "iestot"
  !> for the Navier-Stokes equation:
  !> - 0: estimator not calculated,
  !> - 1: the estimator  \f$ \eta^{*}_{i,1}\f$ is calculated,
  !> without contribution of the volume,
  !> - 2: the estimator \f$ \eta^{*}_{i,2}\f$ is calculated,
  !> with contribution of the volume (norm \f$L^2\f$),
  !> except for \ref paramx::iescor "iescor", for which
  !> \f$|\Omega_i|\ \eta^{corr}_{i,1}\ \f$
  !> is calculated. The names of the estimators appearing
  !> in the listing and the post-processing are made up of
  !> the default name (given before), followed by the value of
  !> \ref iescal}. For instance, EsPre2 is the estimator
  !> \ref paramx::iespre "iespre" calculated with \ref iescal = 2.
  integer, save :: iescal(nestmx)

  !> \ref n_buoyant_scal is the number of buoyant scalar
  !> It will be zero if there is no buoyant scalar
  integer(c_int), pointer, save :: n_buoyant_scal

  !> \}

  !----------------------------------------------------------------------------
  ! Homogeneous mixture modelling
  !----------------------------------------------------------------------------

  !> \defgroup homegeneous_mixture Homogeneous mixture modelling

  !> \addtogroup homegeneous_mixture
  !> \{

  !> cavitation module
  !>    - -1: module not activated
  !>    -  0: no vaporization/condensation model
  !>    -  1: Merkle's model
  integer, save :: icavit

  !> VOF method
  !>    - -1: module not activated
  !>    -  0: method activated
  integer, save :: ivofmt

  !> \}

  !----------------------------------------------------------------------------
  ! Additional source terms
  !----------------------------------------------------------------------------

  !> \defgroup additional_source_terms Additional source terms

  !> \addtogroup additional_source_terms
  !> \{

  !> Indicateur pertes de charge global (ie somme sur les processeurs
  !>   de ncepdc)
  integer, save :: ncpdct

  !> Indicateur termes sources de masse global (ie somme sur les processeurs
  !>   de ncetsm)
  integer, save :: nctsmt

  !> Global indicator of condensation source terms (ie. sum on the processors
  !> of nfbpcd) cells associated to the face with condensation phenomenon
  integer, save :: nftcdt

  !> \anchor iporos
  !> take the porosity fomulation into account
  !>    - 1: Taking porosity into account
  !>    - 0: Standard algorithm (Without porosity)
  integer, save :: iporos

  !TODO move it elsewhere?
  ! Indicateur de passage dans l'initialisation des
  !                       variables par l'utilisateur
  !          iusini = 1 passage dans usiniv ou ppiniv
  !                   0 pas de passage (ni iusini ni ppiniv)

  integer, save :: iusini

  !> \}

  !----------------------------------------------------------------------------
  ! Numerical parameters for the wall distance calculation
  !----------------------------------------------------------------------------

  !> \defgroup num_wall_distance Numerical parameters for the wall distance calculation

  !> \addtogroup num_wall_distance
  !> \{

  !> - 1, the wall distance must be computed,
  !> - 0, the wall distance computation is not necessary.
  integer, save :: ineedy

  !> - 1, the wall distance is up to date,
  !> - 0, the wall distance has not been updated.
  integer, save :: imajdy

  !> Specifies the method used to calculate the distance to the wall y
  !> and the non-dimensional distance \f$ y+ \f$ for all the cells of
  !> the calculation domain (when necessary):
  !> - 1: standard algorithm (based on a Poisson equation for y and
  !> convection equation for \f$ y+ \f$), with reading of the distance
  !> to the wall from the restart file if possible
  !> - -1: standard algorithm (based on a Poisson equation for y and
  !> convection equation for \f$ y+ \f$ ), with systematic recalculation
  !> of the distance to the wall in case of calculation restart
  !> - 2: former algorithm (based on geometrical considerations), with
  !> reading of the distance to the wall from the restart file if possible\n
  !> - -2: former algorithm (based on geometrical considerations) with systematic
  !> recalculation of the distance to the wall in case of calculation restart.\n\n
  !> In case of restart calculation, if the position of the walls haven’t changed,
  !> reading the distance to the wall from the restart file can save a fair amount
  !> of CPU time.\n Useful in \f$ R_{ij}-\epsilon \f$ model with wall echo
  !> (\ref iturb=30 and \ref irijec=1), in LES with van Driest damping
  !> (\ref iturb=40 and \ref idries=1) and in \f$ k-\omega\f$ SST (\ref iturb=60).
  !> By default, \ref icdpar is initialised to -1, in case there has been a change
  !> in the definition of the boundary conditions between two computations (change
  !> in the number or the positions of the walls). Yet, with the \f$k-\omega\f$ SST model,
  !> the distance to the wall is needed to calculate the turbulent viscosity, which is
  !> done before the calculation of the distance to the wall. Hence, when this model
  !> is used (and only in that case), \ref icdpar is set to 1 by default, to ensure
  !> total continuity of the calculation at restart. As a consequence, with the
  !> \f$k-\omega\f$ SST model, if the number and positions of the walls are changed
  !> at a calculation restart, it is mandatory for the user to set \ref icdpar
  !> explicitly to -1, otherwise the distance to the wall used will not correspond
  !> to the actual position of the walls.\n The former algorithm is not compatible
  !> with parallelism nor periodicity. Also, whatever the value chosen for \ref icdpar,
  !> the calculation of the distance to the wall is made at the most once for all at the
  !> beginning of the calculation; it is therefore not compatible with moving walls.
  !> Please contact the development team if you need to override this limitation.
  integer, save :: icdpar

  !> to the wall \f$ y+ \f$.\n
  !> useful when \ref icdpar \f$\neq\f$ 0 for the calculation of \f$ y+ \f$.
  integer, save :: ntcmxy

  !> Target Courant number for the calculation of the non-dimensional distance
  !> to the wall\n
  !> useful when \ref icdpar \f$\neq\f$ 0 for the calculation of \f$ y+ \f$.
  double precision, save :: coumxy

  !> relative precision for the convergence of the pseudo-transient regime
  !> for the calculation of the non-dimensional distance to the wall \n
  !> useful when \ref icdpar \f$\neq\f$ 0 for the calculation of \f$ y+ \f$.
  double precision, save :: epscvy

  !> value of the non-dimensional distance to the wall above which the
  !> calculation of the distance is not necessary (for the damping)
  !> useful when \ref icdpar \f$\neq\f$ 0 for the calculation of \f$ y+ \f$.
  double precision, save :: yplmxy

  !> \}

  !----------------------------------------------------------------------------
  ! Transported scalars parameters
  !----------------------------------------------------------------------------

  !> \defgroup scalar_params Transported scalars parameters

  !> \addtogroup scalar_params
  !> \{

  !> \anchor iscacp
  !> iscacp : 0 : scalar does not behave like a temperature
  !>          1 : scalar behaves like a temperature (use Cp for wall law)
  !>        > 1 : not yet allowed, could be used for multiple Cp definitions
  integer, save ::          iscacp(nscamx)

  !> iclvfl : 0 : clip variances to zero
  !>          1 : clip variances to zero and to f(1-f)
  !>          2 : clip variances to  max(zero,scamin) and scamax
  !> for every scalar iscal representing the average of the square of the
  !> fluctuations of another scalar ii= \ref iscavr (iscal) (noted \$f\$),
  !> indicator of the clipping method:
  !> - -1: no clipping because the scalar does not represent
  !> the average of the square of the fluctuations of another scalar
  !> - 0: clipping to 0 for the lower range of values
  !> - 1: clipping to 0 for the lower range of values and to
  !> \f$(f-f_{min})(f_{max}-f)\f$ for higher values, where \f$f\f$ is
  !> the associated scalar, \f$f_{min}\f$ and \f$f_{max}\f$ its minimum and maximum
  !> values specified by the user (i.e. scamin (ii) and scamax (ii))
  !> - 2: clipping to max(0,scamin(iscal)) for lower values and to
  !>  scamax(iscal) for higher values.scamin and scamax are limits
  !> specified by the user.\n Useful for the scalars iscal for
  !> which \ref iscavr (iscal) \f$>\f$0.
  integer, save ::          iclvfl(nscamx)

  !> iscasp(ii) : index of the ii^th species (0 if not a species)
  integer, save ::          iscasp(nscamx)

  !> reference molecular diffusivity related to the scalar J (\f$kg.m^{-1}.s^{-1}\f$).\n
  !>
  !> Negative value: not initialised\n
  !> Useful if 1\f$\leqslant\f$J\f$\leqslant\f$ \ref dimens::nscal "nscal",
  !> unless the user specifies the molecular diffusivity in the appropriate
  !> user subroutine (\ref cs_user_physical_properties for the standard
  !> physics) (field_get_key_id (ivarfl(isca(iscal)),kivisl,...)
  !> \f$>\f$ -1)\n Warning: \ref visls0 corresponds to the diffusivity.
  !> For the temperature, it is therefore defined as \f$\lambda/C_p\f$
  !> where \f$\lambda\f$ and \f$C_p\f$ are the conductivity and specific
  !> heat. When using the Graphical Interface, \f$\lambda\f$ and \f$C_p\f$
  !> are specified separately, and \ref visls0 is calculated automatically.\n
  !> With the compressible module, \ref visls0 (given in \ref uscfx2) is
  !> directly the thermal conductivity \f$W.m^{-1}.K^{-1}\f$.\n With gas or
  !> coal combustion, the molecular diffusivity of the enthalpy
  !> (\f$kg.m^{-1}.s^{-1}\f$) must be specified by the user in the variable
  !> \ref ppthch::diftl0 "diftl0"(\ref cs_user_combustion).\n
  !> With the electric module, for the Joule effect, the diffusivity is
  !> specified by the user in \ref cs_user_physical_properties.c (even if
  !> it is constant). For the electric arcs, it is calculated from the
  !> thermochemical data file.
  double precision, save :: visls0(nscamx)

  !> When iscavr(iscal)>0, \ref rvarfl is the coefficient \f$R_f\f$ in
  !> the dissipation term \f$\-\frac{\rho}{R_f}\frac{\varepsilon}{k}\f$
  !> of the equation concerning the scalar,
  !> which represents the root mean square of the fluctuations
  !> of the scalar.\n
  !> Useful if and only if there is 1\f$\leqslant\f$ iscal \f$\leqslant\f$
  !> \ref dimens::nscal "nscal" such as iscavr(iscal)>0
  double precision, save :: rvarfl(nscamx)

  !> ctheta : coefficient des modeles de flux turbulents GGDH et AFM
  double precision, save :: ctheta(nscamx)

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

  parameter (DRIFT_SCALAR_ADD_DRIFT_FLUX=1)
  parameter (DRIFT_SCALAR_THERMOPHORESIS=2)
  parameter (DRIFT_SCALAR_TURBOPHORESIS=3)
  parameter (DRIFT_SCALAR_ELECTROPHORESIS=4)
  parameter (DRIFT_SCALAR_CENTRIFUGALFORCE=5)
  parameter (DRIFT_SCALAR_IMPOSED_MASS_FLUX=6)

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

  !----------------------------------------------------------------------------
  ! electric model parameters
  !----------------------------------------------------------------------------

  !> \defgroup electric_model parameters

  !> \addtogroup electric_model_params
  !> \{

  !> ngazge  : number of species for electric arc
  integer(c_int), pointer, save :: ngazge

  !> ielcor : 0 : electric arc scaling desactivate
  !>          1 : electric arc scaling activate
  integer(c_int), pointer, save :: ielcor

  !> pot_diff : potential between electrods
  real(c_double), pointer, save :: pot_diff

  !> coejou : scaling coefficient
  real(c_double), pointer, save :: coejou

  !> elcou : current
  real(c_double), pointer, save :: elcou

  !> pot_diff : imposed value for current
  real(c_double), pointer, save :: couimp

  !> irestrike : 0 : restrike mode off
  !>             1 : restrike mode on
  integer(c_int), pointer, save :: irestrike

  !> restrike_point : coordinate of restrike point
  real(c_double), pointer, save :: restrike_point_x
  real(c_double), pointer, save :: restrike_point_y
  real(c_double), pointer, save :: restrike_point_z

  !> ntdcla : start iteration for restrike
  integer(c_int), pointer, save :: ntdcla

  !> \}

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
                                           t_prev, t_cur, t_max)             &
      bind(C, name='cs_f_time_step_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: nt_prev, nt_cur, nt_max, nt_ini
      type(c_ptr), intent(out) :: t_prev, t_cur, t_max
    end subroutine cs_f_time_step_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global time step options structure

    subroutine cs_f_time_step_options_get_pointers(inpdt0, iptlro, idtvar, &
                                                   dtref, coumax, cflmmx,  &
                                                   foumax, varrdt, dtmin,  &
                                                   dtmax, relxst)          &
      bind(C, name='cs_f_time_step_options_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: inpdt0, iptlro, idtvar, dtref, coumax, cflmmx
      type(c_ptr), intent(out) :: foumax, varrdt, dtmin, dtmax, relxst
    end subroutine cs_f_time_step_options_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global thermal model structure

    subroutine cs_f_thermal_model_get_pointers(itherm, itpscl, iscalt) &
      bind(C, name='cs_f_thermal_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: itherm, itpscl, iscalt
    end subroutine cs_f_thermal_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global turbulence model structure

    subroutine cs_f_turb_model_get_pointers(iturb, itytur) &
      bind(C, name='cs_f_turb_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iturb, itytur
    end subroutine cs_f_turb_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global wall functions structure

    subroutine cs_f_wall_functions_get_pointers(iwallf, iwalfs, iwallt, &
                                                ypluli)                  &
      bind(C, name='cs_f_wall_functions_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iwallf, iwalfs, iwallt, ypluli
    end subroutine cs_f_wall_functions_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! RANS turbulence model structure

    subroutine cs_f_turb_rans_model_get_pointers(irccor, itycor, idirsm, &
                                                 iclkep, igrhok, igrake, &
                                                 igrari, ikecou, reinit_turb, &
                                                 irijco, iddes, irijnu,  &
                                                 irijrb, irijec, idifre, &
                                                 iclsyr, iclptr)         &
      bind(C, name='cs_f_turb_rans_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: irccor, itycor, idirsm, iclkep, igrhok
      type(c_ptr), intent(out) :: igrake, igrari, ikecou, reinit_turb, irijco, irijnu, irijrb, iddes
      type(c_ptr), intent(out) :: irijec, idifre, iclsyr, iclptr
    end subroutine cs_f_turb_rans_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! LES turbulence model structure

    subroutine cs_f_turb_les_model_get_pointers(idries, ivrtex) &
      bind(C, name='cs_f_turb_les_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: idries, ivrtex
    end subroutine cs_f_turb_les_model_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! Stokes options structure

    subroutine cs_f_stokes_options_get_pointers(ivisse, irevmc, iprco,         &
                                                arak  ,ipucou, iccvfg,         &
                                                idilat, epsdp ,itbrrb, iphydr, &
                                                igprij, igpust,                &
                                                iifren, icalhy, irecmf)        &
      bind(C, name='cs_f_stokes_options_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ivisse, irevmc, iprco, arak
      type(c_ptr), intent(out) :: ipucou, iccvfg, idilat, epsdp, itbrrb, iphydr
      type(c_ptr), intent(out) :: igprij, igpust, iifren, icalhy, irecmf
    end subroutine cs_f_stokes_options_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global spatial discretisation options structure

    subroutine cs_f_space_disc_get_pointers(imvisf, imrgra, anomax, iflxmw) &
      bind(C, name='cs_f_space_disc_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: imvisf, imrgra, anomax, iflxmw
    end subroutine cs_f_space_disc_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global PISO options structure

    subroutine cs_f_piso_get_pointers(nterup, epsup, xnrmu, xnrmu0,         &
                                      n_buoyant_scal)                       &
      bind(C, name='cs_f_piso_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: nterup, epsup, xnrmu, xnrmu0, n_buoyant_scal
    end subroutine cs_f_piso_get_pointers

    ! Interface to C function retrieving pointers to members of the
    ! global electric model structure

    subroutine cs_f_elec_model_get_pointers(ngazge, ielcor, pot_diff, coejou,  &
                                            elcou, couimp, irestrike, ntdcla,  &
                                            restrike_point_x,                  &
                                            restrike_point_y,                  &
                                            restrike_point_z)                  &
      bind(C, name='cs_f_elec_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ngazge, ielcor, pot_diff, coejou, elcou
      type(c_ptr), intent(out) :: couimp, irestrike, ntdcla, restrike_point_x
      type(c_ptr), intent(out) :: restrike_point_y, restrike_point_z
    end subroutine cs_f_elec_model_get_pointers

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

  !> \brief Initialize Fortran time step API.
  !> This maps Fortran pointers to global C structure members.

  subroutine time_step_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ntpabs, c_ntcabs, c_ntmabs, c_ntinit
    type(c_ptr) :: c_ttpabs, c_ttcabs, c_ttmabs

    call cs_f_time_step_get_pointers(c_ntpabs, c_ntcabs, c_ntmabs, c_ntinit, &
                                     c_ttpabs, c_ttcabs, c_ttmabs)

    call c_f_pointer(c_ntpabs, ntpabs)
    call c_f_pointer(c_ntcabs, ntcabs)
    call c_f_pointer(c_ntmabs, ntmabs)
    call c_f_pointer(c_ntinit, ntinit)

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

    type(c_ptr) :: c_inpdt0, c_iptlro, c_idtvar
    type(c_ptr) :: c_dtref, c_coumax, c_cflmmx
    type(c_ptr) :: c_foumax, c_varrdt, c_dtmin
    type(c_ptr) :: c_dtmax, c_relxst

    call cs_f_time_step_options_get_pointers(c_inpdt0, c_iptlro, c_idtvar, &
                                             c_dtref, c_coumax, c_cflmmx,  &
                                             c_foumax, c_varrdt, c_dtmin,  &
                                             c_dtmax, c_relxst)

    call c_f_pointer(c_inpdt0, inpdt0)
    call c_f_pointer(c_iptlro, iptlro)
    call c_f_pointer(c_idtvar, idtvar)
    call c_f_pointer(c_dtref,  dtref)
    call c_f_pointer(c_coumax, coumax)
    call c_f_pointer(c_cflmmx, cflmmx)
    call c_f_pointer(c_foumax, foumax)
    call c_f_pointer(c_varrdt, varrdt)
    call c_f_pointer(c_dtmin,  dtmin)
    call c_f_pointer(c_dtmax,  dtmax)
    call c_f_pointer(c_relxst, relxst)

  end subroutine time_step_options_init

  !> \brief Initialize Fortran thermal model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine thermal_model_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_itherm, c_itpscl, c_iscalt

    call cs_f_thermal_model_get_pointers(c_itherm, c_itpscl, c_iscalt)

    call c_f_pointer(c_itherm, itherm)
    call c_f_pointer(c_itpscl, itpscl)
    call c_f_pointer(c_iscalt, iscalt)

  end subroutine thermal_model_init

  !> \brief Initialize Fortran turbulence model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine turb_model_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_iturb, c_itytur

    call cs_f_turb_model_get_pointers(c_iturb, c_itytur)

    call c_f_pointer(c_iturb, iturb)
    call c_f_pointer(c_itytur, itytur)

  end subroutine turb_model_init

  !> \brief Initialize Fortran wall functions API.
  !> This maps Fortran pointers to global C structure members.

  subroutine wall_functions_init

    use, intrinsic :: iso_c_binding
    use cstphy, only: ypluli
    implicit none

    ! Local variables

    type(c_ptr) :: c_iwallf, c_iwalfs, c_iwallt, c_ypluli

    call cs_f_wall_functions_get_pointers(c_iwallf, c_iwalfs, c_iwallt, &
                                          c_ypluli)

    call c_f_pointer(c_iwallf, iwallf)
    call c_f_pointer(c_iwalfs, iwalfs)
    call c_f_pointer(c_iwallt, iwallt)
    call c_f_pointer(c_ypluli, ypluli)

  end subroutine wall_functions_init

  !> \brief Initialize Fortran RANS turbulence model API.
  !> This maps Fortran pointers to global C structure members.

  subroutine turb_rans_model_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_irccor, c_itycor, c_idirsm, c_iclkep, c_igrhok, c_igrake
    type(c_ptr) :: c_igrari, c_ikecou, c_reinit_turb, c_irijco, c_irijnu, c_irijrb, c_irijec, c_idifre, c_iddes
    type(c_ptr) :: c_iclsyr, c_iclptr

    call cs_f_turb_rans_model_get_pointers( c_irccor, c_itycor, c_idirsm, &
                                            c_iclkep, c_igrhok, c_igrake, &
                                            c_igrari, c_ikecou, c_reinit_turb, &
                                            c_irijco, c_iddes, c_irijnu, &
                                            c_irijrb, c_irijec, c_idifre, &
                                            c_iclsyr, c_iclptr)

    call c_f_pointer(c_irccor, irccor)
    call c_f_pointer(c_itycor, itycor)
    call c_f_pointer(c_idirsm, idirsm)
    call c_f_pointer(c_iclkep, iclkep)
    call c_f_pointer(c_igrhok, igrhok)
    call c_f_pointer(c_igrake, igrake)
    call c_f_pointer(c_igrari, igrari)
    call c_f_pointer(c_ikecou, ikecou)
    call c_f_pointer(c_reinit_turb, reinit_turb)
    call c_f_pointer(c_irijco, irijco)
    call c_f_pointer(c_iddes, iddes)
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

    type(c_ptr) :: c_idries, c_ivrtex

    call cs_f_turb_les_model_get_pointers( c_idries, c_ivrtex)

    call c_f_pointer(c_idries, idries)
    call c_f_pointer(c_ivrtex, ivrtex)

  end subroutine turb_les_model_init

  !> \brief Initialize Fortran Stokes options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine stokes_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ivisse, c_irevmc, c_iprco, c_arak
    type(c_ptr) :: c_ipucou, c_iccvfg, c_idilat, c_epsdp, c_itbrrb, c_iphydr
    type(c_ptr) :: c_igprij, c_igpust, c_iifren, c_icalhy, c_irecmf


    call cs_f_stokes_options_get_pointers(c_ivisse, c_irevmc, c_iprco ,  &
                                          c_arak  , c_ipucou, c_iccvfg, &
                                          c_idilat, c_epsdp , c_itbrrb, c_iphydr, c_igprij, &
                                          c_igpust, c_iifren, c_icalhy, c_irecmf)

    call c_f_pointer(c_ivisse, ivisse)
    call c_f_pointer(c_irevmc, irevmc)
    call c_f_pointer(c_iprco , iprco )
    call c_f_pointer(c_arak  , arak  )
    call c_f_pointer(c_ipucou, ipucou)
    call c_f_pointer(c_iccvfg, iccvfg)
    call c_f_pointer(c_idilat, idilat)
    call c_f_pointer(c_epsdp , epsdp )
    call c_f_pointer(c_itbrrb, itbrrb)
    call c_f_pointer(c_iphydr, iphydr)
    call c_f_pointer(c_igprij, igprij)
    call c_f_pointer(c_igpust, igpust)
    call c_f_pointer(c_iifren, iifren)
    call c_f_pointer(c_icalhy, icalhy)
    call c_f_pointer(c_irecmf, irecmf)

  end subroutine stokes_options_init

  !> \brief Initialize Fortran space discretisation options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine space_disc_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_imvisf, c_imrgra, c_anomax, c_iflxmw

    call cs_f_space_disc_get_pointers(c_imvisf, c_imrgra, c_anomax, &
                                      c_iflxmw)

    call c_f_pointer(c_imvisf, imvisf)
    call c_f_pointer(c_imrgra, imrgra)
    call c_f_pointer(c_anomax, anomax)
    call c_f_pointer(c_iflxmw, iflxmw)

  end subroutine space_disc_options_init

  !> \brief Initialize Fortran PISO options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine piso_options_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_nterup, c_epsup, c_xnrmu, c_xnrmu0, c_n_buoyant_scal

    call cs_f_piso_get_pointers(c_nterup, c_epsup, c_xnrmu, c_xnrmu0, &
                                c_n_buoyant_scal)

    call c_f_pointer(c_nterup, nterup)
    call c_f_pointer(c_epsup, epsup)
    call c_f_pointer(c_xnrmu, xnrmu)
    call c_f_pointer(c_xnrmu0, xnrmu0)
    call c_f_pointer(c_n_buoyant_scal, n_buoyant_scal)

  end subroutine piso_options_init

  !> \brief Initialize Fortran ELEC options API.
  !> This maps Fortran pointers to global C structure members.

  subroutine elec_option_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ngazge, c_ielcor, c_pot_diff, c_coejou, c_couimp
    type(c_ptr) :: c_elcou, c_irestrike, c_ntdcla, c_restrike_point_x
    type(c_ptr) :: c_restrike_point_y, c_restrike_point_z

    call cs_f_elec_model_get_pointers(c_ngazge, c_ielcor, c_pot_diff,          &
                                      c_coejou, c_elcou, c_couimp,             &
                                      c_irestrike, c_ntdcla,                   &
                                      c_restrike_point_x, c_restrike_point_y,  &
                                      c_restrike_point_z)

    call c_f_pointer(c_ngazge,           ngazge)
    call c_f_pointer(c_ielcor,           ielcor)
    call c_f_pointer(c_pot_diff,         pot_diff)
    call c_f_pointer(c_coejou,           coejou)
    call c_f_pointer(c_elcou,            elcou)
    call c_f_pointer(c_couimp,           couimp)
    call c_f_pointer(c_irestrike,        irestrike)
    call c_f_pointer(c_ntdcla,           ntdcla)
    call c_f_pointer(c_restrike_point_x, restrike_point_x)
    call c_f_pointer(c_restrike_point_y, restrike_point_y)
    call c_f_pointer(c_restrike_point_z, restrike_point_y)

  end subroutine elec_option_init

  !=============================================================================

end module optcal
