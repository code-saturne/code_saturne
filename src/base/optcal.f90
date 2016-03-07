!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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
  ! Equation types
  !----------------------------------------------------------------------------

  !> \defgroup equation_types Equation types

  !> \addtogroup equation_types
  !> \{

  !> take unsteady term into account:
  !>    - 1 account for unsteady term
  !>    - 0 ignore unsteady term
  integer, save :: istat(nvarmx)

  !> take convection into account:
  !>    - 1 account for convection
  !>    - 0 ignore convection
  integer, save :: iconv(nvarmx)

  !> take diffusion into account:
  !>    - 1: true
  !>    - 0: false
  integer, save :: idiff(nvarmx)

  !> take turbulent diffusion into account:
  !>    - 1: true
  !>    - 0: false
  integer, save :: idifft(nvarmx)

  !> type of diffusivity:
  !>    - 1: scalar diffusivity
  !>    - 3: orthotropic diffusivity
  !>    - 6: symmetric tensor diffusivity
  integer, save :: idften(nvarmx)

  !> \}

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
  !>    - 2: theta scheme with theta > 0 (theta=0.5 means 2nd order)
  !>    - 0: theta scheme with theta = 0 (explicit)
  !>    - 1: implicit scheme (default)
  integer, save ::          istmpf

  !> number of interations on the pressure-velocity coupling on Navier-Stokes
  !> (for the PISO algorithm)
  integer(c_int), pointer, save ::          nterup

  !> extrapolation of source terms in the Navier-Stokes equations
  !>    - 1: true
  !>    - 0: false (default)
  integer, save ::          isno2t

  !> extrapolation of turbulent quantities
  !>    - 1: true
  !>    - 0: false (default)
  integer, save ::          isto2t

  !> extrapolation of source terms in the transport equation of scalars
  !>    - 1: true
  !>    - 0: false (default)
  integer, save ::          isso2t(nscamx)

  !> extrapolation of the density field
  !>    - 1: true
  !>    - 0: false (default)
  integer, save ::          iroext

  !> extrapolation of the total viscosity field
  !>    - 1: true
  !>    - 0: false (default)
  integer, save ::          iviext

  !> extrapolation of the specific heat field \f$ C_p \f$
  !>    - 1: true
  !>    - 0: false (default)
  integer, save ::          icpext

  !> extrapolation of the scalar diffusivity
  !>    - 1: true
  !>    - 0: false (default)
  integer, save ::          ivsext(nscamx)

  !> initvi : =1 if total viscosity read from checkpoint file
  integer, save ::          initvi

  !> initro : =1 if density read from checkpoint file
  integer, save ::          initro

  !> initcp : =1 if specific heat read from checkpoint file
  integer, save ::          initcp

  !> ibdtso : backward differential scheme in time order
  integer, save ::          ibdtso(nvarmx)

  !> initvs : =1 if scalar diffusivity read from checkpoint file
  integer, save ::          initvs(nscamx)

  !> \f$ \theta \f$-scheme for the main variables
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetav(nvarmx)

  !> \f$ \theta_S \f$-scheme for the source terms in the Navier-Stokes equations
  !>    -  0 : second viscosity explicit
  !>    - 1/2: second viscosity extrapolated in n+1/2
  !>    -  1 : second viscosity extrapolated in n+1
  double precision, save :: thetsn

  !> \f$ \theta \f$-scheme for the source terms of turbulent equations
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetst

  !> \f$ \theta \f$-scheme for the source terms of transport equations of scalars
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetss(nscamx)

  !> \f$ \theta \f$-scheme for the mass flux
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetfl

  !> \f$ \theta \f$-scheme for the total viscosity
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetvi

  !> \f$ \theta \f$-scheme for the density field
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetro

  !> \f$ \theta \f$-scheme for the scpecific heat field
  !>    -  0 : explicit
  !>    - 1/2: extrapolated in n+1/2
  !>    -  1 : extrapolated in n+1
  double precision, save :: thetcp

  !> \f$ \theta \f$-scheme for the diffusivity
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

  !> percentage of upwind:
  !>    - 1: no upwind (except if the slope test is activated)
  !>    - 0: total upwind
  double precision, save :: blencv(nvarmx)

  !> type of convective scheme
  !>    - 0: second order linear upwind
  !>    - 1: centered
  !>    - 2: pure upwind gradient in SOLU
  integer, save ::          ischcv(nvarmx)

  !> Slope test, Min/MAx limiter or Roe and Sweby limiters
  !>    - 0: swich on the slope test
  !>    - 1: swich off the slope test (default)
  !>    - 2: continuous limiter ensuring positivness
  !>    - 3: Roe-Sweby limiter
  !>         (ensuring  Decreasing Total Variation)
  integer, save ::          isstpc(nvarmx)

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

  !> anomax : angle de non orthogonalite des faces en radian au dela duquel
  !> on retient dans le support etendu des cellules voisines
  !> de la face les cellules dont un noeud est sur la face
  real(c_double), pointer, save :: anomax

  !> max number of iterations for the iterative gradient
  integer, save ::          nswrgr(nvarmx)

  !> relative precision of the iterative gradient calculation
  double precision, save :: epsrgr(nvarmx)

  !> type of gradient clipping
  !>    - < 0: no clipping
  !>    -   0: first order
  !>    -   1: second order
  integer, save ::          imligr(nvarmx)

  !>   climgr : facteur de limitation (>=1, =1 : forte limitation)
  double precision, save :: climgr(nvarmx)

  !> gradient extrapolation at the boundary
  !>    - 0: false
  !>    - 1: true
  double precision, save :: extrag(nvarmx)

  !> gradient calculation
  !>    - 0: standard
  !>    - 1: weighted
  integer, save :: iwgrec(nvarmx)

  !> \}

  !> \defgroup diffusive_scheme Diffusive scheme
  !> \addtogroup diffusive_scheme
  !> \{

  !> face flux reconstruction:
  !>    - 0: false
  !>    - 1: true
  integer, save ::          ircflu(nvarmx)

  !> face viscosity field interpolation
  !>    - 1: harmonic
  !>    - 0: arithmetic (default)
  integer(c_int), pointer, save :: imvisf

  !> \}

  !> \defgroup iterative_process Iterative process for the convection diffusion
  !>           equation
  !> \addtogroup iterative_process
  !> \{

  !> max number of iteration for the iterative process used to solved
  !> the convection diffusion equations
  integer, save ::          nswrsm(nvarmx)

  !> relative precision of the iterative process used to solved
  !> the convection diffusion equations
  double precision, save :: epsrsm(nvarmx)

  !> dynamic relaxation type:
  !>    - 0 no dynamic relaxation
  !>    - 1 dynamic relaxation depending on \f$ \delta \varia^k \f$
  !>    - 2 dynamic relaxation depending on \f$ \delta \varia^k \f$ and \f$ \delta \varia^{k-1} \f$
  integer, save :: iswdyn(nvarmx)

  !> \}

  !> \defgroup linear_solver Linear solver
  !> \addtogroup linear_solver
  !> \{

  !> relative precision of the linear solver
  double precision, save :: epsilo(nvarmx)

  !> strengthening of the diagonal part of the matrix if no Dirichlet is set
  !>    - 0: false
  !>    - 1: true
  !> \remark
  !> the code computes automatically for each variable the number of Dirichlet
  !> BCs
  integer, save ::          idircl(nvarmx)

  !> number of Dirichlet BCs
  integer, save ::          ndircl(nvarmx)

  !> \}

  !> \}

  !TODO doxygen it
  ! Gestion du calcul
  !   isuite : suite de calcul
  !     = 0 pour sfs
  !     = 1 pour suite de calcul
  !   iecaux : ecriture du suite auxiliaire
  !   ileaux : lecture  du suite auxiliaire
  !   isuit1 : suite du module thermique 1D en paroi
  !   isuict : suite du module aerorefrigerant
  !   isuivo : suite de la methode des vortex
  !   isuisy : suite des methodes d entree LES

  integer, save :: isuite , ileaux, iecaux,                        &
                   isuit1 , isuict, isuivo, isuisy

  !----------------------------------------------------------------------------
  ! Time stepping options
  !----------------------------------------------------------------------------

  !> \defgroup time_step_options Time step options and variables

  !> \addtogroup time_step_options
  !> \{

  !> Absolute time step number for previous calculation.
  integer(c_int), pointer, save :: ntpabs

  !> Current absolute time step number.
  !> In case of restart, this is equal to ntpabs + number of new iterations.
  integer(c_int), pointer, save :: ntcabs

  !> Maximum absolute time step number.
  integer(c_int), pointer, save :: ntmabs

  !> Number of time steps for initalization.
  integer(c_int), pointer, save :: ntinit

  !> Absolute time value for previous calculation.
  real(c_double), pointer, save :: ttpabs

  !> Current absolute time.
  !> In case of restart, this is equal to ttpabs + additional computed time.
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
  !>    - 0: false
  !>    - 1: true
  integer(c_int), pointer, save :: iptlro

  !> option for a variable time step
  !>    - -1: steady algorithm
  !>    -  0: constant time step
  !>    -  1: time step constant in space but variable in time
  !>    -  2: variable time step in space and in time
  integer(c_int), pointer, save :: idtvar

  !> reference time step
  real(c_double), pointer, save :: dtref

  !> maximum Courant number (when idtvar is different from 0)
  real(c_double), pointer, save :: coumax

  !> maximum Courant number for the continuity equation in compressible model
  real(c_double), pointer, save :: cflmmx

  !> maximum Fourier number (when idtvar is different from 0)
  real(c_double), pointer, save :: foumax

  !> relative allowed variation of dt (when idtvar is different from 0)
  real(c_double), pointer, save :: varrdt

  !> minimum value of dt (when idtvar is different from 0).
  !> Take dtmin = min (ld/ud, sqrt(lt/(gdelta rho/rho)), ...)
  real(c_double), pointer, save :: dtmin

  !> maximum value of dt (when idtvar is different from 0).
  !> Take dtmax = max (ld/ud, sqrt(lt/(gdelta rho/rho)), ...)
  real(c_double), pointer, save :: dtmax

  !> multiplicator coefficient for the time step of each variable
  !>    - useless for u,v,w,p
  !>    - for k,e     the same value is taken (value of k)
  !>    - for Rij, e  the same value is taken (value of r11)
  double precision, save :: cdtvar(nvarmx)

  !> relaxation of variables (1 stands fo no relaxation)
  double precision, save :: relaxv(nvarmx)

  !> relaxation coefficient for the steady algorithm
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
  !>    - 3: total energy (only for compressible module)
  integer(c_int), pointer, save :: itherm

  !> temperature scale
  !>    - 0: none
  !>    - 1: Kelvin
  !>    - 2: Celsius
  integer(c_int), pointer, save :: itpscl

  !> index of the thermal scalar (temperature, energy of enthalpy),
  !> the index of the corresponding variable is isca(iscalt)
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
  !>    - 0: no wall functions
  !>    - 1: one scale of friction velocities (power law)
  !>    - 2: one scale of friction velocities (log law)
  !>    - 3: two scales of friction velocities (log law)
  !>    - 4: two scales of friction velocities (log law) - scalable wall functions
  !>    - 5: two scales of friction velocities (mixing
  !>          length based on V. Driest analysis)
  integer(c_int), pointer, save :: iwallf

  !>  Wall functions for scalar
  !>    - 0: three layer wall function of Arpaci and Larsen
  !>    - 1: Van Driest wall function
  integer(c_int), pointer, save :: iwalfs

  !> exchange coefficient correlation
  !>    - 0: not use by default
  !>    - 1: exchange coefficient computed with a correlation
  integer(c_int), pointer, save :: iwallt

  !> clipping of k and epsilon
  !>    - 0 absolute value clipping
  !>    - 1 coupled clipping based on physical relationships
  integer(c_int), pointer, save :: iclkep

  !> take \f$ 2/3 \rho \grad k \f$ in the momentum equation
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: igrhok

  !> buoyant term in \f$ k- \varepsilon \f$
  !>    - 1: true (default if \f$ \rho \f$ is variable)
  !>    - 0: false
  integer(c_int), pointer, save :: igrake

  !> buoyant term in \f$ R_{ij}- \varepsilon \f$
  !>    - 1: true (default if \f$ \rho \f$ is variable)
  !>    - 0: false
  integer(c_int), pointer, save :: igrari

  !> partially coupled version of \f$ k-\varepsilon \f$ (only for iturb=20)
  !>    - 1: true (default)
  !>    - 0: false
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
  integer(c_int), pointer, save :: irijnu

  !> accurate treatment of \f$ \tens{R} \f$ at the boundary (see \ref condli)
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: irijrb

  !> wall echo term of \f$ \tens{R} \f$
  !>    - 1: true
  !>    - 0: false (default)
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

  !> Van Driest smoothing at the wall (only for itytur=4)
  !>    - 1: true
  !>    - 0: false
  integer(c_int), pointer, save :: idries

  !> vortex method (in LES)
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: ivrtex

  !> turbulent flux model for \f$ \overline{\varia^\prime \vect{u}^\prime} \f$
  !> for any scalar \f$ \varia \f$, iturt(isca)
  !>    - 0: SGDH
  !>    - 10: GGDH
  !>    - 20: AFM
  !>    - 30: DFM (Transport equation modelized)
  integer, save :: iturt(nscamx)
  !    - 11: EB-GGDH
  !    - 21: EB-AFM
  !    - 31: EB-DFM

  !> class turbulent flux model (=iturt/10)
  integer, save :: ityturt(nscamx)

  !> index of the turbulent flux for the scalar iscal
  integer, save :: ifltur(nscamx)

  !> number of variable plus number of turbulent fluxes
  !> (used by the boundary conditions)
  integer(c_int), pointer, save :: nvarcl

  !> \}

  !----------------------------------------------------------------------------
  ! Stokes
  !----------------------------------------------------------------------------

  !> \defgroup stokes Stokes options

  !> \addtogroup stokes
  !> \{

  !> take \f$ \divs \left( \mu \transpose{\gradt \, \vect{u}} - 2/3 \mu \trace{\gradt \, \vect{u}} \right) \f$
  !> into account in the momentum equation
  !>    - 1: true (default)
  !>    - 0: false
  integer(c_int), pointer, save :: ivisse

  !> Reconstruction of the velocity field with the updated pressure option
  !>    - 0: default
  integer(c_int), pointer, save ::          irevmc

  !> Compute the pressure step thanks to the continuity equation
  !>    - 1: true (default)
  !>    - 0: false
  integer(c_int), pointer, save ::          iprco

  !> Compute the normed residual for the pressure step in the prediction step
  !>    - 1: true (default)
  !>    - 0: false
  integer(c_int), pointer, save ::          irnpnw

  !> normed residual for the pressure step
  real(c_double), pointer, save :: rnormp

  !> Arakawa multiplicator for the Rhie and Chow filter (1 by default)
  real(c_double), pointer, save :: arak

  !> Pseudo coupled pressure-velocity solver
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: ipucou

  !> \anchor iccvfg
  !> calculation with a fixed velocity field
  !>    - 1: true
  !>    - 0: false (default)
  integer(c_int), pointer, save :: iccvfg

  !> Algorithm to take into account the density variation in time
  !>    - 1: dilatable steady algorithm (default)
  !>    - 2: dilatable unsteady algorithm
  !>    - 3: low-Mach algorithm
  !>    - 4: algorithm for fire
  !    - 0: boussinesq algorithm with constant density
  integer(c_int), pointer, save :: idilat

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

  !> indicates if the scalar isca is coupled with syrthes
  !>    - 1: coupled with syrthes
  !>    - 0: uncoupled
  !>
  !> \remark
  !> only one scalar can be coupled with syrthes
  integer, save :: icpsyr(nscamx)

  !> improve hydrostatic pressure algorithm
  !>    - 1: impose the equilibrium of the hydrostaic part of the pressure with any external force, even head losses
  !>    - 2: compute an hydrostatic pressure due to buoyancy forces before the prediction step
  !>    - 0: no treatment (default)
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


  !> compute error estimators
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: iescal(nestmx)

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

  !> Indicateur module thermique 1d global (ie somme sur les processeurs
  !>   de nfpt1d)
  integer, save :: nfpt1t

  !> Indicateur termes sources de masse global (ie somme sur les processeurs
  !>   de ncetsm)
  integer, save :: nctsmt

  !> Global indicator of condensation source terms (ie. sum on the processors
  !> of nfbpcd) cells associated to the face with condensation phenomenon
  integer, save :: nftcdt

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

  !> ineedy : = 1 distance a la paroi est necessaire pour le calcul
  !>          = 0 distance a la paroi n'est pas necessaire
  integer, save :: ineedy

  !> imajdy : = 1 distance a la paroi a ete mise a jour
  !>          = 0 distance a la paroi n'a pas ete mise a jour
  integer, save :: imajdy

  !> icdpar : = 1 calcul standard (et relecture en suite de calcul)
  !>          = 2 calcul ancien   (et relecture en suite de calcul)
  !>          =-1 forcer le recalcul en suite (par calcul standard)
  !>          =-2 forcer le recalcul en suite (par calcul ancien)
  integer, save :: icdpar

  !> nitmay : nombre max d'iterations pour les resolutions iteratives
  integer, save :: nitmay

  !> nswrsy : nombre de sweep pour reconstruction des s.m.
  integer, save :: nswrsy

  !> nswrgy : nombre de sweep pour reconstruction des gradients
  integer, save :: nswrgy

  !> imligy : methode de limitation du gradient
  integer, save :: imligy

  !> ircfly : indicateur pour reconstruction des flux
  integer, save :: ircfly

  !> ischcy : indicateur du schema en espace
  integer, save :: ischcy

  !> isstpy : indicateur pour test de pente
  integer, save :: isstpy

  !> iwarny : niveau d'impression
  integer, save :: iwarny

  !> ntcmxy : nombre max d'iteration pour la convection de y
  integer, save :: ntcmxy

  ! blency : 1 - proportion d'upwind
  double precision, save :: blency

  ! epsily : precision pour resolution iterative
  double precision, save :: epsily

  ! epsrsy : precision pour la reconstruction du second membre
  double precision, save :: epsrsy

  ! epsrgy : precision pour la reconstruction des gradients
  double precision, save :: epsrgy

  ! climgy : coef gradient*distance/ecart
  double precision, save :: climgy

  ! extray : coef d'extrapolation des gradients
  double precision, save :: extray

  ! coumxy : valeur max   du courant pour equation convection
  double precision, save :: coumxy

  ! epscvy : precision pour convergence equation convection stationnaire
  double precision, save :: epscvy

  ! yplmxy : valeur max   de yplus au dessus de laquelle l'amortissement de
  !          Van Driest est sans effet et donc pour laquelle un calcul de
  !          yplus moins precis est suffisant
  double precision, save :: yplmxy

  !> \}

  !----------------------------------------------------------------------------
  ! Transported scalars parameters
  !----------------------------------------------------------------------------

  !> \defgroup scalar_params Transported scalars parameters

  !> \addtogroup scalar_params
  !> \{

  !> iscacp : 0 : scalar does not behave like a temperature
  !>          1 : scalar behaves like a temperature (use Cp for wall law)
  !>        > 1 : not yet allowed, could be used for multiple Cp definitions
  integer, save ::          iscacp(nscamx)

  !> iclvfl : 0 : clip variances to zero
  !>          1 : clip variances to zero and to f(1-f)
  !>          2 : clip variances to  max(zero,scamin) and scamax
  integer, save ::          iclvfl(nscamx)

  !> iscasp(ii) : index of the ii^th species (0 if not a species)
  integer, save ::          iscasp(nscamx)

  !> visls0 : viscosity of scalars if constant
  double precision, save :: visls0(nscamx)

  !> sigmas : prandtl of scalars
  double precision, save :: sigmas(nscamx)

  !> rvarfl : coeff de dissipation des variances
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

  parameter (DRIFT_SCALAR_ADD_DRIFT_FLUX=1)
  parameter (DRIFT_SCALAR_THERMOPHORESIS=2)
  parameter (DRIFT_SCALAR_TURBOPHORESIS=3)
  parameter (DRIFT_SCALAR_ELECTROPHORESIS=4)
  parameter (DRIFT_SCALAR_CENTRIFUGALFORCE=5)

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

    subroutine cs_f_turb_model_get_pointers(iturb, itytur, nvarcl) &
      bind(C, name='cs_f_turb_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: iturb, itytur, nvarcl
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
                                                 igrari, ikecou, reinit_turb, irijco, irijnu, &
                                                 irijrb, irijec, idifre, &
                                                 iclsyr, iclptr)         &
      bind(C, name='cs_f_turb_rans_model_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: irccor, itycor, idirsm, iclkep, igrhok
      type(c_ptr), intent(out) :: igrake, igrari, ikecou, reinit_turb, irijco, irijnu, irijrb
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

    subroutine cs_f_stokes_options_get_pointers(ivisse, irevmc, iprco, irnpnw, &
                                                rnormp, arak  ,ipucou, iccvfg, &
                                                idilat, epsdp ,itbrrb, iphydr, &
                                                igprij, igpust,                &
                                                iifren, icalhy, irecmf)        &
      bind(C, name='cs_f_stokes_options_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: ivisse, irevmc, iprco, irnpnw, rnormp, arak
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

    subroutine cs_f_piso_get_pointers(nterup, epsup, xnrmu, xnrmu0) &
      bind(C, name='cs_f_piso_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: nterup, epsup, xnrmu, xnrmu0
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

    type(c_ptr) :: c_iturb, c_itytur, c_nvarcl

    call cs_f_turb_model_get_pointers(c_iturb, c_itytur, c_nvarcl)

    call c_f_pointer(c_iturb, iturb)
    call c_f_pointer(c_itytur, itytur)
    call c_f_pointer(c_nvarcl, nvarcl)

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
    type(c_ptr) :: c_igrari, c_ikecou, c_reinit_turb, c_irijco, c_irijnu, c_irijrb, c_irijec, c_idifre
    type(c_ptr) :: c_iclsyr, c_iclptr

    call cs_f_turb_rans_model_get_pointers( c_irccor, c_itycor, c_idirsm, &
                                            c_iclkep, c_igrhok, c_igrake, &
                                            c_igrari, c_ikecou, c_reinit_turb, c_irijco, c_irijnu, &
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

    type(c_ptr) :: c_ivisse, c_irevmc, c_iprco, c_irnpnw, c_rnormp, c_arak
    type(c_ptr) :: c_ipucou, c_iccvfg, c_idilat, c_epsdp, c_itbrrb, c_iphydr
    type(c_ptr) :: c_igprij, c_igpust, c_iifren, c_icalhy, c_irecmf


    call cs_f_stokes_options_get_pointers(c_ivisse, c_irevmc, c_iprco ,  &
                                          c_irnpnw, c_rnormp, c_arak  , c_ipucou, c_iccvfg, &
                                          c_idilat, c_epsdp , c_itbrrb, c_iphydr, c_igprij, &
                                          c_igpust, c_iifren, c_icalhy, c_irecmf)

    call c_f_pointer(c_ivisse, ivisse)
    call c_f_pointer(c_irevmc, irevmc)
    call c_f_pointer(c_iprco , iprco )
    call c_f_pointer(c_irnpnw, irnpnw)
    call c_f_pointer(c_rnormp, rnormp)
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

    type(c_ptr) :: c_nterup, c_epsup, c_xnrmu, c_xnrmu0

    call cs_f_piso_get_pointers(c_nterup, c_epsup, c_xnrmu, c_xnrmu0)

    call c_f_pointer(c_nterup, nterup)
    call c_f_pointer(c_epsup, epsup)
    call c_f_pointer(c_xnrmu, xnrmu)
    call c_f_pointer(c_xnrmu0, xnrmu0)

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
