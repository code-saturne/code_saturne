!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

  !> take non-stationary term into account:
  !>    - 1 prise en compte du terme instationnaire
  !>    - 0 prise en compte du terme instationnaire
  integer, save :: istat(nvarmx)

  !> take convection into account:
  !>    - 1 prise en compte de la convection
  !>    - 0 non prise en compte de la convection
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

  !> variable density field \f$ \rho \f$:
  !>    - 1: true
  !>    - 0: false
  integer, save :: irovar

  !> variable viscosity field \f$ \mu \f$:
  !>    - 1: true
  !>    - 0: false
  integer, save :: ivivar

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
  integer, save ::          nterup

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

  !> extrapolation of the scpecific heat field \f$ C_p \f$
  !>    - 1: true
  !>    - 0: false (default)
  integer, save ::          icpext

  !> extrapolation of the scalar diffusivity
  !>    - 1: true
  !>    - 0: false (default)
  integer, save ::          ivsext(nscamx)

  !> initvi : =1 si viscosite totale relue dans un suite
  integer, save ::          initvi

  !> initro : =1 si masse volumique relue dans un suite
  integer, save ::          initro

  !> initcp : =1 si  chaleur specifique relue dans un suite
  integer, save ::          initcp

  !> initvs : =1 si  diffusivite scalaire relue dans un suite
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
  double precision, save :: epsup

  !> norm  of the increment \f$ \vect{u}^{k+1} - \vect{u}^k \f$
  !> of the iterative process on pressure-velocity coupling (PISO)
  double precision, save :: xnrmu

  !> norm of \f$ \vect{u}^0 \f$ (used by PISO algorithm)
  double precision, save :: xnrmu0

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
  !>    - 1: centre
  !>    - 0: second order
  integer, save ::          ischcv(nvarmx)

  !> switch off the slope test:
  !>    - 1: swich off the slope test
  !>    - 0: swich on the slope test
  integer, save ::          isstpc(nvarmx)

  !> method to compute interior mass flux due to ALE mesh velocity
  !>    - 1: based on cell center mesh velocity
  !>    - 0: based on nodes displacement
  integer, save ::          iflxmw

  !> \}

  !> \defgroup gradient_calculation Gradient calculation
  !> \addtogroup gradient_calculation
  !> \{

  !> type of gradient reconstruction
  !>    - 0: iterative process
  !>    - 1: standard least suqare methode
  !>    - 2: least suqare methode with extended neighbourhood
  !>    - 3: least suqare methode with reduced extended neighbourhood
  !>    - 4: iterative precess initialized by the least square methode
  integer, save ::          imrgra

  !> anomax : angle de non orthogonalite des faces en radian au dela duquel
  !> on retient dans le support etendu des cellules voisines
  !> de la face les cellules dont un noeud est sur la face
  double precision, save :: anomax

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
  integer, save :: imvisf

  !> \}

  !> \defgroup iterative_process Iterative process for the convection diffusion equation
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

  !> maximal number of iteration for the linear solver
  integer, save ::          nitmax(nvarmx)

  !> relative precision of the linear solver
  double precision, save :: epsilo(nvarmx)

  !> type of linear solver
  !>    - (-1): automatic choice
  !>    -    0: conjugate gradient
  !>    -    1: Jacobi
  !>    -    2: bi-CGSTAB
  !> \remark
  !>  we add ipol*1000 to iresol(ivar) where ipol is the degree of the polynome
  !>  of Neumann preconditionning.
  integer, save ::          iresol(nvarmx)

  !> strengthen of the diagonal part of the matrix if no Dirichlet is set
  !>    - 0: false
  !>    - 1: true
  !> \remark
  !> the code computes automatically for each variable the number of Dirichlet
  !> BCs
  integer, save ::          idircl(nvarmx)

  !> number of Dirichlet BCs
  integer, save ::          ndircl(nvarmx)

  !> multigrid algorithm
  !>    - 0: false
  !>    - 1: algebraic multigrid
  integer, save ::          imgr(nvarmx)

  !> maximal number of cycles in the multigrid algorithm
  integer, save ::          ncymax(nvarmx)

  !> number of iterations on the finer mesh
  integer, save ::          nitmgf(nvarmx)

  !> relaxation parameter for the multigrid
  double precision, save :: rlxp1

  !> \}

  !> \}

  !TODO doxygen it
  ! Gestion du calcul
  !   isuite : suite de calcul
  !     = 0 pour sfs
  !     = 1 pour suite de calcul
  !   iscold : correspondance nouveaux-anciens scalaires
  !   iecaux : ecriture du suite auxiliaire
  !   ileaux : lecture  du suite auxiliaire
  !   isuit1 : suite du module thermique 1D en paroi
  !   isuict : suite du module aerorefrigerant
  !   isuivo : suite de la methode des vortex
  !   isuisy : suite des methodes d entree LES

  integer, save :: isuite , ileaux, iecaux, iscold(nscamx),        &
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
  !>         - pour les calculs non suite :
  !>           on saute uniquement les resolutions (Navier-Stokes,
  !>           turbulence, scalaires...)
  !>         - pour les calculs suite :
  !>           on saute les resolutions (navier-stokes,
  !>           turbulence, scalaires...) et le calcul des proprietes
  !>           physiques, les conditions aux limites (les grandeurs
  !>           sont lues dans le fichier suite)
  integer, save ::          inpdt0

  !> Clip the time step with respect to the buoyant effects
  !>    - 0: false
  !>    - 1: true
  integer, save ::          iptlro

  !> option for a variable time step
  !>    - -1: stationary algorithm
  !>    -  0: constant time step
  !>    -  1: time step constant in space but variable in time
  !>    -  2: variable time step in space and in time
  integer, save ::          idtvar

  !> reference time step
  double precision, save :: dtref

  !> maximum Courant number (when idtvar is different from 0)
  double precision, save :: coumax

  !> maximum Courant number for the continuity equation in compressible model
  double precision, save :: cflmmx

  !> maximum Fourier number (when idtvar is different from 0)
  double precision, save :: foumax

  !> relative allowed variation of dt (when idtvar is different from 0)
  double precision, save :: varrdt

  !> minimum value of dt (when idtvar is different from 0).
  !> Take dtmin = min (ld/ud, sqrt(lt/(gdelta rho/rho)), ...)
  double precision, save :: dtmin

  !> maximum value of dt (when idtvar is different from 0).
  !> Take dtmax = max (ld/ud, sqrt(lt/(gdelta rho/rho)), ...)
  double precision, save :: dtmax

  !> multiplicator coefficient for the time step of each variable
  !>    - useless for u,v,w,p
  !>    - for k,e     the same value is taken (value of k)
  !>    - for Rij, e  the same value is taken (value of r11)
  double precision, save :: cdtvar(nvarmx)

  !> relaxation of variables (1 stands fo no relaxation)
  double precision, save :: relaxv(nvarmx)

  !> relaxation coefficient for the stationary algorithm
  double precision, save :: relxst

  !> \}

  !----------------------------------------------------------------------------
  ! turbulence
  !----------------------------------------------------------------------------

  !> \defgroup turbulence turbulence options

  !> \addtogroup turbulence
  !> \{

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
  integer, save :: iturb

  !> Class of turbulence model (integer value iturb/10)
  integer, save :: itytur

  !> Activation of rotation/curvature correction for an eddy viscosity turbulence models
  !>    - 0: false
  !>    - 1: true
  integer, save :: irccor

  !> Type of rotation/curvature correction for an eddy viscosity turbulence models
  !>    - 1 Cazalbou correction (default when irccor=1 and itytur=2 or 5)
  !>    - 2 Spalart-Shur correction (default when irccor=1 and iturb=60 or 70)
  integer, save :: itycor

  !>  Wall functions
  !>    - 0: one scale of friction velocities
  !>    - 1: two scale of friction velocities
  !>    - 2: scalable wall functions
  integer, save :: ideuch

  !> exchange coefficient correlation
  !>    - 0: not use by default
  !>    - 1: exchange coefficient computed with a correlation
  integer, save :: iwallt

  !> wall function with
  !>    - 0 a power lay (deprecated)
  !>    - 1 a log lay
  integer, save :: ilogpo

  !> clipping of k and epsilon
  !>    - 0 absolute value clipping
  !>    - 1 coupled clipping based on physical relationships
  integer, save :: iclkep

  !> take \f$ 2/3 \rho \grad k \f$ in the momentum equation
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: igrhok

  !> buoyant term in \f$ k- \varepsilon \f$
  !>    - 1: true (default if \f$ \rho \f$ is variable)
  !>    - 0: false
  integer, save :: igrake

  !> buoyant term in \f$ R_{ij}- \varepsilon \f$
  !>    - 1: true (default if \f$ \rho \f$ is variable)
  !>    - 0: false
  integer, save :: igrari

  !> index of the thermal scalar (temperature, energy of enthalpy),
  !> the index of the corresponding variable is isca(iscalt)
  integer, save :: iscalt

  !> partially coupled version of \f$ k-\varepsilon \f$ (only for iturb=20)
  !>    - 1: true (default)
  !>    - 0: false
  integer, save :: ikecou

  !> pseudo eddy viscosity in the matrix of momentum equation to partially
  !> implicit \f$ \divv \left( \rho \tens{R} \right) \f$
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: irijnu

  !> accurate treatment of \f$ \tens{R} \f$ at the boundary (see \ref condli)
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: irijrb

  !> wall echo term of \f$ \tens{R} \f$
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: irijec

  !> whole treatment of the diagonal part of the dissusion tensor of
  !> \f$ \tens{R} \f$ and \f$ \varepsilon \f$
  !>    - 1: true (default)
  !>    - 0: simplified treatment
  integer, save :: idifre

  !> partial implicitation of symmetry BCs of \f$ \tens{R} \f$
  !>    - 1: true (default)
  !>    - 0: false
  integer, save :: iclsyr

  !> partial implicitation of wall BCs of \f$ \tens{R} \f$
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: iclptr

  !> Van Driest smoothing at the wall (only for itytur=4)
  !>    - 1: true
  !>    - 0: false
  integer, save :: idries

  !> vortex method (in LES)
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: ivrtex

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
  !> (used by the Boundary Conditions)
  integer, save :: nvarcl

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
  integer, save :: ivisse

  !> Reconstruction of the velocity field with the updated pressure option
  !>    - 0: default
  integer, save ::          irevmc

  !> Compute the pressure step thanks to the continuity equation
  !>    - 1: true (default)
  !>    - 0: false
  integer, save ::          iprco

  !> Compute the normed residual for the pressure step in the prediction step
  !>    - 1: true (default)
  !>    - 0: false
  integer, save ::          irnpnw

  !> normed residual for the pressure step
  double precision, save :: rnormp

  !> Arakawa multiplicator for the Rhie and Chow filter (1 by default)
  double precision, save :: arak

  !> Pseudo coupled pressure-velocity solver
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: ipucou

  !> calculation with a fixed velocity field
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: iccvfg

  !> Algorithm to take into account the density variation in time
  !>    - 1: dilatable steady algorithm (default)
  !>    - 2: dilatable unsteady algorithm
  !>    - 3: low-Mach algorithm
  !>    - 4: algorithm for fire
  !    - 0: boussinesq algorithm with constant density
  integer, save :: idilat

  !> parameter of diagonal pressure strengthening
  double precision, save :: epsdp

  !TODO doxygen
  ! Type des conditions limites et index min et max
  !                 des sous listes defaces de bord
  integer, save :: idebty(ntypmx), ifinty(ntypmx)

  !> accurate treatment of the wall temperature
  !>    - 1: true
  !>    - 0: false (default)
  !> (see \ref condli, usefull in case of coupling with syrthes)
  integer, save :: itbrrb

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
  integer, save :: iphydr

  !> compute the hydrostatic pressure in order to compute the Dirichlet
  !> conditions on the pressure at outlets
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: icalhy

  !TODO doxygen
  ! icond: Handling condensation source terms
  !        1: condensation source terms activated
  !        2: condensation source terms with metal structures activated
  !        0: by default (without condensation source terms)
  integer, save :: icond

  !> compute error estimators
  !>    - 1: true
  !>    - 0: false (default)
  integer, save :: iescal(nestmx)

  !> \}

  !----------------------------------------------------------------------------
  ! Temporal mean and moments calculation
  !----------------------------------------------------------------------------

  !> \defgroup mean_moments Temporal mean and moments calculation

  !> \addtogroup mean_moments
  !> \{

  !> number of moments
  integer, save ::          nbmomt

  !> nombre de tableaux ncel pour le temps cumule
  integer, save ::          nbdtcm

  !> index of the initial time step for computing the moment
  integer, save ::          ntdmom(nbmomx)

  !> numero de l'ancien moment correspondant en cas de suite
  integer, save ::          imoold(nbmomx)

  !> icmome : pointeur pour les moments (donne un numero de propriete)
  !>           s'utilise ainsi propce(iel,ipproc(icmome(imom)))
  integer, save ::          icmome(nbmomx)

  !> numero du temps cumule associe aux moments
  !> ce numero va de 1 a n pour les temps cumules non uniformes
  !> et de -1 a -p pour les temps cumules uniformes
  !> s'utilise ainsi
  !>    - si idtmom(imom) > 0 propce(iel,ipropc(icdtmo(idtmom(imom))))
  !>    - si idtmom(imom) < 0 dtcmom(-idtmom(imom))
  integer, save ::          idtmom(nbmomx)

  !> numero des variables composant le moment idfmom(jj,imom)
  integer, save ::          idfmom(ndgmox,nbmomx)

  !> moment degree
  integer, save ::          idgmom(nbmomx)

  !> numero de propriete du temps cumule (voir idtmom)
  integer, save ::          icdtmo(nbmomx)

  !> valeur du pas de temps cumule quand il est uniforme (see \ref idtmom).
  double precision, save :: dtcmom(nbmomx)

  !> initial time for computing the moment
  double precision, save :: ttdmom(nbmomx)

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

  !> imgrpy : multigrille
  integer, save :: imgrpy

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

  !TODO move it elsewhere?
  ! Parametres numeriques pour le calcul des efforts aux bords

  !> ineedf : = 1 on calcule les efforts aux parois
  !>          = 0 on ne calcule pas les efforts aux parois
  integer, save :: ineedf

  !> \}

  !----------------------------------------------------------------------------
  ! Transported scalars parameters
  !----------------------------------------------------------------------------

  !> \defgroup scalar_params Transported scalars parameters

  !> \addtogroup scalar_params
  !> \{

  !> iscsth
  !>   -1 : de type temperature en C (      Cp pour la loi de paroi)
  !>    0 : scalaire passif      (ie pas de Cp pour la loi de paroi)
  !>    1 : de type temperature en K (      Cp pour la loi de paroi)
  !>    2 : enthalpie            (ie pas de Cp pour la loi de paroi)
  !>    3 : energie (en compressible, pas de Cp pour la loi de paroi)
  !>      la distinction C/K sert en rayonnement
  integer, save ::          iscsth(nscamx)

  !> ivisls : si positif strictement, indique que la viscosite associee
  !>            au scalaire est variable, et la valeur est le numero
  !>            d'ordre de la viscosite dans le tableau des viscosites
  !>            variables
  integer, save ::          ivisls(nscamx)

  !> ivissa : comme ivisls sauf que sert au stockage de la viscosite au
  !>          pas de temps precedent
  integer, save ::          ivissa(nscamx)

  !> iclvfl : 0 : clipping des variances a zero
  !>          1 : clipping des variances a zero et a f(1-f)
  !>          2 : clipping des variances a max(zero,scamin) et scamax
  integer, save ::          iclvfl(nscamx)

  !> iscavr : numero du scalaire associe a la variance ou zero
  !>          si le scalaire n'est pas une variance
  integer, save ::          iscavr(nscamx)

  !> iscasp : 0 : le scalaire associe n est pas une espece
  !>          1 : le scalaire associe est une espece
  integer, save ::          iscasp(nscamx)

  !> scamin, scamax : min et max pour clipping des scalaires
  !>                  on ne clippe que si scamin < scamax
  double precision, save :: scamin(nscamx), scamax(nscamx)

  !> visls0 : viscosite des scalaires si constante
  double precision, save :: visls0(nscamx)

  !> sigmas : prandtl des scalaires
  double precision, save :: sigmas(nscamx)

  !> molar fraction for multi-species scalars
  !> \remarks
  !> wmolsp(0) is associated to the deduced species.
  double precision, save :: wmolsp(0:nscamx)

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

  !> \}

  !=============================================================================

  interface

    !---------------------------------------------------------------------------

    !> \cond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

    ! Interface to C function retrieving pointers to members of the
    ! global time step structure

    subroutine cs_f_time_step_get_pointers(nt_prev, nt_cur, nt_max,  &
                                           t_prev, t_cur, t_max)     &
      bind(C, name='cs_f_time_step_get_pointers')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(out) :: nt_prev, nt_cur, nt_max
      type(c_ptr), intent(out) :: t_prev, t_cur, t_max
    end subroutine cs_f_time_step_get_pointers

    !---------------------------------------------------------------------------

    !> \endcond DOXYGEN_SHOULD_SKIP_THIS

    !---------------------------------------------------------------------------

  end interface

  !=============================================================================

contains

  !=============================================================================

  !> \brief Initialize Fortran time step API.
  !> This maps Fortran pointers to global C structure members.

  subroutine time_step_init

    use, intrinsic :: iso_c_binding
    implicit none

    ! Local variables

    type(c_ptr) :: c_ntpabs, c_ntcabs, c_ntmabs
    type(c_ptr) :: c_ttpabs, c_ttcabs, c_ttmabs

    call cs_f_time_step_get_pointers(c_ntpabs, c_ntcabs, c_ntmabs, &
                                     c_ttpabs, c_ttcabs, c_ttmabs)

    call c_f_pointer(c_ntpabs, ntpabs)
    call c_f_pointer(c_ntcabs, ntcabs)
    call c_f_pointer(c_ntmabs, ntmabs)

    call c_f_pointer(c_ttpabs, ttpabs)
    call c_f_pointer(c_ttcabs, ttcabs)
    call c_f_pointer(c_ttmabs, ttmabs)

  end subroutine time_step_init

  !=============================================================================

end module optcal
