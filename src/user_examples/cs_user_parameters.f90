!-------------------------------------------------------------------------------

!VERS

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_parameters.f90
!>
!> \brief User subroutines for input of calculation parameters (Fortran modules).
!>        These subroutines are called in all cases.
!>
!>  See \ref f_parameters for examples.
!>
!>   If the code_saturne GUI is used, this file is not required (but may be
!>   used to override parameters entered through the GUI, and to set
!>   parameters not accessible through the GUI).
!>
!>   Several routines are present in the file, each destined to defined
!>   specific parameters.
!>
!>   To modify the default value of parameters which do not appear in the
!>   examples provided, code should be placed as follows:
!>   - usipsu   for numerical and physical options
!>   - usipes   for input-output related options
!>
!>   As a convention, "specific physics" defers to the following modules only:
!>   pulverized coal, gas combustion, electric arcs.
!>
!>   In addition, specific routines are provided for the definition of some
!>   "specific physics" options.
!>   These routines are described at the end of this file and will be activated
!>   when the corresponding option is selected in the usppmo routine.
!-------------------------------------------------------------------------------

!===============================================================================

!> \brief User subroutine for selection of specific physics module

!> Define the use of a specific physics amongst the following:
!>   - combustion with gas / coal / heavy fuel oil
!>   - compressible flows
!>   - electric arcs
!>   - atmospheric modelling
!>   - radiative transfer
!>   - cooling towers modelling
!>
!>    Only one specific physics module can be activated at once.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ixmlpu        indicates if an XML file from the GUI is used  !
!______________________________________________________________________________!

subroutine usppmo &
 ( ixmlpu )

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use cstphy
use ppppar
use ppthch
use ppincl
use ppcpfu
use coincl
use radiat
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer ixmlpu

!< [usppmo]

!===============================================================================
! 1.  Choice for a specific physics
!===============================================================================

! --- cod3p: Diffusion flame with complete fast chemistry (3 points)
! ==========

!        if = -1   module not activated
!        if =  0   adiabatic model
!        if =  1   extended model with enthalpy source term

ippmod(icod3p) = -1

! --- coebu: Eddy-Break Up pre-mixed flame
! ==========

!        if = -1   module not activated
!        if =  0   reference Spalding model
!                   (adiabatic, homogeneous mixture fraction)
!        if =  1   extended model with enthalpy source term
!                   (homogeneous mixture fraction : perfect premix)
!        if =  2   extended model with mixture fraction transport
!                   (adiabatic, no variance of mixture fraction)
!        if =  3   extended model with enthalpy and mixture fraction transport
!                   (dilution, thermal losses, etc.)

ippmod(icoebu) = -1

! --- colwc: Libby-Williams pre-mixed flame
! ==========

!        if = -1   module not activated
!        if =  0   reference two-peak model with adiabatic condition
!        if =  1   extended two-peak model with enthapy source terms
!        if =  2   extended three-peak model, adiabatic
!        if =  3   extended three-peak model with enthalpy source terms
!        if =  4   extended four-peak model, adiabatic
!        if =  5   extended four-peak model with enthalpy source terms

ippmod(icolwc) = -1


! --- Soot model
! =================

!        if = -1   module not activated
!        if =  0   constant soot yield
!        if =  1   2 equations model of Moss et al.

isoot = 0

xsoot  = 0.1d0 ! (if isoot = 0 and only if the soot yield is not
               !  defined in the thermochemistry data file)
rosoot = 2000.d0 ! kg/m3

! --- cfuel: Heavy fuel oil combustion
! ==========

!        Progressive evaporation (temperature gap)
!        Char residue
!        Sulphur tracking

!        if = -1   module not activated
!        if = 0    module activated

ippmod(icfuel) = -1

! --- coal :
! ==========
!
!     Pulverized coal combustion
!        Description of granulometry
!        Assumption of diffusion flame around particles
!         (extension of 3-point fast chemistry "D3P")
!        Between a mixture of gaseous fuels (volatiles matters, CO from char
!                                            oxydation)
!            and a mixture of oxidisers (air and water vapor)
!        Enthalpy for both mix and solid phase are solved
!
!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying

ippmod(iccoal) = -1

! Activate the drift: 0 (no activation),
!                     1 (transported particle velocity)
!                     2 (limit drop particle velocity)

i_comb_drift = 1


! --- cpl3c: Pulverized coal with Lagrangian reciprocal approach
! ==========

!        Not recently tested... at least outdated, may be obsolete

!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying (NOT functional)

ippmod(icpl3c) = -1

! --- eljou: Joule effect
! ==========

!        if = -1   module not activated
!        if = 1    real potential
!        if = 2    complex potential
!        if = 3    real potential + Transfo
!        if = 4    complex potential + Transfo

ippmod(ieljou) = -1

! --- elarc: Electric arcs
! ==========

!        if = -1   module not activated
!        if = 1    electric potential
!        if = 2    electric potential and vector potential (hence 3D modelling)

ippmod(ielarc) = -1

! --- aeros: Cooling towers
! ==========

!        if = -1   module not activated
!        if = 0    no model (NOT functional)
!        if = 1    Poppe's model
!        if = 2    Merkel's model

ippmod(iaeros) = -1

! Radiative transfer module (iirayo)
!--------------------------
!        if = 0: not activated (Default)
!        if = 1: DOM
!        if = 2: approximation P1 method

iirayo = 1

! --- richards model
! ==========

!        if = -1   module not activated
!        if =  1   module activated

ippmod(idarcy) = -1

!===============================================================================
! 2.  Specific options related to herebefore modules
!===============================================================================

! These options are defined here at the moment, this might change in the future

! --- Enthalpy-Temperature conversion law (for gas combustion modelling)

!       if = 0   user-specified
!       if = 1   tabulated by JANAF (default)

indjon = 1

!===============================================================================
! 2.  Data file related to modules above
!===============================================================================

! Combustion

if (     ippmod(icod3p).ge.0                                          &
    .or. ippmod(icoebu).ge.0 .or. ippmod(icolwc).ge.0) then

  if (indjon.eq.1) then
    ficfpp = 'dp_C3P'
  else
    ficfpp = 'dp_C3PSJ'
  endif

endif

! Fuel combustion

if (ippmod(icfuel).ge.0) then
  ficfpp = 'dp_FUE'
endif

! Specific condensation modelling (used with ippmod(igmix) >= 0)

! wall condensation model
!   -1: not activated
!    0: condensation source terms
icondb = -1

! internal condensation model
!   -1: not activated
!    0: condensation source terms with metal structures
icondv = -1

!< [usppmo]

!----
! End
!----

return
end subroutine usppmo

!===============================================================================

!> \brief User subroutine for the input of additional user parameters.
!
!>  This subroutine allows setting parameters
!>  which do not already appear in the other subroutines of this file.
!>
!>  It is possible to add or remove parameters.
!>  The number of physical properties and variables is known here.
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp         number of active specific physics models
!______________________________________________________________________________!

subroutine usipsu &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use albase
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation
use post
use rotation
use atincl
use atsoil
use atchem
use atimbr
use sshaerosol
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer       ii, jj
integer       iscal, ifcvsl
integer       kturt, turb_flux_model

!===============================================================================

!< [usipsu]

! Calculation options (optcal)
! ============================

! In case of restart, read auxiliary restart file ileaux (= 1) or not (0).

! By default, this file is read, but it may be useful to deactivate
! its use when restarting after a preprocessing stage possibly leading
! to a different number of faces (such as simply joining meshes on
! a different architecture or optimization level or with different options).

! Writing of auxiliary restart files may also be deactivated using: iecaux = 0

ileaux = 0

! --- Algorithm to take into account the thermodynamical pressure variation in time
!     (not used by default except if idilat = 3)

!     by default:
!     ----------
!      - the thermodynamic pressure (pther) is initialized with p0 = p_atmos
!      - the maximum thermodynamic pressure (pthermax) is initialized with -1
!        (no maximum by default, this term is used to model a venting effect when
!         a positive value is given by the user)
!      - a global leak can be set through a leakage surface sleak with a head
!      loss kleak of 2.9 (Idelcick)

ipthrm = 0

pthermax= -1.d0

sleak = 0.d0
kleak = 2.9d0


! --- Temperature or enthalpy

!   When used without specific physics, if we have chosen to solve in temperature
!     (that is if itherm = 1), the fluid temperature is considered to be in
!     degrees Kelvin by default (be careful for boundary conditions an expression
!     of physical properties depending on temperature)t.

!     If we wish for the fluid solver to work with a temperature in degrees Celsius,
!     we must set itpscl = 2.

!     This is recommended for Syrthes Coupling, but not recommended for the
!     radiative model, as it is a source of user errors in this case:
!     Indeed, the boundary conditions for the fluid temperature will then be
!     in degrees Celsius, while the boundary conditions for radiation in
!     cs_user_radiative_transfer_bcs must still be in Kelvin.

if (nmodpp.eq.0) then
  itpscl = 2
endif

!   If a USER scalar behaves like a temperature (relative to Cp):
!     we set the "is_temperature" keyword to 1.

if (nscaus.gt.0) then
  do ii = 1, nscaus
    call field_set_key_int(ivarfl(isca(ii)), kscacp, 1)
  enddo
endif

! --- Calculation (restart) with frozen velocity field (1 yes, 0 no)

iccvfg = 1


! --- Rotation/curvature correction for eddy-viscosity turbulence models
!      0: deactivated
!      1: activated


irccor = 1

! --- rotation of the reference frame (omega in rad/s)

!       If the rotation is not nul, then
!          icorio = 0: rotation is taken into account by rotating the mesh
!                      (simulation in the absolute frame)
!                 = 1: rotation is taken into account by Coriolis source terms
!                      (simulation in the relative frame)


icorio = 0

call rotation_define(0.d0, 0.d0, 0.d0,  &    ! rotation vector
                     0.d0, 0.d0, 0.d0)       ! invariant point


! --- Variable diffusivity field id (ifcvsl>=0) or constant
!     diffusivity (ifcvsl=-1) for the thermal scalar and USER scalars.

!     With ifcvsl = 0, the field will be added automatically, and later calls to
!       field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
!       will return its id.
!     With ifcvsl > 0, the id of an existing, predifined field is given. This
!       may allow sharing a diffusivity between multiple scalars.

!     For user scalars iscal which represent the variance of another user
!       scalar, the diffusivity of the variance of a scalar is assumed to
!       have the same behavior as the diffusivity of this scalar,
!       so values set here will be ignored.

!     For non-user scalars relative to specific physics (coal, combustion,
!       electric arcs: see usppmo) implicitly defined in the model,
!       the diffusivity should not be modified here.

!     Caution:   complete cs_user_physical_properties with the law defining the diffusivity
!     ========   if and only if ifcvsl = 0 has been set here.


! For thermal scalar
if (ippmod(icompf).ge.0) then
  ifcvsl = -1
  call field_set_key_int(ivarfl(isca(itempk)), kivisl, ifcvsl)
else if (iscalt.gt.0) then
  ifcvsl = -1
  call field_set_key_int(ivarfl(isca(iscalt)), kivisl, ifcvsl)
endif

do iscal = 1, nscaus
  if (iscavr(iscal).le.0) then
    ifcvsl = -1
    call field_set_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
  endif
enddo


! --- Variable density field id (ifcvsl>=0) or bulk
!     density (ifcvsl=-1) for USER scalars.

!     With ifcvsl = 0, the field will be added automatically, and later calls to
!       field_get_key_int(ivarfl(isca(iscal)), kromsl, ifcvsl)
!       will return its id.
!     With ifcvsl > 0, the id of an existing, predifined field is given. This
!       may allow sharing a density between multiple scalars.

!     For user scalars iscal which represent the variance of another user
!       scalar, the density of the variance of a scalar is assumed to
!       have the same behavior as the density of this scalar,
!       so values set here will be ignored.

!     Caution:   complete cs_user_physical_properties with the law defining the density
!     ========   if and only if ifcvsl = 0 has been set here.

do iscal = 1, nscaus
  if (iscavr(iscal).le.0) then
    ifcvsl = -1
    call field_set_key_int(ivarfl(isca(iscal)), kromsl, ifcvsl)
  endif
enddo


! --- Turbulent flux model u'T' for the scalar T
!     Algebraic Model
!      0  SGDH
!      10 GGDH
!      11 EB-GGDH (Elliptic Blending)
!      20 AFM
!      21 EB-AFM (Elliptic Blending)
!     Model with transport equations
!      30 DFM
!      31 EB-DFM (Elliptic Blending)

call field_get_key_id('turbulent_flux_model', kturt)

! GGDH for thermal scalar:
if (iscalt.gt.0) then
  turb_flux_model = 10
  call field_set_key_int(ivarfl(isca(iscalt)), kturt, turb_flux_model)
endif

! GGDH for all the scalars:
do jj = 1, nscaus
  turb_flux_model = 10
  call field_set_key_int(ivarfl(isca(jj)), kturt, turb_flux_model)
enddo

! Error estimators for Navier-Stokes (non-frozen velocity field)

! We recommend running a calculation restart on a few time steps
! with the activation of the most interesting of those.
! (=2 to activate, =0 to deactivate).

iescal(iescor) = 2   ! div(rho u) -Gamma
iescal(iestot) = 2   ! resolution precision for the momentum

! ALE (Arbitrary Lagrangian Eulerian) related options
!====================================================

! Number of iterations for fluid initialization. Contrary to ntmabs,
! nalinf is not an absolute iteration number, meaning that in case of
! restart calculation nalinf corresponds to the number of iterations
! for fuid initialization beginning from the first current iteration of
! the calculation restart. In general nalinf = 0 in that case.

nalinf = 75


! Maximum number of iterations in case of implicit Fluid Structure Coupling
! with structural calculations (internal and/or external
! (i.e. using code_aster)).
! nalimx = 1, in case of explicit FSI algorithm.

nalimx = 15


! Relative precision of sub-cycling Fluid Structure Coupling algorithm.

epalim = 1.d-5


!< [usipsu]

!===============================================================================
! Initialize non-standard calculation options for the atmospheric version.
!===============================================================================

!< [usati1]
!  -----------------------------------------------------------------------------
!  Atmospheric imbrication on large scale meteo (atimbr module)
!  -----------------------------------------------------------------------------
!
! --------------------------------------------------------------
! activation flag
! --------------------------------------------------------------
imbrication_flag    = .false.
imbrication_verbose = .false.

! ------------------------------------------------------------------------------
! flags for activating the cressman interpolation for the boundary conditions
! ------------------------------------------------------------------------------
cressman_u     = .true.
cressman_v     = .true.
cressman_tke   = .true.
cressman_eps   = .true.
cressman_theta = .true.
cressman_qw    = .true.
cressman_nc    = .true.

! --------------------------------------------------------------
! numerical parameters for the cressman interpolation formulas
! --------------------------------------------------------------
horizontal_influence_radius = 8500.d0
vertical_influence_radius = 100.d0

! --------------------------------------------------------------

! ifilechemistry: choice to read (=1,2,3,4, according to the scheme)
! or not (0) a concentration profile file
! if ichemistry>0 ifilechemistry is automaticaly set to ichemistry
ifilechemistry = 0

! Change the name of the chemistry concentration profile
call atmo_set_chem_conc_file_name('chem_01_01_2001')

! Change the name of the aerosol concentration profile
call atmo_set_aero_conc_file_name('aero_01_01_2001')

! isepchemistry: splitted (=1) or semi-coupled (=2, pu-sun)
! resolution of chemistry.
! Splitted (=1) mandatory for aerosols.
! Semi-coupled (=2) by default.
isepchemistry = 1

! dtchemmax: maximal time step (s) for chemistry resolution
dtchemmax = 10.0d0

! computation / storage of downward and upward infrared radiative fluxes
irdu = 1

!< [usati1]

!----
! Formats
!----

return
end subroutine usipsu


!===============================================================================

!> \brief User subroutine for the input of additional user parameters for
!>        input/output.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nmodpp       number of active specific physics models
!______________________________________________________________________________!

subroutine usipes &
 ( nmodpp )

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use field
use parall
use period
use post
use ppppar
use ppthch
use ppincl
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii, f_id

type(var_cal_opt) :: vcopt

!===============================================================================

!>     This subroutine allows setting parameters
!>     which do not already appear in the other subroutines of this file.
!>
!>     It is possible to add or remove parameters.
!>     The number of physical properties and variables is known here.

!===============================================================================

!===============================================================================
! 1. Logging
!===============================================================================

! Frequency of log output

ntlist = 1


! Log verbosity

do ii = 1, nvar
  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)
  vcopt%iwarni = 1
  call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)
enddo

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)
vcopt%iwarni = 2
call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt)
vcopt%iwarni = 2
call field_set_key_struct_var_cal_opt(ivarfl(iu), vcopt)

!===============================================================================
! 2. Definition of deformable structure time plots
!===============================================================================

! structures output step

nthist = 1
frhist = -1.d0

tplfmt = 1 ! time plot format (1: .dat, 2: .csv, 3: both)

!===============================================================================
! 3. Fine control of variables output
!===============================================================================

! Per variable output control.
! More examples are provided in cs_user_parameters-output.f90

! User scalar variables.

if (isca(1).gt.0.and.nscaus.ge.1) then
  f_id = ivarfl(isca(1))
  call field_set_key_str(f_id, keylbl, 'Scalar 1')
  call field_set_key_int(f_id, keyvis, POST_ON_LOCATION + POST_MONITOR)
  call field_set_key_int(f_id, keylog, 1)
endif

if (isca(2).gt.0.and.nscaus.ge.2) then
  f_id = ivarfl(isca(2))
  call field_set_key_str(f_id, keylbl, 'Scalar 2')
  call field_set_key_int(f_id, keyvis, POST_ON_LOCATION + POST_MONITOR)
  call field_set_key_int(f_id, keylog, 1)
endif


return
end subroutine usipes

!===============================================================================
! Purpose:
! -------
!
!> 1. Additional Calculation Options
!>    a. Density Relaxation
!>
!> 2. Physical Constants
!>    a.Dynamic Diffusion Coefficient
!>    b.Constants of the chosen model (EBU, Libby-Williams, ...)
!
!> This routine is called:
!>
!>
!>  - Eddy Break Up pre-mixed flame
!>  - Diffusion flame in the framework of ``3 points'' rapid complete chemistry
!>  - Libby-Williams pre-mixed flame
!>  - Lagrangian module coupled with pulverized coal:
!>    Eulerian combustion of pulverized coal and
!>    Lagrangian transport of coal particles
!>  - Pulverised coal combustion
!>  - Fuel (oil) combustion
!
!===============================================================================

subroutine cs_user_combustion

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use radiat

!===============================================================================

implicit none

!< [cs_user_combustion]

!===============================================================================
! 1. Additional Calculation Options
!===============================================================================

! --- Kinetic model for CO <=> CO2

!     if = 0  unused (maximal conversion in turbulent model)
!     if = 1  transport of CO2 mass fraction
!     if = 2  transport of CO mass fraction (coal and fuel only)

ieqco2 = 0

! --- Density Relaxation
!     RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0

!===============================================================================
! 2. Physical Constants
!===============================================================================

! diftl0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5


! -----------------------------------------------------------------------------
! 2.1 For 3 points combusution model ONLY
! -----------------------------------------------------------------------------

! Reference temperature for fuel and oxydant (K)
tinfue = 436.d0
tinoxy = 353.d0


! -----------------------------------------------------------------------------
! 2.2 For EBU-model ONLY
! -----------------------------------------------------------------------------

! cebu: EBU-model constant
cebu   = 2.5d0


! -----------------------------------------------------------------------------
! 2.3 For Libby-Williams model ONLY
! -----------------------------------------------------------------------------

! Reference velocity
vref = 60.d0
! Reference length scale
lref = 0.1d0
! Activation Temperature
ta   = 0.2d5
! Cross-over Temperature (combustion of propane)
tstar= 0.12d4

!< [cs_user_combustion]

!----
! End
!----

return
end subroutine cs_user_combustion


!===============================================================================

!> \brief User subroutine.

!> Initialize non standard options for the compressible flow scheme such
!> as the variability of the thermal conductivity and the volume viscosity.
!> Their values can be given in the subroutine \ref uscfx2 .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine uscfx1

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use field

!===============================================================================

implicit none

! Arguments


! Local variables

integer :: ifcvsl

!===============================================================================

!===============================================================================

!< [uscfx1]

!===============================================================================
! 1. Properties options
!===============================================================================

! --> Molecular thermal conductivity
!       constant  : ifcvsl = -1
!       variable  : ifcvsl = 0

ifcvsl = -1
call field_set_key_int(ivarfl(isca(itempk)), kivisl, ifcvsl)

! --> Volumetric molecular viscosity
!       iviscv = -1 : uniform  in space and constant in time
!              =  0 : variable in space and time

iviscv = -1

!< [uscfx1]

!----
! End
!----

return
end subroutine uscfx1

!===============================================================================


!> \brief User subroutine.
!>
!> Set values for the reference volumic viscosity, the reference
!> conductivity and the molar mass for compressible flow.
!>
!> Initialize non standard options for the compressible flow scheme such
!> as the hydrostatic equilibrium.
!>
!> In addition to options set in the user subroutine \ref uscfx1 (or in
!> the GUI): this subroutine allows to set a switch to indicate if the
!> molecular viscosity is constant, its values being given in the user
!> subroutine \ref usipsu .

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine uscfx2

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

! Local variables

!< [uscfx2]

!===============================================================================
! 1. Physical properties
!===============================================================================

! --> Molecular viscosity
!       constant  : ivivar = 0
!       variable  : ivivar = 1

ivivar = 0

! --> Reference molecular thermal conductivity
!       visls0 = lambda0 (molecular thermal conductivity, W/(m K))

!     Reference conductivity must be strictly positive
!      (set a realistic value here even if conductivity is variable)

call field_set_key_double(ivarfl(isca(itempk)), kvisl0, 3.d-2)

!       If the molecular thermal conductivity is variable, its values
!         must be provided in the user subroutine 'cs_user_physical_properties'

! --> Volumetric molecular viscosity

!       Reference volumetric molecular viscosity

!       viscv0 = kappa0  (volumetric molecular viscosity, kg/(m s))

viscv0 = 0.d0

!       If the volumetric molecular viscosity is variable, its values
!         must be provided in the user subroutine 'cs_user_physical_properties'

! --> Molar mass of the gas (kg/mol)

!       For example with dry air, xmasml is around 28.8d-3 kg/mol

xmasmr = 0.028966

! --> Hydrostatic equilibrium

!       Specify if the hydrostatic equilibrium must be accounted for
!         (yes = 1 , no = 0)

icfgrp = 1

!< [uscfx2]

!----
! End
!----

return
end subroutine uscfx2

!===============================================================================

!> \brief Definition of cooling tower model and exchange zones

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine cs_user_cooling_towers

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use ctincl

!===============================================================================

implicit none

!===============================================================================

!===============================================================================
! 1. Parameters for prescibed temperature difference
!===============================================================================

! Air and liquid properties
!< [cs_user_cooling_towers]
cp_a    = 1006.0d0
cp_v    = 1831.0d0
cp_l    = 4179.0d0
hv0    = 2501600.0d0
rho_l   = 997.85615d0
viscl0 = 1.765d-5
lambda_l = 0.02493d0
humidity0 = 0.d0
droplet_diam = 0.005d0
!< [cs_user_cooling_towers]

!----
! End
!----

return
end subroutine cs_user_cooling_towers

!===============================================================================

!> \brief User routine for definition of computation parameters dealing
!>        with Darcy module

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine user_darcy_ini1

!===============================================================================
! Module files
!===============================================================================

use entsor
use darcy_module

!===============================================================================

implicit none

!===============================================================================

!< [user_darcy_ini1]

darcy_anisotropic_permeability = 0 ! permeability : 0 isotrop, 1 anisotrop

darcy_anisotropic_dispersion = 0 ! dispersion : 0 isotrop, 1 anisotrop

darcy_unsteady = 0 ! 0 steady flow, 1 unsteady flow

darcy_convergence_criterion = 0 ! convergence criterion of Newton scheme:
                                ! 0, over pressure, 1, over velocity

darcy_gravity = 0 ! gravity is taken into account : 0 no, 1 yes

!< [user_darcy_ini1]

!----
! End
!----

return

end subroutine user_darcy_ini1


!===============================================================================
