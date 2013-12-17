!-------------------------------------------------------------------------------

!VERS

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

! Purpose:
! -------

! User subroutines for input of calculation parameters (Fortran modules).
!   These subroutines are called in all cases.

! If the Code_Saturne GUI is used, this file is not required (but may be
!   used to override parameters entered through the GUI, and to set
!   parameters not accessible through the GUI).

! Several routines are present in the file, each destined to defined
!   specific parameters.

! To modify the default value of parameters which do not appear in the
!   examples provided, code should be placed as follows:
!   - usipsu   for numerical and physical options
!   - usipes   for input-output related options

! As a convention, "specific physics" defers to the following modules only:
!   pulverized coal, gas combustion, electric arcs.

! In addition, specific routines are provided for the definition of some
!   "specific physics" options.
!   These routines are described at the end of this file and will be activated
!   when the corresponding option is selected in the usppmo routine.

!-------------------------------------------------------------------------------


!===============================================================================


subroutine usppmo &
!================
 ( ixmlpu )


!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Define the use of a specific physics amongst the following:
!      - combustion with gas / coal / heavy fuel oil
!      - compressible flows
!      - electric arcs
!      - atmospheric modelling
!      - cooling towers modelling

!    Only one specific physics module can be activated at once.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ixmlpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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

!===============================================================================

implicit none

! Arguments

integer ixmlpu


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!===============================================================================

if (1.eq.1) then
  return
endif

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  Choice for a specific physics
!===============================================================================

! --- cod3p: Diffusion flame with complete fast chemistry (3 points)
! ==========

!        if = -1   module not activated
!        if =  0   adiabatic model
!        if =  1   extended model with enthalpy source term

if (ixmlpu.eq.0) then

  ippmod(icod3p) = -1

endif

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

if (ixmlpu.eq.0) then

  ippmod(icoebu) = -1

endif

! --- colwc: Libby-Williams pre-mixed flame
! ==========

!        if = -1   module not activated
!        if =  0   reference two-peak model with adiabatic condition
!        if =  1   extended two-peak model with enthapy source terms
!        if =  2   extended three-peak model, adiabatic
!        if =  3   extended three-peak model with enthalpy source terms
!        if =  4   extended four-peak model, adiabatic
!        if =  5   extended four-peak model with enthalpy source terms

if (ixmlpu.eq.0) then

  ippmod(icolwc) = -1

endif


! --- Soot model
! =================

!        if = -1   module not activated
!        if =  0   constant fraction of fuel Xsoot
!        if =  1   2 equations model of Moss et al.

if (.false.) then
  isoot = 0

  xsoot  = 0.1d0 ! ( if isoot = 0 )
  rosoot = 2000.d0 ! kg/m3
endif



! --- cfuel: Heavy fuel oil combustion
! ==========

!        Progressive evaporation (temperature gap)
!        Char residue
!        Sulphur tracking

!        if = -1   module not activated
!        if = 0    module activated

if (ixmlpu.eq.0) then

  ippmod(icfuel) = -1

endif

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

!     Activate the drift for the coal scalars
!        i_coal_drift = 0  not activated (default)
!        i_coal_drift = 1  activated

if (ixmlpu.eq.0) then

  ippmod(iccoal) = -1

endif

! --- cpl3c: Pulverized coal with Lagrangian reciprocal approach
! ==========

!        Not recently tested... at least outdated, may be obsolete

!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying (NOT functional)

if (ixmlpu.eq.0) then

  ippmod(icpl3c) = -1

endif

! --- compf: Compressible flows
! ==========

!        if = -1   module not activated
!        if = 0    module activated

if (ixmlpu.eq.0) then

  ippmod(icompf) = -1

endif

! --- eljou: Joule effect
! ==========

!        if = -1   module not activated
!        if = 1    Potentiel reel
!        if = 2    Potentiel complexe
!        if = 3    Potentiel reel     + CDL Transfo
!        if = 4    Potentiel complexe + CDL Transfo

if (ixmlpu.eq.0) then

  ippmod(ieljou) = -1

endif

! --- elarc: Electric arcs
! ==========

!        if = -1   module not activated
!        if = 1    electric potential
!        if = 2    electric potential and vector potential (hence 3D modelling)

if (ixmlpu.eq.0) then

  ippmod(ielarc) = -1

endif

! --- atmos: Atmospheric flows
! ==========

!        if = -1   module not activated
!        if = 0    standard modelling
!        if = 1    dry atmosphere
!        if = 2    humid atmosphere (experimental)

if (ixmlpu.eq.0) then

  ippmod(iatmos) = -1

endif

! --- aeros: Cooling towers
! ==========

!        if = -1   module not activated
!        if = 0    no model (NOT functional)
!        if = 1    Poppe's model
!        if = 2    Merkel's model

if (ixmlpu.eq.0) then

  ippmod(iaeros) = -1

endif

!===============================================================================
! 2.  Specific options related to herebefore modules
!===============================================================================

! These options are defined here at the moment, this might change in the future

! --- Enthalpy-Temperature conversion law (for gas combustion modelling)

!       if = 0   user-specified
!       if = 1   tabulated by JANAF (default)

if (ixmlpu.eq.0) then

  indjon = 1

endif

! --- Kinetic model for NOx formation

!         Only compatible with heavy fuel oil combustion

!         if = 0  unused
!         if = 1  activated

if (ixmlpu.eq.0) then

  ieqnox = 1

endif

! --- NOx Formation Features
!
!     imdnox = 0: - HCN is the only intermediate nitrogen species
!                   liberated during the devolatilisation process.
!                 - HCN is the only intermediate nitrogen species
!                   liberated during char combustion.
!                 - Constant ratio of the nitrogen mass liberated
!                   during devolatilisation and the nitrogen mass
!                   remaining in char.
!            = 1 (only if iccoal = 0(1))
!                 - HCN and NH3 are the intermediate nitrogen
!                   species liberated during the devolatilisation
!                   process.
!                 - HCN and NO are the intermediate nitrogen species
!                   liberated during char combustion.
!                 - Temperature dependent ratios of the nitrogen
!                   mass liberated during devolatilisation and the
!                   nitrogen mass remaining in char.
!                 - Activation of Reburning kinetics is possible.
!

if (ixmlpu.eq.0) then

  imdnox = 0

endif


! --- Reburning

!     Only compatible with ieqnox = 1 and imdnox = 1 and iccoal = 0(1)

!     if = 0  unused
!     if = 1  Model of Chen et al.
!     if = 2  Model of Dimitriou et al.


if (ixmlpu.eq.0) then

  irb = 0

endif


! --- Kinetic model for CO <=> CO2

!         Compatible with coal and heavy fuel oil combustion

!         if = 0  unused (maximal conversion in turbulent model)
!         if = 1  transport of CO2 mass fraction
!         if = 2  transport of CO mass fraction

if (ixmlpu.eq.0) then

  ieqco2 = 0

endif

! --- Heteregoneous combustion by CO2

!         Needs the activation of the CO2 transport equation
!         Account for the reaction between char and CO2: C(s) + CO2 => 2 CO

!         if = 0  unused
!         if = 1  activated

if (ixmlpu.eq.0) then

  ihtco2 = 0

endif

!===============================================================================
! 2.  Data file related to modules above
!===============================================================================

if (ixmlpu.eq.0) then

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

  ! Electric arcs

  if (ippmod(ielarc).ge.1) then

    ficfpp = 'dp_ELE'

  endif

  ! Joule effect

  if (ippmod(ieljou).eq.1 .or. ippmod(ieljou).eq.2) then
    ficfpp = 'dp_ELE'
  else if (ippmod(ieljou).eq.3 .or. ippmod(ieljou).eq.4) then
    ficfpp = 'dp_transfo'
  endif

  ! Atmospheric flows

  if (ippmod(iatmos).ge.0) then
    ficmet = 'meteo'
  endif

endif

!----
! End
!----

return
end subroutine usppmo


!===============================================================================


subroutine usipph &
!================
 ( ixmlpu, nfecra , iturb , irccor , itherm, icp )


!===============================================================================
! Purpose:
! --------

! User subroutine for input of parameters.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ixmlpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! iturb            ! ia ! <-> ! turbulence model                               !
! irccor           ! ia ! <-> ! flag for rotation/curvature correction or not  !
! itherm           ! ia ! <-> ! thermal model                                  !
! icp              ! ia ! <-> ! flag for uniform Cp or not                     !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

! No module should appear here


!===============================================================================

implicit none

! Arguments

integer ixmlpu, nfecra
integer iturb, irccor, itherm, icp

! Local variables

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (ixmlpu.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''usipph'' must be completed',/,       &
'@       in file cs_user_parameters.f90',/,                       &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================


!     In this subroutine, only the parameters which already appear may

!       be set, to the exclusion of any other.
!               ================


!     If we are not using the Code_Saturne GUI:

!       All the parameters which appear in this subroutine must be set.
!       ===


!     If we are using the Code_Saturne GUI:

!       parameters protected by a test of the form:

!       if (ixmlpu.eq.0) then
!          ...
!       endif

!       should already have been defined using the GUI, so only
!       experts should consider removing the test and adapting them here.

!===============================================================================

! --- Turbulence
!       0: Laminar
!      10: Mixing length
!      20: k-epsilon
!      21: k-epsilon (linear production)
!      30: Rij-epsilon, (standard LRR)
!      31: Rij-epsilon (SSG)
!      32: Rij-epsilon (EBRSM)
!      40: LES (Smagorinsky)
!      41: LES (Dynamic)
!      42: LES (WALE)
!      50: v2f (phi-model)
!      51: v2f (BL-v2/k)
!      60: k-omega SST
!      70: Spalart Allmaras
!  For 10, contact the development team before use

if (ixmlpu.eq.0) then

  iturb = 21

endif

! --- Rotation/curvature correction for eddy-viscosity turbulence models
!      0: deactivated
!      1: activated

if (.false.) then

  irccor = 1

endif

! --- Thermal model
!      0: none
!      1: temperature
!      2: enthalpy
!      3: total energy (only for compressible module)
!
!  For temperature, the temperature scale may be set later using itpscl
!  (1 for Kelvin, 2 for Celsius).
!
!  When using specific physics, this value is set automatically by the physics model.

if (ixmlpu.eq.0) then

  itherm = 0

endif

! --- Variable specific heat (ICP=1) or not (ICP=0)

!     Should be set only if specific physics (coal, combustion, electric arcs)
!       ARE NOT activated.

!     For these specific physics, ICP MUST NOT BE MODIFIED here, and the
!       following options are forced:
!          coal and combustion: constant Cp;
!          electric arcs:       variable Cp.

!     Caution:    complete usphyv with the law defining Cp
!     =========   if and only if variable Cp has been selected here
!                 (with icp=1)

if (ixmlpu.eq.0) then

  icp = 0

endif

!----
! Formats
!----


return
end subroutine usipph


!===============================================================================


subroutine usinsc &
!================

 ( ixmlpu, nfecra , nscaus )


!===============================================================================
! Purpose:
! -------

! User subroutine for input of the number of user scalars.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ixmlpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! nscaus           ! i  ! <-> ! number of user scalars                         !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================


! No module should appear here


!===============================================================================

implicit none

! Arguments

integer ixmlpu, nfecra
integer nscaus

! Local variables


!===============================================================================

!     In this subroutine, only the parameters which already appear may

!       be set, to the exclusion of any other.
!               ================


!     If we are not using the Code_Saturne GUI:

!       All the parameters which appear in this subroutine must be set.
!       ===


!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================

! --- Number of USER scalars (thermal or not).
!       These scalars come in addition to the following "basic" scalars
!       (which are naturally included in the model):
!        - pressure
!        - turbulent variables
!        - nscapp scalars introduced by an active combustion, coal,
!          or electric arc module.

!     Thus, for a calculation with no specific physics, the user scalars
!       may for example be:
!        - temperature or enthalpy,
!        - mass fractions of transported scalars
!        - the variance of another user scalar

!     The maximum number of scalars is defined by 'nscamx' in paramx;
!       it is the maximum admissible value for: nscaus + nscapp.


!     Set nscaus = 0 if there is no user scalar.

if (ixmlpu.eq.0) then

  nscaus = 0

endif

!----
! Formats
!----


return
end subroutine usinsc


!===============================================================================


subroutine usipsc &
!================

 ( nscmax, nscaus, ixmlpu, nfecra, iscavr, ivisls )


!===============================================================================
! Purpose:
! -------

! User subroutine for the input of parameters depending on the
!   number of user scalars.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nscmax           ! i  ! <-- ! maximum number of scalars                      !
! nscaus           ! i  ! <-- ! number of user scalars                         !
! ixmlpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! iscavr(nscmax)   ! ia ! <-- ! associated scalar number for variance scalars  !
! ivisls(nscmax)   ! ia ! <-> ! uniform scalar diffusivity flag                !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

! No module should appear here


!===============================================================================

implicit none

! Arguments

integer nscmax, nscaus, ixmlpu, nfecra
integer iscavr(nscmax), ivisls(nscmax)

! Local variables

integer iscal

!===============================================================================

!     In this subroutine, only the parameters which already appear may

!       be set, to the exclusion of any other.
!               ================


!     If we are not using the Code_Saturne GUI:

!       All the parameters which appear in this subroutine must be set.
!       ===


!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================

! --- Variance of a USER scalar:
!     If we wish a user scalar j to represent the variance of a
!       user scalar k, we set
!       iscavr(j) = k.
!     The values taken by iscavr are thus naturally greater or equal to 1
!       and less than or equal to the total number of scalars.
!       So, if we set iscavr(j) = k, we must have
!       0 < j < nscaus+1, 0< k < nscaus+1 and j different from k.

!     For example for user scalar 3 to be the variance of user scalar 2,
!       we set:
!       iscavr(3) = 2
!       with nscaus at least equal to 3.

!     Do not intervene if you do not wish to explicitly include the
!       variance of a user scalar in the simulation.

!     For non-user scalars relative to specific physics (coal, combustion,
!       electric arcs: see usppmo) implicitly defined in the model,
!       the corresponding information is given automatically, and
!       iscavr should not be modified.

if (.false.) then
  iscavr(3) = 2
endif

! --- Variable diffusivity (ivisls=1) or constant diffusivity (ivisls=0) for
!       each USER scalar, EXCEPT those which represent the variance
!       of another.

!     For user scalars iscal which represent the variance of another user
!       scalar, we do not set ivisls(iscal) here.
!       This is the purpose of the test on iscavr(ISCAL) in the example below.
!       Indeed, the diffusivity of the variance of a scalar is assumed to
!       have the same behavior as the diffusivity of this scalar.

!     For non-user scalars relative to specific physics (coal, combustion,
!       electric arcs: see usppmo) implicitly defined in the model,
!       the corresponding information is given automatically, and
!       ivisls should not be modified here.

!     Caution:    complete usphyv with the law defining the diffusivity
!     =========   if and only if ivisls = 1 has been set here.

if (.false.) then

  do iscal = 1, nscaus

    ! For user scalars which do not represent the variance of another scalar
    if (iscavr(iscal).le.0) then

      ivisls(iscal) = 0

    endif

  enddo

endif

!----
! Formats
!----


return
end subroutine usipsc


!===============================================================================


subroutine usipgl &
!================

 ( nesmax,                                                        &
   iespre, iesder, iescor, iestot,                                &
   ixmlpu, nfecra,                                                &
   idtvar, ipucou, idilat, iphydr, ialgce , iescal )


!===============================================================================
! Purpose:
! -------

! User subroutine for the setting of global parameters.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nesmax           ! i  ! <-- ! maximum number of error estimators             !
! iespre           ! i  ! <-- ! number of the prediction error estimator       !
! iesder           ! i  ! <-- ! number of the derivative error estimator       !
! iescor           ! i  ! <-- ! number of the correction error estimator       !
! iestot           ! i  ! <-- ! number of the total error estimator            !
! ixmlpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! idtvar           ! i  ! --> ! variable time step flag                        !
! ipucou           ! i  ! --> ! reinforced u-p coupling flag                   !
! idilat           ! i  ! --> ! algorithm with density variation in time       !
! iphydr           ! i  ! --> ! flag for handling of the equilibrium between   !
!                  !    !     ! the pressure gradient and the gravity and      !
!                  !    !     ! head-loss terms                                !
! ialgce           ! i  ! <-- ! option for the method of calculation of        !
!                  !    !     ! cell centers                                   !
! iescal(nesmax)   ! ia ! <-- ! flag for activation of error estimators for    !
!                  !    !     ! Navier-Stokes                                  !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================


! No module should appear here


!===============================================================================

implicit none

! Arguments

integer nesmax
integer iespre, iesder, iescor, iestot
integer ixmlpu, nfecra
integer idtvar, ipucou, idilat, iphydr, ialgce
integer iescal(nesmax)

! Local variables

!===============================================================================

!     In this subroutine, only the parameters which already appear may

!       be set, to the exclusion of any other.
!               ================


!     If we are not using the Code_Saturne GUI:

!       All the parameters which appear in this subroutine must be set.
!       ===


!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================

! --- Time step  (0 : uniform and constant
!                 1 : variable in time, uniform in space
!                 2 : variable in time and space
!                -1 : steady algorithm)

if (.false.) then
  idtvar = 0
endif

! --- Velocity/pressure coupling (0 : classical algorithm,
!                                 1 : transient coupling)

if (.false.) then
  ipucou = 0
endif

! Algorithm to take into account the density variation in time
!
!     idilat = 0 : boussinesq algorithm with constant density (not available)
!              1 : dilatable steady algorithm (default)
!              2 : dilatable unsteady algorithm
!              3 : low-Mach algorithm
!

if (.false.) then
  idilat = 1
endif

! --- Handling of hydrostatic pressure
!     iphydr = 0 : ignore hydrostatic pressure (by default)
!              1 : with hydrotatic pressure computation to handle the balance
!                  between the pressure gradient and source terms (gravity and
!                  head losses)
!              2 : with hydrostatic pressure computation to handle the imbalance
!                  between the pressure gradient and gravity source term

if (.false.) then
  iphydr = 1
endif

! --- Estimators for Navier-Stokes (non-frozen velocity field)
!     We recommend running a calculation restart on a few time steps
!       with the activation of the most interesting of those.
!        (=2 to activate, =0 to deactivate).

if (.false.) then
  iescal(iescor) = 2   ! div(rho u) -Gamma
  iescal(iestot) = 2   ! resolution precision for the momentum
endif

!----
! Formats
!----

return
end subroutine usipgl


!===============================================================================


subroutine usipsu &
!================

 ( nmodpp )


!===============================================================================
! Purpose:
! -------

! User subroutine for the input of additional user parameters.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nmodpp           ! i  ! <-- ! number of active specific physics models       !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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
use mltgrd
use parall
use period
use ihmpre
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use elincl
use field

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii, jj, imom, kscmin, kscmax

!===============================================================================

!     This subroutine allows setting parameters

!       which do not already appear in the other subroutines of this file.


!     It is possible to add or remove parameters.


!     The number of physical properties and variables is known here.


!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================


! Calculation options (optcal)
! ============================

! In case of restart, read auxiliary restart file ileaux (= 1) or not (0).

! By default, this file is read, but it may be useful to deactivate
! its use when restarting after a preprocessing stage possibly leading
! to a different number of faces (such as simply joining meshes on
! a different architecture or optimization level or with different options).

! Writing of auxiliary restart files may also be deactivated using: iecaux = 0

if (.false.) then
  ileaux = 0
endif

! --- Duration
!       ntmabs = absolute number of the last time step required
!         if we have already run 10 time steps and want to
!         run 10 more, ntmabs must be set to 10 + 10 = 20

if (.false.) then
  ntmabs = 10
endif

! --- Reference time step
!     The example given below is probably not adapted to your case.

if (.false.) then
  dtref  = 0.01d0
endif

! --- Maximum time step: dtmax
!     Set a value base on characteristic values of your case.
!      otherwise, the code will use a multiple of dtref by default.
!     Example with
!        Ld: "dynamic" length (for example, the domain length)
!        Ud: characteristic flow velocity
!        Lt: thermal length (for example, the domain height gravity-wise)
!        Delta_rho/rho: relative density difference
!        g: gravity acceleration

!     dtmax = min(Ld/Ud, sqrt(Lt/(g.Delta_rho/rho)))


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
!     usray2 must still be in Kelvin.

if (.false.) then

  if (nmodpp.eq.0) then
    itpscl = 2
  endif

endif

!   If a USER scalar behaves like a temperature (relative to Cp):
!     we set iscacp(isca) = 1.
!
!   Otherwise, we do not modify iscacp(isca)

if (.false.) then

  if (nscaus.gt.0) then
    do ii = 1, nscaus
      iscacp(isca(ii)) = 1
    enddo
  endif

endif

! --- Solver taking a pscalar porosity into account:
!       0 No porosity taken into account (Standard)
!       1 Porosity taken into account
!

if (.false.) then
  iporos = 1
endif

! --- Calculation (restart) with frozen velocity field (1 yes, 0 no)

if (.false.) then
  iccvfg = 1
endif

! --- Vortex method for inlet conditions in L.E.S.
!       (0: not activated,  1: activated)
!     The vortex method only regards the L.E.S. models
!     To use the vortex method, edit the 'usvort.f90' user file.

if (.false.) then

  if (itytur.eq.4) then
    ivrtex = 1
  endif

endif

! --- Convective scheme

!     blencv = 0 for upwind (order 1 in space, "stable but diffusive")
!            = 1 for centered/second order (order 2 in space)
!       we may use intermediate real values.
!       Here we choose:
!         for the velocity and user scalars:
!           an upwind-centered scheme with 100% centering (blencv=1)
!         for other variables
!           the default code value (upwind standard, centered in LES)

!     Specifically, for user scalars
!       if we suspect an excessive level of numerical diffusion on
!         a variable ivar representing a user scalar
!         iscal (with ivar=isca(iscal)), it may be useful to set
!         blencv(ivar) = 1.0d0 to use a second-order scheme in space for
!         convection. For temperature or enthalpy in particular, we
!         may thus choose in this case:
!          blencv(isca(iscalt)) = 1.0d0

!       For non-user scalars relative to specific physics (coal, combustion,
!         electric arcs: see usppmo) implicitly defined by the model,
!         the corresponding information is set automatically elsewhere:
!         we do not modify blencv here.

if (.false.) then

  blencv(iu) = 1.0d0
  blencv(iv) = 1.0d0
  blencv(iw) = 1.0d0
  if (nscaus.ge.1) then
    do ii = 1, nscaus
      blencv(isca(ii)) = 1.0d0
    enddo
  endif

endif

! --- Linear solver parameters (for each unknown)

!     iresol = -1:           default
!     iresol = 1000*ipol +j: ipol is the degree of the Neumann polynomial
!                            used for preconditioning,
!                            j = 0: conjugate gradient,
!                            j = 200: conjugate gradient, single reduction
!                            j = 1: Jacobi
!                            j = 2: bi-CgStab
!                            j = 3: GMRES

!     nitmax: maximum number of iterations for each unknown ivar
!     epsilo: relative precision for the solution of the linear system.

if (.false.) then

  iresol(iu) = 2
  iresol(iv) = 2
  iresol(iw) = 2
  if (nscaus.ge.1) then
    do ii = 1, nscaus
      iresol(isca(ii)) = 2
      nitmax(isca(ii)) = 5000
      epsilo(isca(ii)) = 1.d-6
    enddo
  endif

endif

! --- Algebraic multigrid parameters

! imgr = 0: no multigrid
! imgr = 1: algebraic multigrid

! Only available for pressure and purely diffusive variables.

! mltmmn = 300  ! mean number of cells under which merging takes place
! mltmgl = 500  ! global number of cells under which merging takes place
! mltmmr = 1    ! number of active ranks under which no merging is done
! mltmst = 4    ! number of ranks over which merging takes place
! mlttyp = 0    ! 0: loop over faces to coarsen in natural order
!               ! 1: loop over faces to coarsen in criteria order
!               ! 3: loop over faces to coarsen in Hilbert order

if (.false.) then
  imgr(ipr) = 1
endif

! --- Dynamic reconstruction sweeps to handle non-orthogonlaities
!     This parameter computes automatically a dynamic relax factor,
!     and can be activated for any variable.
!      - iswdyn(ipr) = 1: means that the last increment is relaxed
!      - iswdyn(ipr) = 2: means that the last two increments are used to
!                         relax
!     NB: when iswdyn is greater than 1, then the number of
!         non-orthogonality sweeps is increased to 20.

if (.false.) then
  iswdyn(ipr) = 1
endif

!=========================================================================

! --- Stabilization in turbulent regime

!     For difficult cases, a stabilization may be obtained by not
!     reconstructing the convective and diffusive flux for variables
!     of the turbulence model, that is
!       in k-epsilon: if (itytur.eq.2) then
!          ircflu(ik)   = 0 and ircflu(iep)  = 0
!       in Rij-epsilon: if (itytur.eq.3) then
!          ircflu(ir11) = 0,    ircflu(ir22) = 0,
!          ircflu(ir33) = 0,
!          ircflu(ir12) = 0,    ircflu(ir23) = 0,
!          ircflu(ir23) = 0,
!                                  and ircflu(iep)  = 0
!     (note that variable itytur is equal to iturb/10)

if (.false.) then

  if (itytur.eq.2) then
    ircflu(ik)   = 0
    ircflu(iep)  = 0
  endif

endif


! Physical constants (cstphy)
! ===========================

! --- gravity (g in m/s2, with the sign in the calculation coordinate axes).

if (.false.) then

  gx = 0.d0
  gy = 0.d0
  gz = 0.d0

endif

! --- rotation vector of the reference frame (omega in s-1)

!       If the rotation is not nul, then
!          icorio = 0: rotation is taken into account by rotating the mesh
!                      (simulation in the absolute frame)
!                 = 1: rotation is taken into account by Coriolis source terms
!                      (simulation in the relative frame)

if (.false.) then

  icorio = 0

  omegax = 0.d0
  omegay = 0.d0
  omegaz = 0.d0

endif

! --- Reference fluid properties

!       ro0        : density in kg/m3
!       viscl0     : dynamic viscosity in kg/(m s)
!       cp0        : specific heat in J/(Kelvin kg)
!       t0         : reference temperature in Kelvin
!       p0         : total reference pressure in Pascal
!                    the calculation is based on a
!                    reduced pressure P*=Ptot-ro0*g.(x-xref)
!                    (except in compressible case)
!       xyzp0(3)   : coordinates of the reference point for
!                    the total pressure (where it is equal to p0)

!     In general, it is not necessary to furnish a reference point xyz0.
!       If there are outlets, the code will take the center of the
!       reference outlet face.
!       On the other hand, if we plan to explicitly fix Dirichlet conditions
!       for pressure, it is better to indicate to which reference the
!       values relate (for a better resolution of reduced pressure).


!     Other properties are given by default in all cases.

!     Nonetheless, we may note that:

!       In the standard case (no gas combustion, coal, electric arcs,
!                             compressibility):
!       ---------------------
!         ro0, viscl0 and cp0
!             are useful and represent either the fluid properties if they
!             are constant, either simple mean values for the initialization
!             if properties are variable and defined in usphyv.
!         t0  is not useful
!         p0  is useful but is not used in an equation of state. p0
!             is a reference value for the incompressible solver
!             which will serve to set the (possible) domain outlet pressure.
!             We may also take it as 0 or as a physical value in Pascals.

!       With the electric module:
!       ------------------------
!         ro0, viscl0 and cp0
!             are useful but simply represent mean initial values;
!             the density, molecular dynamic viscosity, and specific
!             heat are necessarily given in propce (whether they are
!             physically variable or not): see uselph for the Joule effect
!             module and the electric arcs dp_ELE data file.
!         t0  is useful an must be in Kelvin (> 0) but represents a simple
!             initialization value.
!         p0  is useful bu is not used in the equation of state. p0
!             is a reference value for the incompressible solver which
!             will be used to calibrate the (possible) outlet pressure
!             of the domain. We may take it as zero or as a physical
!             value in Pascals.

!       With gas combustion:
!       --------------------
!         ro0 is not useful (it is automatically recalculated by the
!             law of ideal gases from t0 and p0).
!         viscl0 is indispensable: it is the molecular dynamic viscosity,
!             assumed constant for the fluid.
!         cp0 is indispensable: it is the heat capacity, assumed constant,
!             (modelization of source terms involving a local Nusselt in
!             the Lagrangian module, reference value allowing the
!             calculation of a radiative
!             (temperature, exchange coefficient) couple).
!         t0  is indispensible and must be in Kelvin (> 0).
!         p0  is indispensable and must be in Pascal (> 0).

!       With pulverized coal:
!       ---------------------
!         ro0 is not useful (it is automatically recalculated by the
!             law of ideal gases from t0 and p0).
!         viscl0 is indispensable: it is the molecular dynamic viscosity,
!             assumed constant for the fluid (its effect is expected to
!             be small compared to turbulent effects).
!         cp0 is indispensable: it is the heat capacity, assumed constant,
!             (modelization of source terms involving a local Nusselt in
!             the coal or Lagrangian module, reference value allowing the
!             calculation of a radiative
!             (temperature, exchange coefficient) couple).
!         t0  is indispensable and must be in Kelvin (> 0).
!         p0  is indispensable and must be in Pascal (> 0).

!       With compressibility:
!       ---------------------
!         ro0 is not useful, stricto sensu; nonetheless, as experience
!             shows that users often use this variable, it is required
!             to assign to it a strictly positive value (for example,
!             an initial value).
!         viscl0 is useful and represents the molecular dynamic viscosity,
!             when it is constant, or a value which will be used during
!             initializations (or in inlet turbulence conditions,
!             depending on the user choice.
!         cp0 is indispensable: it is the heat capacity, assumed constant
!             in the thermodynamics available by default
!         t0  is indispensable and must be in Kelvin (> 0).
!         p0  is indispensable and must be in Pascal (> 0).
!             With the thermodynamic law available by default,
!             t0 and p0 are used for the initialization of the density.
!         xyzp0 is not useful because the pressure variable directly
!             represents the total pressure.

if (.false.) then
  ro0    = 1.17862d0
  viscl0 = 1.83337d-5
  cp0    = 1017.24d0
endif

if (.false.) then
  t0 = 20.d0 + 273.15d0
  p0 = 1.01325d5
endif

! We only specify XYZ0 if we explicitely fix Dirichlet conditions
! for the pressure.

if (.false.) then
  xyzp0(1) = 0.d0
  xyzp0(2) = 0.d0
  xyzp0(3) = 0.d0
endif

! --- irovar, ivivar: density and viscosity constant or not ?

!     When a specific physics module is active
!       (coal, combustion, electric arcs, compressible: see usppmo)
!       we MUST NOT set variables 'irovar' and 'ivivar' here, as
!       they are defined automatically.
!     Nonetheless, for the compressible case, ivivar may be modified
!       in the uscfx2 user subroutine.

!     When no specific physics module is active, we may specify if the
!         density and the molecular viscosity
!         are constant (irovar=0, ivivar=0), which is the default
!          or variable (irovar=1, ivivar=1)

!       if they are variable, the law must be defined in usphyv
!         (incs_user_physical_properties.f90);
!       if they are constant, they take values ro0 and viscl0.

if (.false.) then
  irovar = 1
  ivivar = 1
endif

! --- Minimum and maximum admissible values for each USER scalar:

!      Results are clipped at the end of each time step.

!      If min > max, we do not clip.

!      For a scalar jj representing the variance of another, we may
!        abstain from defining these values
!        (a default clipping is set in place).
!        This is the purpose of the test on iscavr(jj) in the example below.

!      For non-user scalars relative to specific physics (coal, combustion,
!        electric arcs: see usppmo) implicitly defined according to the
!        model, the information is automatically set elsewhere: we
!        do not set min or max values here.

if (.false.) then

  call field_get_key_id("min_scalar_clipping", kscmin)
  call field_get_key_id("max_scalar_clipping", kscmax)

  ! Loop on user scalars:
  do jj = 1, nscaus
    ! For scalars which are not variances
    if (iscavr(jj).le.0) then
      ! We define the min and max bounds
      call field_set_key_double(ivarfl(isca(jj)), kscmin, -grand)
      call field_set_key_double(ivarfl(isca(jj)), kscmax, +grand)
    endif
  enddo

endif

! --- Reference diffusivity visls0 in kg/(m s) for each
!        USER scalar except those which represent the variance of another.

!     For non-user scalars relative to specific physics (coal, combustion,
!       electric arcs: see usppmo) implicitly defined in the model,
!       the information is given automatically elsewhere:
!       we do not modify visls0 here.

!     For user scalars JJ which represent the variance of another user
!       scalar, we do not define visls0(jj) here.
!       This is the purpose of the test on iscavr(jj) in the example below.
!       Indeed the diffusivity of the variance of a scalar is assumed
!       identical to that scalar's diffusivity.

!     When no specific physics has been activated
!       (coal, combustion, electric arcs) and if a user scalar represents:
!       - the temperature:
!       visls0(iscalt) = Lambda
!        because the Cp is outside of the diffusion term in the temperature equation
!       - the enthalpy:
!       visls0(iscalt) = Lambda/Cp

!     Here, as an example, we assign to viscl0 the viscosity of the fluid
!       phase, which is fitting for passive tracers which follow the fluid
!       (this is also the default used if not modified here or using the GUI).

if (.false.) then

  ! We loop on user scalars:
  do jj = 1, nscaus
    ! For scalars which are not variances
    if (iscavr(jj).le.0) then
      ! We define the diffusivity
      visls0(jj) = viscl0
    endif
  enddo

endif

! --- Turbulent flux model u'T' for the scalar T
!     Algebraic Model
!      0  SGDH
!      10 GGDH
!      20 AFM
!     Model with transport equations
!      30 DFM

if (.false.) then

  ! GGDH for all the scalars:
  do jj = 1, nscaus
    iturt(jj) = 10
  enddo

endif


! --- Define scalar (among nscaus) which are species:
!     If a user scalar isca represents the species Yk,
!     iscasp(isca) is set to 1. By default, iscasp(isca)= 0.
!
!     To use the Low-Mach algorithm with a multi-species state law (idilat = 3),
!     we also have to specify the molar mass associated (wmolsp(isca))
!     to the species Yk.
!
!     WARNING: This algorithm assumes that the last species is deduced from
!              the others thank to the relation Ym = 1 - Sum_k Yk.
!              The molar mass associated to this species has to be
!              specified in wmolsp(0).

! The example set 4 species, the molar mass associated to the last one (not
! computed) is stored in  wmolsp(0).

if (.false.) then

  iscasp(2) =  1
  wmolsp(2) =  0.032d0

  iscasp(3) =  1
  wmolsp(3) =  0.016d0

  iscasp(4) =  1
  wmolsp(4) =  0.016d0

  wmolsp(0) =  0.028d0

endif

! --- Reference velocity for turbulence initialization (m2/s)
!       (useful only with turbulence)

if (.false.) then
  uref = 1.d0
endif

! --- Reference length scale in meters for initialization
!       of epsilon (and specific clipping of turbulence, but
!       this is not the default option)
!       Assign a value of the order of the largest dimension of the
!       physical domain in which the flow may develop.
!       If a negative value is set here, or no value set and the GUI not
!       used, the cubic root of the domain will be used.
!       (useful only for turbulence).

if (.false.) then
  almax = 0.5
endif

! --- Definition of moments
!     (at the most nbmomx moments, correlations of maximum order ndgmox)

!     We calculate temporal means of the type <f1*f2*f3*...*fn>
!     The fi's are cell-defined variables (arrays rtp and propce).

!        idfmom(i,imom) ientifies the variable fi of moment imom
!          if idfmom > 0 it is a resolved variable (rtp)
!          if idfmom < 0 it is an auxiliary variable (propce)
!        imoold(imom) defined in the case of a restart the number, in the
!          previous calculation, of the moment to use to initialize moment
!          imom of the new calculation (by default imoold(imom)=imom).
!            Value -1 indicates the we must reinitialize moment imom.
!        ntdmom(imom) defined the time step at which the moment calculation
!          is started.
!        ttdmom(imom) defined the time at which the moment calculation is started.

!     We give below the example of the calculation of moments <u> and <rho u v>
!       the moment <u> is reread in the restart file if we are restarting,
!         the moment <rho u v> is reinitialized to zero.
!       Moment <u> is calculated starting from time step 1000
!         Moment <rho u v> is calculated from time step 10000.

if (.false.) then

  ! First moment: <u>
  imom  = 1
  idfmom(1,imom) =  iu
  ntdmom(imom)   =  1000
  ttdmom(imom)   =  0.d0

  ! TODO: allow fields instead of properties for moments
  ! Second moment: <rho u v>
  ! imom  = 2
  ! idfmom(1,imom) = -irom
  ! idfmom(2,imom) =  iu
  ! idfmom(3,imom) =  iv
  ! imoold(imom)   = -1
  ! ntdmom(imom)   =  10000
  ! ttdmom(imom)   =  10.d0

endif

!----
! Formats
!----

return
end subroutine usipsu


!===============================================================================


subroutine usipes &
!================

 ( nmodpp )


!===============================================================================
! Purpose:
! --------

! User subroutine for the input of additional user parameters for
! input/output.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nmodpp           ! i  ! <-- ! number of active specific physics models       !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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
use ihmpre
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii, ipp, imom, idirac, icla, icha
integer idimve, iesp

!===============================================================================

!     This subroutine allows setting parameters

!       which do not already appear in the other subroutines of this file.


!     It is possible to add or remove parameters.


!     The number of physical properties and variables is known here.


!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================

!===============================================================================
! 1. Input-output (entsor)
!===============================================================================

! Frequency of log output

if (.false.) then

  ntlist = 1

endif

! Log (listing) verbosity

if (.false.) then

  do ii = 1, nvar
    iwarni(ii) = 1
  enddo

  iwarni(ipr) = 2
  iwarni(iu) = 2
  iwarni(iv) = 2
  iwarni(iw) = 2

endif

! --- probes output step

if (.false.) then

  nthist = 1
  frhist = -1.d0

endif

! --- Number of monitoring points (probes) and their positions
!     (limited to ncaptm=100)

if (.false.) then

  ncapt  = 4
  tplfmt = 1 ! time plot format (1: .dat, 2: .csv, 3: both)

  xyzcap(1,1) = 0.30d0
  xyzcap(2,1) = 0.15d0
  xyzcap(3,1) = 0.01d0

  xyzcap(1,2) = 0.30d0
  xyzcap(2,2) = 0.00d0
  xyzcap(3,2) = 0.01d0

  xyzcap(1,3) = 0.30d0
  xyzcap(2,3) =-0.08d0
  xyzcap(3,3) = 0.01d0

  xyzcap(1,4) = 0.60d0
  xyzcap(2,4) =-0.05d0
  xyzcap(3,4) = 0.01d0

endif

! Per variable output control
! Many more examples are provided in cs_user_parameters-output.f90

! User scalar variables.

! We may modify here the arrays relative to user scalars, but scalars
!   reserved for specific physics are handled automatically. This explains
!   the tests on 'nscaus', which ensure that the targeted scalars are
!   truly user scalars.
! By specific physics, we mean only those which are handled in specific
!   modules of the code, such as coal, combustion, electric arcs (see usppmo).

if (.false.) then

  if (isca(1).gt.0.and.nscaus.ge.1) then
    call field_set_key_str(ivarfl(isca(1)), keylbl, 'Scalar 1')
    ipp = ipprtp(isca(1))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

  if (isca(2).gt.0.and.nscaus.ge.2) then
    call field_set_key_str(ivarfl(isca(2)), keylbl, 'Scalar 1')
    ipp = ipprtp(isca(2))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

endif

!----
! Formats
!----


return
end subroutine usipes


!===============================================================================

subroutine user_field_parameters
!===============================

!===============================================================================
! Purpose:
! --------

! Define (redefine) key-value pairs on calculation fields.

! This subroutine is called at the end of the parameters initialization
! stage, after all other routines from this file have been called.

! Note that to determine which fields are defined in a computation, you
! may check the 'config.log' file after a first execution.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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
use ihmpre
use ppppar
use ppthch
use ppincl
use field

!===============================================================================

implicit none

! Local variables

logical       ilved, inoprv
integer       f_id, idim1, iflpst, itycat, ityloc, iscdri, iscal
integer       keydri


!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================

! --- Scalar with a drift (key work "drift_scalar_model">0) or without drift
!       ((key work "drift_scalar_model"=0, default option) for each USER scalar.
!       - to specify that a scalar have a drift and need the drift computation:
!       iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
!
! --- Then, for each scalar with a drift, add a flag to specify if
!     specific terms have to be taken into account:
!       - thermophoresis terms:
!       iscdri = ibset(iscdri, DRIFT_SCALAR_THERMOPHORESIS)
!       - turbophoresis terms:
!       iscdri = ibset(iscdri, DRIFT_SCALAR_TURBOPHORESIS)
!       - centrifugal force terms:
!       iscdri = ibset(iscdri, DRIFT_SCALAR_CENTRIFUGALFORCE)

if (.false.) then

  ! Key id for drift scalar
  call field_get_key_id("drift_scalar_model", keydri)

  if (nscaus.ge.1) then

    iscdri = 1
    iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)

    if (.false.) then
      iscdri = ibset(iscdri, DRIFT_SCALAR_THERMOPHORESIS)
    endif

    if (.false.) then
      iscdri = ibset(iscdri, DRIFT_SCALAR_TURBOPHORESIS)
    endif

    if (.false.) then
      iscdri = ibset(iscdri, DRIFT_SCALAR_CENTRIFUGALFORCE)
    endif

    iscal = 1
    f_id = ivarfl(isca(iscal))

    ! Set the key word "drift_scalar_model" into the field structure
    call field_set_key_int(f_id, keydri, iscdri)

  endif

endif

! Example: force postprocessing of projection of some variables at boundary
!          with no reconstruction.
!          This is handled automatically if the second bit of a field's
!          'post_vis' key value is set to 1 (which amounts to adding 2
!          to that key value).
!
!          field_get_id returns -1 if field does not exist

if (.false.) then

  f_id = ivarfl(iu)
  if (iand(iflpst, 2) .eq. 0) then
    iflpst = ior(iflpst, 2)
    call field_set_key_int(f_id, keyvis, iflpst)
  endif

  f_id = ivarfl(ipr)
  if (iand(iflpst, 2) .eq. 0) then
    iflpst = ior(iflpst, 2)
    call field_set_key_int(f_id, keyvis, iflpst)
  endif

endif

!-------------------------------------------------------------------------------

! Example: enforce existence of 'tplus' and 'tstar' fields, so that
!          a boundary temperature or Nusselt number may be computed using the
!          post_boundary_temperature or post_boundary_nusselt subroutines.
!          When postprocessing of these quantities is activated, those fields
!          are present, but if we need to compute them in the
!          cs_user_extra_operations user subroutine without postprocessing them,
!          forcing the definition of these fields to save the values computed
!          for the boundary layer is necessary.

if (.false.) then

  itycat = FIELD_INTENSIVE + FIELD_PROPERTY
  ityloc = 3 ! boundary faces
  idim1 = 1 ! dimension
  ilved = .true. ! interleaved
  inoprv = .false. ! no previous time step values needed

  call field_get_id('tplus', f_id)
  if (f_id.lt.0) then
    call field_create('tplus', itycat, ityloc, idim1, ilved, inoprv, f_id)
  endif

  call field_get_id('tstar', f_id)
  if (f_id.lt.0) then
    call field_create('tstar', itycat, ityloc, idim1, ilved, inoprv, f_id)
  endif

endif

!===============================================================================

!----
! Formats
!----

return
end subroutine user_field_parameters

!===============================================================================


subroutine usalin
!================

!===============================================================================
!  Purpose :
! --------

! --- User subroutine dedicated to the use of ALE (Arbitrary Lagrangian Eulerian)
!     method :
!
!          Here one defines parameters and input data dedicated to the use ALE
!          method
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use optcal
use albase

!===============================================================================

implicit none

! Arguments

! Local variables

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
!
! Here are some examples that can be adapted and changed by Code Saturne users.
!
!
! --- Activation of ALE (Arbitrary Lagrangian Eulerian) method

if (.false.) then
  iale = 1
endif

! --- Number of iterations for fluid initialization. Contrary to ntmabs (for example)
!     nalinf is not an absolute iteration number, meaning that in case of
!     restart calculation nalinf corresponds to the number of iterations
!     for fuid initialization beginning from the first current iteration of
!     the calculation restart. In general nalinf = 0 in that case.

if (.false.) then
  nalinf = 75
endif

! --- Maximum number of iterations in case of implicit Fluid Structure Coupling
!     with structural calculations (internal and/or external
!     (i.e. using Code_Aster)).
!     NALIMX = 1, in case of explicit FSI algorithm.

if (.false.) then
  nalimx = 15
endif

! --- Relative precision of sub-cycling Fluid Structure Coupling algorithm.

if (.false.) then
  epalim = 1.d-5
endif

! --- Mesh viscosity modeling (cf. usvima)
!     0 : isotropic
!     1 : orthotropic

if (.false.) then
  iortvm = 0
endif

!----
! Formats
!----

!----
! End
!----

return
end subroutine usalin


!===============================================================================


subroutine usati1
!================


!===============================================================================
! Purpose:
! --------

! Initialize non-standard calculation options for the atmospheric version.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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
use atincl
use atsoil
use atchem
use siream

!===============================================================================

implicit none

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. Example of calculation options to modify
!===============================================================================

!!! Reading the meteo file

imeteo = 1

!!! For radiative model or chemistry

! Time of the simulation
!  syear  --> starting year
!  squant --> starting quantile
!  shour  --> starting hour (UTC)
!  smin   --> starting minute
!  ssec   --> starting second

syear = 1994
squant = 1
shour = 1
smin = 0
ssec = 0.d0

! Geographic position
! xlon --> longitude of the domain origin
! xlat --> latitude of the domain origin

xlon = 0.d0
xlat = 45.d0

!!! Gaseous chemistry

! ichemistry: choice of chemistry resolution scheme
!0 --> no atmospheric chemistry
!1 --> quasi stationary equilibrium NOx scheme with 4 species and 5 reactions
!2 --> scheme with 20 species and 34 reactions
!3 --> scheme CB05 with 52 species and 155 reactions
!4 --> user defined schema
ichemistry = 1

! ificchemistry: choice to read (=1,2,3,4, according to the scheme)
! or not (0) a concentration profile file
! if ichemistry>0 ifilechemistry is automaticaly set to ichemistry
ifilechemistry = 0

! isepchemistry: splitted (=1) or semi-coupled (=2, pu-sun)
! resolution of chemistry
isepchemistry = 1

! iphotolysis: inclusion (=1) or not (=2) of photolysis reactions
iphotolysis = 1

! dtchemmax: maximal time step (s) for chemistry resolution
dtchemmax = 10.0d0

!!! Aerosol chemistry
! iaerosol: flag to activate aerosol chemistry
! if iaerosol = 1, ichemistry is automatically set to 3 (scheme 3)
iaerosol = 1

! inogaseouschemistry: flag to prevent automatic resolution (=1)
! of gaseous chemistry (scheme 3)
inogaseouschemistry = 0

! ncycle_aer: number of iterations for time splitting
ncycle_aer = 1

! icoag_siream: flag to activate (=1) or not (=0) coagulation
icoag_siream = 1

! icond_siream: flag to activate (=1) or not (=0) condensation/evaporation
icond_siream = 1

! inucl_siream: flag to activate (=1) or not (=0) nucleation
inucl_siream = 1

! icut_siream: cutting bin between equilibrium (1 to icut_siream)
! and dynamic bins (icut_siream to nbin_aer)
icut_siream = nbin_aer

!----
! End
!----

return
end subroutine usati1


!===============================================================================


subroutine usd3p1
!================


!===============================================================================
!  Features of this subroutine:
!  ----------------------------
!  1. Additional Calculation Options
!     a. Density Relaxation
!
!  2. Physical Constants
!     a.Dynamic Diffusion Coefficient
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat

!===============================================================================

implicit none

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Additional Calculation Options
!===============================================================================

! -->  Density Relaxation
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0

!===============================================================================
! 2. Physical Constants
!===============================================================================

!       DIFTL0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5

! Reference temperature for fuel and oxydant (K)
tinfue = 436.d0
tinoxy = 353.d0

!----
! End
!----

return
end subroutine usd3p1


!===============================================================================


subroutine usebu1

!===============================================================================
!  PURPOSE:
!  --------
!  1. Additional Calculation Options
!     a. Density Relaxation
!
!  2. Physical Constants
!     a.Dynamic Diffusion Coefficient
!===============================================================================

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
use radiat

!===============================================================================

implicit none

!===============================================================================

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Additional Calculation Options
!===============================================================================

! -->  Density Relaxation
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0


!===============================================================================
! 2. Physical Constants
!===============================================================================

! DIFTL0: Dynamic Diffusion Coefficient (kg/(m s))

diftl0 = 4.25d-5

! cebu: EBU-model constant

cebu   = 2.5d0

!----
! End
!----

return
end subroutine usebu1


!===============================================================================


subroutine uslwc1


!===============================================================================
!  PURPOSE:
!  --------
!  1. Additional Calculation Options
!     a. Density Relaxation
!
!  2. Physical Constants
!     a.Dynamic Diffusion Coefficient
!     b.Constants of the Libby-Williams Model
!===============================================================================

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
use radiat

!===============================================================================

implicit none

integer          ipp, idirac

!===============================================================================

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Additional Calculation Options
!===============================================================================

! -->  Density Relaxation
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1))

srrom = 0.95d0


!===============================================================================
! 2. Physical Constants
!===============================================================================

! --> DIFTL0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5

! --> Constants of the Libby-Williams Model

! Reference velocity
vref = 60.d0
! Reference length scale
lref = 0.1d0
! Activation Temperature
ta   = 0.2d5
! Cross-over Temperature (combustion of propane)
tstar= 0.12d4

!----
! End
!----

return
end subroutine uslwc1


!===============================================================================


subroutine uscfx1
!================


!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Initialize non standard options for the compressible flow scheme such
!    what kind of equation of state must be used.

!    In addition to options set in the user subroutine 'uscfx2' (or in
!    the GUI): this subroutine allows to set switches to indicate if the
!    volumetric viscosity and the conductivity are constants, their
!    values being given in the subroutine 'uscfx2'.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use ihmpre
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

!===============================================================================

implicit none

! Arguments


! Local variables

!===============================================================================

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     This subroutine is  mandatory for compressible flow,
!       thus the default (library reference) version stops immediately.
!===============================================================================

if (iihmpr.eq.0) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in compressible flow options definition',/,&
'@    =======',/,                                                 &
'@     The user subroutine ''uscfx1'' must be completed.',/,      &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Scheme options
!===============================================================================

if (iihmpr.eq.0) then   !  Remove test to set values here when also using GUI.

! Equation of state choice
! --> ieos = 1: Perfect gas with constant Gamma
! --> ieos = 2: Perfect gas with variable Gamma (please fill-in the source code
!               cfther in this case, otherwise it won't work!)
  ieos = 1

! --> Molecular thermal conductivity
!       constant  : ivisls = 0
!       variable  : ivisls = 1

  ivisls(itempk) = 0

! --> Volumetric molecular viscosity
!       iviscv = 0 : uniform  in space and constant in time
!              = 1 : variable in space and time

  iviscv = 0

endif

!----
! End
!----

return
end subroutine uscfx1


!===============================================================================


subroutine uscfx2
!================


! Purpose:
! -------

!    User subroutine.

!    Set values for the reference volumic viscosity, the reference
!    conductivity and the molar mass for compressible flow.

!    Initialize non standard options for the compressible flow scheme such
!    as the hydrostatic equilibrium.

!    In addition to options set in the user subroutine 'uscfx1' (or in
!    the GUI): this subroutine allows to set a switch to indicate if the
!    molecular viscosity is constant, its values being given in the user
!    subroutine 'usipsu'.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!


!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use ihmpre
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

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     This subroutine is  mandatory for compressible flow,
!       thus the default (library reference) version stops immediately.
!===============================================================================

if (iihmpr.eq.0) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input for compressible flow',/,    &
'@    =======',/,                                                 &
'@     The user subroutine ''uscfx2'' must be completed',/,       &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Physical properties
!===============================================================================

if (iihmpr.eq.0) then   !  Remove test to set values here when also using GUI.

! --> Molecular viscosity
!       constant  : ivivar = 0
!       variable  : ivivar = 1

  ivivar = 0

! --> Reference molecular thermal conductivity
!       visls0 = lambda0  (molecular thermal conductivity, W/(m K))

!       WARNING: visls0 must be strictly positive
!         (set a realistic value here even if conductivity is variable)

  visls0(itempk) = 3.d-2

!       If the molecular thermal conductivity is variable, its values
!         must be provided in the user subroutine 'uscfpv'

! --> Volumetric molecular viscosity

!       Reference volumetric molecular viscosity

!       viscv0 = kappa0  (volumetric molecular viscosity, kg/(m s))

  viscv0 = 0.d0

!       If the volumetric molecular viscosity is variable, its values
!         must be provided in the user subroutine 'uscfpv'

! --> Molar mass of the gas (kg/mol)

!       For example with dry air, xmasml is around 28.8d-3 kg/mol
  if (ieos.eq.1) then
    xmasmr = 0.028966
  endif

! --> Hydrostatic equilibrium

!       Specify if the hydrostatic equilibrium must be accounted for
!         (yes = 1 , no = 0)
  icfgrp = 1

endif

!----
! End
!----

return
end subroutine uscfx2


!===============================================================================


subroutine uscpl1
!================


!===============================================================================
!  Purpose:
!  -------

!   Lagrangian module coupled with pulverized coal:
!   -----------------------------------------------

!      Eulerian combustion of pulverized coal and
!      Lagrangian transport of coal particles

!    User subroutine for calculation parameter definitions (modules)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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
use ihmpre
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0. This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''uscpl1'' must be completed',/,       &
'@     for pulverized coal combustion coupled with',/,            &
'@     lagrangian transport of coal particles',/,                 &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Calculation options
!===============================================================================

! Relaxation coefficient for density
! rho(n+1) = srrom * rho(n) + (1-srrom) * rho(n+1)

srrom = 0.8d0


!===============================================================================
! 2. Physical constants
!===============================================================================

! Laminar viscosity associated t Enthalpy scalar
! DIFTL0 (dynamic diffusivity in kg/(m s))
diftl0 = 4.25d-5


!----
! End
!----

return

end subroutine uscpl1


!===============================================================================


subroutine user_coal_ini1
!========================


!===============================================================================
!  Purpose:
!  ---------

!  User's routine to control outing of variables for pulverised coal combustion
!  (these parameters are in a module)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!
!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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
use ihmpre
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl

!===============================================================================

implicit none

!===============================================================================

!===============================================================================
! 1. Computation OPTION
!===============================================================================

! --- Relaxation for density (Advisable when starting combustion computation)
!                            (Forbidden for unstationnary computation)
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.95d0

!===============================================================================
! 2. Physical constants
!===============================================================================

! ---  Laminar viscosity for enthalpy (dynamical diffusivity) kg/(m.s)
diftl0 = 4.25d-5

!----
! End
!----


end subroutine user_coal_ini1


!===============================================================================


subroutine user_fuel_ini1
!========================

!===============================================================================
!  Purpose:
!  --------

!  User routine for allocate computation parameters dealing with fuel

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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
use cs_fuel_incl
use ppincl
use ppcpfu

!===============================================================================

implicit none

!===============================================================================

!===============================================================================
! 1. Computation OPTION
!===============================================================================

! --- Relaxation for density (Advisable when starting combustion computation)
!                            (Forbidden for unstationnary computation)
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.7d0


!===============================================================================
! 2. Physical constants
!===============================================================================

! ---  Laminar viscosity for enthalpy (dynamical diffusivity) kg/(m.s)
diftl0 = 4.25d-5

!----
! End
!----

return

end subroutine user_fuel_ini1


!===============================================================================


subroutine useli1 &
!================

 ( iihmpu )


!===============================================================================
!  Purpose  :
!  -------
!          User subroutines for input of calculation parameters,
!       and to initialize variables used for specific electric models,
!
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iihmpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================
!

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
use elincl
use mesh

!===============================================================================

implicit none

! Arguments

integer iihmpu

integer          ipp, iesp , idimve
integer          ilelt, nlelt, izone, iel
integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0. This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================


if (iihmpu.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''useli1'' must be completed',/,       &
'@     for electric module',/,                                    &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


! Allocate a temporary array for cells selection
allocate(lstelt(ncel))

!===============================================================================
! 1. Calculation options
!===============================================================================

! --> Relaxation coefficient for mass density
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom = 0.d0

! --> "Electric variables" scaling (Joule effect or electric arc version)
!      IELCOR = 0 : NO Correction
!      IELCOR = 1 : CORRECTION
ielcor = 0

!     Imposed current intensity (electric arc) in Amp
!        and Imposed Power (Joule effect for glass melting applications) in Watt
!       These values have to be positive
!
couimp = 0.d0
puisim = 0.d0

!     Initial Potential Difference (positive value)
dpot = 0.d0

! ---> Model for scaling intensity (electric arcs)
!       MODREC = 0 : user defined
!       MODREC = 1 : standard model
!       MODREC = 2 : resetting plane model for electromagnetic quantities
modrec = 1

! ---> Define current density component used to calculate current when MODREC = 2
!       IDRECA (1, 2 or 3) for component (x, y or z)
idreca = 3

! Example : plane z = 3 with epsilon 0.0002

crit_reca(1) = 0.
crit_reca(2) = 0.
crit_reca(3) = 1.
crit_reca(4) = -3.
crit_reca(5) = 0.0002

! Deallocate the temporary array
deallocate(lstelt)

!----
! End
!----

return
end subroutine useli1


!===============================================================================


subroutine uscti1
!================


!===============================================================================
! Purpose:
! -------

! Definition of cooling tower model and exchange zones

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0. This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''uscti1'' must be completed',/,       &
'@     for the cooling tower module',/,                           &
'@                                                            ',/,&
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Parameters for prescibed temperature difference
!===============================================================================

! Activation
iaeeri = 0

! Temperature difference (cooling) to prescribe
vaeeri = 13.d0

! Temperature modification frequency
iaeerp = 5

! Temperature step to compute difference slope tref(teau)
paseri = 0.015d0

! Maximum average hot water temperature
aetemx = 80.d0

! Minimum average cooled water temperature
aetemn = 10.d0

! Number of excange zones with a water inlet boundary
nbzsup = 2

! List of the nbzsup exchange zones at water inlet boundary
lizsup(1) = 1
lizsup(2) = 2

! Number of excange zones with a water outlet boundary
nbzinf = 2

! List of the nbzinf exchange zones at water outlet boundary
lizinf(1) = 1
lizinf(2) = 2

! Prescribed difference activation start time

inbaei = 1000.d0

!===============================================================================
! 2. Post-processing of exchange zones
!===============================================================================

ichrze = 1

!===============================================================================
! 3. Cooling tower restart
!===============================================================================

isuict = isuite

!----
! End
!----

return
end subroutine uscti1
