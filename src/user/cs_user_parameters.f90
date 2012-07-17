!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

! User subroutines for input of calculation parameters (Fortran commons).
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


subroutine usppmo
!================


!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Define the use of a specific physics amongst the following:
!      - combustion with gaz / coal / heavy fuel oil
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

!===============================================================================

implicit none


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!===============================================================================

if(1.eq.1) then
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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(icod3p) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(icoebu) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- colwc: Libby-Williams pre-mixed flame
! ==========

!        if = -1   module not activated
!        if =  0   reference two-peak model with adiabatic condition
!        if =  1   extended two-peak model with enthapy source terms
!        if =  2   extended three-peak model, adiabatic
!        if =  3   extended three-peak model with enthalpy source terms
!        if =  4   extended four-peak model, adiabatic
!        if =  5   extended four-peak model with enthalpy source terms

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(icolwc) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- cp3pl: Pulverized coal combustion
! ==========

!        Description of granulometry
!        Assumption of diffusion flame around particles
!         (extension of 3-point fast chemistry "D3P")
!        Between a mixture of gaseous fuels (volatiles matters, CO from char
!                                            oxydation)
!            and a mixture of oxidisers (air and water vapor)
!        Enthalpy for both mix and solid phase are solved

!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(icp3pl) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- cpl3c: Pulverized coal with Lagrangian reciprocal approach
! ==========

!        Not recently tested... at least outdated, may be obsolete

!        if = -1   module not activated
!        if = 0    module activated
!        if = 1    with drying (NOT functional)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(icpl3c) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- cfuel: Heavy fuel oil combustion
! ==========

!        Progressive evaporation (temperature gap)
!        Char residue
!        Sulphur tracking

!        if = -1   module not activated
!        if = 0    module activated

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(icfuel) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(iccoal) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- compf: Compressible flows
! ==========

!        if = -1   module not activated
!        if = 0    module activated

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(icompf) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- eljou: Joule effect
! ==========

!        if = -1   module not activated
!        if = 1    Potentiel reel
!        if = 2    Potentiel complexe
!        if = 3    Potentiel reel     + CDL Transfo
!        if = 4    Potentiel complexe + CDL Transfo

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(ieljou) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- elarc: Electric arcs
! ==========

!        if = -1   module not activated
!        if = 1    electric potential
!        if = 2    electric potential and vector potential (hence 3D modelling)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(ielarc) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- atmos: Atmospheric flows
! ==========

!        if = -1   module not activated
!        if = 0    standard modelling
!        if = 1    dry atmosphere
!        if = 2    humid atmosphere (NOT functional)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(iatmos) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- aeros: Cooling towers
! ==========

!        if = -1   module not activated
!        if = 0    no model (NOT functional)
!        if = 1    Poppe's model
!        if = 2    Merkel's model

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ippmod(iaeros) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


!===============================================================================
! 2.  Specific physics module not available at the moment
!===============================================================================

! WARNING: The following modules ARE NOT functional!
! =======

! --- cobml: Premix model of Bray - Moss - Libby
! ==========

!        if = -1   module not activated

ippmod(icobml) = -1

! --- codeq: Diffusion flame with fast equilibrium chemistry
! ==========

!        if = -1   module not activated

ippmod(icodeq) = -1

! --- elion: Ionic mobility
! ==========

!        if = -1   module not activated

ippmod(ielion) = -1


!===============================================================================
! 3.  Specific options related to herebefore modules
!===============================================================================

! These options are defined here at the moment, this might change in the future

! --- Enthalpy-Temperature conversion law (for gas combustion modelling)

!       if = 0   user-specified
!       if = 1   tabulated by JANAF (default)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

indjon = 1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Kinetic model for NOx formation

!         Only compatible with heavy fuel oil combustion

!         if = 0  unused
!         if = 1  activated

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ieqnox = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Kinetic model for CO <=> CO2

!         Compatible with coal and heavy fuel oil combustion

!         if = 0  unused (maximal conversion in turbulent model)
!         if = 1  transport of CO2 mass fraction
!         if = 2  transport of CO mass fraction

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ieqco2 = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Heteregoneous combustion by CO2

!         Needs the activation of the CO2 transport equation
!         Account for the reaction between char and CO2: C(s) + CO2 => 2 CO

!         if = 0  unused
!         if = 1  activated

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ihtco2 = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!----
! Formats
!----

!----
! End
!----

return
end subroutine


!===============================================================================


subroutine usipph &
!================

 ( iihmpu, nfecra , iturb , icp , iverif )


!===============================================================================
! Purpose:
! --------

! User subroutine for input of parameters.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iihmpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! iturb            ! ia ! <-> ! turbulence model                               !
! icp              ! ia ! <-> ! flag for uniform Cp or not                     !
! iverif           ! i  ! <-- ! flag for elementary tests                      !
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

integer iihmpu, nfecra
integer iturb, icp
integer iverif

! Local variables

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (iverif.eq.0) then
  if (iihmpu.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
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

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================

! --- Turbulence
!       0...Laminar
!      10...Mixing length
!      20...k-epsilon
!      21...k-epsilon (linear production)
!      30...Rij-epsilon, (standard LRR)
!      31...Rij-epsilon (SSG)
!      32...Rij-epsilon (EBRSM)
!      40...LES (Smagorinsky)
!      41...LES (Dynamic)
!      42...LES (WALE)
!      50...v2f (phi-model)
!      51...v2f (BL-v2/k)
!      60...k-omega SST
!      70...Spalart Allmaras
!  For 10, contact the development team before use

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iturb = 20

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

icp = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!----
! Formats
!----


return
end subroutine


!===============================================================================


subroutine usinsc &
!================

 ( iihmpu, nfecra , nscaus , iverif )


!===============================================================================
! Purpose:
! -------

! User subroutine for input of the number of user scalars.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iihmpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! nscaus           ! i  ! <-> ! number of user scalars                         !
! iverif           ! i  ! <-- ! flag for elementary tests                      !
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

integer iihmpu, nfecra
integer nscaus
integer iverif

! Local variables


!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (iverif.eq.0) then
  if (iihmpu.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''usinsc'' must be completed',/,       &
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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

nscaus = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!----
! Formats
!----


return
end subroutine


!===============================================================================


subroutine usipsc &
!================

 ( nscmax, nscaus, iihmpu, nfecra, iscavr, ivisls , iverif )


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
! iihmpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! iscavr(nscmax)   ! ia ! <-- ! associated scalar number for variance scalars  !
! ivisls(nscmax)   ! ia ! <-> ! uniform scalar diffusivity flag                !
! iverif           ! i  ! <-- ! flag for elementary tests                      !
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

integer nscmax, nscaus, iihmpu, nfecra
integer iscavr(nscmax), ivisls(nscmax)
integer iverif

! Local variables

integer iutile, iscal

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (iverif.eq.0) then
  if (iihmpu.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''usipsc'' must be completed',/,       &
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


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     The test on iutile allows deactivation of the instructions
!       (which are only given as an example).

iutile = 0
if (iutile.eq.1) then
  iscavr(3) = 2
endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

do iscal = 1, nscaus

  ! For user scalars which do not represent the variance of another scalar
  if (iscavr(iscal).le.0) then

    ivisls(iscal) = 0

  endif

enddo

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!----
! Formats
!----



return
end subroutine


!===============================================================================


subroutine usipgl &
!================

 ( nesmax,                                                        &
   iespre, iesder, iescor, iestot,                                &
   iihmpu, nfecra,                                                &
   idtvar, ipucou, idilat, iphydr, ialgce , iescal , iverif ,     &
   icwfps, cwfthr )


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
! iihmpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! idtvar           ! i  ! --> ! variable time step flag                        !
! ipucou           ! i  ! --> ! reinforced u-p coupling flag                   !
! idilat           ! i  ! --> ! algorithm with density variation in time       !
! iphydr           ! i  ! --> ! flag for handling of the equilibrium between   !
!                  !    !     ! the pressure gradient and the gravity and      !
!                  !    !     ! head-loss terms                                !
! ialgce           ! i  ! <-- ! option for the method of calculation of        !
!                  !    !     !  cell centers                                  !
! iescal(nesmax)   ! ia ! <-- ! flag for activation of error estimators for    !
!                  !    !     ! Navier-Stokes                                  !
! iverif           ! i  ! <-- ! flag for elementary tests                      !
! cwfthr           ! i  ! <-- ! Treshold angle to cut warped faces (do not     !
!                  !    !     !  cut warped faces if value is negative)        !
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
integer iihmpu, nfecra
integer idtvar, ipucou, idilat, iphydr, ialgce
integer iescal(nesmax)
integer iverif, icwfps

double precision cwfthr

! Local variables

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (iverif.eq.0) then
  if (iihmpu.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''usipgl'' must be completed',/,       &
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

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

!===============================================================================

! --- Time step  (0 : uniform and constant
!                 1 : variable in time, uniform in space
!                 2 : variable in time and space
!                -1 : steady algorithm)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

idtvar = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Velocity/pressure coupling (0 : classical algorithm,
!                                 1 : transient coupling)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ipucou = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! Algorithm to take into account the density variation in time
!
!     idilat = 0 : boussinesq algorithm with constant density (not available)
!              1 : dilatable steady algorithm (default)
!              2 : dilatable unsteady algorithm
!              3 : low-Mach algorithm
!

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

idilat = 1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Handling of hydrostatic pressure
!                               (0 : usual algorithm
!                                1 : specific handling)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iphydr = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Estimators for Navier-Stokes (non-frozen velocity field)
!     We recommend running a calculation restart on a few time steps
!       with the activation of the most interesting of those.
!        (=2 to activate, =0 to deactivate).

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!       div(rho u) -Gamma
iescal(iescor) = 0
!       resolution precision for the momentum
iescal(iestot) = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Triangulate warped faces:
!       If cwfthr is positive, faces whose warping angle are greater than
!         the given value (in degrees) are subdivided into triangles;
!       if cwfthr negative, faces are not subdivided.
!       If icwfps = 1, additional postprocessing will be activated to
!         show faces before and after cutting.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

icwfps = 0
cwfthr= -1.d0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


!----
! Formats
!----


return
end subroutine


!===============================================================================


subroutine usipsu &
!================

 ( nmodpp , iverif )


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
! iverif           ! i  ! <-- ! flag for elementary tests                      !
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

!===============================================================================

implicit none

! Arguments

integer nmodpp
integer iverif

! Local variables

integer iutile, ii, jj, imom

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (iverif.eq.0) then
  if (iihmpr.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''usipsu'' must be completed',/,       &
'@       in file cs_user_parameters.f90',/,                       &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

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

!     In case of restart, read auxiliary restart file ileaux (= 1) or not (0).

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ileaux = 1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Duration
!       ntmabs = absolute number of the last time step required
!         if we have already run 10 time steps and want to
!         run 10 more, ntmabs must be set to 10 + 10 = 20

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ntmabs = 10

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Reference time step
!     The example given below is probably not adapted to your case.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

dtref  = 0.01d0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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



!   When specific physics are activated (coal, combustion, electric arcs)
!     we DO NOT edit this section: we DO NOT modify 'iscalt' nor 'iscsth'
!    (the test: if (nmodpp.eq.0) is used for this).


!   On the other hand, if specific physics are NOT activated:

!     If a USER scalar represents the temperature or enthalpy:
!       we define the number of this scalar in iscalt and
!       we set iscsth(iscalt) = 1 if it is the temperature
!          or  iscsth(iscalt) = 2 if it is the enthalpy.

!     If no scalar represents the temperature or enthalpy
!       we set iscalt = -1
!       and we do not define iscsth(iscalt).


!     For the radiative module when used without specific physics, if we
!      have chosen to solve in temperature (that is if
!      iscsth(iscalt) = 1), the fluid temperature is considered to
!      be in degrees KELVIN (be careful for boundary conditions an expression
!      of physical properties depending on temperature).
!      Nonetheless, even though it is not recommended, if we wish for the
!      fluid solver to work with a temperature in degrees Celsius, we must set
!      iscsth(iscalt) = -1.
!      This choice is a source of user errors. Indeed, the boundary conditions
!      for the fluid temperature will then be in degrees Celsius, while the
!      boundary conditions for radiation in usray2 must still be in Kelvin.


!    If specific physics are not activated
!       (coal, combustion, electric arcs: see usppmo):

! --- Segregated or coupled solver for the velocity components:
!       0 for the segregated solver
!       1 for the coupled solver
!
!     The coupled solver may improve the accuracy and the robustness of the
!     simulation in case of periodicity of rotation, Corriolis source terms.
!     It implicits the wall shear stress.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ivelco = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Solver taking a pscalar porosity into account:
!       0 No porosity taken into account (Standard)
!       1 Porosity taken into account
!

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iporos = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


if (nmodpp.eq.0) then

  ! Number of the scalar representing temperature or enthalpy,
  !   or -1 if there is none.
  ! When the choice is done by the Code_Saturne GUI, the scalar representing
  !   the temperature or enthalpy is always the first.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

  iscalt = -1

! If there is a temperature or enthalpy variable:
  if (iscalt.gt.0) then
    ! we indicate if it is the temperature (=1) or the enthalpy (=2).
    iscsth(iscalt) = 1
  endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Calculation (restart) with frozen velocity field (1 yes, 0 no)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iccvfg = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Vortex method for inlet conditions in L.E.S.
!       (0: not activated,  1: activated)
!     The vortex method only regards the L.E.S. models
!     To use the vortex method, edit the 'usvort.f90' user file.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

if (itytur.eq.4) then
  ivrtex = 0
endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

blencv(iu) = 1.0d0
blencv(iv) = 1.0d0
blencv(iw) = 1.0d0
if (nscaus.ge.1) then
  do ii = 1, nscaus
    blencv(isca(ii)) = 1.0d0
  enddo
endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Linear solver parameters (for each unknown)

!     iresol = -1:           default
!     iresol = 1000*ipol +j: ipol is the degree of the Neumann polynomial
!                            used for preconditioning,
!                            j = 0: conjugate gradient,
!                            j = 1: Jacobi
!                            j = 2: bi-CgStab
!                            j = 3: GMRES

!     nitmax: maximum number of iterations for each unknown ivar
!     epsilo: relative precision for the solution of the linear system.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iutile = 0
if (iutile.eq.1) then

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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Algebraic multigrid parameters

!     imgr = 0: no multigrid
!     imgr = 1: algebraic multigrid

!     Only available for pressure and purely diffusive variables.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! mltmmn = 300  ! mean number of cells under which merging takes place
! mltmgl = 500  ! global number of cells under which merging takes place
! mltmmr = 1    ! number of active ranks under which no merging is done
! mltmst = 4    ! number of ranks over which merging takes place
! mlttyp = 0    ! 0: loop over faces to coarsen in natural order
!               ! 1: loop over faces to coarsen in criteria order
!               ! 3: loop over faces to coarsen in Hilbert order

imgr(ipr) = 1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     The test on iutile allows deactivation of the instructions
!       (which are only given as an example).

iutile = 0
if (iutile.eq.1) then

  if (iturb.eq.20) then
    ircflu(ik)   = 0
    ircflu(iep)  = 0
  endif

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! Physical constants (cstphy)
! ===========================

! --- gravity (g in m/s2, with the sign in the calculation coordinate axes).

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

gx = 0.d0
gy = 0.d0
gz = 0.d0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- rotation vector of the reference frame (omega in s-1)

!       If the rotation is not nul, then
!          icorio = 0: rotation is taken into account by rotating the mesh
!                      (simulation in the absolute frame)
!                 = 1: rotation is taken into account by Coriolis source terms
!                      (simulation in the relative frame)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

icorio = 0

omegax = 0.d0
omegay = 0.d0
omegaz = 0.d0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ro0    = 0.235d0
viscl0 = 0.84d-6
cp0    = 1219.d0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

t0 = 1000.d0 + 273.15d0
p0 = 1.01325d5
! We only specify XYZ0 if we explicitely fix Dirichlet conditions
! for the pressure.
! xyzp0(1) = 0.d0
! xyzp0(2) = 0.d0
! xyzp0(3) = 0.d0


! --- irovar, ivivar: density and viscosity constant or not ?

!     When a specific physics module is active
!       (coal, combustion, electric arcs, compressible: see usppmo)
!       we DO NOT set variables 'irovar' and 'ivivar' here, as
!       they are defined automatically.
!     Nonetheless, for the compressible case, ivivar may be modified
!       in the uscfx1 user subroutine.

!     When no specific physics module is active, it is necessary to
!       specify if the density and the molecular viscosity
!         are constant (irovar=0, ivivar=0)
!          or variable (irovar=1, ivivar=1)

!       if they are variable, the law must be defined in usphyv;
!       if they are constant, they take values ro0 and viscl0.

!       as an example, we assume below that they are constant.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

if (nmodpp.eq.0) then
  irovar = 0
  ivivar = 0
endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Minimum (scamin) and maximum (scamax) admissible values for
!        each USER scalar:

!      Results are clipped at the end of each time step.

!      If scamin > scamax, we do not clip.

!      For a scalar jj representing the variance of another, we may
!        abstain from defining these values
!        (a default clipping is set in place).
!        This is the purpose of the test on iscavr(jj) in the example below.

!      For non-user scalars relative to specific physics (coal, combustion,
!        electric arcs: see usppmo) implicitly defined according to the
!        model, the information is automatically set elsewhere: we
!        do not set scamin or scamax.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! If there are user scalars
if (nscaus.gt.0) then

  ! Loop on user scalars:
  do jj = 1, nscaus
    ! For scalars which are not variances
    if (iscavr(jj).le.0) then
      ! We define the min and max bounds
      scamin(jj) =-grand
      scamax(jj) =+grand
    endif
  enddo

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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
!       (coal, combustion, electric arcs) and if a user scalar represents
!       the temperature or enthalpy:
!       visls0(iscalt) = Lambda/Cp

!     Here, as an example, we assign to viscl0 the viscosity of the
!       carrier phase, which is fitting for passive tracers which
!       follow the fluid.


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! If there are user scalars
if (nscaus.gt.0) then

  ! We loop on user scalars:
  do jj = 1, nscaus
    ! For scalars which are not variances
    if (iscavr(jj).le.0) then
      ! We define the diffusivity
      visls0(jj) = viscl0
    endif
  enddo

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! The example set 4 species, the molar mass associated to the last one (not
! computed) is stored in  wmolsp(0).

if (nscaus.gt.0) then
  iscasp(2) =  1
  wmolsp(2) =  0.032d0

  iscasp(3) =  1
  wmolsp(3) =  0.016d0

  iscasp(4) =  1
  wmolsp(4) =  0.016d0


  wmolsp(0) =  0.028d0
endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Reference velocity for turbulence initialization (m2/s)
!       (useful only with turbulence)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

uref    = 1.d0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Reference length scale in meters for initialization
!       of epsilon (and specific clipping of turbulence, but
!       this is not the default option)
!       Assign a value of the order of the largest dimension of the
!       physical domain in which the flow may develop.
!       (useful only for turbulence).

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

almax = -grand

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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

!     We give below the example of the calculation of moments <u> and <rho u v>
!       the moment <u> is reread in the restart file if we are restarting,
!         the moment <rho u v> is reinitialized to zero.
!       Moment <u> is calculated starting from time step 1000
!         Moment <rho u v> is calculated from time step 10000.


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     The test on iutile allows deactivation of the instructions
!       (which are only given as an example).

iutile = 0
if (iutile.eq.1) then

  ! First moment: <u>
  imom  = 1
  idfmom(1,imom) =  iu
  ntdmom(imom)   =  1000
  ! Second moment: <rho u v>
  imom  = 2
  idfmom(1,imom) = -irom
  idfmom(2,imom) =  iu
  idfmom(3,imom) =  iv
  imoold(imom)   = -1
  ntdmom(imom)   =  10000

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!----
! Formats
!----


return
end subroutine


!===============================================================================


subroutine usipes &
!================

 ( nmodpp , iverif )


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
! iverif           ! i  ! <-- ! flag for elementary tests                      !
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
use period
use ihmpre
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

integer nmodpp
integer iverif

! Local variables

integer ii, ipp, imom, iutile

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if (iverif.eq.0) then
  if (iihmpr.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
else
  if(iihmpr.eq.1) then
    return
  endif
endif

 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''usipes'' must be completed',/,       &
'@       in file cs_user_parameters.f90',/,                       &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

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

! --- write auxiliary restart file iecaux = 1 yes, 0 no

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iecaux = 1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! Frequency of log output

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ntlist = 1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! Log (listing) verbosity

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iutile = 0
if (iutile.eq.1) then

  do ii = 1, nvar
    iwarni(ii) = 1
  enddo

  iwarni(ipr) = 2
  iwarni(iu) = 2
  iwarni(iv) = 2
  iwarni(iw) = 2

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- history output step

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

nthist = 1
frhist = -1.d0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Number of monitoring points (probes) and their positions
!     (limited to ncaptm=100)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

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

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- current variable

!     As for other variables,
!       if we do not assign the following array values,
!       default values will be used

!     nomvar( ) = variable name
!     ichrvr( ) = chonological output (yes 1/no 0)
!     ilisvr( ) = logging in listing (yes 1/no 0)
!     ihisvr( ) = history output (number of probes and their numbers)
!     if ihisvr(.,1)  = -1, output for all probes

!     Note: Only the fist 8 characters of a name will be used in the most
!           detailed log.



! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! Current dynamic variables

! pressure variable
ipp = ipprtp(ipr)
nomvar(ipp)   = 'Pressure'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
if (icorio.eq.1) then
  nomvar(ipp)   = 'Rel Pressure'
endif

! variable v1x
ipp = ipprtp(iu)
nomvar(ipp)   = 'VelocityX'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
if (icorio.eq.1) then
  nomvar(ipp)   = 'Rel VelocityX'
endif

! v1y variable
ipp = ipprtp(iv)
nomvar(ipp)   = 'VelocityY'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
if (icorio.eq.1) then
  nomvar(ipp)   = 'Rel VelocityY'
endif

! v1z variable
ipp = ipprtp(iw)
nomvar(ipp)   = 'VelocityZ'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
if (icorio.eq.1) then
  nomvar(ipp)   = 'Rel VelocityZ'
endif

if (itytur.eq.2) then

  ! turbulent kinetic energy
  ipp = ipprtp(ik)
  nomvar(ipp)   = 'Turb Kinetic Energy'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! turbulent dissipation
  ipp = ipprtp(iep)
  nomvar(ipp)   = 'Turb Dissipation'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif (itytur.eq.3) then

  ! Reynolds stresses
  ipp = ipprtp(ir11)
  nomvar(ipp)   = 'R11'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir22)
  nomvar(ipp)   = 'R22'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir33)
  nomvar(ipp)   = 'R33'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir12)
  nomvar(ipp)   = 'R12'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir13)
  nomvar(ipp)   = 'R13'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir23)
  nomvar(ipp)   = 'R23'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! turbulent dissipation
  ipp = ipprtp(iep)
  nomvar(ipp)   = 'Turb Dissipation'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif (iturb.eq.50) then

  ! turbulent kinetic energy
  ipp = ipprtp(ik)
  nomvar(ipp)   = 'Turb Kinetic Energy'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! turbulent dissipation
  ipp = ipprtp(iep)
  nomvar(ipp)   = 'Turb Dissipation'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! phi
  ipp = ipprtp(iphi)
  nomvar(ipp)   = 'Phi'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! f_bar
  ipp = ipprtp(ifb)
  nomvar(ipp)   = 'f_bar'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif (iturb.eq.51) then

  ! turbulent kinetic energy
  ipp = ipprtp(ik)
  nomvar(ipp)   = 'Turb Kinetic Energy'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! turbulent dissipation
  ipp = ipprtp(iep)
  nomvar(ipp)   = 'Turb Dissipation'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! phi
  ipp = ipprtp(iphi)
  nomvar(ipp)   = 'Phi'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! alpha
  ipp = ipprtp(ial)
  nomvar(ipp)   = 'Alpha'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif (iturb.eq.60) then

  ! turbulent kinetic energy
  ipp = ipprtp(ik)
  nomvar(ipp)   = 'Turb Kinetic Energy'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! omega
  ipp = ipprtp(iomg)
  nomvar(ipp)   = 'Omega'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif (iturb.eq.70) then

  ! Spalart-Allmaras variable (viscosity-like)
  ipp = ipprtp(inusa)
  nomvar(ipp)   = 'NuTilda'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! User scalar variables.

! We may modify here the arrays relative to user scalars, but scalars
!   reserved for specific physics are handled automatically. This explains
!   the tests on 'nscaus', which ensure that the targeted scalars are
!   truly user scalars.
! By specific physics, we mean only those which are handled in specific
!   modules of the code, such as coal, combustion, electric arcs (see usppmo).

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

if (isca(1).gt.0.and.nscaus.ge.1) then
  ipp = ipprtp(isca  (1))
  nomvar(ipp)  = 'Scalar 1'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

if (isca(2).gt.0.and.nscaus.ge.2) then
  ipp = ipprtp(isca  (2))
  nomvar(ipp)  = 'Scalar 2'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! Other variables

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! Density variable (output for post-processing only if variable or
!                   in the case of specific physics)
ipp = ipppro(ipproc(irom))
nomvar(ipp)   = 'Density'
ichrvr(ipp)   = max(irovar,nmodpp)
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! specific heat
if (icp .gt. 0) then
  ipp = ipppro(ipproc(icp))
  nomvar(ipp)   = 'Specific Heat'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = 0
endif

! laminar viscosity
ipp = ipppro(ipproc(iviscl))
nomvar(ipp)   = 'Laminar Viscosity'
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = 0

! turbulent viscosity
ipp = ipppro(ipproc(ivisct))
nomvar(ipp)   = 'Turb Viscosity'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! Courant number
ipp = ipppro(ipproc(icour))
nomvar(ipp)   = 'CFL'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

! Fourier number
ipp = ipppro(ipproc(ifour))
nomvar(ipp)   = 'Fourier Number'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

! 'csmago' variable for dynamic L.E.S. models
!    (square of the Samgorinsky "constant")
if (ismago.gt.0) then
  ipp = ipppro(ipproc(ismago))
  nomvar(ipp)   = 'Csdyn2'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
endif

! temporal means (example for moment 1)
if (nbmomt.gt.0) then
  imom = 1
  ipp = ipppro(ipproc(icmome(imom)))
  nomvar(ipp) = 'Time Average 01'
  ichrvr(ipp) = 1
  ilisvr(ipp) = 1
  ihisvr(ipp,1) = -1
endif

! total pressure (not defined in compressible case)
if (ippmod(icompf).lt.0) then
  ipp = ipppro(ipproc(iprtot))
  nomvar(ipp)   = 'Total Pressure'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
endif

! local time step
ipp = ippdt
nomvar(ipp)   = 'Local Time Step'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! characteristic time of transient velocity/pressure coupling
ipp = ipptx
nomvar(ipp)   = 'Tx'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

ipp = ippty
nomvar(ipp)   = 'Ty'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

ipp = ipptz
nomvar(ipp)   = 'Tz'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


!----
! Formats
!----



return
end subroutine


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

integer iutile

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
!
!     Here are some examples that can be adapted and changed by Code Saturne
!     users.
!
!
! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Activation of ALE (Arbitrary Lagrangian Eulerian) method
iutile = 0
if (iutile.eq.1) then
  iale = 1
endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Number of iterations for fluid initialization. Contrary to ntmabs (for example)
!     nalinf is not an absolute iteration number, meaning that in case of
!     restart calculation nalinf corresponds to the number of iterations
!     for fuid initialization beginning from the first current iteration of
!     the calculation restart. In general nalinf = 0 in that case.

nalinf = 75

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Maximum number of iterations in case of implicit Fluid Structure Coupling with structural
!     calculations (internal and/or external(i.e. using Code_Aster)). NALIMX = 1, in case of
!     explicit FSI algorithm.

nalimx = 15

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Relative precision of sub-cycling Fluid Structure Coupling algorithm.
epalim = 1.d-5

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END
! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! --- Mesh viscosity modeling (cf. usvima)
!     0 : isotropic
!     1 : orthotropic
iortvm = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!----
! Formats
!----

!----
! End
!----

return
end subroutine


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

!===============================================================================

implicit none

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. Example of calculation options to modify
!===============================================================================

! Reading the meteo file

imeteo = 1

! Option for variable density and viscosity

irovar = 0
ivivar = 0

!----
! End
!----

return
end subroutine


!===============================================================================


subroutine usd3p1
!================


!===============================================================================
!  Features of this subroutine:
!  ----------------------------
!  1. Variable Output
!     a. Transported Variables
!     b. Variables of State; User definied Variables
!
!  2. Additional Calculation Options
!     a. Density Relaxation
!
!  3. Physical Constants
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

integer          ipp

!===============================================================================
!===============================================================================
! 1. Variable Output
!===============================================================================
!    Function                             |  Key Word |   Indicator
!    ---------------------------------------------------------------
!    Variable Output in the result file   | ICHRVR()  | yes= 1  ; no=0
!    Variable Output in the listing file  | ILISVR()  | yes= 1  ; no=0
!    Output of the temporal evolution of  | IHISVR()  | yes=-1* ; no=0
!    the variable at monitoring points    |           |
!    -----------------------------------------------------------------
!    *: Output for all monitoring points
!
!===============================================================================
! a. Transported Variables
!===============================================================================

! ---- Mean mixture fraction
ipp = ipprtp(isca(ifm))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Variance of mixture fraction
ipp = ipprtp(isca(ifp2m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Enthalpy
 if (ippmod(icod3p).eq.1) then
   ipp = ipprtp(isca(ihm))
   ichrvr(ipp)  = 1
   ilisvr(ipp)  = 1
   ihisvr(ipp,1)= -1
  endif


!===============================================================================
! b. Variables of State; User definied Variables
!===============================================================================

! ---- Temperature
ipp = ipppro(ipproc(itemp))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Fuel Mass fraction :    YM_Fuel
ipp = ipppro(ipproc(iym(1)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Oxydizer Mass fraction : YM_Oxy
ipp = ipppro(ipproc(iym(2)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Product Mass fraction : YM_Prod
ipp = ipppro(ipproc(iym(3)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Diffusion flame including gas radiation

if (iirayo.gt.0) then

! ---- Absorption Coefficient
  ipp = ipppro(ipproc(ickabs))
  nomvar(ipp)   = 'KABS'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

! ---- Term T^4
  ipp = ipppro(ipproc(it4m))
  nomvar(ipp)   = 'TEMP4'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

! ---- Term T^3
  ipp = ipppro(ipproc(it3m))
  nomvar(ipp)   = 'TEMP3'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

endif


!===============================================================================
! 2. Additional Calculation Options
!===============================================================================

! -->  Density Relaxation
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0


!===============================================================================
! 3. Physical Constants
!===============================================================================

!       DIFTL0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5


!----
! End
!----

return
end subroutine


!===============================================================================


subroutine usebu1

!===============================================================================
!  PURPOSE:
!  --------
!  1. Variable Output
!     a. Transported Variables
!     b. Variables of State; User definied Variables
!
!  2. Additional Calculation Options
!     a. Density Relaxation
!
!  3. Physical Constants
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

integer          ipp

!===============================================================================
!===============================================================================
! 1. Variable Output
!===============================================================================
!    Function                             |  Key Word |   Indicator
!    ---------------------------------------------------------------
!    Variable Output in the result file   | ICHRVR()  | yes= 1  ; no=0
!    Variable Output in the listing file  | ILISVR()  | yes= 1  ; no=0
!    Output of the temporal evolution of  | IHISVR()  | yes=-1* ; no=0
!    the variable at monitoring points    |           |
!    -----------------------------------------------------------------
!    *: Output for all monitoring points
!
!===============================================================================
! a. Transported Variables
!===============================================================================
! ---- Mass fraction of unburned (or fresh)  gas
if (ippmod(icoebu).ge.0) then
  ipp = ipprtp(isca(iygfm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

! ---- Mean Mixture Fraction
if (ippmod(icoebu).ge.2) then
  ipp = ipprtp(isca(ifm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


! ---- Enthalpy
if (ippmod(icoebu).eq.1 .or. ippmod(icoebu).eq.3) then
  ipp = ipprtp(isca(ihm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


!===============================================================================
! b. Variables of State; User definied Variables
!===============================================================================

! ---- Temperature
ipp = ipppro(ipproc(itemp))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Mean mass fraction of Fuel:    YM_Fuel
ipp = ipppro(ipproc(iym(1)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Mean mass fraction of Oxidizer : YM_Oxy
ipp = ipppro(ipproc(iym(2)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Mean mass fraction of Product: YM_Prod
ipp = ipppro(ipproc(iym(3)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! ---- Premixed flame including gas radiation

if (iirayo.gt.0) then

! ---- Absorption Coefficient
  ipp = ipppro(ipproc(ickabs))
  nomvar(ipp)   = 'KABS'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

! ---- Term T^4
  ipp = ipppro(ipproc(it4m))
  nomvar(ipp)   = 'TEMP4'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

! ---- Term T^3
  ipp = ipppro(ipproc(it3m))
  nomvar(ipp)   = 'TEMP3'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

endif


!===============================================================================
! 2. Additional Calculation Options
!===============================================================================

! -->  Density Relaxation
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0


!===============================================================================
! 3. Physical Constants
!===============================================================================

!       DIFTL0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5

!       cebu: EBU-model constant

 cebu   = 2.5d0


!----
! End
!----

return
end subroutine


!===============================================================================


subroutine uslwc1


!===============================================================================
!  PURPOSE:
!  --------
!  1. Variable Output
!     a. Transported Variables
!     b. Variables of State; User definied Variables
!
!  2. Additional Calculation Options
!     a. Density Relaxation
!
!  3. Physical Constants
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
! 1. Variable Output
!===============================================================================
!    Function                             |  Key Word |   Indicator
!    ---------------------------------------------------------------
!    Variable Output in the result file   | ICHRVR()  | yes= 1  ; no=0
!    Variable Output in the listing file  | ILISVR()  | yes= 1  ; no=0
!    Output of the temporal evolution of  | IHISVR()  | yes=-1* ; no=0
!    the variable at monitoring points    |           |
!    -----------------------------------------------------------------
!    *: Output for all monitoring points
!
!===============================================================================
! a. Transported Variables
!===============================================================================

! ---- Mean Mixture Fraction
if (ippmod(icolwc).ge.0) then
  ipp = ipprtp(isca(ifm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Variance of Mixture Fraction
  ipp = ipprtp(isca(ifp2m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Fuel Mass fraction
  ipp = ipprtp(isca(iyfm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

! ---- Variance of Fuel Mass fraction
  ipp = ipprtp(isca(iyfp2m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

if (ippmod(icolwc).ge.2) then
    ipp = ipprtp(isca(icoyfp))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
endif

! ---- Enthalpy
if (ippmod(icolwc).eq.1 .or.                                     &
    ippmod(icolwc).eq.3 .or.                                     &
    ippmod(icolwc).eq.5) then
  ipp = ipprtp(isca(ihm))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


!===============================================================================
! b. Variables of State; User definied Variables
!===============================================================================

! --- Source term
  ipp = ipppro(ipproc(itsc))
  nomvar(ipp)   = 'T.SOURCE'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- Temperature in K
  ipp = ipppro(ipproc(itemp))
  nomvar(ipp)   = 'Temperature'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- Fuel Mass fraction
  ipp = ipppro(ipproc(iym(1)))
  nomvar(ipp)   = 'YM_Fuel'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- Oxidizer Mass fraction
  ipp = ipppro(ipproc(iym(2)))
  nomvar(ipp)   = 'YM_Oxyd'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
! --- Products Mass fraction
  ipp = ipppro(ipproc(iym(3)))
  nomvar(ipp)   = 'YM_Prod'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  do idirac = 1, ndirac
    ipp = ipppro(ipproc(irhol(idirac)))
    write(nomvar(ipp),'(A4,I1)') 'RHOL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iteml(idirac)))
    write(nomvar(ipp),'(A4,I1)') 'TEML', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmel(idirac)))
    write(nomvar(ipp),'(A4,I1)') 'FMEL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(ifmal(idirac)))
    write(nomvar(ipp),'(A4,I1)') 'FMAL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(iampl(idirac)))
    write(nomvar(ipp),'(A4,I1)') 'AMPL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(itscl(idirac)))
    write(nomvar(ipp),'(A4,I1)') 'TSCL', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1

    ipp = ipppro(ipproc(imaml(idirac)))
    write(nomvar(ipp),'(A4,I1)') 'MAML', idirac
    ichrvr(ipp)   = 1
    ilisvr(ipp)   = 1
    ihisvr(ipp,1) = -1
  enddo

! ---- Premixed flame including gas radiation

if (iirayo.gt.0) then

! ---- Absorption Coefficient
  ipp = ipppro(ipproc(ickabs))
  nomvar(ipp)   = 'KABS'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

! ---- Term T^4
  ipp = ipppro(ipproc(it4m))
  nomvar(ipp)   = 'TEMP4'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

! ---- Term T^3
  ipp = ipppro(ipproc(it3m))
  nomvar(ipp)   = 'TEMP3'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

endif


!===============================================================================
! 2. Additional Calculation Options
!===============================================================================

! -->  Density Relaxation
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1))

srrom = 0.95d0


!===============================================================================
! 3. Physical Constants
!===============================================================================

! --> DIFTL0: Dynamic Diffusion Coefficient (kg/(m s))
diftl0 = 4.25d-5

! --> Constants of the Libby-Williams Model

! --- Reference velocity
 vref = 60.d0
! --- Reference length scale
 lref = 0.1d0
! --- Activation Temperature
 ta   = 0.2d5
! --- Cross-over Temperature (combustion of propane)
 tstar= 0.12d4

!----
! End
!----

return
end subroutine


!===============================================================================


subroutine uscfx1
!================


!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Initialize non standard options for the compressible flow scheme.


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

if(1.eq.1) then
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

! Specify if the hydrostatic equilibrium must be accounted for
!     (yes = 1 , no = 0)

icfgrp = 1

!----
! End
!----

return
end subroutine


!===============================================================================


subroutine uscfx2
!================


! Purpose:
! -------

!    User subroutine.

!    Set options for viscosity and conductivity for compressible flow.

!    In addition to options set in the user subroutine 'uscfx1' (or in
!    the GUI): this subroutine allows to set switches to indicate if the
!    volumetric viscosity and the conductivity are constants. If they are,
!    the subroutines allows to set their values.


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

if(1.eq.1) then
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

! --> Molecular thermal conductivity

!       constant  : ivisls = 0
!       variable  : ivisls = 1

ivisls(itempk) = 0

!       Reference molecular thermal conductivity
!       visls0 = lambda0  (molecular thermal conductivity, W/(m K))

!       WARNING: visls0 must be strictly positive
!         (set a realistic value here even if conductivity is variable)

visls0(itempk) = 3.d-2

!       If the molecular thermal conductivity is variable, its values
!         must be provided in the user subroutine 'uscfpv'


! --> Volumetric molecular viscosity

!       Reference volumetric molecular viscosity

!       viscv0 = kappa0  (volumetric molecular viscosity, kg/(m s))
!       iviscv = 0 : uniform  in space and constant in time
!              = 1 : variable in space and time

iviscv = 0
viscv0 = 0.d0

!       If the volumetric molecular viscosity is variable, its values
!         must be provided in the user subroutine 'uscfpv'


!----
! End
!----

return
end subroutine


!===============================================================================


subroutine uscpi1
!================


!===============================================================================
!  PURPOSE   :
!  ---------

!  User's routine to control outing of variables for pulverised coal combustion
!  (these parameters are in COMMON)

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

!===============================================================================

implicit none

integer          ipp , icla , icha

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  THIS TEST CERTIFY THIS VERY ROUTINE IS USED
!     IN PLACE OF LIBRARY'S ONE
!===============================================================================

if (iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ BEWARE : STOP during data inlet for pulverised coal     ',/,&
'@    =========                                               ',/,&
'@     THE USER SUBROUTINE uscpi1 have to be modified         ',/,&
'@                                                            ',/,&
'@  The computation will not start                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. TRANSPORTED VARIABLES
!===============================================================================

! OUTLET chrono, listing, and histo
!     if below vector are not allocated, default values will be used

!       ICHRVR( ) =  chono outlet (Yes 1/No  0)
!       ILISVR( ) =  listing outlet (Yes 1/No  0)
!       IHISVR( ) =  histo outlet (number of roiqu and number)
!       if IHISVR(.,1)  = -1 every monitoring point


! --> Variables for the mix (carrying gas and coal particles)

!      - Enthalpy
ipp = ipprtp(isca(ihm))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! --> Variables for coal particles

do icla = 1, nclacp

!       - Char mass fraction (in class ICLA)
  ipp = ipprtp(isca(ixck(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Coal mass fraction (in class ICLA)
  ipp = ipprtp(isca(ixch(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Number of particles for 1 kg mix (from class ICLA)
  ipp = ipprtp(isca(inp(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Enthalpy J/kg (for class ICLA)
  ipp = ipprtp(isca(ih2(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Water mass fraction (in class ICLA)
  if (ippmod(icp3pl) .eq. 1) then
    ipp = ipprtp(isca(ixwt(icla)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
enddo

! --> Variables for the carrier phase

do icha = 1, ncharb

!       - Mean of 1 mixture fraction
!         (from light volatiles of char ICHA)
  ipp = ipprtp(isca(if1m(icha)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Mean of 2 mixture fraction
!         (from heavy volatiles of char ICHA)
  ipp = ipprtp(isca(if2m(icha)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

enddo

!     - Mean of 3 mixture fraction
!       (C from heterogeneoux oxidation, of char, by O2)
ipp = ipprtp(isca(if3m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - Meam of (6 ?) mixture fraction
!       (C from heterogeneous reaction between char and CO2)
if ( ihtco2 .eq. 1) then
  ipp = ipprtp(isca(if3mc2))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!     - Variance of 4 mixture fraction
!       (oxidisers)
ipp = ipprtp(isca(if4p2m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - Mean of 5 mixture fraction
!       (water vapor from drying)
if (ippmod(icp3pl) .eq. 1) then
  ipp = ipprtp(isca(if5m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!     - Mass fraction of CO2 or CO (relaxation to equilibrium)

if (ieqco2 .ge. 1) then
  ipp = ipprtp(isca(iyco2))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!===============================================================================
! 2. Sate variables
!===============================================================================

! OUTLET chrono, listing, and histo
!     if below vector are not allocated, default values will be used

!       ICHRVR( ) =  chono outlet (Yes 1/No  0)
!       ILISVR( ) =  listing outlet (Yes 1/No  0)
!       IHISVR( ) =  histo outlet (number of roiqu and number)
!       if IHISVR(.,1)  = -1 every monitoring point

! --> State varables for the mix

!     - Mean Molar Mass
ipp = ipppro(ipproc(immel))
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

! --> State variables for coal particles

do icla = 1, nclacp

!       - Particles' Temperature K (of class ICLA)
  ipp = ipppro(ipproc(itemp2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Particles' Density kg/m3 (of class ICLA)
  ipp = ipppro(ipproc(irom2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Particles' Diameter m (of class ICLA)
  ipp = ipppro(ipproc(idiam2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Rate of coal consumption  (s-1) < 0
!         (for class ICLA)
  ipp = ipppro(ipproc(igmdch(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of light volatiles exhaust (s-1) < 0
!         (for class ICLA)
  ipp = ipppro(ipproc(igmdv1(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of heavy volatile exhaust (s-1) < 0
!         (de la classe ICLA)
  ipp = ipppro(ipproc(igmdv2(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of char oxidation by O2 (s-1) < 0
!         (from class ICLA)
  ipp = ipppro(ipproc(igmhet(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of char gazeification by CO2 (s-1) < 0
!         (from class ICLA)
  if (ihtco2 .eq. 1) then
    ipp = ipppro(ipproc(ighco2(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1
  endif

!       - Rate of drying (s-1) < 0
!         (from class ICLA)
  if (ippmod(icp3pl) .eq. 1) then
    ipp = ipppro(ipproc(igmsec(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1
  endif

!       - Mass fraction (of class ICLA) in mix
  ipp = ipppro(ipproc(ix2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

enddo

! --> State variables for carrier gas phase

!     - Temperature of gas mixture
ipp = ipppro(ipproc(itemp1))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Mass fraction (among gases) of  CHx1m
ipp = ipppro(ipproc(iym1(1)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CHx2m
ipp = ipppro(ipproc(iym1(2)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CO
ipp = ipppro(ipproc(iym1(3)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of O2
ipp = ipppro(ipproc(iym1(4)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CO2
ipp = ipppro(ipproc(iym1(5)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2O
ipp = ipppro(ipproc(iym1(6)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of N2
ipp = ipppro(ipproc(iym1(7)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1


!===============================================================================
! 3. Computation OPTION
!===============================================================================

! --- Relaxation for density (Advisable when starting combustion computation)
!                            (Forbidden for unstationnary computation)
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.95d0


!===============================================================================
! 4. Physical constants
!===============================================================================

! ---  Laminar viscosity for enthalpy (dynamical diffusivity) kg/(m.s)
diftl0 = 4.25d-5


!----
! End
!----

return

end subroutine


!===============================================================================


subroutine uscpl1
!================


!===============================================================================
!  FONCTION  :
!  ---------


!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!    ROUTINE UTILISATEUR POUR ENTREE DES PARAMETRES DE CALCUL
!      (COMMONS)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
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

!===============================================================================

implicit none

integer          ipp , icha

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     COMBUSTION CHARBON PULVERISE COUPLE AU                 ',/,&
'@     TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON :       ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uscpl1 DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

!  Sortie chrono, suivi listing, sortie histo
!     Si l'on n'affecte pas les tableaux suivants,
!     les valeurs par defaut seront utilisees

!       ICHRVR( ) = sortie chono (oui 1/non 0)
!       ILISVR( ) = suivi listing (oui 1/non 0)
!       IHISVR( ) = sortie historique (nombre de sondes et numeros)
!       si IHISVR(.,1)  = -1 sortie sur toutes les sondes


! --> Variables propres a la phase gaz continue

!      - Enthalpie de la phase gaz continue
ipp = ipprtp(isca(ihm))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! --> Variables propres a la phase continue

do icha = 1, ncharb

!       - Moyenne du traceur 1
!         (representatif des MV legeres du charbon ICHA)
  ipp = ipprtp(isca(if1m(icha)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Moyenne du traceur 2
!         (representatif des MV lourdes du charbon ICHA)
  ipp = ipprtp(isca(if2m(icha)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

enddo

!     - Moyenne du traceur 3 (representatif du C libere sous forme de CO
!       lors de la combustion heterogene)
ipp = ipprtp(isca(if3m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - Variance associe au traceur 4 (representatif de l'air)
ipp = ipprtp(isca(if4p2m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!===============================================================================
! 2. VARIABLES ALGEBRIQUES OU D'ETAT
!===============================================================================

!  Sortie chrono, suivi listing, sortie histo
!     Si l'on n'affecte pas les tableaux suivants,
!     les valeurs par defaut seront utilisees

!       ICHRVR( ) = sortie chono (oui 1/non 0)
!       ILISVR( ) = suivi listing (oui 1/non 0)
!       IHISVR( ) = sortie historique (nombre de sondes et numeros)
!       si IHISVR(.,1)  = -1 sortie sur toutes les sondes

! --> Variables algebriques propres a la suspension gaz - particules

!     - Masse molaire du melange gazeux
ipp = ipppro(ipproc(immel))
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

!     - Temperature du melange gazeux
ipp = ipppro(ipproc(itemp1))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Fraction massique (dans le melange gazeux) du CHx1m
ipp = ipppro(ipproc(iym1(1)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Fraction massique (dans le melange gazeux) du CHx2m
ipp = ipppro(ipproc(iym1(2)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Fraction massique (dans le melange gazeux) du CO
ipp = ipppro(ipproc(iym1(3)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Fraction massique (dans le melange gazeux) du O2
ipp = ipppro(ipproc(iym1(4)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Fraction massique (dans le melange gazeux) du CO2
ipp = ipppro(ipproc(iym1(5)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Fraction massique (dans le melange gazeux) du H2O
ipp = ipppro(ipproc(iym1(6)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Fraction massique (dans le melange gazeux) du N2
ipp = ipppro(ipproc(iym1(7)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1


!===============================================================================
! 3. OPTIONS DE CALCUL
!===============================================================================

! --- Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.8d0


!===============================================================================
! 4. CONSTANTES PHYSIQUES
!===============================================================================

! ---> Viscosite laminaire associee au scalaire enthalpie
!       DIFTL0 (diffusivite dynamique en kg/(m s))
diftl0 = 4.25d-5


!----
! End
!----

return

end subroutine


!===============================================================================


subroutine user_coal_ini1
!========================


!===============================================================================
!  PURPOSE   :
!  ---------
!  User's routine to control outing of variables for pulverised coal combustion
!  (these parameters are in COMMON)
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

integer          ipp , icla , icha

!===============================================================================


!===============================================================================
! 1. TRANSPORTED VARIABLES
!===============================================================================

! OUTLET chrono, listing, and histo
!     if below vector are not allocated, default values will be used

!       ICHRVR( ) =  chono outlet (Yes 1/No  0)
!       ILISVR( ) =  listing outlet (Yes 1/No  0)
!       IHISVR( ) =  histo outlet (number of roiqu and number)
!       if IHISVR(.,1)  = -1 every probe


! --> Variables for the mix (carrying gas and coal particles)

!      - Enthalpy
ipp = ipprtp(isca(ihm))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! --> Variables for coal particles

do icla = 1, nclacp

!       - Char mass fraction (in class ICLA)
  ipp = ipprtp(isca(ixck(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Coal mass fraction (in class ICLA)
  ipp = ipprtp(isca(ixch(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Number of particles for 1 kg mix (from class ICLA)
  ipp = ipprtp(isca(inp(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Enthalpy J/kg (for class ICLA)
  ipp = ipprtp(isca(ih2(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Water mass fraction (in class ICLA)
  if (ippmod(icp3pl) .eq. 1) then
    ipp = ipprtp(isca(ixwt(icla)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
enddo

! --> Variables for the carrier phase

do icha = 1, ncharb

!       - Mean of 1 mixture fraction
!         (from light volatiles of char ICHA)
  ipp = ipprtp(isca(if1m(icha)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Mean of 2 mixture fraction
!         (from heavy volatiles of char ICHA)
  ipp = ipprtp(isca(if2m(icha)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

enddo

! ---- Variables propres a la phase continue
  if (noxyd .ge. 2) then
    ipp = ipprtp(isca(if4m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if (noxyd .eq. 3) then
    ipp = ipprtp(isca(if5m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if (ippmod(iccoal) .ge. 1) then
    ipp = ipprtp(isca(if6m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  ipp = ipprtp(isca(if7m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  if (ihtco2 .eq. 1) then
    ipp = ipprtp(isca(if8m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if (ihth2o .eq. 1) then
    ipp = ipprtp(isca(if9m))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
!

  ipp = ipprtp(isca(ifvp2m))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!
  if (ieqco2 .ge. 1) then
    ipp = ipprtp(isca(iyco2))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
  if (ieqnox .ge. 1) then
    ipp = ipprtp(isca(iyhcn))
    nomvar(ipp)  = 'FR_HCN'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(iyno))
    nomvar(ipp)  = 'FR_NO'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
    ipp = ipprtp(isca(ihox))
    nomvar(ipp)  = 'Enth_Ox'
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif

!===============================================================================
! 2. Sate variables
!===============================================================================

! OUTLET chrono, listing, and histo
!     if below vector are not allocated, default values will be used

!       ICHRVR( ) =  chono outlet (Yes 1/No  0)
!       ILISVR( ) =  listing outlet (Yes 1/No  0)
!       IHISVR( ) =  histo outlet (number of roiqu and number)
!       if IHISVR(.,1)  = -1 every probe

! --> State varables for the mix

!     - Mean Molar Mass
ipp = ipppro(ipproc(immel))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! --> State variables for coal particles

do icla = 1, nclacp

!       - Particles' Temperature K (of class ICLA)
  ipp = ipppro(ipproc(itemp2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Particles' Density kg/m3 (of class ICLA)
  ipp = ipppro(ipproc(irom2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Particles' Diameter m (of class ICLA)
  ipp = ipppro(ipproc(idiam2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Rate of coal consumption  (s-1) < 0
!         (for class ICLA)
  ipp = ipppro(ipproc(igmdch(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of light volatiles exhaust (s-1) < 0
!         (for class ICLA)
  ipp = ipppro(ipproc(igmdv1(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of heavy volatile exhaust (s-1) < 0
!         (de la classe ICLA)
  ipp = ipppro(ipproc(igmdv2(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of char oxidation by O2 (s-1) < 0
!         (from class ICLA)
  ipp = ipppro(ipproc(igmhet(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Rate of char gazeification by CO2 (s-1) < 0
!         (from class ICLA)
  if (ihtco2 .eq. 1) then
    ipp = ipppro(ipproc(ighco2(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1
  endif

!       - Rate of char gazeification by H2O (s-1) < 0
!         (from class ICLA)
  if (ihth2o .eq. 1) then
    ipp = ipppro(ipproc(ighh2o(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1
  endif

!       - Rate of drying (s-1) < 0
!         (from class ICLA)
  if (ippmod(icp3pl) .eq. 1) then
    ipp = ipppro(ipproc(igmsec(icla)))
    ichrvr(ipp)   = 0
    ilisvr(ipp)   = 0
    ihisvr(ipp,1) = -1
  endif

!       - Mass fraction (of class ICLA) in mix
  ipp = ipppro(ipproc(ix2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

enddo

! --> State variables for carrier gas phase

!     - Temperature of gas mixture
ipp = ipppro(ipproc(itemp1))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Mass fraction (among gases) of  CHx1m
ipp = ipppro(ipproc(iym1(1)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CHx2m
ipp = ipppro(ipproc(iym1(2)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CO
ipp = ipppro(ipproc(iym1(3)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2S
ipp = ipppro(ipproc(iym1(4)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2
ipp = ipppro(ipproc(iym1(5)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of HCN
ipp = ipppro(ipproc(iym1(6)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of NH3
ipp = ipppro(ipproc(iym1(7)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of O2
ipp = ipppro(ipproc(iym1(8)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CO2
ipp = ipppro(ipproc(iym1(9)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2O
ipp = ipppro(ipproc(iym1(10)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of SO2
ipp = ipppro(ipproc(iym1(11)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of N2
ipp = ipppro(ipproc(iym1(12)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!===============================================================================
! 3. Computation OPTION
!===============================================================================

! --- Relaxation for density (Advisable when starting combustion computation)
!                            (Forbidden for unstationnary computation)
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.95d0

!===============================================================================
! 4. Physical constants
!===============================================================================

! ---  Laminar viscosity for enthalpy (dynamical diffusivity) kg/(m.s)
diftl0 = 4.25d-5

!----
! End
!----


end subroutine


!===============================================================================


subroutine user_fuel_ini1
!========================

!===============================================================================
!  PURPOSE   :
!  ---------

!  USER ROUTINE FOR ALLOCATE COMPUTATION PARAMETERS DEALING WITH FUEL
!    (COMMONS)

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

integer          ipp , icla

!===============================================================================

!===============================================================================
! 1. TRANSPORTED VARIABLES
!===============================================================================

! OUTLET chrono, listing, and histo
!     if below vector are not allocated, default values will be used

!       ICHRVR( ) =  chono outlet (Yes 1/No  0)
!       ILISVR( ) =  listing outlet (Yes 1/No  0)
!       IHISVR( ) =  histo outlet (number of roiqu and number)
!       if IHISVR(.,1)  = -1 every probe


! --> Variables for the mix (carrying gas and coal particles)

!      - Enthalpy

ipp = ipprtp(isca(ihm))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! --> Variables for droplets

do icla = 1, nclafu
!       - Fuel mass fraction
  ipp = ipprtp(isca(iyfol(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Number of droplets in mix (1/kg)
  ipp = ipprtp(isca(ing(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

!       - Fuel enthalpy (J/kg)
  ipp = ipprtp(isca(ih2(icla)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
enddo


! --> Variables for carrying gas

!       - Mean of 1 mixture fraction (fuel vapor)
ipp = ipprtp(isca(ifvap))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - Mean of 3 mixture fraction
!       (carbon from heterogeneous oxidation of char)
ipp = ipprtp(isca(if7m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - Variance of 4 mixture fraction (air)
ipp = ipprtp(isca(ifvp2m))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!     - YCO2

if (ieqco2 .ge. 1) then
  ipp = ipprtp(isca(iyco2))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!     - HCN and NO

if (ieqnox .eq. 1) then
  ipp = ipprtp(isca(iyhcn))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(iyno))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
  ipp = ipprtp(isca(ihox))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

!===============================================================================
! 2. State variables
!===============================================================================

! OUTLET chrono, listing, and histo
!     if below vector are not allocated, default values will be used

!       ICHRVR( ) =  chono outlet (Yes 1/No  0)
!       ILISVR( ) =  listing outlet (Yes 1/No  0)
!       IHISVR( ) =  histo outlet (number of roiqu and number)
!       if IHISVR(.,1)  = -1 every monitoring point


! --> Variables for the mix (carrying gas and coal particles)

!     - Mean Molar Mass of gases in kg
ipp = ipppro(ipproc(immel))
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

! --> Variables for droplets

do icla = 1, nclafu
!       - Droplets' Temperature in K
  ipp = ipppro(ipproc(itemp2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Droplet's Density in kg/m3
  ipp = ipppro(ipproc(irom2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Droplet's Diameter
  ipp = ipppro(ipproc(idiam2(icla)))
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!       - Heat flux (between gases and ICLA class droplets)
  ipp = ipppro(ipproc(ih1hlf(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Evaporation mass flow rate (s-1) < 0
  ipp = ipppro(ipproc(igmeva(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1

!       - Char combsution mass flow rate
  ipp = ipppro(ipproc(igmhtf(icla)))
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = -1
enddo

! --> State variables for carrying gas

!     - Temperature for gases only (not mixed with droplets)
ipp = ipppro(ipproc(itemp1))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!      - Nothing
ipp = ipppro(ipproc(iym1(1)))
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

!     -  Mass fraction of fuel vapor
!          (relative to pure gases : not mixed with droplets ..)
ipp = ipppro(ipproc(iym1(2)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CO
ipp = ipppro(ipproc(iym1(3)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2S
ipp = ipppro(ipproc(iym1(4)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2
ipp = ipppro(ipproc(iym1(5)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of HCN
ipp = ipppro(ipproc(iym1(6)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of NH3
ipp = ipppro(ipproc(iym1(7)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of O2
ipp = ipppro(ipproc(iym1(8)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of CO2
ipp = ipppro(ipproc(iym1(9)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of H2O
ipp = ipppro(ipproc(iym1(10)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of SO2
ipp = ipppro(ipproc(iym1(11)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - mass fraction (among gases) of N2
ipp = ipppro(ipproc(iym1(12)))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Carbone Bilan
ipp = ipppro(ipproc(ibcarbone))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Oxygen bilan
ipp = ipppro(ipproc(iboxygen))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     - Hydrogen bilan
ipp = ipppro(ipproc(ibhydrogen))
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!===============================================================================
! 3. Computation OPTION
!===============================================================================

! --- Relaxation for density (Advisable when starting combustion computation)
!                            (Forbidden for unstationnary computation)
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)

srrom = 0.7d0


!===============================================================================
! 4. Physical constants
!===============================================================================

! ---  Laminar viscosity for enthalpy (dynamical diffusivity) kg/(m.s)
diftl0 = 4.25d-5

!----
! End
!----

return

end subroutine


!===============================================================================


subroutine useli1
!================


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

!===============================================================================

implicit none

integer          ipp, iesp , idimve

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0. This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif
!
 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in data input',/,                          &
'@    =======',/,                                                 &
'@     The user subroutine ''useli1'' must be completed',/, &
'@     for electric module',/,                                    &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. SOLVED VARIABLES
!===============================================================================

!  Chronological output, logging in listing, history output
!       if we do not assign the following array values,
!       default values will be used!
!
!     ichrvr( ) = chonological output (yes 1/no 0)
!     ilisvr( ) = logging in listing (yes 1/no 0)
!     ihisvr( ) = history output (number of probes and their numbers)
!     if ihisvr(.,1)  = -1, output for all probes
!
! --> Current variables for electric modules

! ---- Enthalpy
ipp = ipprtp(isca(ihm))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Real component of the electrical potential
ipp = ipprtp(isca(ipotr))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

!---- Mass fraction of the different constituants of the phase
if (ngazg .gt. 1) then
  do iesp = 1, ngazg-1
    ipp = ipprtp(isca(iycoel(iesp)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
endif

! --> Specific variables for Joule effect for direct conduction
!     Imaginary component of electrical potential
if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
  ipp = ipprtp(isca(ipoti))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

! --> Specific variables for electric arc in 3D
!     vector potential components
if (ippmod(ielarc) .ge. 2) then
  do idimve = 1, ndimve
    ipp = ipprtp(isca(ipotva(idimve)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
endif

! --> Ionic conduction module
!     Not available in the present version of the code

!===============================================================================
! 2. Algebric or state variables
!===============================================================================

! ---- Temperature
ipp = ipppro(ipproc(itemp))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Electric conductivity
ipp = ipppro(ipproc(ivisls(ipotr)))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Joule effect Power
ipp = ipppro(ipproc(iefjou))
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! ---- Real component of the current density
do idimve = 1, ndimve
  ipp = ipppro(ipproc(idjr(idimve)))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
enddo

! ---- Imaginary component of the current density
if (ippmod(ieljou).eq.4) then
  do idimve = 1, ndimve
    ipp = ipppro(ipproc(idji(idimve)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo
endif

if (ippmod(ielarc).ge.1) then

! ---- Electromagnetic Forces (Laplace forces)
  do idimve = 1, ndimve
    ipp = ipppro(ipproc(ilapla(idimve)))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  enddo

! ---- Absorption oefficient  or Radiative sources term
  if (ixkabe.gt.0) then
    ipp = ipppro(ipproc(idrad))
    ichrvr(ipp)  = 1
    ilisvr(ipp)  = 1
    ihisvr(ipp,1)= -1
  endif
endif

! ---- Electric charge (volumic)
if (ippmod(ielion).ge.1) then
  ipp = ipppro(ipproc(iqelec))
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif


!===============================================================================
! 3. Calculation options
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


!----
! End
!----

return
end subroutine


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

double precision cpa,cpe,cpv,hv0,rhoe,visc,conduc

!===============================================================================

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@                      MODULE AEROREFRIGERANTS               ',/,&
'@                                                            ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usctin DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@     Ce sous-programme utilisateur permet de definir les    ',/,&
'@       options generales. Il est indispensable.             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  PARAMETRES POUR L'ECART DE TEMPERATURE IMPOSE
!===============================================================================

!     ACTIVATION
iaeeri = 0

!     ECART DE REFRIGERATION A IMPOSER
vaeeri = 13.d0

!     FREQUENCE DE MODIFICATION DE LA TEMPERATURE
iaeerp = 5

!     PAS DE TEMPERATURE POUR LE CALCUL DE LA PENTE DE ECARTREF(TEAU)
paseri = 0.015d0

!     MAXIMUM DE LA TEMPERATURE D'EAU CHAUDE MOYENNE PONDEREE
aetemx = 80.d0

!     MINIMUM DE LA TEMPERATURE D'EAU REFROIDIE MOYENNE PONDEREE
aetemn = 10.d0

!     NOMBRE DE ZONES D'ECHANGES AYANT UNE FRONTIERE ENTREE EAU

nbzsup = 2

!     LISTE DES NBZSUP ZONES D'ECHANGES EN BORD DE L'ENTREE EAU

lizsup(1) = 1
lizsup(2) = 2

!     NOMBRE DE ZONES D'ECHANGES AYANT UNE FRONTIERE SORTIE EAU
nbzinf = 2

!     LISTE DES NBZINF ZONES D'ECHANGES EN BORD DE LA SORTIE EAU

lizinf(1) = 1
lizinf(2) = 2

!     INSTANT ACTIVATION ECART IMPOSE

inbaei = 1000.D0

!===============================================================================
! 2.  POST-PROCESSING DES ZONES D'ECHANGES
!===============================================================================

ichrze = 1

!===============================================================================
! 3.  CALCUL SUITE AEROREFRIGERANT
!===============================================================================

isuict = isuite

!===============================================================================
! 4.  PROPRIETES DE L'AIR
!===============================================================================

! Il est deconseille de modifier ici ces proprietes

cpa    = 1006.0d0
cpv    = 1831.0d0
cpe    = 4179.0d0
hv0    = 2501600.0d0
rhoe   = 997.85615d0
visc   = 1.765d-5
conduc = 0.02493d0

call ctprof &
!==========
( cpa, cpv, cpe, hv0, rhoe, visc, conduc, gx, gy, gz )

!----
! End
!----

return
end subroutine
