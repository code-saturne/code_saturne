!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

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

!-------------------------------------------------------------------------------


!===============================================================================


subroutine usipph &
!================

 ( nphmax, nphas , iihmpu, nfecra , iturb , icp , iverif )


!===============================================================================
! Purpose:
! --------

! User subroutine for input of parameters depending on the number of phases.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nphmax           ! i  ! <-- ! maximum number of phases                       !
! nphas            ! i  ! <-- ! number of active phases                        !
! iihmpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! iturb(nphmax)    ! ia ! <-> ! turbulence model                               !
! icp(nphmax)      ! ia ! <-> ! flag for uniform Cp or not                     !
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

integer nphmax, nphas, iihmpu, nfecra
integer iturb(nphmax), icp(nphmax)
integer iverif

! Local variables

integer iphas

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
'@     The user subroutine ''usipph'' must be completed',/, &
'@       in file usini1.f90',/,                                   &
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


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!===============================================================================

! --- Turbulence (for each phase)
!       0...Laminar
!      10...Mixing length
!      20...k-epsilon
!      21...k-epsilon (linear production)
!      30...Rij-epsilon, (standard LRR)
!      31...Rij-epsilon (SSG)
!      40...LES (Smagorinsky)
!      41...LES (Dynamic)
!      42...LES (WALE)
!      50...v2f (phi-model)
!      60...k-omega SST
!  For 10, contact the development team before use

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iphas = 1
iturb(iphas) = 20

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Variable specific heat (ICP=1) or not (ICP=0)
!       for each phase IPHAS

!     Should be set only if specific physics (coal, combustion, electric arcs)
!       ARE NOT activated.

!     For these specific physics, ICP MUST NOT BE MODIFIED here, and the
!       following options are forced:
!          coal and combustion: constant CP constant;
!          electric arcs:       variable CP.

!     Caution:    complete usphyv with the law defining Cp
!     =========   if and only if variable Cp has been selected here
!                 (with icp(iphas)=1)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iphas = 1
icp(iphas) = 0

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
'@     The user subroutine ''usinsc'' must be completed',/, &
'@       in file usini1.f90',/,                                   &
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


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!===============================================================================

! --- Number of USER scalars (thermal or not, and whatever their carrier phase).
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

!     The maximum number of scalars is defined by 'nscamx' in paramx.h;
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
'@     The user subroutine ''usipsc'' must be completed',/, &
'@       in file usini1.f90',/,                                   &
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


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!===============================================================================

! --- Variance of a USER scalar:
!     If we wish a user scalar j to represent the variance of a
!       user scalar k, we set
!       iscavr(j) = k.
!     The values taken by iscavr are thus naturally greater or equal to 1
!       and less than or equal to the total number of scalars.
!       So, if we set iscavr(j) = k, we must have
!       0 < j < nscaus+1, 0< k < nscaus+1 and j different from k.

!     For example for user scalar 3 to be the variance of user scalar 3,
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

 ( nphmax, nesmax,                                                &
   iespre, iesder, iescor, iestot,                                &
   nphas , iihmpu, nfecra,                                        &
   idtvar, ipucou, iphydr, ialgce , iescal , iverif ,             &
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
! nphmax           ! i  ! <-- ! maximum number of phases                       !
! nesmax           ! i  ! <-- ! maximum number of error estimators per phase   !
! iespre           ! i  ! <-- ! number of the prediction error estimator       !
! iesder           ! i  ! <-- ! number of the derivative error estimator       !
! iescor           ! i  ! <-- ! number of the correction error estimator       !
! iestot           ! i  ! <-- ! number of the total error estimator            !
! nphas            ! i  ! <-- ! number of active phases                        !
! iihmpu           ! i  ! <-- ! indicates if the XML file from the GUI is      !
!                  !    !     ! used (1: yes, 0: no)                           !
! nfecra           ! i  ! <-- ! Fortran unit number for standard output        !
! idtvar           ! i  ! --> ! variable time step flag                        !
! ipucou           ! i  ! --> ! reinforced u-p coupling flag                   !
! iphydr           ! i  ! --> ! flag for handling of the equilibrium between   !
!                  !    !     ! the pressure gradient and the gravity and      !
!                  !    !     ! head-loss terms                                !
! ialgce           ! i  ! <-- ! option for the method of calculation of        !
!                  !    !     !  cell centers                                  !
! iescal           ! ia ! <-- ! flag for activation of error estimators for    !
!  (nesmax,nphmax) !    !     ! Navier-Stokes                                  !
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

integer nphmax, nesmax
integer iespre, iesder, iescor, iestot
integer nphas , iihmpu, nfecra
integer idtvar, ipucou, iphydr, ialgce
integer iescal(nesmax,nphmax)
integer iverif, icwfps

double precision cwfthr

! Local variables

integer iphas

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
'@     The user subroutine ''usipgl'' must be completed',/, &
'@       in file usini1.f90',/,                                   &
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


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

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
!     Only in single-phase

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ipucou = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Handling of hydrostatic pressure
!                               (0 : usual algorithm
!                                1 : specific handling)
!     Only in single-phase

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iphydr = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Estimators for Navier-Stokes (non-frozen velocity field)
!     We recommend running a calculation restart on a few time steps
!       with the activation of the most interesting of those.
!        (=2 to activate, =0 to deactivate).

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iphas = 1
!       div(rho u) -Gamma
iescal(iescor,iphas) = 0
!       resolution precision for the momentum
iescal(iestot,iphas) = 0

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

integer iphas, iutile, ii, jj, imom

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
'@     The user subroutine ''usipsu'' must be completed',/, &
'@       in file usini1.f90',/,                                   &
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


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!===============================================================================


! Calculation options (optcal.h)
! ==============================

! --- Calculation restart: isuite (= 1) or not (0)
!     In case of restart, read auxiliary restart file ileaux (= 1) or not (0).

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

isuite = 0
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

!     If a USER scalar represents the temperature or enthalpy (of phase iphas):
!       we define the number of this scalar in iscalt(iphas) and
!       we set iscsth(iscalt(iphas)) = 1 if it is the temperature
!          or  iscsth(iscalt(iphas)) = 2 if it is the enthalpy.

!     If no scalar represents the temperature or enthalpy (of phase iphas)
!       we set iscalt(iphas) = -1
!       and we do not define iscsth(iscalt(iphas)).


!     For the radiative module when used without specific physics, if we
!      have chosen to solve in temperature (that is if
!      iscsth(iscalt(iphas)) = 1), the fluid temperature is considered to
!      be in degrees KELVIN (be careful for boundary conditions an expression
!      of physical properties depending on temperature).
!      Nonetheless, even though it is not recommended, if we wish for the
!      fluid solver to work with a temperature in degrees Celsius, we must set
!      iscsth(iscalt(iphas)) = -1.
!      This choice is a source of user errors. Indeed, the boundary conditions
!      for the fluid temperature will then be in degrees Celsius, while the
!      boundary conditions for radiation in usray2 must still be in Kelvin.


!    If specific physics are not activated
!       (coal, combustion, electric arcs: see usppmo):

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

if (nmodpp.eq.0) then

  iphas = 1

  ! Number of the scalar representing temperature or enthalpy,
  !   or -1 if there is none.
  ! When the choice is done by the Code_Saturne GUI, the scalar representing
  !   the temperature or enthalpy is always the first.
  iscalt(iphas) = -1

! If there is a temperature or enthalpy variable:
  if (iscalt(iphas).gt.0) then
    ! we indicate if it is the temperature (=1) or the enthalpy (=2).
    iscsth(iscalt(iphas)) = 1
  endif

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Calculation (restart) with frozen velocity field (1 yes, 0 no)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iccvfg = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Vortex method for inlet conditions in L.E.S.
!       (0: not activated,  1: activated)
!     The vortex method only regards the L.E.S. models
!       and is only valid with one phase.
!     To use the vortex method, edit the 'usvort.f90' user file.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iphas = 1
if (itytur(iphas).eq.4) then
  ivrtex = 0
endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Convective scheme

!     blencv = 0 for upwind (order 1 in space, "stable but diffusive")
!            = 1 for centered/second order (order 2 in space)
!       we may use intermediate real values.
!       Here we choose:
!         for the velocity of phase 1 and user scalars:
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
!          blencv(isca(iscalt(iphas))) = 1.0d0

!       For non-user scalars relative to specific physics (coal, combustion,
!         electric arcs: see usppmo) implicitly defined by the model,
!         the corresponding information is set automatically elsewhere:
!         we do not modify blencv here.


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iphas = 1

blencv(iu(iphas)) = 1.0d0
blencv(iv(iphas)) = 1.0d0
blencv(iw(iphas)) = 1.0d0
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

!     nitmax: maximum number of iterations for each unknown ivar
!     epsilo: relative precision for the solution of the linear system.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iutile = 0
if (iutile.eq.1) then

  iphas = 1
  iresol(iu(iphas)) = 2
  iresol(iv(iphas)) = 2
  iresol(iw(iphas)) = 2
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

iphas = 1
imgr(ipr(iphas)) = 1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


!=========================================================================

! --- Stabilization in turbulent regime

!     For difficult cases, a stabilization may be obtained by not
!     reconstructing the convective and diffusive flux for variables
!     of the turbulence model, that is
!       in k-epsilon: if (itytur(iphas).eq.2) then
!          ircflu(ik(iphas))   = 0 and ircflu(iep(iphas))  = 0
!       in Rij-epsilon: if (itytur(iphas).eq.3) then
!          ircflu(ir11(iphas)) = 0,    ircflu(ir22(iphas)) = 0,
!          ircflu(ir33(iphas)) = 0,
!          ircflu(ir12(iphas)) = 0,    ircflu(ir23(iphas)) = 0,
!          ircflu(ir23(iphas)) = 0,
!                                  and ircflu(iep(iphas))  = 0
!     (note that variable itytur is equal to iturb/10)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     The test on iutile allows deactivation of the instructions
!       (which are only given as an example).

iutile = 0
if (iutile.eq.1) then

  iphas = 1
  if (iturb(iphas).eq.20) then
    ircflu(ik(iphas))   = 0
    ircflu(iep(iphas))  = 0
  endif

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! Physical constants (cstphy.h)
! =============================

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

! --- Reference fluid properties (for each phase)

!       ro0        : density in kg/m3
!       viscl0     : dynamic viscosity in kg/(m s)
!       cp0        : specific heat in J/(degres kg)
!       t0         : reference temperature in Kelvin
!       p0         : total reference pressure in Pascal
!                    the calculation is based on a
!                    reduced pressure P*=Ptot-ro0*g.(x-xref)
!                    (except in compressible case)
!       xyzp0(3,.) : coordinates of the reference point for
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

iphas = 1

ro0(iphas)    = 0.235d0
viscl0(iphas) = 0.84d-6
cp0(iphas)    = 1219.d0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

iphas = 1

t0(iphas) = 1000.d0 + 273.15d0
p0(iphas) = 1.01325d5
! We only specify XYZ0 if we explicitely fix Dirichlet conditions
! for the pressure.
! xyzp0(1,iphas) = 0.d0
! xyzp0(2,iphas) = 0.d0
! xyzp0(3,iphas) = 0.d0


! --- irovar, ivivar: density and viscosity constant or not ?

!     When a specific physics module is active
!       (coal, combustion, electric arcs, compressible: see usppmo)
!       we DO NOT set variables 'irovar' and 'ivivar' here, as
!       they are defined automatically.
!     Nonetheless, for the compressible case, ivivar may be modified
!       in the uscfx1 user subroutine.

!     When no specific physics module is active, it is necessary to
!       specify is the density and the molecular viscosity
!         are constant (irovar=0, ivivar=0)
!          or variable (irovar=1, ivivar=1)

!       if they are variable, the law must be defined in usphyv;
!       if they are constant, they take values ro0 and viscl0.

!       as an example, we assume below that they are constant.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

if (nmodpp.eq.0) then
  iphas = 1
  irovar(iphas) = 0
  ivivar(iphas) = 0
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
!       visls0(iscalt(iphas)) = Lambda/Cp

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
      visls0(jj) = viscl0(iphsca(jj))
    endif
  enddo

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- Reference velocity for turbulence initialization (m2/s)
!       (useful only with turbulence)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iphas = 1
uref(iphas)    = 1.d0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- Reference length scale in meters for initialization
!       of epsilon (and specific clipping of turbulence, but
!       this is not the default option)
!       Assign a value of the order of the largest dimension of the
!       physical domain in which the flow may develop.
!       (useful only for turbulence).

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

iphas = 1
almax(iphas) = -grand

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


!     The test on iutile allows deactivation of the instructions
!       (which are only given as an example).

iutile = 0
if (iutile.eq.1) then

  ! First moment: <u>
  imom  = 1
  iphas = 1
  idfmom(1,imom) =  iu(iphas)
  ntdmom(imom)   =  1000
  ! Second moment: <rho u v>
  imom  = 2
  iphas = 1
  idfmom(1,imom) = -irom(iphas)
  idfmom(2,imom) =  iu(iphas)
  idfmom(3,imom) =  iv(iphas)
  imoold(imom)   = -1
  ntdmom(imom)   =  10000

endif

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

integer ii, iphas, ipp, imom, iutile

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
'@     The user subroutine ''usipes'' must be completed',/, &
'@       in file usini1.f90',/,                                   &
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


! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

!     If we are using the Code_Saturne GUI:

!       we will find in the user subroutines commented examples
!       on the model of the present section.

!       If necessary, the user may uncomment them and adapt them to
!       his needs.

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!===============================================================================

!===============================================================================
! 1. Input-output (entsor.h)
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

  iphas = 1

  do ii = 1, nvar
    iwarni(ii) = 1
  enddo

  iwarni(ipr(iphas)) = 2
  iwarni(iu(iphas)) = 2
  iwarni(iv(iphas)) = 2
  iwarni(iw(iphas)) = 2

endif

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

! --- post-processing output

!     ichrvl: post-processing of the fluid domain (yes 1/no 0)
!     ichrbo: post-processing of the domain boundary (yes 1/no 0)
!     ichrsy: post-processing of zones coupled with SYRTHES (yes 1/ no 0)
!     ichrmd: indicates if the meshes output are:
!               0: fixed,
!               1: deformable with constant connectivity,
!               2: modifyable (may be completely redefined during the
!                  calculation using the usmpst subroutine).
!              10: as indmod = 0, with a displacement field
!              11: as indmod = 1, with a displacement field
!              11: as indmod = 2, with a displacement field

!     fmtchr: output format, amid
!               'EnSight Gold', 'MED', or 'CGNS'
!     optchr: options associated with the output format, separated by
!             commas, from the following list:
!               'text'              (text format, for EnSight)
!               'binary'            (binary format, default choice)
!               'big_endian'        (forces binary EnSight output to
!                                   'big-endian' mode)
!               'discard_polygons'  (ignore polygon-type faces)
!               'discard_polyhedra' (ignore polyhedron-type cells)
!               'divide_polygons'   (subdivides polygon-type faces)
!               'divide_polyhedra'  (subdivides polyhedron-type cells)
!               'split_tensors'     (writes tensors as separate scalars)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ichrvl = 1
ichrbo = 0
ichrsy = 0

ichrmd = 0

fmtchr = 'EnSight Gold'
optchr = 'binary'

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END


! --- chronological output step
!       (-1: only one value at calculation end)
!       (strictly positive value: output periodicity)

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

ntchr = -1
frchr = -1.d0

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

! Examples for phase 1
iphas = 1

! pressure variable
ipp = ipprtp(ipr   (iphas))
nomvar(ipp)   = 'Pressure'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
if (icorio.eq.1) then
  nomvar(ipp)   = 'Rel Pressure'
endif

! variable v1x
ipp = ipprtp(iu    (iphas))
nomvar(ipp)   = 'VelocityX'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
if (icorio.eq.1) then
  nomvar(ipp)   = 'Rel VelocityX'
endif

! v1y variable
ipp = ipprtp(iv    (iphas))
nomvar(ipp)   = 'VelocityY'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
if (icorio.eq.1) then
  nomvar(ipp)   = 'Rel VelocityY'
endif

! v1z variable
ipp = ipprtp(iw    (iphas))
nomvar(ipp)   = 'VelocityZ'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1
if (icorio.eq.1) then
  nomvar(ipp)   = 'Rel VelocityZ'
endif

if (itytur(iphas).eq.2) then

  ! turbulent kinetic energy
  ipp = ipprtp(ik    (iphas))
  nomvar(ipp)   = 'Turb Kinetic Energy'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! turbulent dissipation
  ipp = ipprtp(iep   (iphas))
  nomvar(ipp)   = 'Turb Dissipation'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif (itytur(iphas).eq.3) then

  ! Reynolds stresses
  ipp = ipprtp(ir11  (iphas))
  nomvar(ipp)   = 'R11'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir22  (iphas))
  nomvar(ipp)   = 'R22'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir33  (iphas))
  nomvar(ipp)   = 'R33'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir12  (iphas))
  nomvar(ipp)   = 'R12'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir13  (iphas))
  nomvar(ipp)   = 'R13'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! Reynolds stresses
  ipp = ipprtp(ir23  (iphas))
  nomvar(ipp)   = 'R23'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! turbulent dissipation
  ipp = ipprtp(iep   (iphas))
  nomvar(ipp)   = 'Turb Dissipation'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif (iturb(iphas).eq.50) then

  ! turbulent kinetic energy
  ipp = ipprtp(ik    (iphas))
  nomvar(ipp)   = 'Turb Kinetic Energy'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! turbulent dissipation
  ipp = ipprtp(iep   (iphas))
  nomvar(ipp)   = 'Turb Dissipation'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! phi
  ipp = ipprtp(iphi  (iphas))
  nomvar(ipp)   = 'Phi'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! f_bar
  ipp = ipprtp(ifb   (iphas))
  nomvar(ipp)   = 'f_bar'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif (iturb(iphas).eq.60) then

  ! turbulent kinetic energy
  ipp = ipprtp(ik    (iphas))
  nomvar(ipp)   = 'Turb Kinetic Energy'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

  ! omega
  ipp = ipprtp(iomg  (iphas))
  nomvar(ipp)   = 'Omega'
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

iphas = 1

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! Density variable (output for post-processing only if variable or
!                   in the case of specific physics)
ipp = ipppro(ipproc(irom  (iphas)))
nomvar(ipp)   = 'Density'
ichrvr(ipp)   = max(irovar(iphas),nmodpp)
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! specific heat
if (icp   (iphas).gt.0) then
  ipp = ipppro(ipproc(icp   (iphas)))
  nomvar(ipp)   = 'Specific Heat'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = 0
endif

! laminar viscosity
ipp = ipppro(ipproc(iviscl(iphas)))
nomvar(ipp)   = 'Laminar Viscosity'
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = 0

! turbulent viscosity
ipp = ipppro(ipproc(ivisct(iphas)))
nomvar(ipp)   = 'Turb Viscosity'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! Courant number
ipp = ipppro(ipproc(icour(iphas)))
nomvar(ipp)   = 'CFL'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

! Fourier number
ipp = ipppro(ipproc(ifour(iphas)))
nomvar(ipp)   = 'Fourier Number'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

! 'csmago' variable for dynamic L.E.S. models
!    (square of the Samgorinsky "constant")
if (ismago(iphas).gt.0) then
  ipp = ipppro(ipproc(ismago(iphas)))
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
  ipp = ipppro(ipproc(iprtot(iphas)))
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


subroutine ustbtr &
!================

 ( ncel   , ncelet , nfac   , nfabor , nnod   , longia , longra )

!===============================================================================
! Purpose:
! -------

! User subroutine to define the sizes of macro-arrays ia and ra.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncel             ! i  ! <-- ! number of cells                                !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nnod             ! i  ! <-- ! number of vertices                             !
! longia           ! i  ! --> ! size of array ia                               !
! longra           ! i  ! --> ! size of array ra                               !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

integer          ncel  , ncelet, nfac  , nfabor, nnod
integer          longia, longra

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Size of macro arrays ia and ra:

!  The user may need to modify the size of integer and real work
!    arrays here: longia and longra respectively.

!  The number of integers 'longia' and the number of reals 'longra' depend
!    on calculation options, on the element type and mesh characteristics
!    (2d, 3d, hybrid, non-conforming, ...) and on the number of variables.
!  In k-epsilon, if we note 'ncel' the local number of cells in the mesh,
!    we ay usually use the following coarse overestimation:
!    longia = 45*ncel and longra = 220*ncel. In Rij-epsilon, an additional
!    20% may be applied.
!  These values are relatively high so as to account for 2D meshes which
!    have many boundary faces. A more precise but complex formula would be
!    necessary. For large 3D cases, a more precise estimation is given by:
!    longia = 25*ncel  and longra = 120*ncel.

!  If longia and longra are left at 0, then these values are estimated
!    and set automatically.

!===============================================================================

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_START

! Size of main integer work array 'ia'

longia = 0

! Size of main real work array 'ra'

longra = 0

! EXAMPLE_CODE_TO_BE_ADAPTED_BY_THE_USER_END

!----
! Formats
!----


return
end subroutine
