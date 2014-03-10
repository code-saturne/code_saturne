!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

use coincl
use cs_coal_incl
use cs_fuel_incl
use cpincl
use elincl
use ppcpfu
use radiat

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
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

! Frequency of log output

if (.false.) then
!< [init_01]
  ntlist = 1
!< [init_01]


endif

! Log (listing) verbosity

if (.false.) then

!< [init_02]
  do ii = 1, nvar
    iwarni(ii) = 1
  enddo

  iwarni(ipr) = 2
  iwarni(iu) = 2
  iwarni(iv) = 2
  iwarni(iw) = 2
!< [init_02]

endif

! --- probes output step

if (.false.) then

!< [init_03]
  nthist = 1
  frhist = -1.d0
!< [init_03]

endif

! --- Number of monitoring points (probes) and their positions
!     (limited to ncaptm=100)

if (.false.) then

!< [init_04]
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
!< [init_04]

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
integer       fldid, idim1, iflpst, itycat, ityloc

!===============================================================================

! Example: force postprocessing of projection of some variables at boundary
!          with no reconstruction.
!          This is handled automatically if the second bit of a field's
!          'post_vis' key value is set to 1 (which amounts to adding 2
!          to that key value).
!
!          field_get_id returns -1 if field does not exist

!< [example_1]
fldid = ivarfl(iu)
if (iand(iflpst, 2) .eq. 0) then
  call field_get_key_int(fldid, keyvis, iflpst)
  iflpst = ior(iflpst, 2)
  call field_set_key_int(fldid, keyvis, iflpst)
endif

fldid = ivarfl(ipr)
if (iand(iflpst, 2) .eq. 0) then
  call field_get_key_int(fldid, keyvis, iflpst)
  iflpst = ior(iflpst, 2)
  call field_set_key_int(fldid, keyvis, iflpst)
endif
!< [example_1]

!-------------------------------------------------------------------------------

! Example: enforce existence of 'tplus' and 'tstar' fields, so that
!          a boundary temperature or Nusselt number may be computed using the
!          post_boundary_temperature or post_boundary_nusselt subroutines.
!          When postprocessing of these quantities is activated, those fields
!          are present, but if we need to compute them in the
!          cs_user_extra_operations user subroutine without postprocessing them,
!          forcing the definition of these fields to save the values computed
!          for the boundary layer is necessary.

!< [example_2]
itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 3 ! boundary faces
ilved = .true. ! interleaved
inoprv = .false. ! no previous time step values needed

call field_get_id('tplus', fldid)
if (fldid.lt.0) then
  call field_create('tplus', itycat, ityloc, idim1, ilved, inoprv, fldid)
endif

call field_get_id('tstar', fldid)
if (fldid.lt.0) then
  call field_create('tstar', itycat, ityloc, idim1, ilved, inoprv, fldid)
endif
!< [example_2]

return

!===============================================================================

!----
! Formats
!----

return
end subroutine user_field_parameters

!===============================================================================

