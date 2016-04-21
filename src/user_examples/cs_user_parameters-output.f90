!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_parameters-ouput.f90
!>
!> \brief Output parameters example
!>
!>  See \subpage f_parameters for examples.
!>
!>   If the Code_Saturne GUI is used, this file is not required (but may be
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
use parall
use period
use ihmpre
use ppppar
use ppthch
use ppincl
use coincl
use cpincl
use field
use cavitation

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

logical       inoprv
integer       f_id, idim1, itycat, ityloc, iscdri, iscal, ivar, kislts

!===============================================================================

!     This subroutine allows setting parameters

!       which do not already appear in the other subroutines of this file.


!     It is possible to add or remove parameters.


!     The number of physical properties and variables is known here.

!===============================================================================

! Enforce existence of 'tplus' and 'tstar' fields, so that
! a boundary temperature or Nusselt number may be computed using the
! post_boundary_temperature or post_boundary_nusselt subroutines.
! When postprocessing of these quantities is activated, those fields
! are present, but if we need to compute them in the
! cs_user_extra_operations user subroutine without postprocessing them,
! forcing the definition of these fields to save the values computed
! for the boundary layer is necessary.

!< [usipsu_ex_1]
itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 3 ! boundary faces
inoprv = .false. ! no previous time step values needed

call field_get_id('tplus', f_id)
if (f_id.lt.0) then
  call field_create('tplus', itycat, ityloc, idim1, inoprv, f_id)
endif

call field_get_id('tstar', f_id)
if (f_id.lt.0) then
  call field_create('tstar', itycat, ityloc, idim1, inoprv, f_id)
endif
!< [usipsu_ex_1]

!< [usipsu_ex_2]
call field_get_key_id("slope_test_upwind_id", kislts)

do ivar = 1, nvar
  call field_set_key_int(ivarfl(ivar), kislts, 0)
enddo
!< [usipsu_ex_2]

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

use coincl
use cs_coal_incl
use cs_fuel_incl
use cpincl
use ppcpfu
use radiat

!===============================================================================

implicit none

! Arguments

integer nmodpp

! Local variables

integer ii
integer f_id, idim1, ipp, iflpst, ifllog

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

!< [usipes_ex_01]
ntlist = 1
!< [usipes_ex_01]

! Log (listing) verbosity

!< [usipes_ex_02]
do ii = 1, nvar
  iwarni(ii) = 1
enddo

iwarni(ipr) = 2
iwarni(iu) = 2
iwarni(iv) = 2
iwarni(iw) = 2
!< [usipes_ex_02]

! Logging of variables and properties
! (example for velocity: 1 to activate logging output, 0 to deactivate)

!< [usipes_ex_03]
f_id = ivarfl(iu)
ifllog = 0
call field_set_key_int(f_id, keyvis, ifllog)
!< [usipes_ex_03]

! Change a property's label
! (here for specific heat, first checking if it is variable)

!< [usipes_ex_04]
if (icp.ne.0) then
  f_id = iprpfl(icp)
  call field_set_key_str (f_id, keylbl, 'Cp')
endif
!< [usipes_ex_04]

! --- probes output step

!< [usipes_ex_05]
nthist = 1
frhist = -1.d0
!< [usipes_ex_05]

! --- Number of monitoring points (probes) and their positions
!     (limited to ncaptm=200)

!< [usipes_ex_06]
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
!< [usipes_ex_06]

! Postprocessing of variables and properties
! (example for velocity: 1 to activate postprocessing output, 0 to deactivate)

!< [usipes_ex_07]
f_id = ivarfl(iu)
iflpst = 0
call field_set_key_int(f_id, keyvis, iflpst)
!< [usipes_ex_07]

! Probes for variables and properties
! (example for velocity, probes on u component only)

!< [usipes_ex_08]
f_id = ivarfl(iu)
ipp = field_post_id(f_id)
ihisvr(ipp,1) = -1
ihisvr(ipp+1,1) = 0
ihisvr(ipp+2,1) = 0
!< [usipes_ex_08]

! Force postprocessing of projection of some variables at boundary
! with no reconstruction.
! This is handled automatically if the second bit of a field's
! 'post_vis' key value is set to 1 (which amounts to adding 2
! to that key value).

!< [usipes_ex_09]
f_id = ivarfl(iu)
if (iand(iflpst, 2) .eq. 0) then
  call field_get_key_int(f_id, keyvis, iflpst)
  iflpst = ior(iflpst, 2)
  call field_set_key_int(f_id, keyvis, iflpst)
endif

f_id = ivarfl(ipr)
if (iand(iflpst, 2) .eq. 0) then
  call field_get_key_int(f_id, keyvis, iflpst)
  iflpst = ior(iflpst, 2)
  call field_set_key_int(f_id, keyvis, iflpst)
endif
!< [usipes_ex_09]

! Probes for Radiative Transfer (Luminance and radiative density flux vector)
! for all probes (ihisvr = -1)

!< [usipes_ex_10]
call field_get_id_try('luminance', f_id)
ipp = field_post_id(f_id)
ihisvr(ipp,1) = -1

call field_get_id_try('radiative_flux', f_id)
! X component
ipp = field_post_id(f_id)
ihisvr(ipp,1) = -1
! Y component
ihisvr(ipp+1,1) = -1
! Z component
ihisvr(ipp+2,1) = -1
!< [usipes_ex_10]

!--------
! Formats
!--------

return
end subroutine usipes
