!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine defptr &
!================

 ( numper, crit, fract, plane, iwarnj, tx, ty, tz )

!===============================================================================
! Purpose:
! -------

! Definition of a periodicity of translation.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! numper           ! i  ! --> ! periodicity number                             !
! crit             ! a  ! <-- ! selection criteria for the border faces to     !
!                  !    !     ! transform                                      !
! fract            ! r  ! <-- ! fraction parameter                             !
! plane            ! r  ! <-- ! face coplanarity parameter                     !
! iwarnj           ! i  ! <-- ! verbosity level                                !
! tx               ! r  ! <-- ! X coordinate of the translation vector         !
! ty               ! r  ! <-- ! Y coordinate of the translation vector         !
! tz               ! r  ! <-- ! Z coordinate of the translation vector         !
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

character*(*)    crit
integer          numper, iwarnj
double precision fract, plane
double precision tx, ty, tz

! Local variables

integer       lcrit

!===============================================================================

lcrit = len(crit)

call defpt1(numper, crit, fract, plane, iwarnj, tx, ty, tz, lcrit)
!==========

return

end subroutine

!===============================================================================

subroutine defpro &
!================

 ( numper, crit, fract, plane, iwarnj, ax, ay, az, theta, ix, iy, iz )

!===============================================================================
! Purpose:
! -------

! Definition of a periodicity of rotation.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! numper           ! i  ! --> ! periodicity number                             !
! crit             ! a  ! <-- ! selection criteria for the border faces to     !
!                  !    !     ! transform                                      !
! fract            ! r  ! <-- ! fraction parameter                             !
! plane            ! r  ! <-- ! face coplanarity parameter                     !
! iwarnj           ! i  ! <-- ! verbosity level                                !
! ax               ! r  ! <-- ! X coordinate of the rotation axis              !
! ay               ! r  ! <-- ! Y coordinate of the rotation axis              !
! az               ! r  ! <-- ! Z coordinate of the rotation axis              !
! theta            ! r  ! <-- ! angle of the rotation (radian)                 !
! ix               ! r  ! <-- ! X coordinate of the invariant point            !
! iy               ! r  ! <-- ! Y coordinate of the invariant point            !
! iz               ! r  ! <-- ! Z coordinate of the invariant point            !
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

character*(*)    crit
integer          numper, iwarnj
double precision fract, plane
double precision ax, ay, az, theta, ix, iy, iz

! Local variables

integer       lcrit

!===============================================================================

lcrit = len(crit)

call defpr1(numper, crit, fract, plane, iwarnj, &
!==========
            ax, ay, az, theta, ix, iy, iz, lcrit)

return

end subroutine

!===============================================================================

subroutine defpge &
!================

 ( numper, crit, fract, plane, iwarnj,   &
   r11, r12, r13, tx,                    &
   r21, r22, r23, ty,                    &
   r31, r32, r33, tz )

!===============================================================================
! Purpose:
! -------

! Definition of a general periodicity.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! numper           ! i  ! --> ! periodicity number                             !
! crit             ! a  ! <-- ! selection criteria for the border faces to     !
!                  !    !     ! transform                                      !
! fract            ! r  ! <-- ! fraction parameter                             !
! plane            ! r  ! <-- ! face coplanarity parameter                     !
! iwarnj           ! i  ! <-- ! verbosity level                                !
! r11              ! r  ! <-- ! coef. (1,1) of the homogeneous matrix          !
! r12              ! r  ! <-- ! coef. (1,2) of the homogeneous matrix          !
! r13              ! r  ! <-- ! coef. (1,3) of the homogeneous matrix          !
! tx               ! r  ! <-- ! coef. (1,4) of the homogeneous matrix          !
! r21              ! r  ! <-- ! coef. (2,1) of the homogeneous matrix          !
! r22              ! r  ! <-- ! coef. (2,2) of the homogeneous matrix          !
! r23              ! r  ! <-- ! coef. (2,3) of the homogeneous matrix          !
! ty               ! r  ! <-- ! coef. (2,4) of the homogeneous matrix          !
! r31              ! r  ! <-- ! coef. (3,1) of the homogeneous matrix          !
! r32              ! r  ! <-- ! coef. (3,2) of the homogeneous matrix          !
! r33              ! r  ! <-- ! coef. (3,3) of the homogeneous matrix          !
! tz               ! r  ! <-- ! coef. (3,4) of the homogeneous matrix          !
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

character*(*)    crit
integer          numper, iwarnj
double precision fract, plane
double precision r11, r12, r13, r21, r22, r23, r31, r32, r33, tx, ty, tz

! Local variables

integer       lcrit

!===============================================================================

lcrit = len(crit)

call defpg1(numper, crit, fract, plane, iwarnj, &
!==========
            r11, r12, r13, tx, r21, r22, r23, ty, r31, r32, r33, tz, lcrit)

return

end subroutine
