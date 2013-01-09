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

subroutine getfac &
!=================

 ( fstr , facnb, faces)


!===============================================================================
! Purpose:
! -------

! Build the list of interior faces matching a criteria string.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! fstr             ! a  ! <-- ! criteria string                                !
! facnb            ! i  ! --> ! number of selected faces                       !
! faces            ! ia ! --> ! selected faces                                 !
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

character*(*) fstr
integer       faces(*), facnb

! Local variables

integer          lenstr

!===============================================================================

lenstr=len(fstr)
call csgfac(fstr, lenstr, facnb, faces)

return

end subroutine

!===============================================================================

subroutine getfbr &
!=================

 ( fstr , facnb, faces)

!===============================================================================
! Purpose:
! -------

! Build the list of boundary faces matching a criteria string.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! fstr             ! a  ! <-- ! criteria string                                !
! facnb            ! i  ! --> ! number of selected faces                       !
! faces            ! ia ! --> ! selected faces                                 !
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

character*(*)    fstr
integer      faces(*), facnb

! Local variables

integer          lenstr

!===============================================================================

lenstr=len(fstr)
call csgfbr(fstr, lenstr, facnb, faces)

return

end subroutine

!===============================================================================

subroutine getcel &
!=================

 ( fstr , cellnb, cells)

!===============================================================================
! Purpose:
! -------

! Build the list of cells matching a criteria string.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! fstr             ! a  ! <-- ! criteria string                                !
! cellnb           ! i  ! --> ! number fo selected cells                       !
! cells            ! ia ! --> ! selected cells                                 !
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

character*(*) fstr
integer       cells(*), cellnb

! Local variables

integer          lenstr

!===============================================================================

lenstr=len(fstr)
call csgcel(fstr, lenstr, cellnb, cells)

return

end subroutine

!===============================================================================

subroutine getceb &
!=================

 ( fstr , ifacnb, bfacnb, ifaces, bfaces )

!===============================================================================
! Purpose:
! -------

! Build the lists of faces at the boundary of cells matching a criteria string.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! fstr             ! a  ! <-- ! criteria string                                !
! ifaces           ! i  ! --> ! selected interior faces                        !
! bfaces           ! i  ! --> ! selected boundary faces                        !
! ifacnb           ! ia ! --> ! number fo selected interior faces              !
! bfacnb           ! ia ! --> ! number fo selected boundary faces              !
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

character*(*) fstr
integer       ifaces(*), bfaces(*), ifacnb, bfacnb

! Local variables

integer       lenstr

!===============================================================================

lenstr=len(fstr)
call csgceb(fstr, lenstr, ifacnb, bfacnb, ifaces, bfaces)

return

end subroutine

!===============================================================================

subroutine getfam &
!=================

 ( fstr , famnb, families)

!===============================================================================
! Purpose:
! -------

! Build the list of mesh element families matching a criteria string.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! fstr             ! a  ! <-- ! criteria string                                !
! families         ! i  ! <-- ! selected families                              !
! famnb            ! i  ! <-- ! number of selected families                    !
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

character*(*) fstr
integer       families(*), famnb

! Local variables

integer          lenstr

!===============================================================================

lenstr=len(fstr)
call csgfac(fstr, lenstr, famnb, families)

return

end subroutine
