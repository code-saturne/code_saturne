!-------------------------------------------------------------------------------

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

!===============================================================================
! Function:
! ---------

!> \file cs_selector_f2c.f90

!> \brief Build the list of interior faces matching a criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    facnb         number of selected faces
!> \param[out]    faces         selected faces
!_______________________________________________________________________________

subroutine getfac &
 ( fstr , facnb, faces)

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

!===============================================================================
! Function:
! ---------

!> \brief Build the list of boundary faces matching a criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    facnb         number of selected faces
!> \param[out]    faces         selected faces
!_______________________________________________________________________________

subroutine getfbr &
 ( fstr , facnb, faces)

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

!===============================================================================
! Function:
! ---------

!> \brief Build the list of cells matching a criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    cellnb        number of selected cells
!> \param[out]    cells         selected cells
!_______________________________________________________________________________


subroutine getcel &
 ( fstr , cellnb, cells)

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

!===============================================================================
! Function:
! ---------

!> \brief Build the lists of faces at the boundary of cells matching a
!> criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    ifaces        selected interior faces
!> \param[out]    bfaces        selected boundary faces
!> \param[out]    ifacnb        number fo selected interior faces
!> \param[out]    bfacnb        number fo selected boundary faces
!_______________________________________________________________________________

subroutine getceb &
 ( fstr , ifacnb, bfacnb, ifaces, bfaces )

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

!===============================================================================
! Function:
! ---------

!> \brief Build the list of mesh element families matching a criteria string.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     fstr          criteria string
!> \param[out]    families      selected families
!> \param[out]    famnb         number of selected families
!_______________________________________________________________________________

subroutine getfam &
 ( fstr , famnb, families)

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
