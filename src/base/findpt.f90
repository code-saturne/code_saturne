!-------------------------------------------------------------------------------

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
! Function:
! ---------

!> \file findpt.f90
!>
!> \brief This subroutine looks for the nearest element to the position
!>  (xx, yy, zz) among the element of xyzcen array.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ncelet        number of extended (real + ghost) elements
!> \param[in]     ncel          number of element
!> \param[in]     xyzcen        interleaved array containing vector positions
!> \param[in]     xx            X component of the point position to be compared
!>                              with
!> \param[in]     yy            Y component of the point position to be compared
!>                              with
!> \param[in]     zz            Z component of the point position to be compared
!>                              with
!> \param[in]     node          element of xyzcen which is the nearest to
!>                              (xx,yy,zz)
!> \param[in]     ndrang        index of the Rank containing the nearest element
!_______________________________________________________________________________


subroutine findpt &
 ( ncelet , ncel   , xyzcen ,                                     &
   xx     , yy     , zz     , node   , ndrang)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use parall

!===============================================================================

implicit none

integer          ncelet, ncel, node, ndrang
double precision xyzcen(3,ncelet)
double precision xx, yy, zz

integer          ii
double precision xx1, yy1, zz1, dis2, dis2mn

!===============================================================================
! 1. Initialization
!===============================================================================

node = int((ncel+1)/2)

xx1 = xyzcen(1,node)
yy1 = xyzcen(2,node)
zz1 = xyzcen(3,node)
dis2mn = (xx-xx1)**2+(yy-yy1)**2+(zz-zz1)**2

do ii = 1, ncel
  xx1 = xyzcen(1,ii)
  yy1 = xyzcen(2,ii)
  zz1 = xyzcen(3,ii)
  dis2 = (xx-xx1)**2+(yy-yy1)**2+(zz-zz1)**2
  if (dis2.lt.dis2mn) then
    node = ii
    dis2mn = dis2
  endif
enddo

if (irangp.ge.0) then
  call parfpt (node, ndrang, dis2mn)
else
  ndrang = -1
endif

return
end subroutine findpt
