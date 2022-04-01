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

!> \file initi2.f90
!> \brief End of commons initialization.
!>
!------------------------------------------------------------------------------

subroutine initi2

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use optcal
use entsor
use cstphy
use mesh

!===============================================================================

implicit none

! Arguments

! Local variables

!===============================================================================

! If almax < 0, recompute it

write(nfecra, 1000)

if (almax.le.0.d0) then
  almax = voltot**(1.d0/3.d0)
  write(nfecra, 1100) almax
  write(nfecra, 1102)
  if (itytur.eq.2.or.itytur.eq.3                  &
       .or. itytur.eq.5 .or. iturb.eq.60          &
       .or. iturb.eq.70) then
    write(nfecra, 1101)
  endif
endif

 1000 format('')
 1100 format(                                                           &
'       ALMAX  = ', e14.5, ' (Characteristic length)')
 1101 format(                                                           &
'       ALMAX is the length used to initialize the turbulence.')
 1102 format(                                                           &
'       ALMAX is the cubic root of the domain volume.',/)

!===============================================================================
! End
!===============================================================================

return
end subroutine
