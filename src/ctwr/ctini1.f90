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

subroutine ctini1
!================


!===============================================================================
! Purpose:
! --------

! Initialize global settings for cooling towers module.

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
use ppthch
use ppincl
use ctincl

!===============================================================================

implicit none

! Local variables

integer ii, jj, isc, ipp
double precision cpa,cpe,cpv,hv0,rhoe,visc,conduc

!===============================================================================

!===============================================================================
! 1. Transported variables
!===============================================================================

itherm = 1
iscacp(itemp4) = 1
iscacp(ihumid) = 0

iscalt = itemp4

irovar = 1
ivivar = 0

! --> Physical or numerical properties specific to scalars

do isc = 1, nscapp

  jj = iscapp(isc)

  if (iscavr(jj).le.0) then
    visls0(jj) = viscl0
  endif

  blencv(isca(jj)) = 1.d0

enddo

ipp = ipprtp(isca(itemp4))
nomvar(ipp)  = 'Temperature'
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

ipp = ipprtp(isca(ihumid))
nomvar(ipp)  = 'Humidity'
ichrvr(ipp)  = 1
ilisvr(ipp)  = 1
ihisvr(ipp,1)= -1

! Postprocessing output

ichrze = 1

! Air properties

cpa    = 1006.0d0
cpv    = 1831.0d0
cpe    = 4179.0d0
hv0    = 2501600.0d0
rhoe   = 997.85615d0
visc   = 1.765d-5
conduc = 0.02493d0

call ctprof(cpa, cpv, cpe, hv0, rhoe, visc, conduc, gx, gy, gz)
!==========

!===============================================================================
! 2. Define user settings
!===============================================================================

call uscti1
!==========

return
end subroutine
