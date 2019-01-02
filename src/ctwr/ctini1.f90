!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer jj, isc

type(var_cal_opt) :: vcopt

!===============================================================================

!===============================================================================
! 1. Transported variables
!===============================================================================

irovar = 1 ! Variable density
ivivar = 0 ! Constant molecular viscosity

do isc = 1, nscapp

  jj = iscapp(isc)

  if (iscavr(jj).le.0) then  ! iscavr = 0 for scalars which are not mean square errors of other scalars
    visls0(jj) = viscl0
  endif

  call field_get_key_struct_var_cal_opt(ivarfl(isca(jj)), vcopt)

  ! Upwind for scalars in the packing zones
  if (jj.eq.iyml .or. jj.eq.ihml) then
    vcopt%blencv = 0.d0
    vcopt%idiff  = 0
    vcopt%idifft = 0
  else
    vcopt%blencv = 1.d0
  endif

  ! Set beta limiter to maintain y_p in the limits
  if (jj.eq.iy_p_l) then
    vcopt%isstpc = 2
  endif

  call field_set_key_struct_var_cal_opt(ivarfl(isca(jj)), vcopt)

enddo

!===============================================================================
! 2. Define user settings
!===============================================================================

call cs_user_cooling_towers()

return
end subroutine
