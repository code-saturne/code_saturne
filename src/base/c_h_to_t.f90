!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
! Function :
! --------

!> \file c_h_to_t.f90
!>
!> \brief Convert enthalpy to temperature at cells
!>
!> This handles both user and model enthalpy conversions, so can be used
!> safely whenever conversion is needed.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     h             enthalphy
!> \param[out]    t             temperature
!_______________________________________________________________________________


subroutine c_h_to_t            &
 ( h, t )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use radiat
use mesh
use field

!===============================================================================

implicit none

! Arguments

double precision h(ncelet), t(ncelet)

! Local variables

integer          iel, mode
double precision hl

double precision, dimension(:), pointer :: cpro_t

!===============================================================================

mode = 1

! Non-specific physics

if (ippmod(iphpar).le.1) then

  do iel = 1,ncel
    hl = h(iel)
    call usthht(mode, hl, t(iel))
  enddo

  return

endif

! Gas combustion: premix or diffusion flame

if (ippmod(icoebu).ge.0 .or. ippmod(icod3p).ge.0) then

  call field_get_val_s(itemp, cpro_t)

  do iel = 1, ncel
    t(iel) = cpro_t(iel)
  enddo

! Pulverized coal or fuel combustion

else if (ippmod(iccoal).ge.0 .or. ippmod(icfuel).ge.0) then

  call field_get_val_s(itemp1, cpro_t)

  do iel = 1, ncel
    t(iel) = cpro_t(iel)
  enddo

! Electric arcs

else if (ippmod(ieljou).ge.1 .or.                              &
         ippmod(ielarc).ge.1) then

  call field_get_val_s(itemp, cpro_t)

  do iel = 1, ncel
    t(iel) = cpro_t(iel)
  enddo

! Other cases ?

else

  do iel = 1, ncel
    hl = (iel)
    call usthht(mode, h(iel), t(iel))
  enddo

endif

return
end subroutine
