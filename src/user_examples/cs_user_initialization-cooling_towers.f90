!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Purpose:
! -------

!> \file cs_user_initialization-cooling_towers.f90
!>
!> \brief Cooling towers example
!>
!> See \subpage cs_user_initialization for examples.
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________


subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel, iutile
integer          ilelt, nlelt

double precision d2s3

double precision, dimension(:,:), pointer :: vel

integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: cvar_temp, cpro_humid
!< [loc_var_dec]

!===============================================================================

!---------------
! Initialization
!---------------

!< [init]
allocate(lstelt(ncel)) ! temporary array for cells selection

d2s3 = 2.d0/3.d0

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================
! Variables initialization:
!
!   ONLY done if there is no restart computation
!===============================================================================

if (isuite.eq.0) then

!   --- Initialize tmperature of humid air at 11 deg Celsius
!       and of humidity at 0.0063

  call field_get_val_s(ivarfl(isca(iscalt)), cvar_temp)
  call field_get_val_s(ihumid, cpro_humid)!FIXME

  do iel = 1, ncel

    cvar_temp(iel) = 11.d0
    cpro_humid(iel) = 0.0063d0

  enddo

!   --- Initialize tmperature of humid air at 20 deg Celsius
!       and of humidity at 0.012
!       and of velocity at 0.5 m/s
!       for cells of color 6

  call getcel('6', nlelt, lstelt)

  do ilelt = 1, nlelt

    iel = lstelt(ilelt)

    vel(1,iel) = -0.5d0

    cvar_temp(iel) = 20.d0
    cpro_humid(iel) = 0.012d0

  enddo

endif
!< [init]

!--------
! Formats
!--------

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine cs_user_f_initialization
