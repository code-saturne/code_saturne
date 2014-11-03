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

subroutine daini1
!================


!===============================================================================
! Purpose:
! --------

! Initialize global settings for darcy module.

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
use ihmpre
use darcy_module

!===============================================================================

implicit none

! Local variables

!===============================================================================

!===============================================================================
! 1. Transported variables
!===============================================================================

iwgrec(ipr) = 0
imrgra = -11
imvisf = 1

ro0 = 1.
p0 = 0

darcy_anisotropic_permeability = 0
darcy_anisotropic_diffusion = 0
darcy_unsteady = 0
darcy_convergence_criterion = 0
darcy_gravity = 0
darcy_gravity_x = 0.
darcy_gravity_y = 0.
darcy_gravity_z = 0.

if (iihmpr.eq.1) then
  call uidai1(ippmod(idarcy),                     &
              darcy_anisotropic_permeability,     &
              darcy_anisotropic_diffusion,        &
              darcy_unsteady,                     &
              darcy_convergence_criterion,        &
              darcy_gravity,                      &
              darcy_gravity_x,                    &
              darcy_gravity_y,                    &
              darcy_gravity_z)
endif
!===============================================================================
! 2. Define user settings
!===============================================================================

call user_darcy_ini1
!===================

return
end subroutine
