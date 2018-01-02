!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
!> \file atprop.f90
!> \brief Add if needed the variables fields for temperature and liquid water
!>
!> \brief Add if needed the variable field for temperature and liquid water \n
!>    Nota : ippmod(iatmos) = 1 --> Dry atmosphere, = 2 --> Humid atmosphere
subroutine atprop


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

use atincl

!===============================================================================

implicit none

! Local variables

!===============================================================================

!===============================================================================
! 1. POSITIONNEMENT DES PROPRIETES
!    Atmospheric modules:  dry and humid atmosphere
!===============================================================================

! Temperature (IPPMOD(IATMOS) = 1 or 2)
!--------------------------------------

if (ippmod(iatmos).ge.1) then
  call add_property_field_1d('real_temperature', 'RealTemp', itempc)
endif

! Liquid water content (IPPMOD(IATMOS) = 2)
!------------------------------------------

if (ippmod(iatmos).eq.2) then
  call add_property_field_1d('liquid_water', 'LiqWater', iliqwt)
endif

return
end subroutine atprop
