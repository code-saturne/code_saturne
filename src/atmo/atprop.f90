!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine atprop
!================

!===============================================================================
!  FONCTION  :
!  ---------

!     INIT DES POSITIONS DES VARIABLES D'ETAT
!            POUR LE MODULE ATMOSPHERIQUE

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

use atincl

!===============================================================================

implicit none

! Local variables

integer       nprini

!===============================================================================

nprini = nproce

!===============================================================================
! 1. POSITIONNEMENT DES PROPRIETES
!    Atmospheric modules:  dry and humid atmosphere
!===============================================================================

! Temperature (IPPMOD(IATMOS) = 1 or 2)
!--------------------------------------

if (ippmod(iatmos).ge.1) then
  call add_property_field('real_temperature', 'RealTemp', itempc)
endif

! Liquid water content (IPPMOD(IATMOS) = 2)
!------------------------------------------

if (ippmod(iatmos).eq.2) then
  call add_property_field('liquid_water', 'LiqWater', iliqwt)
endif

! Nb algebraic (or state) variables
!   specific to specific physic: nsalpp
!   total: nsalto

nsalpp = nproce - nprini
nsalto = nproce

return
end subroutine atprop
