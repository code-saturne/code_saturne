!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine cfprop
!================

!===============================================================================
!  FONCTION  :
!  ---------

!     INIT DES POSITIONS DES VARIABLES D'ETAT POUR
!              POUR LE COMPRESSIBLE SANS CHOC
!         (DANS VECTEURS PROPCE)

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

!===============================================================================

implicit none

! Local variables

integer       nprini

!===============================================================================
! 1. POSITIONNEMENT DES PROPRIETES : PROPCE
!    Physique particuliere : Compressible sans choc
!===============================================================================

nprini = nproce

if (ippmod(icompf).ge.0) then

  if (icv.gt.0) then
    call add_property_field('specific_heat_const_vol', &
                            'Specific_Heat_Const_Vol', &
                            icv)
    call hide_property(icv)
    ihisvr(field_post_id(iprpfl(icv)),1) = 0
  endif

  if (iviscv.ne.0) then
    call add_property_field('volume_viscosity', &
                            'Volume_Viscosity', &
                            iviscv)
    call hide_property(iviscv)
    ihisvr(field_post_id(iprpfl(iviscv)),1) = 0
  endif

! Nb algebraic (or state) variables
!   specific to specific physic: nsalpp
!   total: nsalto

nsalpp = nproce - nprini
nsalto = nproce

endif


return
end subroutine
