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

subroutine cplpro
!================

!===============================================================================
!  FONCTION  :
!  ---------


!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!      INIT DES POSITIONS DES VARIABLES D'ETAT

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
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Local variables

!===============================================================================

! Continuous phase (gas mixture)
call add_property_field_1d('t_gas', 'T_Gas', itemp1)

! State variables

! ---> Variables algebriques propres a la phase continue

call add_property_field_1d('ym_chx1m', 'Ym_CHx1m', iym1(1))
call add_property_field_1d('ym_chx2m', 'Ym_CHx2m', iym1(2))
call add_property_field_1d('ym_co',    'Ym_CO',    iym1(3))
call add_property_field_1d('ym_o2',    'Ym_O2',    iym1(4))
call add_property_field_1d('ym_co2',   'Ym_CO2',   iym1(5))
call add_property_field_1d('ym_h2o',   'Ym_H2O',   iym1(6))
call add_property_field_1d('ym_n2',    'Ym_N2',    iym1(7))
call add_property_field_1d('xm',       'Xm',       immel)
call hide_property(immel)

return
end subroutine cplpro
