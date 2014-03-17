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

subroutine cs_fuel_prop
!======================

!===============================================================================
! Purpose:
! --------

! Define state variables for fuel combustion.

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
use ppcpfu
use cpincl
use cs_fuel_incl
use ppincl

!===============================================================================

implicit none

! Local variables

integer       icla , nprini

character(len=80) :: f_name, f_label

!===============================================================================

! Initialization

nprini = nproce

! Continuous phase (gaseous mix)
call add_property_field('t_gas', 'T_Gas', itemp1)
call add_property_field('rho_gas', 'Rho_Gas', irom1)


! Gas mixture fractions
call add_property_field('ym_fo0',   'Ym_FO0',   iym1(1))
call add_property_field('ym_fov',   'Ym_FOV',   iym1(2))
call add_property_field('ym_co',    'Ym_CO',    iym1(3))
call add_property_field('ym_h2s',   'Ym_H2S',   iym1(4))
call add_property_field('ym_h2',    'Ym_H2',    iym1(5))
call add_property_field('ym_hcn',   'Ym_HCN',   iym1(6))
call add_property_field('ym_nh3',   'Ym_NH3',   iym1(7))
call add_property_field('ym_o2',    'Ym_O2',    iym1(8))
call add_property_field('ym_co2',   'Ym_CO2',   iym1(9))
call add_property_field('ym_h2o',   'Ym_H2O',   iym1(10))
call add_property_field('ym_so2',   'Ym_SO2',   iym1(11))
call add_property_field('ym_n2',    'Ym_N2',    iym1(12))

! Algebraic variables specific to gas - particles suspension
call add_property_field('xm',    'Xm',    immel)
call hide_property(immel)

if (ieqnox .eq. 1) then
  call add_property_field('exp1',      'EXP1',      ighcn1)
  call add_property_field('exp2',      'EXP1',      ighcn2)
  call add_property_field('exp3',      'EXP3',      ignoth)
endif

! Dispersed phase (particle classes)
do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 't_fuel_', icla
  write(f_label, '(a,i2.2)') 'T_Fuel_', icla
  call add_property_field(f_name, f_label, itemp2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'rho_fuel_', icla
  write(f_label, '(a,i2.2)') 'Rho_Fuel_', icla
  call add_property_field(f_name, f_label, irom2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'diameter_fuel_', icla
  write(f_label, '(a,i2.2)') 'Diam_Drop_', icla
  call add_property_field(f_name, f_label, idiam2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'h1_hlf_', icla
  write(f_label, '(a,i2.2)') 'H1-Hlf_', icla
  call add_property_field(f_name, f_label, ih1hlf(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'eva_fuel_', icla
  write(f_label, '(a,i2.2)') 'EVA_Fuel_', icla
  call add_property_field(f_name, f_label, igmeva(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'het_ts_fuel_', icla
  write(f_label, '(a,i2.2)') 'Het_TS_Fuel_', icla
  call add_property_field(f_name, f_label, igmhtf(icla))
enddo

! Balance: C,  O,  H

call add_property_field('balance_c', 'Balance_C', ibcarbone)
call add_property_field('balance_o', 'Balance_O', iboxygen)
call add_property_field('balance_h', 'Balance_H', ibhydrogen)

! Nb algebraic (or state) variables
!   specific to specific physic: nsalpp
!   total: nsalto

nsalpp = nproce - nprini
nsalto = nproce

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
