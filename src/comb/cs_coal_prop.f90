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

subroutine cs_coal_prop
!======================

!===============================================================================
! Purpose:
! --------

! Define state variables for pulverized coal combustion.

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
use ppcpfu
use coincl
use cpincl
use ppincl
use ihmpre
use cs_coal_incl
use field

!===============================================================================

implicit none

! Local variables

integer          icla
integer          f_id, itycat, ityloc, idim1, nprini
integer          keyccl
integer          iopchr

logical          ilved, iprev, inoprv

character(len=80) :: f_name, f_label

!===============================================================================

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 1 ! variables defined on cells
idim1  = 1
ilved  = .false.   ! not interleaved by default
iprev  = .true.    ! variables have previous value
inoprv = .false.   ! variables have no previous value
iopchr = 1         ! Postprocessing level for variables

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! ---> Definition des pointeurs relatifs aux variables d'etat

nprini = nproce

! Continuous phase (gaseous mix)
call add_property_field('t_gas', 'T_Gas', itemp1)
call add_property_field('rho_gas', 'Rho_Gas', irom1)

! Gas mixture fractions
call add_property_field('ym_chx1m', 'Ym_CHx1m', iym1(1))
call add_property_field('ym_chx2m', 'Ym_CHx2m', iym1(2))
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

! Algebraic variables specific to continuous phase
if (ieqnox .eq. 1) then
  call add_property_field('exp1',      'EXP1',      ighcn1)
  call add_property_field('exp2',      'EXP1',      ighcn2)
  call add_property_field('exp3',      'EXP3',      ignoth)
  call add_property_field('exp4',      'EXP4',      ignh31)
  call add_property_field('exp5',      'EXP5',      ignh32)
  call add_property_field('f_hcn_dev', 'F_HCN_DEV', ifhcnd)
  call add_property_field('f_hcn_het', 'F_HCN_HET', ifhcnc)
  call add_property_field('f_nh3_dev', 'F_NH3_DEV', ifnh3d)
  call add_property_field('f_nh3_het', 'F_NH3_HET', ifnh3c)
  call add_property_field('f_no_hcn',  'F_NO_HCN',  ifnohc)
  call add_property_field('f_no_nh3',  'F_NO_NH3',  ifnonh)
  call add_property_field('f_no_het',  'F_NO_HET',  ifnoch)
  call add_property_field('f_no_the',  'F_NO_THE',  ifnoth)
  call add_property_field('c_no_hcn',  'C_NO_HCN',  icnohc)
  call add_property_field('c_no_nh3',  'C_NO_NH3',  icnonh)
  call add_property_field('f_hcn_rb',  'F_HCN_RB',  ifhcnr)
  call add_property_field('c_no_rb',   'C_NO_RB',   icnorb)
  call add_property_field('exp_rb',    'Exp_RB',    igrb)
endif

! Dispersed phase (particle classes)
do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 't_coal', icla
  write(f_label, '(a,i2.2)') 'T_Coal', icla
  call add_property_field(f_name, f_label, itemp2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'w_solid_coal', icla
  write(f_label, '(a,i2.2)') 'w_solid_coal', icla
  call add_property_field(f_name, f_label, ix2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'rho_coal', icla
  write(f_label, '(a,i2.2)') 'Rho_Coal', icla
  call add_property_field(f_name, f_label, irom2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'diameter_coal', icla
  write(f_label, '(a,i2.2)') 'Diam_Coal', icla
  call add_property_field(f_name, f_label, idiam2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'dissapear_rate_coal', icla
  write(f_label, '(a,i2.2)') 'D_Rate_Coal', icla
  call add_property_field(f_name, f_label, igmdch(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'm_transfer_v1_coal', icla
  write(f_label, '(a,i2.2)') 'D_V1_Coal', icla
  call add_property_field(f_name, f_label, igmdv1(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'm_transfer_v2_coal', icla
  write(f_label, '(a,i2.2)') 'D_V2_Coal', icla
  call add_property_field(f_name, f_label, igmdv2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'het_ts_o2_coal', icla
  write(f_label, '(a,i2.2)') 'Het_TS_O2_Coal', icla
  call add_property_field(f_name, f_label, igmhet(icla))
enddo

if (i_coal_drift.eq.1) then
  do icla = 1, nclacp
    write(f_name,'(a,i2.2)') 'age_coal', icla
    write(f_name,'(a,i2.2)') 'Age_Coal', icla
    call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! For post-processing
    call field_set_key_int(f_id, keyvis, iopchr)
    ! For log in the listing
    call field_set_key_int(f_id, keylog, 1)
  enddo
endif

if (ihtco2 .eq. 1) then
  do icla = 1, nclacp
    write(f_name,  '(a,i2.2)') 'het_ts_co2_coal', icla
    write(f_label, '(a,i2.2)') 'Het_TS_CO2_Coal', icla
    call add_property_field(f_name, f_label, ighco2(icla))
  enddo
endif

if (ihth2o .eq. 1) then
  do icla = 1, nclacp
    write(f_name,  '(a,i2.2)') 'het_ts_h2o_coal',  icla
    write(f_label, '(a,i2.2)') 'Het_TS_H2O_Coal', icla
    call add_property_field(f_name, f_label, ighh2o(icla))
  enddo
endif

if (ippmod(iccoal) .ge. 1) then
  do icla = 1, nclacp
    write(f_name,  '(a,i2.2)') 'dry_ts_coal',  icla
    write(f_label, '(a,i2.2)') 'Dry_TS_Coal', icla
    call add_property_field(f_name, f_label, igmsec(icla))
  enddo
endif

if (i_coal_drift.eq.1) then
  icla = -1
  f_name = 'Age_Gas'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
  ! For log in the listing
  call field_set_key_int(f_id, keylog, 1)
endif

! Balance: C,  O,  H

call add_property_field('balance_c', 'Balance_C', ibcarbone)
call add_property_field('balance_o', 'Balance_O', iboxygen)
call add_property_field('balance_h', 'Balance_H', ibhydrogen)

! Nb algebraic (or state) variables
!   specific to specific physic: nsalpp
!   total: nsalto

nsalpp = nproce - nprini
nsalto = nproce

!   - Interface Code_Saturne
!     ======================
!     Construction de l'indirection entre la numerotation du noyau et XML
if (iihmpr.eq.1) then
  call uicppr &
  !==========
 ( nclacp, nsalpp, ippmod, iccoal, ipppro,            &
   ipproc, ieqnox, ihtco2, ihth2o, itemp1,            &
   irom1, iym1, ighcn1, ighcn2, ignoth,               &
   ignh31, ignh32, ifhcnd, ifhcnc, ifnh3d, ifnh3c,    &
   ifnohc, ifnonh, ifnoch, ifnoth, icnohc, icnonh,    &
   ifhcnr, icnorb, igrb,   immel,                     &
   itemp2, ix2, irom2, idiam2, igmdch, igmdv1,        &
   igmdv2, igmhet, ighco2, ighh2o, igmsec,            &
   ibcarbone, iboxygen, ibhydrogen)
endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
