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


!===============================================================================
!> Purpose:
!> --------
!> \file cs_coal_prop.f90
!> \brief Define state variables for pulverized coal combustion.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!>
!______________________________________________________________________________!

subroutine cs_coal_prop

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
use post

!===============================================================================

implicit none

! Local variables

integer          icla
integer          f_id, itycat, ityloc, idim1, idim3
integer          keyccl
integer          iopchr

logical          iprev, inoprv

character(len=80) :: f_name, f_label

!===============================================================================

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 1 ! variables defined on cells
idim1  = 1
idim3  = 3
iprev  = .true.    ! variables have previous value
inoprv = .false.   ! variables have no previous value
iopchr = POST_ON_LOCATION + POST_MONITOR ! postprocessing level for variables

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! ---> Definition of pointers related to state variables

! Continuous phase (gaseous mix)
call add_property_field_1d('t_gas', 'T_Gas', itemp1)
call add_property_field_1d('rho_gas', 'Rho_Gas', irom1)

! Gas mixture fractions
call add_property_field_1d('ym_chx1m', 'Ym_CHx1m', iym1(1))
call add_property_field_1d('ym_chx2m', 'Ym_CHx2m', iym1(2))
call add_property_field_1d('ym_co',    'Ym_CO',    iym1(3))
call add_property_field_1d('ym_h2s',   'Ym_H2S',   iym1(4))
call add_property_field_1d('ym_h2',    'Ym_H2',    iym1(5))
call add_property_field_1d('ym_hcn',   'Ym_HCN',   iym1(6))
call add_property_field_1d('ym_nh3',   'Ym_NH3',   iym1(7))
call add_property_field_1d('ym_o2',    'Ym_O2',    iym1(8))
call add_property_field_1d('ym_co2',   'Ym_CO2',   iym1(9))
call add_property_field_1d('ym_h2o',   'Ym_H2O',   iym1(10))
call add_property_field_1d('ym_so2',   'Ym_SO2',   iym1(11))
call add_property_field_1d('ym_n2',    'Ym_N2',    iym1(12))

! Algebraic variables specific to gas - particles suspension
call add_property_field_1d('xm',    'Xm',    immel)
call hide_property(immel)

! Algebraic variables specific to continuous phase
if (ieqnox .eq. 1) then
  call add_property_field_1d('exp1',      'EXP1',      ighcn1)
  call add_property_field_1d('exp2',      'EXP1',      ighcn2)
  call add_property_field_1d('exp3',      'EXP3',      ignoth)
  call add_property_field_1d('exp4',      'EXP4',      ignh31)
  call add_property_field_1d('exp5',      'EXP5',      ignh32)
  call add_property_field_1d('f_hcn_dev', 'F_HCN_DEV', ifhcnd)
  call add_property_field_1d('f_hcn_het', 'F_HCN_HET', ifhcnc)
  call add_property_field_1d('f_nh3_dev', 'F_NH3_DEV', ifnh3d)
  call add_property_field_1d('f_nh3_het', 'F_NH3_HET', ifnh3c)
  call add_property_field_1d('f_no_hcn',  'F_NO_HCN',  ifnohc)
  call add_property_field_1d('f_no_nh3',  'F_NO_NH3',  ifnonh)
  call add_property_field_1d('f_no_het',  'F_NO_HET',  ifnoch)
  call add_property_field_1d('f_no_the',  'F_NO_THE',  ifnoth)
  call add_property_field_1d('c_no_hcn',  'C_NO_HCN',  icnohc)
  call add_property_field_1d('c_no_nh3',  'C_NO_NH3',  icnonh)
  call add_property_field_1d('f_hcn_rb',  'F_HCN_RB',  ifhcnr)
  call add_property_field_1d('c_no_rb',   'C_NO_RB',   icnorb)
  call add_property_field_1d('exp_rb',    'Exp_RB',    igrb)
endif

! Dispersed phase (particle classes)

! NB: 'c' stands for continuous <> 'p' stands for particles

! Temperature of particle class icla
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 't_p_', icla
  write(f_label, '(a,i2.2)') 'Tp_', icla
  call add_property_field_1d(f_name, f_label, itemp2(icla))
enddo

! Temperature of particle class icla
! NB: mixture fraction (fr) (unreactive) <> mass fraction (x) (reactive)
do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'x_p_', icla
  write(f_label, '(a,i2.2)') 'Xp_', icla
  call add_property_field_1d(f_name, f_label, ix2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'rho_p_', icla
  write(f_label, '(a,i2.2)') 'Rhop_', icla
  call add_property_field_1d(f_name, f_label, irom2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'diam_p_', icla
  write(f_label, '(a,i2.2)') 'Diamp_', icla
  call add_property_field_1d(f_name, f_label, idiam2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'dissapear_rate_p_', icla
  write(f_label, '(a,i2.2)') 'D_Rate_Coal', icla
  call add_property_field_1d(f_name, f_label, igmdch(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'm_transfer_v1_p_', icla
  write(f_label, '(a,i2.2)') 'D_V1_Coal', icla
  call add_property_field_1d(f_name, f_label, igmdv1(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'm_transfer_v2_p_', icla
  write(f_label, '(a,i2.2)') 'D_V2_Coal', icla
  call add_property_field_1d(f_name, f_label, igmdv2(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'het_ts_o2_p_', icla
  write(f_label, '(a,i2.2)') 'Het_TS_O2_Coal', icla
  call add_property_field_1d(f_name, f_label, igmhet(icla))
enddo

do icla = 1, nclacp
  write(f_name,  '(a,i2.2)') 'imp_m_transfer_to_g_p_', icla
  write(f_label, '(a,i2.2)') 'Implicit_Mass_transfer', icla
  call add_property_field_1d(f_name, f_label, igmtr(icla))
enddo


if (i_comb_drift.ge.1) then
  do icla = 1, nclacp

    ! Age of the particle class
    write(f_name,'(a,i2.2)') 'age_p_', icla
    write(f_label,'(a,i2.2)') 'Agep_', icla
    call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! For post-processing
    call field_set_key_int(f_id, keyvis, iopchr)

    ! Limit velocity
    write(f_name,'(a,i2.2)')'vg_lim_p_' ,icla
    call field_create(f_name, itycat, ityloc, idim3, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! For post-processing
    call field_set_key_int(f_id, keyvis, iopchr)
    ! For log in the log
    call field_set_key_int(f_id, keylog, 1)

    write(f_name,'(a,i2.2)')'vg_p_' ,icla
    call field_create(f_name, itycat, ityloc, idim3, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! For post-processing
    call field_set_key_int(f_id, keyvis, iopchr)
    ! For log in the log
    call field_set_key_int(f_id, keylog, 1)

    ! Additional drift velocity for the particle class
    write(f_name,'(a,i2.2)')'vd_p_' ,icla
    call field_create(f_name, itycat, ityloc, idim3, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! For post-processing
    call field_set_key_int(f_id, keyvis, iopchr)
    ! For log in the log
    call field_set_key_int(f_id, keylog, 1)

  enddo
endif

if (ihtco2 .eq. 1) then
  do icla = 1, nclacp
    write(f_name,  '(a,i2.2)') 'het_ts_co2_p', icla
    write(f_label, '(a,i2.2)') 'Het_TS_CO2_p', icla
    call add_property_field_1d(f_name, f_label, ighco2(icla))
  enddo
endif

if (ihth2o .eq. 1) then
  do icla = 1, nclacp
    write(f_name,  '(a,i2.2)') 'het_ts_h2o_p',  icla
    write(f_label, '(a,i2.2)') 'Het_TS_H2O_p', icla
    call add_property_field_1d(f_name, f_label, ighh2o(icla))
  enddo
endif

if (ippmod(iccoal) .ge. 1) then
  do icla = 1, nclacp
    write(f_name,  '(a,i2.2)') 'dry_ts_p',  icla!FIXME is it a Source term?
    write(f_label, '(a,i2.2)') 'Dry_TS_p', icla
    call add_property_field_1d(f_name, f_label, igmsec(icla))
  enddo
endif

! Continuous phase variables
!---------------------------

! NB: 'c' stands for continuous <> 'p' stands for particles

if (i_comb_drift.ge.1) then

  ! Additional fields for drift velocity for the gas

  f_name= 'vd_c'
  call field_create(f_name, itycat, ityloc, idim3, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
  ! For log in the log
  call field_set_key_int(f_id, keylog, 1)

endif

! Mass fraction of the continuous phase (X1)
f_name= 'x_c'
call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)
call field_set_key_str(f_id, keylbl, f_name)

! Mass fraction of the continuous phase (X1) BOUNDARY VALUE
f_name= 'b_x_c'
call field_create(f_name, itycat, 3, idim1, inoprv, f_id)
call field_set_key_str(f_id, keylbl, f_name)

! Explicit interfacial source terms for x1 h1 (deduced from thoses of x2 h2)
f_name= 'x_h_c_exp_st'
call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)

! Implicit interfacial source terms for x1 h1 (deduced from thoses of x2 h2)
f_name= 'x_h_c_imp_st'
call field_create(f_name, itycat, ityloc, idim1, inoprv, f_id)

! Bulk
!-----

! Mass fraction of elements: C,  O,  H (used in balances)
call add_property_field_1d('x_carbone', 'Z_Carbone', ibcarbone)
call add_property_field_1d('x_oxygen', 'Z_Oxygen', iboxygen)
call add_property_field_1d('x_hydrogen', 'Z_Hydrogen', ibhydrogen)

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
