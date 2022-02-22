!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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
! --------
!> \file cs_fuel_prop.f90
!> \brief Define state variables for fuel combustion.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!______________________________________________________________________________!

subroutine cs_fuel_prop

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
use cs_fuel_incl
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

! Initialization

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
call add_property_field_1d('temperature', 'T_Gas', itemp)
call add_property_field_1d('rho_gas', 'Rho_Gas', irom1)

! Gas mixture fractions
call add_property_field_1d('ym_fo0',   'Ym_FO0',   iym1(1))
call add_property_field_1d('ym_fov',   'Ym_FOV',   iym1(2))
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
endif

! Dispersed phase (particle classes)
do icla = 1, nclafu
  write(f_name,  '(a,i2.2)') 't_fuel_', icla
  write(f_label, '(a,i2.2)') 'T_Fuel_', icla
  call add_property_field_1d(f_name, f_label, itemp2(icla))
enddo

do icla = 1, nclafu
  write(f_name,  '(a,i2.2)') 'rho_fuel_', icla
  write(f_label, '(a,i2.2)') 'Rho_Fuel_', icla
  call add_property_field_1d(f_name, f_label, irom2(icla))
enddo

do icla = 1, nclafu
  write(f_name,  '(a,i2.2)') 'diameter_fuel_', icla
  write(f_label, '(a,i2.2)') 'Diam_Drop_', icla
  call add_property_field_1d(f_name, f_label, idiam2(icla))
enddo

do icla = 1, nclafu
  write(f_name,  '(a,i2.2)') 'h1_hlf_', icla
  write(f_label, '(a,i2.2)') 'H1-Hlf_', icla
  call add_property_field_1d(f_name, f_label, ih1hlf(icla))
enddo

do icla = 1, nclafu
  write(f_name,  '(a,i2.2)') 'eva_fuel_', icla
  write(f_label, '(a,i2.2)') 'EVA_Fuel_', icla
  call add_property_field_1d(f_name, f_label, igmeva(icla))
enddo

do icla = 1, nclafu
  write(f_name,  '(a,i2.2)') 'het_ts_fuel_', icla
  write(f_label, '(a,i2.2)') 'Het_TS_Fuel_', icla
  call add_property_field_1d(f_name, f_label, igmhtf(icla))
enddo

if (i_comb_drift.ge.1) then
  do icla = 1, nclafu

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
call add_property_field_1d('balance_c', 'Balance_C', ibcarbone)
call add_property_field_1d('balance_o', 'Balance_O', iboxygen)
call add_property_field_1d('balance_h', 'Balance_H', ibhydrogen)

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
