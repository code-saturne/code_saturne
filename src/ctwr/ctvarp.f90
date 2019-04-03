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
! Function:
! ---------

!> \file ctvarp.f90
!>
!> \brief Declare additional transported variables for cooling towers module.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     use_atmo      Also use atmospherique module if greater than 0
!_______________________________________________________________________________


subroutine ctvarp(use_atmo)

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
use ctincl
use field
use atincl, only : iymw

!===============================================================================

implicit none

! Arguments

integer          use_atmo

! Local variables

integer          keyccl, keydri
integer          kscmin, kscmax
integer          icla, ifcvsl, iscdri, f_id

!===============================================================================

! Key id of the scalar class
call field_get_key_id("scalar_class", keyccl)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

!===============================================================================
! 1. Definition of fields
!===============================================================================

! Bulk definition - For cooling towers, the bulk is the humid air
! By definition, humid air is composed of two species: dry air and
! water vapour (whether in gas or condensate form)
! ---------------------------------------------------------------

! Thermal model - Set parameters of calculations (module optcal)

itherm = 1  ! Solve for temperature of the bulk (humid air)

itpscl = 2  ! Temperature in Celsius

icp = 0     ! Cp is variable (>=0 means variable, -1 means constant)
            ! It has to vary with humidity
            ! Needs to be specified here because the automated creation and initialisation
            ! of the cell array for Cp in 'iniva0' depends on its value
            ! (unlike the cell arrays for density and viscosity which are initialised
            ! irrespective of the values of irovar and ivivar)

! The thermal transported scalar is the temperature of the bulk.
! If the atmospheric module in switch off, we create the field
if (ippmod(iatmos).ne.2) then

  call add_model_scalar_field('temperature', 'Temperature humid air', iscalt)

endif

f_id = ivarfl(isca(iscalt))

ifcvsl = 0 ! Set variable diffusivity for the humid air enthalpy
           ! The diffusivity used in the transport equation will be
           ! the cell value of the viscls array for f_id
           ! This value is updated at the top of each time step in 'ctphyv'
           ! along with the other variable properties
call field_set_key_int(f_id, kivisl, ifcvsl)


! Rain zone phase variables
!--------------------------

! Associate the injected liquid water with class 1
icla = 1

! NB: 'c' stands for continuous <> 'p' stands for particles

! Activate the drift for all scalars with key "drift" > 0
iscdri = 1

!TODO make it optionnal
! Mass fraction of liquid
call add_model_scalar_field('y_p', 'Yp liq', iy_p_l)
f_id = ivarfl(isca(iy_p_l))

! Set min and max clipping
call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 1.d0)

call field_set_key_int(f_id, keyccl, icla) ! Set the class index for the field

! Scalar with drift: Create additional mass flux
! This flux will then be reused for all scalars associated with this class

! GNU function to return the value of iscdri
! with the bit value of iscdri at position
! 'DRIFT_SCALAR_ADD_DRIFT_FLUX' set to one
iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)

call field_set_key_int(f_id, keydri, iscdri)

ifcvsl = -1 ! Set constant diffusivity for the injected liquid mass fraction
            ! The diffusivity used in the transport equation will be
            ! the value of visls0(iy_p_l)
call field_set_key_int(f_id, kivisl, ifcvsl)

! Transport and solve for the temperature of the liquid - with the same drift
! as the mass fraction Y_l in rain zones
! NB: Temperature of the liquidus must be transported after the bulk enthalpy
call add_model_scalar_field('y_p_t_l', 'Tp liq', it_p_l)
f_id = ivarfl(isca(it_p_l))

call field_set_key_int(f_id, keyccl, icla)

! Scalar with drift, but do not create an additional mass flux (use 'ibclr' instead of 'ibset')
! for the enthalpy.  It reuses the mass flux of already identified with the mass fraction
iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)

call field_set_key_int(f_id, keydri, iscdri)

ifcvsl = 0   ! Set variable diffusivity for the injected liquid enthalpy transport
             ! The diffusivity used in the transport equation will be
             ! the cell value of the viscls array for f_id
call field_set_key_int(f_id, kivisl, ifcvsl)


! Packing zones variables
!-------------------------

! Associate the injected liquid water with class 2
icla = 2

! Injected liquid water definition - This is the separate phase
! which is injected in the packing zones.  Not to be confused with
! the water content in the humid air.

! Activate the drift for all scalars with key "drift" > 0
iscdri = 1

! Mass fraction of liquid
call add_model_scalar_field('y_l_packing', 'Yl packing', iyml)
f_id = ivarfl(isca(iyml))

! Set min clipping
call field_set_key_double(f_id, kscmin, 0.d0)

call field_set_key_int(f_id, keyccl, icla) ! Set the class index for the field

! Scalar with drift: Create additional mass flux
! This flux will then be reused for all scalars associated with this class

! GNU function to return the value of iscdri
! with the bit value of iscdri at position
! 'DRIFT_SCALAR_ADD_DRIFT_FLUX' set to one
iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)

! GNU function to return the value of iscdri
! with the bit value of iscdri at position
! 'DRIFT_SCALAR_IMPOSED_MASS_FLUX' set to one
iscdri = ibset(iscdri, DRIFT_SCALAR_IMPOSED_MASS_FLUX)

call field_set_key_int(f_id, keydri, iscdri)

ifcvsl = -1 ! Set constant diffusivity for the injected liquid mass fraction
            ! The diffusivity used in the transport equation will be
            ! the value of f_id
call field_set_key_int(f_id, kivisl, ifcvsl)

! Transport and solve for the enthalpy of the liquid - with the same drift
! as the mass fraction Y_l
! NB: Enthalpy of the liquidus must be transported after the bulk enthalpy
call add_model_scalar_field('enthalpy_liquid', 'Enthalpy liq', ihml) ! TODO x_p_h_l or y_p_h_2
f_id = ivarfl(isca(ihml))

call field_set_key_int(f_id, keyccl, icla)

! Scalar with drift, but do not create an additional mass flux (use 'ibclr' instead of 'ibset')
! for the enthalpy.  It reuses the mass flux of already identified with the mass fraction
iscdri = ibclr(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)

iscdri = ibset(iscdri, DRIFT_SCALAR_IMPOSED_MASS_FLUX)

call field_set_key_int(f_id, keydri, iscdri)

ifcvsl = 0   ! Set variable diffusivity for the injected liquid enthalpy transport
             ! The diffusivity used in the transport equation will be
             ! the cell value of the viscls array for f_id
call field_set_key_int(f_id, kivisl, ifcvsl)


! Continuous phase variables
!---------------------------

icla = -1

! NB: 'c' stands for continuous <> 'p' stands for particles

! If not using the atmospheric module, we create the fields
if (ippmod(iatmos).ne.2) then

  ! Total mass fraction water in the bulk, humid air
  call add_model_scalar_field('ym_water', 'Ym water', iymw)
  f_id = ivarfl(isca(iymw))

  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)
endif

call field_set_key_int(f_id, keyccl, icla)

ifcvsl = -1 ! Set constant diffusivity for the dry air mass fraction
            ! The diffusivity used in the transport equation will be
            ! the value of visls0 of f_id
call field_set_key_int(f_id, kivisl, ifcvsl)

! Activate the drift for all scalars with key "drift" > 0
iscdri = 1
! Activated drift. As it is the continuous phase class (==-1)
! the convective flux is deduced for classes > 0 and bulk class (==0)
iscdri = ibset(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)
call field_set_key_int(f_id, keydri, iscdri)

end subroutine
