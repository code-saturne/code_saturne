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

!> \file addfld.f90
!>
!> \brief Add additional property fields for dedicated modules
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________


subroutine ppprop

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
use atincl
use ppincl
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

character(len=80) :: f_name

integer          f_id
integer          itycat

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_elec_add_property_fields()  &
    bind(C, name='cs_elec_add_property_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_elec_add_property_fields

  subroutine cs_gas_mix_add_property_fields()  &
    bind(C, name='cs_gas_mix_add_property_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gas_mix_add_property_fields

end interface

!===============================================================================

! ---> Physique particuliere : Combustion Gaz

if (ippmod(icod3p).ge.0 .or. ippmod(islfm ).ge.0 .or.            &
    ippmod(icoebu).ge.0 .or. ippmod(icolwc).ge.0) then
  call coprop
endif

! ---> Physique particuliere :  Combustion Charbon Pulverise

if (ippmod(iccoal).ge.0) then
  call cs_coal_prop
endif

! ---> Physique particuliere :  Combustion Charbon Pulverise
!      Couplee Transport Lagrangien des particules de charbon

if (ippmod(icpl3c).ge.0) then
  call cplpro
endif

! ---> Physique particuliere : Combustion Fuel

if (ippmod(icfuel).ge.0) then
  call cs_fuel_prop
endif

! ---> Physique particuliere : Compressible

if (ippmod(icompf).ge.0) then
  call cfprop
endif

! ---> Physique particuliere : Versions electriques

if (ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1) then
  call cs_elec_add_property_fields
endif

! ---> Atmospheric modules:
if (ippmod(iatmos).ge.0) then

  ! Momentum source terms
  if (iatmst.ge.1) then
    call add_property_field('momentum_source_terms', 'MomentumSourceTerms', 3, .false., imomst)
    call field_set_key_int(imomst, keylog, 1)
    call field_set_key_int(imomst, keyvis, 1)
  endif

  call atprop
endif

! ---> Cooling towers model
if (ippmod(iaeros).ge.0) then
  call add_property_field_1d('humidity', 'Humidity', ihumid)
  call add_property_field_1d('x_s', 'Humidity sat', f_id)
  call add_property_field_1d('enthalpy', 'Enthalpy humid air', ihm)
  call add_property_field_1d('temperature_liquid', 'Temp liq', itml)
  call add_property_field_1d('vertvel_l', 'Vertical vel liq', ivertvel)

  ! Continuous phase properties
  !----------------------------

  ! NB: 'c' stands for continuous <> 'p' stands for particles

  ! Mass fraction of the continuous phase (X1)
  f_name= 'x_c'
  call add_property_field_1d('x_c', 'Gas mass fraction', f_id)

  ! Mass fraction of the continuous phase (X1) BOUNDARY VALUE
  f_name= 'b_x_c'
  itycat = FIELD_INTENSIVE + FIELD_PROPERTY
  call field_create(f_name,  &
                    itycat,  &
                    3,       & ! location (boundary faces)
                    1,       &! dimension
                    .false., & ! Has previous ?
                    f_id)
  call field_set_key_str(f_id, keylbl, f_name)

endif

! Add the mixture molar mass fraction field
if (ippmod(igmix).ge.0) then

  call cs_gas_mix_add_property_fields

  call field_get_id('mix_mol_mas', igmxml)

  if (ippmod(igmix).ge.0 .and. ippmod(igmix).le.5) then
    iddgas = cs_gas_mix_species_to_field_id(nscasp)
  endif

endif

end subroutine ppprop
