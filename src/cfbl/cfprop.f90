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

!> \file cfprop.f90
!> \brief Properties definition initialization for the compressible module,
!> according to calculation type selected by the user.
!>
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!   mode          name          role
!-------------------------------------------------------------------------------
!______________________________________________________________________________!

subroutine cfprop

!===============================================================================
! Module files
!===============================================================================

use paramx
use ihmpre
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppincl
use field
use cs_cf_bindings

!===============================================================================

implicit none

! Local variables

integer       ifcvsl

!===============================================================================

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_field_pointer_map_compressible()  &
    bind(C, name='cs_field_pointer_map_compressible')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_field_pointer_map_compressible

end interface

!===============================================================================

! Variability of specific heat at constant volume Cv (constant by default)
icv = -1
cv0 = 0.d0
call cs_cf_set_thermo_options

! Variability of volumetric molecular viscosity (gui setting)
if (iihmpr.eq.1) then
  call csvvva(iviscv)
endif

! User settings: variability of molecular thermal conductivity and volume
! viscosity
call uscfx1

! Dynamic viscosity of reference of the scalar total energy (ienerg).
call field_get_key_int(ivarfl(isca(itempk)), kivisl, ifcvsl)
if (ifcvsl.ge.0 .or. icv.ge.0) then
  call field_set_key_int(ivarfl(isca(ienerg)), kivisl, 0)
else
  call field_set_key_int(ivarfl(isca(ienerg)), kivisl, -1)
endif

! Properties definition initialization according to their variability
if (icv.ge.0) then
  call add_property_field_1d('specific_heat_const_vol', &
                          'Specific_Heat_Const_Vol', &
                          icv)
  call hide_property(icv)
endif

if (iviscv.ge.0) then
  call add_property_field_1d('volume_viscosity', &
                          'Volume_Viscosity', &
                          iviscv)
  call hide_property(iviscv)
endif

! MAP to C API
call cs_field_pointer_map_compressible

return
end subroutine
