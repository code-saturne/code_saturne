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

!> \file cfvarp.f90
!> \brief Variables location initialization for the compressible module,
!> according to calculation type selected by the user.
!>
!------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!-------------------------------------------------------------------------------
!   mode          name          role
!-------------------------------------------------------------------------------
!______________________________________________________________________________!

subroutine cfvarp

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
use field
use cs_cf_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

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

  subroutine cs_gui_labels_compressible()  &
    bind(C, name='cs_gui_labels_compressible')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_labels_compressible

end interface

!===============================================================================

if (ippmod(icompf).ge.0) then

  ! Pointers and reference values definition

  ! Total energy
  itherm = 3
  call add_model_scalar_field('total_energy', 'TotEner', ienerg)
  iscalt = ienerg

  ! Alias for B.C.
  irunh = ienerg

  ! Temperature (post)
  call add_model_scalar_field('temperature', 'TempK', itempk)

  ! Pointer and reference value for conductivity of temperature scalar
  ! TODO itempk should be a property
  call field_set_key_int (ivarfl(isca(itempk)), kivisl, -1)
  visls0(itempk) = epzero

  ! Pointer and reference value for diffusivity of total energy scalar
  call field_set_key_int (ivarfl(isca(ienerg)), kivisl, -1)
  visls0(ienerg) = epzero

  ! Pointer and reference value for volumetric molecular viscosity
  ! (constant by default)
  iviscv = 0
  viscv0 = 0.d0

  ! MAP to C API
  call cs_field_pointer_map_compressible

  ! Mapping for GUI
  if (iihmpr.eq.1) then
    call cs_gui_labels_compressible
  endif

  ! Computation parameters

  ! Variability of specific heat at constant volume Cv (constant by default)
  cv0 = 0.d0
  call cs_cf_set_thermo_options

  ! Variability of volumetric molecular viscosity (gui setting)
  if (iihmpr.eq.1) then
    call csvvva(iviscv)
  endif

endif

!--------
! Formats
!--------


return
end subroutine
