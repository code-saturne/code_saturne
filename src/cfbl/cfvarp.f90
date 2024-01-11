!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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
!> \brief Variables definition initialization for the compressible module,
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
use field
use cs_c_bindings
use cs_cf_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer keyrf

type(var_cal_opt) :: vcopt

!===============================================================================

if (ippmod(icompf).ge.0) then

  ! Pointers and reference values definition

  ! Total energy
  itherm = 3
  call add_model_scalar_field('total_energy', 'TotEner', ienerg)
  iscalt = ienerg

  ! Alias for B.C.
  irunh = isca(ienerg)

  ! Temperature (post)
  call add_model_scalar_field('temperature', 'TempK', itempk)

  ! Pointer and reference value for conductivity of temperature scalar
  ! TODO itempk should be a property
  call field_set_key_int (ivarfl(isca(itempk)), kivisl, -1)
  call field_set_key_double(ivarfl(isca(itempk)), kvisl0, epzero)

  ! Pointer and reference value for diffusivity of total energy scalar
  call field_set_key_int (ivarfl(isca(ienerg)), kivisl, -1)
  call field_set_key_double(ivarfl(isca(ienerg)), kvisl0, epzero)

  ! Mixture fractions (two-phase homogeneous flows)
  if (ippmod(icompf).eq.2) then
    ! Volume fraction of phase 1 (with respect to the EOS parameters)
    call add_model_scalar_field('volume_fraction', 'Volume Fraction', ifracv)

    ! Mass fraction of phase 1
    call add_model_scalar_field('mass_fraction', 'Mass Fraction', ifracm)

    ! Energy fraction of phase 1
    call add_model_scalar_field('energy_fraction', 'Energy Fraction', ifrace)

    ! Pointer and reference value for diffusivity of three fractions
    call field_set_key_int (ivarfl(ifracv), kivisl, -1)
    call field_set_key_int (ivarfl(ifracm), kivisl, -1)
    call field_set_key_int (ivarfl(ifrace), kivisl, -1)
    call field_set_key_double(ivarfl(isca(ifracv)), kvisl0, epzero)
    call field_set_key_double(ivarfl(isca(ifracm)), kvisl0, epzero)
    call field_set_key_double(ivarfl(isca(ifrace)), kvisl0, epzero)

    ! Pure convection equation for three fractions
    call field_get_key_struct_var_cal_opt(ivarfl(ifracv), vcopt)
    vcopt%idifft = 0
    call field_set_key_struct_var_cal_opt(ivarfl(ifracv), vcopt)

    call field_get_key_struct_var_cal_opt(ivarfl(ifracm), vcopt)
    vcopt%idifft = 0
    call field_set_key_struct_var_cal_opt(ivarfl(ifracm), vcopt)

    call field_get_key_struct_var_cal_opt(ivarfl(ifrace), vcopt)
    vcopt%idifft = 0
    call field_set_key_struct_var_cal_opt(ivarfl(ifrace), vcopt)

    ! Set restart file for fractions
    call field_get_key_id('restart_file', keyrf)
    call field_set_key_int (ivarfl(ifracv), keyrf, RESTART_MAIN)
    call field_set_key_int (ivarfl(ifracm), keyrf, RESTART_MAIN)
    call field_set_key_int (ivarfl(ifrace), keyrf, RESTART_MAIN)
  endif

endif
!--------
! Formats
!--------

return
end subroutine
