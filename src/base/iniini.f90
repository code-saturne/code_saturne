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

!> \file iniini.f90
!> \brief Commons default initialization before handing over the user.
!>
!------------------------------------------------------------------------------

subroutine iniini () &
  bind(C, name='cs_f_iniini')

!===============================================================================
! Module files
!===============================================================================

use atincl
use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use ppincl
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer        ii, iscal

!===============================================================================
! 0. Global field ids
!===============================================================================

icrom = -1

!===============================================================================
! Map Fortran pointers to C global data
!===============================================================================

call atmo_init
call atmo_init_chem
call time_step_init
call time_step_options_init
call physical_constants_init
call fluid_properties_init

!===============================================================================
! Position of variables in numvar.f90
!===============================================================================

! Initialize mappings of field ids

do ii = 1, nvarmx
  ivarfl(ii) = -1
enddo

! Scalar to variable mappings

do iscal = 1, nscamx
  isca(iscal) = 0
enddo

return
end subroutine
