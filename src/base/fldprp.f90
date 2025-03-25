!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2025 EDF S.A.
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

!> \file fldprp.f90
!> \brief Properties definition initialization, according to calculation type
!> selected by the user.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine fldprp () &
 bind(C, name='cs_f_fldprp')

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use entsor
use atincl
use ppincl
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

!===============================================================================
! 1. PROPRIETES PRINCIPALES
!===============================================================================

! Base properties, always present

call field_get_id_try('density', icrom)

! ---> Atmospheric modules:
if (ippmod(iatmos).ge.0) then
  call field_get_id_try("real_temperature", itempc)
  if (ippmod(iatmos).eq.2) then
    call field_get_id("liquid_water", iliqwt)
  end if
endif

!====
! End
!====

return
end subroutine fldprp
