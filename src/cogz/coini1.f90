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

!> \brief Initialize variables for combustion model.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!

subroutine coini1 () &
  bind(C, name='cs_f_coini1')

!===============================================================================
! Module files
!===============================================================================

use cstphy
use coincl
use ppincl
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

! Compute RO0 from T0 and P0

if (ippmod(islfm).ne.-1) then
  ro0 = flamelet_library(flamelet_rho, 1, 1, 1, 1)
endif

return
end subroutine

!---------------------------------------------------------------------------
