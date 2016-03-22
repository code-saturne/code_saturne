!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Purpose:
! -------
!>
!> \file ushist.f90
!>
!> \brief Non-standard monitoring point definition.
!>
!> See \subpage us_hist for examples.
!

subroutine ushist &
!================

 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================
! Purpose:
! -------
!>
!> \brief Non-standard monitoring point definition.
!>
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode          name          role                                            !
!______________________________________________________________________________!
!> \param[in]    nvar          total number of variables
!> \param[in]    nscal         total number of scalars
!> \param[in]    dt            time step (per cell)
!______________________________________________________________________________!


!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use entsor
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          c_id, f_id0, f_dim

double precision dt(ncelet)

!===============================================================================

!----
! End
!----

return
end subroutine ushist
