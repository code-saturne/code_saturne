!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file cs_user_initialization-time_step.f90
!>
!> \brief Time step example
!>
!> See \subpage cs_user_initialization for examples.
!>
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine cs_user_f_initialization &
 ( nvar   , nscal  ,                                              &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel

integer, allocatable, dimension(:) :: lstelt
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstelt(ncel)) ! temporary array for cells selection

!===============================================================================
! Time step modification
!
! We do a computation restart with a variable in time and constant in space
!  time step or with a variable in time and space time step.
!  We want to modify the time step given by the reading of the restart file
!  (in order to overcome a too slow evolution for instance).
!===============================================================================

if (isuite .eq. 1 .and. (idtvar.eq.1 .or. idtvar.eq.2)) then
  do iel = 1, ncel
    dt (iel) = 10.d0*dt(iel)
  enddo
endif
!< [init]

!--------
! Formats
!--------

!----
! End
!----

deallocate(lstelt)  ! temporary array for cells selection

return
end subroutine cs_user_f_initialization
