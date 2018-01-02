!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file cs_user_initialization-base.f90
!>
!> \brief Basic examples
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
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          iel, iutile
integer, allocatable, dimension(:) :: lstelt
double precision, dimension(:), pointer :: cvar_scal
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!< [alloc]
allocate(lstelt(ncel)) ! temporary array for cells selection
!< [alloc]

!===============================================================================
! Variables initialization:
!
!   isca(1) is the number related to the first user-defined scalar variable.
!   cvar_scal(iel) is the value of this variable in cell number iel.
!
!   ONLY done if there is no restart computation
!===============================================================================

!< [init]
if (isuite.eq.0) then
  call field_get_val_s(ivarfl(isca(1)), cvar_scal)

  do iel = 1, ncel
    cvar_scal(iel) = 25.d0
  enddo

endif
!< [init]

!--------
! Formats
!--------

!----
! End
!----

!< [finalize]
deallocate(lstelt)  ! temporary array for cells selection
!< [finalize]

return
end subroutine cs_user_f_initialization
