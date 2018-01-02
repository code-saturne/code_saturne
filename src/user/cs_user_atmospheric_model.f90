!-------------------------------------------------------------------------------

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
! Function:
! ---------

!> \file cs_user_atmospheric_model.f90
!>
!> \brief User subroutine dedicated to the atmospheric model.
!>
!> See \subpage cs_user_atmospheric_model for examples.

subroutine usatdv &
     ( imode )

!===============================================================================
!  Purpose:
!  -------
!> \brief Atmospheric module subroutine
!>
!> User definition of the vertical 1D arrays
!> User initialization of corresponding 1D ground model
!>
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     imode        number of calls of usatdv
!______________________________________________________________________________!

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl
use atsoil
use mesh

implicit none

!===============================================================================

! Arguments

integer           imode

return
end subroutine usatdv


!===============================================================================


subroutine usatsoil &
     !==================
     ( iappel )

!===============================================================================
! Purpose:
! -------
!
!> \brief Data Entry for the atmospheric ground model.
!>
!>
!> Introduction:
!>
!> Define the different values which can be taken by iappel:
!>
!> iappel = 1 (only one call on initialization):
!>            Computation of the cells number where we impose a
!>            Ground Model
!>
!> iappel = 2 (only one call on initialization):
!>            users may defined the ground face composition
!>            Warning : be coherent with the dimension of the array \c pourcent_sol
!>            It's also possible to modify the \c tab_sol array of the ground
!>            type constants
!
!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use atincl
use atsoil
use mesh

!===============================================================================

implicit none

! Arguments
!-------------------------------------------------------------------
integer          iappel

! Local variables

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END
!===============================================================================

allocate(lstelt(nfabor))

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine usatsoil