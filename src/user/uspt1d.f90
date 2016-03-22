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
!> \file uspt1d.f90
!>
!> \brief Data Entry of the thermal module in 1-Dimension Wall.
!>
!> See \subpage us_pt1d for examples.
!-------------------------------------------------------------------------------

subroutine uspt1d &
!================

 ( nvar   , nscal  , nfpt1d , iappel ,                            &
   ifpt1d , izft1d , nppt1d , iclt1d ,                            &
   tppt1d , rgpt1d , eppt1d ,                                     &
   tept1d , hept1d , fept1d ,                                     &
   xlmt1d , rcpt1d , dtpt1d ,                                     &
   dt     )

!===============================================================================
!> \brief Data Entry of the thermal module in 1-Dimension Wall.
!>
!> Introduction:
!>
!> Define the different values which can be taken by iappel:
!>
!>
!> iappel = 1 (only one call on initialization):
!>            Computation of the cells number where we prescribe a wall
!>
!> iappel = 2 (only one call on initialization):
!>            Locating cells where we prescribe a wall
!>            Data linked to the meshing.
!>
!> iappel = 3 (call on each time step):
!>            Value of the physical computational coefficients and
!>            boundary condition type on the exterior wall:
!>
!>
!>            iclt1d = 1 -> constant temperature prescribed
!>            iclt1d = 3 -> heat flux prescribed
!>
!>            Initialization of the temperature on the wall.
!>
!> Boundary faces identification
!>
!>
!> Boundary faces may be identified using the \ref getfbr subroutine.
!> The syntax of this subroutine is described in the
!> \ref cs_user_boundary_conditions subroutine,
!> but a more thorough description can be found in the user guide.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode         name          role                                             !
!______________________________________________________________________________!
!> \param[in]   nvar          total number of variables
!> \param[in]   nscal         total number of scalars
!> \param[in]   nfpt1d        number of faces with the 1-D thermal module
!> \param[in]   iappel        data type to send
!> \param[in]   ifpt1d        number of the face treated
!> \param[in]   izft1d        boundary faces zone for 1d-module definition
!> \param[in]   nppt1d        number of discretized points
!> \param[in]   iclt1d        boundary condition type
!> \param[in]   eppt1d        wall thickness
!> \param[in]   rgpt1d        geometric ratio of the meshing refinement
!> \param[in]   tppt1d        wall temperature initialization
!> \param[in]   tept1d        exterior temperature
!> \param[in]   hept1d        exterior exchange coefficient
!> \param[in]   fept1d        flux applied to the exterior
!> \param[in]   xlmt1d        lambda wall conductivity coefficient
!> \param[in]   rcpt1d        rhoCp wall coefficient
!> \param[in]   dtpt1d        wall time step
!> \param[in]   dt            time step (per cell)
!______________________________________________________________________________!


!===============================================================================
! Data in common
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , nfpt1d
integer          iappel

integer          ifpt1d(nfpt1d), nppt1d(nfpt1d), iclt1d(nfpt1d)
integer          izft1d(nfabor)

double precision dt(ncelet)
double precision eppt1d(nfpt1d) , rgpt1d(nfpt1d) , tppt1d(nfpt1d)
double precision tept1d(nfpt1d) , hept1d(nfpt1d) , fept1d(nfpt1d)
double precision xlmt1d(nfpt1d) , rcpt1d(nfpt1d) , dtpt1d(nfpt1d)

! Local variables

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
if (1.eq.1) return
!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))

!===============================================================================
! Rereading of the restart file:
!----------------------------------

!     isuit1 = 0        --> No rereading
!                           (meshing and wall temperature reinitialization)
!     isuit1 = 1        --> Rereading of the restart file for the 1-Dimension
!                           thermal module
!     isuit1 = isuite   --> Rereading only if the computational fluid dynamic is
!                           a continuation of the computation.

!     The initialization of isuit1 is mandatory.
!===============================================================================

isuit1 = isuite

!===============================================================================
! End of the uspt1d subroutine
!===============================================================================

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine uspt1d
