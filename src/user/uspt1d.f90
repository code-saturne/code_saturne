!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine uspt1d &
!================

 ( nvar   , nscal  , nfpt1d , iappel ,                            &
   ifpt1d , izft1d , nppt1d , iclt1d ,                            &
   tppt1d , rgpt1d , eppt1d ,                                     &
   tept1d , hept1d , fept1d ,                                     &
   xlmt1d , rcpt1d , dtpt1d ,                                     &
   dt     , rtpa   ,                                              &
   propce , propfb )

!===============================================================================
! Purpose:
! -------

!     User subroutine.

!     Data Entry ot the thermal module in 1-Dimension Wall.


! Introduction:
!=============

! Define the different values which can be taken by iappel:
!--------------------------------------------------------

! iappel = 1 (only one call on initialization):
!            Computation of the cells number where we prescribe a wall

! iappel = 2 (only one call on initialization):
!            Locating cells where we prescribe a wall
!            Data linked to the meshing.

! iappel = 3 (call on each time step):
!            Value of the physical computational coefficients and
!            boundary condition type on the exterior wall:
!            --------------------------------------------
!
!             iclt1d = 1 -> constant temperature prescribed
!             iclt1d = 3 -> heat flux prescribed

!            Initialization of the temperature on the wall.


! Boundary faces identification
! =============================

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the
! 'cs_user_boundary_conditions' subroutine,
! but a more thorough description can be found in the user guide.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nfpt1d           ! i  ! <-- ! number of faces with the 1-D thermal module    !
! iappel           ! i  ! <-- ! data type to send                              !
! ifpt1d           ! ia ! <-- ! number of the face treated                     !
! izft1d           ! ia ! <-- ! boundary faces zone for 1d-module definition   !
! nppt1d           ! ia ! <-- ! number of discretized points                   !
! iclt1d           ! ia ! <-- ! boundary condition type                        !
! eppt1d           ! ra ! <-- ! wall thickness                                 !
! rgpt1d           ! ra ! <-- ! geometric ratio of the meshing refinement      !
! tppt1d           ! ra ! <-- ! wall temperature initialization                !
! tept1d           ! ra ! <-- ! exterior temperature                           !
! hept1d           ! ra ! <-- ! exterior exchange coefficient                  !
! fept1d           ! ra ! <-- ! flux applied to the exterior                   !
! xlmt1d           ! ra ! <-- ! lambda wall conductivity coefficient           !
! rcpt1d           ! ra ! <-- ! rhoCp wall coefficient                         !
! dtpt1d           ! ra ! <-- ! wall time step                                 !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

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

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfb(nfabor,*)
double precision eppt1d(nfpt1d) , rgpt1d(nfpt1d) , tppt1d(nfpt1d)
double precision tept1d(nfpt1d) , hept1d(nfpt1d) , fept1d(nfpt1d)
double precision xlmt1d(nfpt1d) , rcpt1d(nfpt1d) , dtpt1d(nfpt1d)

! Local variables

integer          ifbt1d , ii , ifac
integer          ilelt, nlelt
integer          izone

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

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

izone = 0
ifbt1d = 0

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
if (1.eq.1) return
!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

if (iappel.eq.1.or.iappel.eq.2) then

!===============================================================================
! Faces determining with the 1-D thermal module:
!----------------------------------------------
!
!     nfpt1d    : Total number of faces with the 1-D thermal module
!     ifpt1d(ii): Number of the (ii)th face with the 1-D thermal module

! Remarks:
!--------
!     During the rereading of the restart file, nfpt1d and ifpt1d are
!     compared with the other values from the restart file being the result of
!     the start or restarting computation.
!
!     A total similarity is required to continue with the previous computation.
!     Regarding the test case on ifpt1d, it is necessary that the array be
!     arranged in increasing order (ifpt1d(jj) > ifpt1d(ii) if jj > ii).
!
!     If it is impossible, contact the developer team to deactivate this test.
!===============================================================================

  call getfbr('3', nlelt, lstelt)
  !==========

  izone = izone + 1

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)
    izft1d(ifac) = izone
    ifbt1d = ifbt1d + 1
    if (iappel.eq.2) ifpt1d(ifbt1d) = ifac

  enddo

endif

if (iappel.eq.1) then
  nfpt1d = ifbt1d
endif

!===============================================================================
! Parameters padding of the mesh and initialization:
!--------------------------------------------------
!
!     (Only one pass during the beginning of the computation)

!     nppt1d(ii): number of discretized points associated to the (ii)th face
!                 with the 1-D thermal module.
!     eppt1d(ii): wall thickness associated to the (ii)th face
!                 with the 1-D thermal module.
!     rgpt1d(ii): geometric progression ratio of the meshing refinement
!                 associated to the (ii)th face with the 1-D thermal module.
!                 (with : rgpt1d(ii) > 1 => small meshes  on the fluid side)
!     tppt1d(ii): wall temperature initialization associated to the (ii)th face
!                 with the 1-D thermal module.

! Remarks:
!--------
!     During the rereading of the restart file for the 1-D thermal module,
!     the tppt1d variable is not used.
!
!     The nfpt1d, eppt1d and rgpt1d variables are compared to the previous
!     values being the result of the restart file.
!
!     An exact similarity is necessary to continue with the previous computation.
!===============================================================================
if (iappel.eq.2) then
  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    nppt1d(ii) = 8
    eppt1d(ii) = 0.01144d0
    rgpt1d(ii) = 1.d0
    tppt1d(ii) = 25.d0
  enddo
endif
!===============================================================================
! Padding of the wall exterior boundary conditions:
!-------------------------------------------------
!
!     iclt1d(ii): boundary condition type
!     ----------
!                  iclt1d(ii) = 1: dirichlet condition,
!                                  with exchange coefficient
!                  iclt1d(ii) = 3: flux condition
!
!     tept1d(ii): exterior temperature
!     hept1d(ii): exterior exchange coefficient
!     fept1d(ii): flux applied to the exterior (flux<0 = coming flux)
!     xlmt1d(ii): lambda wall conductivity coefficient (W/m/C)
!     rcpt1d(ii): wall coefficient rho*Cp (J/m3/C)
!     dtpt1d(ii): time step resolution of the thermal equation to the
!                 (ii)th border face with the 1-D thermal module (s)
!===============================================================================
if (iappel.eq.3) then
  do ii = 1, nfpt1d
    iclt1d(ii) = 1
    ! Physical parameters
    ifac = ifpt1d(ii)
    if (cdgfbo(2,ifac).le.0.025d0) then
      iclt1d(ii) = 3
      fept1d(ii) = -1.d4
    else
      iclt1d(ii) = 3
      fept1d(ii) =  1.d4
    endif
    xlmt1d(ii) = 31.5d0
    rcpt1d(ii) = 3.5d6
    dtpt1d(ii) = 0.3d0
  enddo
endif

!===============================================================================
! End of the uspt1d subroutine
!===============================================================================

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine uspt1d

