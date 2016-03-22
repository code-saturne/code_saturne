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

!< [arg]

! Arguments

integer          nvar   , nscal  , nfpt1d
integer          iappel

integer          ifpt1d(nfpt1d), nppt1d(nfpt1d), iclt1d(nfpt1d)
integer          izft1d(nfabor)

double precision dt(ncelet)
double precision eppt1d(nfpt1d) , rgpt1d(nfpt1d) , tppt1d(nfpt1d)
double precision tept1d(nfpt1d) , hept1d(nfpt1d) , fept1d(nfpt1d)
double precision xlmt1d(nfpt1d) , rcpt1d(nfpt1d) , dtpt1d(nfpt1d)

!< [arg]

!< [loc_var_dec]

! Local variables

integer          ifbt1d , ii , ifac
integer          ilelt, nlelt
integer          izone

integer, allocatable, dimension(:) :: lstelt

!< [loc_var_dec]

!===============================================================================

!< [allocate]

! Allocate a temporary array for boundary faces selection
allocate(lstelt(nfabor))

!< [allocate]

!< [restart]

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

!< [restart]

!< [iappel_12]

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

!< [iappel_12]

!< [iappel_2]

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
    tppt1d(ii) = 25.d0 + tkelvi ! FIXME gerer les K et °C
  enddo
endif

!< [iappel_2]

!< [iappel_3]

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
!                 (ii)th boundary face with the 1-D thermal module (s)
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

!< [iappel_3]

!===============================================================================
! End of the uspt1d subroutine
!===============================================================================

!< [deallocate]

! Deallocate the temporary array
deallocate(lstelt)

!< [deallocate]

return
end subroutine uspt1d