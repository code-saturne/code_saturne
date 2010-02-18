!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine uspt1d &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nfpt1d , iphas  , iappel ,          &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ifpt1d , nppt1d , iclt1d ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo ,                   &
   xyznod , volume ,                                              &
   tppt1d , rgpt1d , eppt1d ,                                     &
   tept1d , hept1d , fept1d ,                                     &
   xlmt1d , rcpt1d , dtpt1d ,                                     &
   dt     , rtpa   ,                                              &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! -------

!     User subroutine.

!     Data Entry ot the thermic module in 1-Dimension Wall.


! Introduction:
!=============

! Define the different values which can be taken by iappel:
!--------------------------------------------------------

! iappel = 1 (only one call on initialization):
!            Computation of the cells number where we impose a wall

! iappel = 2 (only one call on initialization):
!            Locating cells where we impose a wall
!            Data linked to the meshing.

! iappel = 3 (call on each time step):
!            Value of the physical computational coefficients and
!            boundary condition type on the exterior wall:
!            --------------------------------------------
!
!             iclt1d = 1 -> constant temperature imposed
!             iclt1d = 3 -> heat flux imposed

!            Initialization of the temperature on the wall.


! Boundary faces identification
! =============================

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nfpt1d           ! i  ! <-- ! number of faces with the 1-D thermic module    !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iappel           ! i  ! <-- ! data type to send                              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! maxelt           !  i ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ifpt1d           ! ia ! <-- ! number of the face treated                     !
! nppt1d           ! ia ! <-- ! number of discretized points                   !
! iclt1d           ! ia ! <-- ! boundary condition type                        !
!--begin. obsolesence ---------------------------------------------------------!
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfac(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
!--end obsolesence ------------------------------------------------------------!
! idevel(nideve)   ! ia ! <-- ! integer work array for temporary developpement !
! ituser(nituse)   ! ia ! <-- ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
!--new ------------------------------------------------------------------------!
! eppt1d           ! ra ! <-- ! wall thickness                                 !<--new
! rgpt1d           ! ra ! <-- ! geometric ratio of the meshing refinement      !<--new
! tppt1d           ! ra ! <-- ! wall temperature initialization                !<--new
! tept1d           ! ra ! <-- ! exterior temperature                           !<--new
! hept1d           ! ra ! <-- ! exterior exchange coefficient                  !<--new
! fept1d           ! ra ! <-- ! flux applied to the exterior                   !<--new
! xlmt1d           ! ra ! <-- ! lambda wall conductivity coefficient           !<--new
! rcpt1d           ! ra ! <-- ! rhoCp wall coefficient                         !<--new
! dtpt1d           ! ra ! <-- ! wall time step                                 !<--new
!--begin. obsolesence ---------------------------------------------------------!--!
! xyzcen           ! ra ! <-- ! cell centers                                   !  !
!  (ndim, ncelet)  !    !     !                                                !  !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !  !
!  (ndim, nfac)    !    !     !                                                !  !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !  !
!  (ndim, nfavor)  !    !     !                                                !  !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !  !
!  (ndim, nfac)    !    !     !                                                !  !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !  !
!  (ndim, nfabor)  !    !     !                                                !  !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !  !
!  (ndim, nnod)    !    !     !                                                !  !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !  !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !  !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !  !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !  !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !  !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !  !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !  !
!  (nfabor, *)     !    !     !                                                !  !
! coefu            ! ra ! --- ! work array                                     !  !
!  (nfabor, 3)     !    !     !  (computation of pressure gradient)            !  !
!-end obsolesence--------------------------------------------------------------!--!
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary developpement    !
! rtuser(nituse    ! ra ! <-- ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================

implicit none

!===============================================================================
! Data in common
!===============================================================================

include "paramx.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments
!-------------------------------------------------------------------
integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , nfpt1d
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas  , iappel

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          ifpt1d(nfpt1d), nppt1d(nfpt1d), iclt1d(nfpt1d)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision eppt1d(nfpt1d) , rgpt1d(nfpt1d) , tppt1d(nfpt1d)
double precision tept1d(nfpt1d) , hept1d(nfpt1d) , fept1d(nfpt1d)
double precision xlmt1d(nfpt1d) , rcpt1d(nfpt1d) , dtpt1d(nfpt1d)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables
!-------------------------------------------------------------------
integer          idebia , idebra
integer          ifbt1d , ii , ifac
integer          ilelt, nlelt

!===============================================================================

idebia = idbia0
idebra = idbra0


!===============================================================================
! Rereading of the restart file:
!----------------------------------

!     isuit1 = 0        --> No rereading
!                           (meshing and wall temperature reinitialization)
!     isuit1 = 1        --> Rereading of the restart file for the 1-Dimension
!                           thermic module
!     isuit1 = isuite   --> Rereading only if the computational fluid dynamic is
!                           a continuation of the computation.

!     The initialization of isuit1 is mandatory.
!===============================================================================

isuit1 = isuite

ifbt1d = 0

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
if(1.eq.1) return
!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

if(iappel.eq.1.or.iappel.eq.2) then

!===============================================================================
! Faces determining with the 1-D thermic module:
!----------------------------------------------
!
!     nfpt1d    : Total number of faces with the 1-D thermic module
!     ifpt1d(ii): Number of the (ii)th face with the 1-D thermic module

! Remarks:
!--------
!     During the rereading of the restart file, nfpt1d and ifpt1d are
!     compared with the other values from the restart file being the result of
!     the start or restarting computation.
!
!     A total similarity is required to continue with the previous computation.
!     Regarding the test case on ifpt1d, it is necessary that the array will be
!     arranged in increasing order
!               (as following : ifpt1d(jj) > ifpt1d(ii) si jj > ii).
!
!     If it is impossible, contact the developer team to deactivate this test.
!===============================================================================

  CALL GETFBR('3',NLELT,LSTELT)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    ifbt1d =ifbt1d + 1
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
!                 with the 1-D thermic module.
!     eppt1d(ii): wall thickness associated to the (ii)th face
!                 with the 1-D thermic module.
!     rgpt1d(ii): geometric progression ratio of the meshing refinement
!                 associated to the (ii)th face with the 1-D thermic module.
!                 (with : rgpt1d(ii) > 1 => small meshes  on the fluid side)
!     tppt1d(ii): wall temperature initialization associated to the (ii)th face
!                 with the 1-D thermic module.

! Remarks:
!--------
!     During the rereading of the restart file for the 1-D thermic module,
!     the tppt1d variable is not used.
!
!     The nfpt1d, eppt1d and rgpt1d variables are compared to the previous
!     values being the result of the restart file.
!
!     An exact similarity is necessary to continue with the previous computation.
!===============================================================================
if (iappel.eq.2) then
   if(iphas.eq.1) then
      do ii = 1, nfpt1d
        ifac = ifpt1d(ii)
        nppt1d(ii) = 8
        eppt1d(ii) = 0.01144d0
        rgpt1d(ii) = 1.d0
        tppt1d(ii) = 25.d0
      enddo
   endif
endif
!===============================================================================
! Padding of the wall exterior boundary conditions:
!-------------------------------------------------
!
!     iclt1d(ii): boundary condition type
!     ----------
!                  iclt1d(ii) = 1: dirichlet's condition ,  with exchange coefficient
!                  iclt1d(ii) = 3: flux condition
!
!     tept1d(ii): exterior temperature
!     hept1d(ii): exterior exchange coefficient
!     fept1d(ii): flux applied to the exterior (flux<0 = coming flux)
!     xlmt1d(ii): lambda wall conductivity coefficient (W/m/°C)
!     rcpt1d(ii): wall coefficient rho*Cp (J/m3/°C)
!     dtpt1d(ii): time step resolution of the thermic equation to the
!                 (ii)th border face with the 1-D thermic module (s)
!===============================================================================
if (iappel.eq.3) then
   if(iphas.eq.1) then
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
endif

!===============================================================================
! END of the uspt1d subroutine =====================================================
!===============================================================================
return
end subroutine

