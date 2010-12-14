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

subroutine usmpst &
!================

 ( idbia0 , idbra0 , ipart  ,                                     &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   ncelps , nfacps , nfbrps ,                                     &
   nideve , nrdeve , nituse , nrtuse , imodif ,                   &
   itypps ,                                                       &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis ,                                     &
   tracel , trafac , trafbr , rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

! Modify list of cells or faces defining an existing post-processing
! output mesh; this subroutine is called for true (non-alias) user meshes,
! for each time step at which output on this mesh is active, and only if
! all writers associated with this mesh allow mesh modification
! (i.e. were defined with 'indmod' = 2 or 12).

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ipart            ! i  ! <-- ! number of the post-processing mesh (< 0 or > 0)!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nvlsta           ! i  ! <-- ! number of Lagrangian statistical variables     !
! ncelps           ! i  ! <-- ! number of cells in post-processing mesh        !
! nfacps           ! i  ! <-- ! number of interior faces in post-process. mesh !
! nfbrps           ! i  ! <-- ! number of boundary faces in post-process. mesh !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! imodif           ! i  ! --> ! 0 if the mesh was not modified by this call,   !
!                  !    !     ! 1 if it has been modified.                     !
! itypps(3)        ! ia ! <-- ! global presence flag (0 or 1) for cells (1),   !
!                  !    !     ! interior faces (2), or boundary faces (3) in   !
!                  !    !     ! post-processing mesh                           !
! lstcel(ncelps)   ! ia ! --> ! list of cells in post-processing mesh          !
! lstfac(nfacps)   ! ia ! --> ! list of interior faces in post-processing mesh !
! lstfbr(nfbrps)   ! ia ! --> ! list of boundary faces in post-processing mesh !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-- ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! statis           ! ra ! <-- ! statistic means                                !
!  (ncelet, nvlsta)!    !     !                                                !
! tracel(*)        ! ra ! --- ! work array for post-processed cell values      !
! trafac(*)        ! ra ! --- ! work array for post-processed face values      !
! trafbr(*)        ! ra ! --- ! work array for post-processed boundary face v. !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use optcal
use numvar
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ipart
integer          nvar   , nscal  , nphas  , nvlsta
integer          ncelps , nfacps , nfbrps
integer          nideve , nrdeve , nituse , nrtuse, imodif

integer          itypps(3)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision dt(ncelet), rtpa(ncelet,*), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision statis(ncelet,nvlsta)
double precision tracel(ncelps*3)
double precision trafac(nfacps*3), trafbr(nfbrps*3)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          ifac  , iphas
integer          ii   , jj
double precision vmin2, v2, w2

!===============================================================================

! Note:

! The 'itypps" array allows determining if the mesh contains at first cells,
! interior faces, or boundary faces (in a global sense when in parallel).

! This enables using "generic" selection criteria, which may function on any
! post-processing mesh, but if such a mesh is empty for a given call to this
! function, we will not know at the next call if it contained cells of faces.
! In this case, it may be preferable to use its number to decide if it should
! contain cells or faces.

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
!     1. TRAITEMENT DES MAILLAGES POST A REDEFINIR
!         A RENSEIGNER PAR L'UTILISATEUR aux endroits indiques
!===============================================================================

! Example: for user meshes 3 and 4, we only keep the cells of faces
!          at which the velocity is greater than a given threshold.

if (ipart.eq.3) then

  imodif = 1

  ncelps = 0
  nfacps = 0
  nfbrps = 0

  vmin2 = (0.5d0)**2

  ! If the mesh contains cells
  ! --------------------------

  if (itypps(1) .eq. 1) then

    do ii = 1, ncel

      iphas = 1

      v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2        &
           + rtp(ii, iw(iphas))**2
      if (v2 .ge. vmin2) then
        ncelps = ncelps + 1
        lstcel(ncelps) = ii
      endif

    enddo

  ! If the mesh contains interior faces
  ! -----------------------------------

  else if (itypps(2) .eq. 1) then

    do ifac = 1, nfac

      iphas = 1

      ii = ifacel(1, ifac)
      jj = ifacel(2, ifac)

      v2 =   rtp(ii, iu(iphas))**2   &
           + rtp(ii, iv(iphas))**2   &
           + rtp(ii, iw(iphas))**2
      w2 =   rtp(jj, iu(iphas))**2   &
           + rtp(jj, iv(iphas))**2   &
           + rtp(jj, iw(iphas))**2

      if (v2 .ge. vmin2 .or. w2 .ge. vmin2) then
        nfacps = nfacps + 1
        lstfac(nfacps) = ifac
      endif

    enddo

  ! If the mesh contains boundary faces
  ! -----------------------------------

  else if (itypps(3) .eq. 1) then

    do ifac = 1, nfabor

      iphas = 1

      ii = ifabor(ifac)

      v2 =   rtp(ii, iu(iphas))**2   &
           + rtp(ii, iv(iphas))**2   &
           + rtp(ii, iw(iphas))**2

      if (v2 .ge. vmin2) then
        nfbrps = nfbrps + 1
        lstfbr(nfbrps) = ifac
      endif

    enddo

  endif ! End of test on pre-existing mesh element types

else if (ipart.eq.4) then

  imodif = 1

  ncelps = 0
  nfacps = 0
  nfbrps = 0

  vmin2 = (0.5d0)**2

  ! Select interior faces
  ! ---------------------

  do ifac = 1, nfac

    iphas = 1

    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)

    v2 =   rtp(ii, iu(iphas))**2   &
         + rtp(ii, iv(iphas))**2   &
         + rtp(ii, iw(iphas))**2
    w2 =   rtp(jj, iu(iphas))**2   &
         + rtp(jj, iv(iphas))**2   &
         + rtp(jj, iw(iphas))**2

    if (     (v2 .ge. vmin2 .and. w2 .lt. vmin2)         &
        .or. (v2 .lt. vmin2 .and. w2 .ge. vmin2)) then
      nfacps = nfacps + 1
      lstfac(nfacps) = ifac
    endif

  enddo

  ! Select boundary faces
  ! ---------------------

  do ifac = 1, nfabor

    iphas = 1

    ii = ifabor(ifac)

    v2 =   rtp(ii, iu(iphas))**2   &
         + rtp(ii, iv(iphas))**2   &
         + rtp(ii, iw(iphas))**2

    if (v2 .ge. vmin2) then
      nfbrps = nfbrps + 1
      lstfbr(nfbrps) = ifac
    endif

  enddo

endif ! end of test on post-processing mesh number

return

end subroutine
