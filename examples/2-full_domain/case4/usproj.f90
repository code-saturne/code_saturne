!-------------------------------------------------------------------------------

!                      Code_Saturne version 2.0.0-rc1
!                      --------------------------

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

subroutine usproj &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvep   , nivep  , ntersl , nvlsta , nvisbr , &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr , itepa  ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , statis , stativ , tslagr , parbor , &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Called at end of each time step, very general purpose
!    (i.e. anything that does not have another dedicated user subroutine)


! Several examples are given here:

!  - compute a thermal balance
!    (if needed, see note  below on adapting this to any scalar)

!  - compute global efforts on a subset of faces

!  - arbitrarily modify a calculation variable

!  - extract a 1 d profile

!  - print a moment

!  - examples on using parallel utility functions

! These examples are valid when using periodicity (iperio .gt. 0)
! and in parallel (irangp .ge. 0).

! The thermal balance compution also illustates a few other features,
! including the required precautions in parallel or with periodicity):
! - gradient calculation
! - computation of a value depending on cells adjacent to a face
!   (see synchronization of Dt and Cp)
! - computation of a global sum in parallel (parsom)


! Cells, boundary faces and interior faces identification
! =======================================================

! Cells, boundary faces and interior faces may be identified using
! the subroutines 'getcel', 'getfbr' and 'getfac' (respectively).

!  getfbr(string, nelts, eltlst):
!  - string is a user-supplied character string containing selection criteria;
!  - nelts is set by the subroutine. It is an integer value corresponding to
!    the number of boundary faces verifying the selection criteria;
!  - lstelt is set by the subroutine. It is an integer array of size nelts
!    containing the list of boundary faces verifying the selection criteria.

!  string may contain:
!  - references to colors (ex.: 1, 8, 26, ...)
!  - references to groups (ex.: inlet, group1, ...)
!  - geometric criteria (ex. x < 0.1, y >= 0.25, ...)
!  These criteria may be combined using logical operators ('and', 'or') and
!  parentheses.
!  Example: '1 and (group2 or group3) and y < 1' will select boundary faces
!  of color 1, belonging to groups 'group2' or 'group3' and with face center
!  coordinate y less than 1.

! Similarly, interior faces and cells can be identified using the 'getfac'
! and 'getcel' subroutines (respectively). Their syntax are identical to
! 'getfbr' syntax.

! For a more thorough description of the criteria syntax, it can be referred
! to the user guide.


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
! nbpmax           ! i  ! <-- ! max. number of particles allowed               !
! nvp              ! i  ! <-- ! number of particle-defined variables           !
! nvep             ! i  ! <-- ! number of real particle properties             !
! nivep            ! i  ! <-- ! number of integer particle properties          !
! ntersl           ! i  ! <-- ! number of return coupling source terms         !
! nvlsta           ! i  ! <-- ! number of Lagrangian statistical variables     !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! maxelt           ! i  ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! itepa            ! ia ! <-- ! integer particle attributes                    !
!  (nbpmax, nivep) !    !     !   (containing cell, ...)                       !
! idevel(nideve)   ! ia ! <-- ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-- ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp, ettpa      ! ra ! <-- ! particle-defined variables                     !
!  (nbpmax, nvp)   !    !     !  (at current and previous time steps)          !
! tepa             ! ra ! <-- ! real particle properties                       !
!  (nbpmax, nvep)  !    !     !  (statistical weight, ...                      !
! statis           ! ra ! <-- ! statistic means                                !
!  (ncelet, nvlsta)!    !     !                                                !
! stativ(ncelet,   ! ra ! <-- ! accumulator for variance of volume statisitics !
!        nvlsta -1)!    !     !                                                !
! tslagr           ! ra ! <-- ! Lagrangian return coupling term                !
!  (ncelet, ntersl)!    !     !  on carrier phase                              !
! parbor           ! ra ! <-- ! particle interaction properties                !
!  (nfabor, nvisbr)!    !     !  on boundary faces                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-- ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "dimfbr.h"
include "paramx.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"
include "parall.h"
include "period.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          itepa(nbpmax,nivep)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,nvlsta), stativ(ncelet,nvlsta-1)
double precision tslagr(ncelet,ntersl)
double precision parbor(nfabor,nvisbr)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iel

double precision sum, sumvol

!===============================================================================



!===============================================================================
! 1. Initialization
!===============================================================================

! ---> Extra memory handling

idebia = idbia0
idebra = idbra0


! Opening of file moy.dat at the first time step of the calculation
!   ntcabs = current time step
!   ntpabs = last time step of former calculation (if restart) or 0
! Only the first processor opens the file
!   irangp = rank of current processor (=-1 if single-processor calculation)
if (ntcabs.eq.ntpabs+1 .and. irangp.le.0) then
  open(99,file="moy.dat")
endif

sum = 0.d0
sumvol = 0.d0
do iel = 1, ncel
  sum = sum + rtp(iel,isca(1))*volume(iel)
  sumvol = sumvol + volume(iel)
enddo

! If the computation is done on more than one processor,
!   the "local" sum calculated above must be cumulated on all the processors.
!   parsom = replaces the argument by its sum over all the processors
!            (Code_Saturne routine encapsulating MPI commands).
if (irangp.ge.0) then
  call parsom(sum)
  call parsom(sumvol)
endif
sum = sum/sumvol

! Only the first processor writes
if (irangp.le.0) write(99,99) ntcabs,sum

 99   format(i6,g15.8)

! Close file moy.dat at last time step (ntmabs).
if (ntcabs.eq.ntmabs .and. irangp.le.0) then
  close(99)
endif

return
end subroutine
