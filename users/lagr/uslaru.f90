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

subroutine uslaru &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   itypfb , itrifb , itepa  ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   surfbn , dt     , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , croule , auxl  ,                    &
   distpa , distyp ,                                              &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! --------
!
! User subroutine of the Lagrangian particle-tracking module:
! -----------------------------------------
!
! User subroutine (non-mandatory intervention)

! Calculation of the function of significance for the Russian roulette


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
!                  !    !     !                                                !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! ntersl           ! i  ! <-- ! number of source terms of return coupling      !
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! nideve nrdeve    ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse nrtuse    ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel           ! ia ! <-- ! interior faces -> cells connectivity           !
! (2, nfac)        !    !     !                                                !
! ifabor           ! ia ! <-- ! boundary faces -> cells connectivity           !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! ia ! <-- ! boundary face family numbers                   !
! (nfabor)         !    !     !                                                !
! ifmcel           ! ia ! <-- ! cell family numbers                            !
! (ncelet)         !    !     !                                                !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml,nprfml    !    !     !                                                !
! ipnfac           ! ia ! <-- ! interior faces -> vertices index (optional)    !
!   (lndfac)       !    !     !                                                !
! nodfac           ! ia ! <-- ! interior faces -> vertices list (optional)     !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! ia ! <-- ! boundary faces -> vertices index (optional)    !
!   (lndfbr)       !    !     !                                                !
! nodfbr           ! ia ! <-- ! boundary faces -> vertices list  (optional)    !
!   (nfabor+1)     !    !     !                                                !
! itypfb(nfabor    ! ia ! <-- ! type of the boundary faces                     !
!  nphas)          !    !     !                                                !
! itrifb(nfabor    ! ia ! --> ! indirection for the sorting of the             !
!  nphas)          !    !     ! boundary faces                                 !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! idevel(nideve    ! ia ! <-- ! complementary dev. array of integers           !
! ituser(nituse    ! ia ! <-- ! complementary user array of integers           !
! ia(*)            ! ia ! --- ! macro array of integers                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! surfbn(nfabor    ! ra ! <-- !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! transported variables at cell centers for      !
! (ncelet,*)       !    !     ! the previous timestep                          !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! coefa, coefb     ! ra ! <-- ! boundary conditions at the boundary faces      !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! vagaus           ! ra ! <-- ! Gaussian random variables                      !
!(nbpmax,nvgaus    !    !     !                                                !
! croule(ncelet    ! ra ! --> ! function of significance for                   !
!                  !    !     ! the Russian roulette                           !
! auxl(nbpmax,3    ! ra ! --- !                                                !
! rdevel(nrdeve    ! ra ! <-- ! dev. complementary array of reals              !
! rtuser(nrtuse    ! ra ! <-- ! user complementary array of reals              !
! ra(*)            ! ra ! --- ! macro array of reals                           !
! distpa(ncelet    ! ra ! <-- ! wall-normal distance arrays                    !
! disty(ncelet)    ! ra ! <-- ! y+ distance                                    !
! w1...w3(ncel)    ! ra ! --- ! work arrays                                    !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array

!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "numvar.f90"
include "optcal.f90"
include "entsor.f90"
include "cstphy.f90"
include "pointe.f90"
include "period.f90"
include "parall.f90"
include "lagpar.f90"
include "lagran.f90"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          itepa(nbpmax,nivep)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision surfbn(nfabor)
double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision vagaus(nbpmax,*) , croule(ncelet)
double precision auxl(nbpmax,3)
double precision distpa(ncelet) , distyp(ncelet)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iel
double precision zref

!===============================================================================


!===============================================================================
! 0.  Memory management
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. Default initialization
!---------------------------

!     Caution : the croule parameter is only initialized in this subroutine.
!               Make sure that it is prescribed for every cell.


!===============================================================================

do iel = 1,ncel
  croule(iel) = 1.d0
enddo

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! -1.  If the user does not intervene, croule = 1 everywhere
!===============================================================================

if(1.eq.1) then
  return
endif

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 2. Calculation of a user-defined function of significance
!===============================================================================

!   CAUTION:   the croule array must be filled with positive
!   ^^^^^^^^^  real numbers enabling to weight the importance
!              of some zones with respect to others.
!
!              (the greater croule, the more important the zone)

!              For instance, we can decide that the zone is as important
!              as it is close to a position z=zref; with an importance equal
!              to 1.e-3 near zref, and with an importance equal to 1.e-6
!              far from zref.


zref = 0

do iel = 1,ncel
  croule(iel) = 1.d0/(max( abs(xyzcen(3,iel)-zref),1.d-3 ))
enddo

do iel = 1,ncel
  croule(iel) = max(croule(iel),1.d-6 )
enddo

!===============================================================================

!----
! End
!----

end subroutine
