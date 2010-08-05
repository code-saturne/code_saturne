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

subroutine uslaen &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivarl  , ivarl1 , ivarlm , iflu   , ilpd1  , icla   ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis , stativ , tracel ,                   &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! --------

!   Subroutine of the Lagrangian particle-tracking module :
!   -------------------------------------------------------

!    User subroutine (non-mandatory intervention)

!    For the writing of the listing and the post-processing:
!    Average of the Lagrangian volume statistical variables
!    Possible intervention of the user.

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
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
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
! idevel(nideve    ! ia ! <-- ! complementary dev. array of integers           !
! ituser(nituse    ! ia ! <-- ! complementary user array of integers           !
! ivarl            !  i ! <-- ! number of the stat (between 1 and nvlsta)      !
! ivarl1           !  i ! <-- ! number of the global stat + group              !
!                  !    !     ! (average or variance)                          !
! ivarlm           !  i ! <-- ! number of the stat mean + group                !
! iflu             !  i ! <-- ! 0: mean of the stat ivarl/ivarl1               !
!                  !    !     ! 1: variance of the stat ivarl/ivarl1           !
! ilpd1            !  i ! <-- ! "pointer" to global statistical wqeight        !
!                  !    !     !                                                !
! icla             !  i ! <-- ! 0: global statistic                            !
                   !    ! <-- ! !=0: stat for the icla group                   !
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
! volume           ! ra ! <-- ! cell volumes                                   !
! (ncelet          !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! transported variables at cell centers          !
! (ncelet,*)       !    !     ! at the current and previous time step          !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! coefa, coefb     ! ra ! <-- ! boundary conditions at the boundary faces      !
!  (nfabor,*)      !    !     !                                                !
! statis(ncelet    ! ra ! <-- ! cumulation of the volume statistics            !
!   nvlsta)        !    !     !                                                !
! stativ           ! ra ! <-- ! cumulation for the variances of the            !
!(ncelet,          !    !     ! volume statistics                              !
!   nvlsta-1)      !    !     !                                                !
! tracel(ncelet    ! ra ! <-- ! real array, values cells post                  !
! rdevel(nrdeve    ! ra ! <-- ! dev. complementary array of reals              !
! rtuser(nrtuse    ! ra ! <-- ! user complementary array of reals              !
!                  !    !     !                                                !
! ra(*)            ! ra ! --- ! macro array of reals                           !
!__________________!____!_____!________________________________________________!


!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!==============================================================================


implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "pointe.f90"
include "numvar.f90"
include "optcal.f90"
include "cstphy.f90"
include "entsor.f90"
include "cstnum.f90"
include "lagpar.f90"
include "lagran.f90"
include "parall.f90"
include "period.f90"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , nvlsta
integer          nideve , nrdeve , nituse , nrtuse
integer          ivarl , ivarl1 , ivarlm , iflu , ilpd1 , icla

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision tracel(ncelet)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia , idebra , iel
double precision aa

!===============================================================================
!===============================================================================
! 0. By default, we consider that the subroutine below fits the user's needs;
!    which means running the Lagrangian module triggers the production of
!    standard statistics.
!
!    The user does not need to modify this subroutine in standard-use conditions.
!    In the case where he wishes to produce non-standard statistics, he must
!    intervene in section number 2.
!
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1 . Zone of standard statistics
!===============================================================================

!--> General case:
!     Component X of the particle velocity: ivarl=ilvx
!     Component Y of the particle velocity: ivarl=ilvy
!     Component Z of the particle velocity: ivarl=ilvz
!     Particle temperature: ivarl=iltp
!     Particle diameter: ivarl=ildp
!     Particle mass: ivarl= ilmp
!     Temperature of the coal particles: ivarl=ilhp
!     Mass of reactive coal of the coal particles: ivarl= ilmch
!     Mass of coke of the coal particles: ivarl=ilmck
!     Diameter of the shrinking core of the coal particles: ivarl=ilmck
!    except volume fraction (ivarl=ilfv) and sum of the statistical weights (ivarl=ilpd)
!


if (ivarl.ne.ilfv .and. ivarl.ne.ilpd) then


!-----> Average


  if (iflu.eq.0) then

    do iel = 1, ncel
      if (statis(iel,ilpd1).gt.seuil ) then
        tracel(iel) = statis(iel,ivarl1) / statis(iel,ilpd1)
      else
        tracel(iel) = zero
      endif
    enddo

!-----> Variance

  else

    do iel = 1, ncel
      if ( statis(iel,ilpd1).gt.seuil ) then
        aa = statis(iel,ivarlm)/statis(iel,ilpd1)
        tracel(iel) =  stativ(iel,ivarl1)/statis(iel,ilpd1)       &
                    -( aa * aa )
        tracel(iel) = sqrt( max(zero,tracel(iel)))
      else
        tracel(iel) = zero
      endif
    enddo

  endif

!--> Volume fraction (ilfv)

else if (ivarl.eq.ilfv) then

!-----> Average

  if (iflu.eq.0) then

    do iel = 1, ncel
      if (statis(iel,ilpd1).gt.seuil .and. npst.gt.0) then
        tracel(iel) = statis(iel,ilfv)                            &
                      / (dble(npst) * volume(iel))
      else if (statis(iel,ilpd1).gt.seuil .and.                   &
                iplas.ge.idstnt                  ) then
        tracel(iel) = statis(iel,ilfv) / volume(iel)
      else
        tracel(iel) = zero
      endif
    enddo

  else

!-----> Variance

    do iel = 1, ncel

      if (statis(iel,ilpd1).gt.seuil .and. npst.gt.0) then

        aa = statis(iel,ivarlm) / (dble(npst) * volume(iel))
        tracel(iel) = stativ(iel,ivarl1)                          &
                 / ( dble(npst) * volume(iel) * volume(iel))      &
                 - aa*aa
        tracel(iel) = sqrt( max(zero,tracel(iel)) )

      else if ( statis(iel,ilpd1).gt.seuil .and.                  &
                iplas.ge.idstnt                  ) then

        aa =  statis(iel,ivarlm) / volume(iel)
        tracel(iel) = stativ(iel,ivarl1) / volume(iel)            &
                         - aa*aa
        tracel(iel) = sqrt( max(zero,tracel(iel)))

      else
        tracel(iel) = zero
      endif

    enddo
  endif

!--> Sum of the statistical weights

else if (ivarl.eq.ilpd) then

  if (iflu .eq.0) then
    do iel = 1, ncel
      tracel(iel) = statis(iel,ivarl1)
    enddo
  else
    write(nfecra,9000) iflu
    do iel = 1, ncel
      tracel(iel) = zero
    enddo
  endif

endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION: ERROR IN THE LAGRANGIAN MODULE (uslaen)        ',/,&
'@    =========                                               ',/,&
'@  IT IS NOT POSSIBLE TO COMPUTE THE VARIANCE OF THE         ',/,&
'@     STATISTICAL WEIGHTS                                    ',/,&
'@                                                            ',/,&
'@  The variance of the statistical weight has been asked     ',/,&
'@    in uslaen (ivarl=',   I10,' et iflu=',  I10,').         ',/,&
'@                                                            ',/,&
'@  The calling of the subroutine uslaen must be checked      ',/, &
'@                                                            ',/,&
'@  The calculation continues.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 2. Zone of user intervention
!===============================================================================

!    --------------------------------------------------
!    Example 1: Statistic calculated in uslast.f90 and
!               stored in the array statis
!    --------------------------------------------------

  if (nvlsts.gt.0) then

    if (ivarl.eq.ilvu(1)) then

!-----> Average for the mass concentration

      if (iflu.eq.0) then

        do iel = 1, ncel
          if ( statis(iel,ilpd1).gt.seuil .and. npst.gt.0 ) then
            tracel(iel) = statis(iel,ivarl1)                      &
                        / ( dble(npst) *ro0(ilphas) *volume(iel) )
          else if ( statis(iel,ilpd1).gt.seuil .and.              &
                  iplas.ge.idstnt                  ) then
            tracel(iel) = statis(iel,ivarl1)                      &
                        / ( ro0(ilphas) *volume(iel) )
          else
            tracel(iel) = zero
          endif
        enddo

      else

!-----> Variance of the mass concentration

        do iel = 1, ncel
          if (statis(iel,ilpd1).gt.seuil) then
            aa = statis(iel,ivarlm)/statis(iel,ilpd1)
            tracel(iel) = stativ(iel,ivarl1)/statis(iel,ilpd1)    &
                        -( aa * aa)
            tracel(iel) = sqrt( max(zero,tracel(iel)))
          else
            tracel(iel) = zero
          endif

        enddo

      endif

    endif

  endif

!===============================================================================
! End
!===============================================================================
end subroutine
