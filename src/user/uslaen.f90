!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine uslaen &
!================

 ( nvlsta ,                                                       &
   ivarl  , ivarl1 , ivarlm , iflu   , ilpd1  , icla   ,          &
   tracel )

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
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
! ivarl            ! i  ! <-- ! number of the stat (between 1 and nvlsta)      !
! ivarl1           ! i  ! <-- ! number of the global stat + group              !
!                  !    !     ! (average or variance)                          !
! ivarlm           ! i  ! <-- ! number of the stat mean + group                !
! iflu             ! i  ! <-- ! 0: mean of the stat ivarl/ivarl1               !
!                  !    !     ! 1: variance of the stat ivarl/ivarl1           !
! ilpd1            ! i  ! <-- ! "pointer" to global statistical weight         !
!                  !    !     !                                                !
! icla             ! i  ! <-- ! 0: global statistic                            !
!                  !    ! <-- ! !=0: stat for the icla group                   !
! tracel(ncelet)   ! ra ! <-- ! real array, values cells post                  !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!==============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use entsor
use cstnum
use lagpar
use lagran
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvlsta
integer          ivarl , ivarl1 , ivarlm , iflu , ilpd1 , icla

double precision tracel(ncelet)

! Local variables

integer          iel
double precision aa

integer ii, jj, dimtab

integer, allocatable, dimension(:) :: tabstat


!===============================================================================

!===============================================================================
! 0. By default, we consider that the subroutine below fits the user's needs;
!    which means running the Lagrangian module triggers the production of
!    standard statistics.
!
!    The user does not need to modify this subroutine in standard-use conditions.
!    In the case where he wishes to produce non-standard statistics, he must
!    intervene in section number 2.
!===============================================================================

!===============================================================================
! 1 . Zone of standard statistics
!===============================================================================

allocate(tabstat(ncel))

! Pinpoint the cells where stats are to be calculated

ii = 0
do iel = 1, ncel
  if (statis(iel,ilpd1).gt.seuil) then
    ii = ii + 1
    tabstat(ii) = iel
  endif
  tracel(iel) = 0.d0
enddo
dimtab = ii

! General case:
!
!   Component X of the particle velocity: ivarl=ilvx
!   Component Y of the particle velocity: ivarl=ilvy
!   Component Z of the particle velocity: ivarl=ilvz
!   Particle temperature: ivarl=iltp
!   Particle diameter: ivarl=ildp
!   Particle mass: ivarl= ilmp
!   Temperature of the coal particles: ivarl=ilhp(ilayer)
!   Mass of moisture of the coal particles: ivarl= ilmwat
!   Mass of reactive coal of the coal particles: ivarl= ilmch
!   Mass of coke of the coal particles: ivarl=ilmck
!   Diameter of the shrinking core of the coal particles: ivarl=ilmck
!     except volume fraction (ivarl=ilfv) and sum of the statistical weights
!     (ivarl=ilpd)

if (ivarl.ne.ilfv .and. ivarl.ne.ilpd) then

  !-----> Average

  if (iflu.eq.0) then

    do jj = 1, dimtab
      tracel(tabstat(jj)) =    statis(tabstat(jj),ivarl1)                     &
                             / statis(tabstat(jj),ilpd1)
    enddo

  !-----> Variance

  else

    do jj = 1, dimtab
      aa = statis(tabstat(jj),ivarlm)/statis(tabstat(jj),ilpd1)
      tracel(tabstat(jj)) =   stativ(tabstat(jj),ivarl1)                      &
                            / statis(tabstat(jj),ilpd1) - (aa * aa)
    enddo

  endif

  !--> Volume fraction (ilfv)

else if (ivarl.eq.ilfv) then

  !-----> Average

  if (iflu.eq.0) then

    do jj = 1, dimtab
      tracel(tabstat(jj)) =   statis(tabstat(jj),ilfv)                        &
                            / (dble(npst) * volume(tabstat(jj)))
    enddo

  else

  !-----> Variance

    do jj = 1, dimtab

      if (npst.gt.1) then
        aa = statis(tabstat(jj),ivarlm) / (dble(npst) * volume(tabstat(jj)))
        tracel(tabstat(jj)) =   stativ(tabstat(jj),ivarl1)                    &
                              / (dble(npst) * volume(tabstat(jj)))**2         &
                              - aa*aa
      else
        tracel(tabstat(jj)) = zero
      endif

    enddo
  endif

  !--> Sum of the statistical weights

else if (ivarl.eq.ilpd) then

  if (iflu .eq.0) then
    do jj = 1, dimtab
      tracel(tabstat(jj)) = statis(tabstat(jj),ivarl1)
    enddo
  else
    write(nfecra,9000) iflu
    do jj = 1, dimtab
      tracel(tabstat(jj)) = zero
    enddo
  endif

endif

 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION: ERROR IN THE LAGRANGIAN MODULE (uslaen)        ',/,&
'@    =========                                               ',/,&
'@  It is not possible to compute the variance of the         ',/,&
'@     statistical weights                                    ',/,&
'@                                                            ',/,&
'@  The variance of the statistical weight has been asked     ',/,&
'@    in uslaen (ivarl=',   I10,' et iflu=',  I10,').         ',/,&
'@                                                            ',/,&
'@  The call to subroutine uslaen must be checked             ',/,&
'@                                                            ',/,&
'@  The calculation continues.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 2. Zone of user intervention
!===============================================================================

! --------------------------------------------------
! Example 1: Statistic calculated in uslast.f90 and
!            stored in the array statis
! --------------------------------------------------

if (nvlsts.gt.0) then

  if (ivarl.eq.ilvu(1)) then

    !-----> Average for the mass concentration

    if (iflu.eq.0) then

      do jj = 1, dimtab
        if (npst.gt.0) then
          tracel(tabstat(jj)) =   statis(tabstat(jj),ivarl1)                  &
                                / (dble(npst) *ro0 *volume(tabstat(jj)))
        else if (iplas.ge.idstnt) then
          tracel(tabstat(jj)) =   statis(tabstat(jj),ivarl1)                  &
                                / (ro0 *volume(tabstat(jj)))
        else
          tracel(tabstat(jj)) = zero
        endif
      enddo

    else

    !-----> Variance of the mass concentration

      do jj = 1, dimtab

        aa = statis(tabstat(jj),ivarlm)/statis(tabstat(jj),ilpd1)
        tracel(tabstat(jj)) =   stativ(tabstat(jj),ivarl1)                    &
                              / statis(tabstat(jj),ilpd1)                     &
                              - (aa * aa)

      enddo

    endif

  endif

endif

! Free memory

deallocate(tabstat)

!===============================================================================
! End
!===============================================================================

end subroutine uslaen
