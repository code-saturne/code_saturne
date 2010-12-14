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

subroutine uslast &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itepa  ,                                                       &
   idevel , ituser , ia     ,                                     &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ettp   , ettpa  , tepa   , taup   , tlag   , tempct ,          &
   statis , stativ ,                                              &
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

! User-defined modifications on the variables at the end of the
! Lagrangian iteration and calculation of user-defined
! additional statistics on the particles.
!
! About the user-defined additional statistics, we recall that:
!

!   isttio = 0 : unsteady Lagrangian calculation
!          = 1 : stationary Lagrangian calculation

!   istala : calculation of the statistics if >= 1, else no stats

!   isuist : Restart of statistics calculation if >= 1, else no stats

!   idstnt : Number of the time step for the start of the statistics calculation

!   nstist : Number of the Lagrangian iteration of the start of the stationary computation

!   npst   : Number of iterations of the computation of the stationary statistics

!   npstt  : Total number of iterations of the statistics calculation since the
!            beginning of the calculation, including the unsteady part

!   tstat  : Physical time of the recording of the stationary volume statistics
!            (for the unsteady part, tstat = dtp the Lagrangian time step)
!

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
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
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! idevel(nideve    ! ia ! <-- ! complementary dev. array of integers           !
! ituser(nituse    ! ia ! <-- ! complementary user array of integers           !
! ia(*)            ! ia ! --- ! macro array of integers                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! transported variables at cell centers at       !
! (ncelet,*)       !    !     ! the current and previous time step             !
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
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the previous time step        !
! tepa(nbpmax,     ! ra ! <-- ! properties of the particles (weight..)         !
!       nvep)      !    !     !                                                !
! taup(nbpmax)     ! ra ! <-- ! particle relaxation time                       !
! tlag(nbpmax)     ! ra ! <-- ! relaxation time for the flow                   !
! tempct           ! ra ! <-- ! thermal relaxation time                        !
!  (nbpmax,2)      !    !     !                                                !
! statis           ! ra ! <-- ! cumul. for the averages of the volume stats.   !
!(ncelet,nvlsta    !    !     !                                                !
! stativ           ! ra ! <-- ! cumul. for the variance of the volume stats.   !
!(ncelet,          !    !     !                                                !
!   nvlsta-1)      !    !     !                                                !
! w1..w3(ncelet    ! ra ! --- ! work arrays                                    !
! rdevel(nrdeve    ! ra ! <-- ! dev. complementary array of reals              !
! rtuser(nrtuse    ! ra ! <-- ! user complementary array of reals              !
! ra(*)            ! ra ! --- ! macro array of reals                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstnum
use optcal
use pointe
use entsor
use lagpar
use lagran
use cstphy
use ppppar
use ppthch
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse

integer          itepa(nbpmax,nivep)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision taup(nbpmax) , tlag(nbpmax,3) , tempct(nbpmax,2)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia , idebra
integer          ifinia, ifinra
integer          npt ,  iel , iphas

integer          ivf , ivff , itabvr , iflu , icla

! User-defined local variables

integer          nxlist
parameter       (nxlist=100)

integer          iplan
integer          ii, ind, il
integer          inoeud, irang0, indic
integer          ist(6)

double precision zz(4), zzz(8), xlist(nxlist,8), xyzpt(3)

character        name(8)*4

double precision debm(4)
save             debm

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START

if(istala.eq.1 .and. iplas.ge.idstnt .and. nvlsts.gt.0) then

!
! if needed, the user must fill and adapt the following example
!

  if(1.eq.1) then
    write(nfecra,9000)nvlsts
    call csexit (1)
  endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ CAUTION: STOP IN THE LAGRANGIAN MODULE                  ',/,&
'@    =========                                               ',/,&
'@    THE USER SUBROUTINER uslast MUST BE MODIFIED            ',/,&
'@                                                            ',/,&
'@  The calculation will not be run                           ',/,&
'@                                                            ',/,&
'@  Additional statistics variables have been asked           ',/,&
'@   in uslag1 (nvlsts =',   I10,')                           ',/,&
'@  The subroutine uslast must be adapted to                  ',/, &
'@  precise the computation of their cumulation.              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

else

! During a Lagrangian calculation, we always enter this subroutine
! if we wish to do nothing, we exit immediately
!
  return

endif

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 0.  Memory management
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. Initialization
!===============================================================================

iphas = ilphas

!===============================================================================
! 2 - Computation of user-defined particle statistics
!===============================================================================

!   From a general point of view, we carry out in this subroutine the cumulations of
!   the variables about which we wish to perform statistics. The mean and the
!   variance are calculated in the routine uslaen. This computation is most often
!   carried out by dividing the cumulations by either the stationary cumulation time
!   in the variable tstat, either by the number of particles in statistical weight.
!   This division is applied in each writing in the listing and in
!   the post-processing files.

if (1.eq.0) then

 if(istala.eq.1 .and. iplas.ge.idstnt .and. nvlsts.gt.0) then

  do npt = 1,nbpart

    if( itepa(npt,jisor).gt.0 ) then

      iel = itepa(npt,jisor)

! -------------------------------------------------
! EXAMPLE 1: Cumulation for mass concentration
! -------------------------------------------------

      statis(iel,ilvu(1)) = statis(iel,ilvu(1))                   &
        + tepa(npt,jrpoi) *ettp(npt,jmp)

      stativ(iel,ilvu(1)) = stativ(iel,ilvu(1))                   &
        + tepa(npt,jrpoi) *ettp(npt,jmp) *ettp(npt,jmp)

    endif

  enddo

 endif

endif

!===============================================================================
! 3 - User-defined computation of the particle mass flow rate on 4 plans
!===============================================================================

!  This example is unactivated and must be adapted to the case

if (1.eq.0) then

  zz(1) = 0.1d0
  zz(2) = 0.15d0
  zz(3) = 0.20d0
  zz(4) = 0.25d0

! If we are in an unsteady case, or if the beginning of the stationary stats
! is not reached yet, all statistics are reset to zero at each time step before entering
! this subroutine.

  if(isttio.eq.0 .or. npstt.le.nstist) then
    do iplan = 1,4
      debm(iplan) = 0.d0
    enddo
  endif

  do iplan = 1,4

    do npt = 1,nbpart

      if(itepa(npt,jisor).gt.0) then

        iel = itepa(npt,jisor)

        if( ettp(npt,jxp).gt.zz(iplan) .and.                      &
            ettpa(npt,jxp).le.zz(iplan)      ) then
          debm(iplan) = debm(iplan) +tepa(npt,jrpoi)*ettp(npt,jmp)
        endif

      endif

    enddo
  enddo

  do iplan = 1,4
    write(nfecra,1001)iplan,debm(iplan)/tstat
  enddo

 1001   format(' Debit massique particulaire en Z(',I10,') : ',E14.5)

endif


!===============================================================================
! 4 - Extraction of volume statistics at the end of the calculation
!===============================================================================

!  This example is unactivated and must be adapted to the case

if (1.eq.0) then

  if(ntcabs.eq.ntmabs) then

    zzz(1) = 0.005d0
    zzz(2) = 0.025d0
    zzz(3) = 0.050d0
    zzz(4) = 0.075d0
    zzz(5) = 0.100d0
    zzz(6) = 0.150d0
    zzz(7) = 0.200d0
    zzz(8) = 0.250d0

    NAME(1) = 'XB01'
    NAME(2) = 'XB05'
    NAME(3) = 'XB10'
    NAME(4) = 'XB15'
    NAME(5) = 'XB20'
    NAME(6) = 'XB30'
    NAME(7) = 'XB40'
    NAME(8) = 'XB50'

    ist(1) = ilvx
    ist(2) = ilvz
    ist(3) = ilfv
    ist(4) = ilpd

    npts = nxlist

    ifinia = idebia
    itabvr = idebra
    ifinra = itabvr + ncelet
    CALL RASIZE('USLAST',IFINRA)
    !==========

    do iplan = 1,8

!  Concerning the following file:
!  the user will check if he has not let the unit
!  impusr(1) opened in another user subroutine.
!
      OPEN(FILE=NAME(IPLAN),UNIT=IMPUSR(1),FORM='formatted')

      xyzpt(1) = zzz(iplan)

      do ivf = 1,4

        ivff = ist(ivf)
        icla = 0
        iflu = 0

        call uslaen                                               &
        !==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivff   , ivff   , ivff   , iflu   , ilpd   , icla   ,          &
   idevel , ituser , ia     ,                                     &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis , stativ , ra(itabvr) ,               &
   rdevel , rtuser , ra     )

        ind = 0
        do ii = 1, npts

          xyzpt(2) = 0.d0
          xyzpt(3) = float(ii-1)/float(npts-1)*150.d-3

          call findpt                                             &
          !==========
          (ncelet, ncel, xyzcen,                                  &
           xyzpt(1), xyzpt(2), xyzpt(3), inoeud, irang0)

          indic = ituser(inoeud)
          ituser(inoeud) = 1
          if (indic.eq.1) then
            ind = ind +1
            xlist(ind,1) = xyzcen(1,inoeud)
            xlist(ind,2) = xyzcen(3,inoeud) * (1.d3 / 5.d0)
            xlist(ind,ivf+2) = ra(itabvr+inoeud-1)
          endif
        enddo
      enddo

      do il = 1, ind
        WRITE (IMPUSR(1),'(8E13.5)') (XLIST(IL,II), II=1,6)
      enddo

      close(impusr(1))

    enddo

  endif

endif



!===============================================================================

!====
! End
!====

return

end subroutine
