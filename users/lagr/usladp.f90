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

subroutine usladp &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itepa  ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   ettp   , ettpa  , tepa   , statis ,                            &
   taup   , tlag   , piil   ,                                     &
   vagaus , gradpr , gradvf ,                                     &
   romp   ,                                                       &
   dppar  , dnxpar , dnypar , dnzpar ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module :
!   -------------------------------------

!    User subroutine (non-mandatory intervention)

!    For each particle :
!      - Computation of the wall-normal distance
!      - Computation of the normal to the wall

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
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep    !    !     !                                                !
! ibord            ! ia ! --> ! if nordre=2, contains the number               !
!   (nbpmax)       !    !     ! of the boundary face of part./wall interaction !
! idevel(nideve    ! ia ! <-- ! complementary dev. array of integers           !
! ituser(nituse    ! ia ! <-- ! complementary user array of integers           !
! ia(*)            ! ia ! <-- ! macro array of integers                        !
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
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! transported variables at the current           !
! (ncelet,*)       !    !     ! and previous time step                         !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the current time step         !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     ! the particles at the previous time step        !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
! (nbpmax,nvep)    !    !     !                                                !
! statis           ! ra ! <-- ! cumul for the averages of the volume stats.    !
!(ncelet,nvlsta    !    !     !                                                !
! taup(nbpmax)     ! ra ! <-- ! particle relaxation time                       !
! tlag(nbpmax)     ! ra ! <-- ! relaxation time for the flow                   !
! piil(nbpmax,3    ! ra ! <-- ! term in the sede integration                   !
! tsup(nbpmax,3    ! ra ! <-- ! prediction 1st substep                         !
!                  !    !     ! for the particle velocity                      !
! tsuf(nbpmax,3    ! ra ! <-- ! prediction 1st substep                         !
!                  !    !     ! for the velocity of the flow seen              !
! bx(nbpmax,3,2    ! ra ! <-- ! characteristics of the turbulence              !
! tsfext(nbpmax    ! ra ! <-- ! infos for the return coupling                  !
! vagaus           ! ra ! <-- ! Gaussian random variables                      !
!(nbpmax,nvgaus    !    !     !                                                !
! gradpr(ncel,3    ! ra ! <-- ! pressure gradient                              !
! gradvf(ncel,3    ! ra ! <-- ! flow-velocity gradient                         !
! romp             ! ra ! --- ! particle densite                               !
! fextla           ! ra ! --> ! exterior force field                           !
!(ncelet,3)        !    !     !                                                !
! rdevel(nrdeve    ! ra ! <-- ! dev. complementary array of reals              !
! rtuser(nrtuse    ! ra ! <-- ! user complementary array of reals              !
! ra(*)            ! ra ! --- ! macro array of reals                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array

!===============================================================================


implicit none

!===============================================================================
!     Common blocks
!===============================================================================

include "paramx.f90"
include "numvar.f90"
include "cstnum.f90"
include "cstphy.f90"
include "optcal.f90"
include "entsor.f90"
include "lagpar.f90"
include "lagran.f90"
include "ppppar.f90"
include "ppthch.f90"
include "ppincl.f90"
include "cpincl.f90"

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
integer          itepa(nbpmax,nivep)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp),ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3)
double precision vagaus(nbpmax,*)
double precision dppar(nbpart)  , dnxpar(nbpart)
double precision dnypar(nbpart) , dnzpar(nbpart)

double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision romp(nbpmax)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ip

! User local variables

double precision xnorm

!===============================================================================

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END
!===============================================================================
! 0.  Memory management
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. Example
!===============================================================================

!   This example is unactivated

if (1.eq.0) then

  do ip = 1,nbpart

! Example : for every paricle, we take as
!      - Wall-normal distance: the radius of the duct
!      - Normal: with respect to the radius, by supposing
!                the duct along the Z-axis. Be careful, the
!                convention of Code_Saturne is respected, with
!                the normal oriented from the fluid towards the outside of the domain


    dppar(ip)  = 0.00631942286d0                                  &
                -sqrt( ettp(ip,jxp)*ettp(ip,jxp)                  &
                      +ettp(ip,jyp)*ettp(ip,jyp) )

    xnorm = sqrt( ettp(ip,jxp)*ettp(ip,jxp)                       &
                 +ettp(ip,jyp)*ettp(ip,jyp) )
    dnxpar(ip) = ettp(ip,jxp)/xnorm
    dnypar(ip) = ettp(ip,jyp)/xnorm
    dnzpar(ip) = 0.d0

  enddo

!==============================================================================
! Control: do not modify
!==============================================================================

  do ip = 1,nbpart

    if ( dppar(ip) .le. dparmn ) then
      dppar(ip) = dparmn-dparmn/100.d0
    endif

  enddo

endif

!--------
! Formats
!--------


!----
! End
!----

end subroutine
