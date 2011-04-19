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

subroutine uslafe &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itepa  , ibord  ,                                              &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   ettp   , ettpa  , tepa   , statis , stativ ,                   &
   taup   , tlag   , piil   ,                                     &
   tsuf   , tsup   , bx     , tsfext ,                            &
   vagaus , gradpr , gradvf ,                                     &
   romp   , fextla ,                                              &
   ra     )

!===============================================================================
! Purpose:
! ----------

!   Subroutine of the Lagrangian particle-tracking module:
!   -------------------------------------

!    User subroutine (non-mandatory intervention)

!    Management of an external force field acting on the particles
!    It must be prescribed in every cell and be homogeneous to gravity (m/s^2)
!
!    By default gravity and drag force are the only forces acting on the particles
!    (the gravity components gx gy gz are assigned in usini1)
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
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep)   !    !     !                                                !
! ibord(nbpmax)    ! ia ! --> ! if nordre=2, contains the number of the        !
!                  !    !     ! boundary face of particle/wall interaction     !
! ia(*)            ! ia ! --- ! macro array of integers                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! transported variables at cell centers for      !
! (ncelet,*)       !    !     ! the current and the previous time step         !
! propce           ! ra ! <-- ! physical properties at cell centers            !
! (ncelet,*)       !    !     !                                                !
! propfa           ! ra ! <-- ! physical properties at interior face centers   !
!  (nfac,*)        !    !     !                                                !
! propfb           ! ra ! <-- ! physical properties at boundary face centers   !
!  (nfabor,*)      !    !     !                                                !
! ettp             ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     !                                                !
! ettpa            ! ra ! <-- ! array of the variables associated to           !
!  (nbpmax,nvp)    !    !     !                                                !
! tepa             ! ra ! <-- ! particle information (real) (statis. weight..) !
!  (nbpmax,nvep)   !    !     !                                                !
! statis           ! ra ! <-- ! cumul for the averages of the volume stats.    !
!  (ncelet,nvlsta) !    !     !                                                !
! stativ           ! ra ! <-- ! cumulation for the variance of the volume      !
!  (ncelet,        !    !     ! statistics                                     !
!   nvlsta-1)      !    !     !                                                !
! taup(nbpmax)     ! ra ! <-- ! particle relaxation time                       !
! tlag(nbpmax)     ! ra ! <-- ! relaxation time for the flow                   !
! piil(nbpmax,3)   ! ra ! <-- ! term in the integration of the sde             !
! tsup(nbpmax,3)   ! ra ! <-- ! prediction 1st substep for                     !
!                  !    !     ! the velocity of the particles                  !
! tsuf(nbpmax,3)   ! ra ! <-- ! prediction 1st substep for                     !
!                  !    !     ! the velocity of the flow seen                  !
! bx(nbpmax,3,2)   ! ra ! <-- ! characteristics of the turbulence              !
! tsfext(nbpmax    ! ra ! <-- ! infos for the return coupling                  !
! vagaus           ! ra ! <-- ! Gaussian random variables                      !
!  (nbpmax,nvgaus) !    !     !                                                !
! gradpr(ncel,3)   ! ra ! <-- ! pressure gradient                              !
! gradvf(ncel,3)   ! ra ! <-- ! gradient of the flow velocity                  !
! romp             ! ra ! --- ! particle density                               !
! fextla(ncelet,3) ! ra ! --> ! user external force field (m/s^2)              !
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
use cstphy
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cpincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itepa(nbpmax,nivep) , ibord(nbpmax)
integer          ia(*)

double precision dt(ncelet) , rtp(ncelet,*) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision ettp(nbpmax,nvp) , ettpa(nbpmax,nvp)
double precision tepa(nbpmax,nvep)
double precision statis(ncelet,*),stativ(ncelet,*)
double precision taup(nbpmax) , tlag(nbpmax,3)
double precision piil(nbpmax,3) , bx(nbpmax,3,2)
double precision tsuf(nbpmax,3) , tsup(nbpmax,3)
double precision tsfext(nbpmax)
double precision vagaus(nbpmax,*)
double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision romp(nbpmax)
double precision fextla(nbpmax,3)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ip


!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     If a file from the GUI is used, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================
!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 0. Memory management
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. Example
!===============================================================================

!   This example is unactivated

if (1.eq.0) then


  do ip = 1,nbpart

    fextla(ip,1) = 0.d0
    fextla(ip,2) = 0.d0
    fextla(ip,3) = 0.d0

  enddo


endif

!==============================================================================

!--------
! Formats
!--------


!----
! End
!----

end subroutine
