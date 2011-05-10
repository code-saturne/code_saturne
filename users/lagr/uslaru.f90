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
   nvar   , nscal  ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   itypfb , itrifb , itepa  ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , croule , auxl  ,                    &
   distpa , distyp ,                                              &
   w1     , w2     , w3     ,                                     &
   ra     )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nbpmax           ! i  ! <-- ! maximum number of particles allowed            !
! nvp              ! i  ! <-- ! number of particle variables                   !
! nvp1             ! i  ! <-- ! nvp minus position, fluid and part. velocities !
! nvep             ! i  ! <-- ! number of particle properties (integer)        !
! nivep            ! i  ! <-- ! number of particle properties (integer)        !
! ntersl           ! i  ! <-- ! number of source terms of return coupling      !
! nvlsta           ! i  ! <-- ! nb of Lagrangian statistical variables         !
! nvisbr           ! i  ! <-- ! number of boundary statistics                  !
! itypfb(nfabor)   ! ia ! <-- ! type of the boundary faces                     !
! itrifb(nfabor)   ! ia ! --> ! indirection for the sorting of the             !
! itepa            ! ia ! <-- ! particle information (integers)                !
! (nbpmax,nivep)   !    !     !                                                !
! ia(*)            ! ia ! --- ! macro array of integers                        !
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
! ra(*)            ! ra ! --- ! macro array of reals                           !
! distpa(ncelet    ! ra ! <-- ! wall-normal distance arrays                    !
! disty(ncelet)    ! ra ! <-- ! y+ distance                                    !
! w1...w3(ncel)    ! ra ! --- ! work arrays                                    !
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
use optcal
use entsor
use cstphy
use pointe
use parall
use period
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr

integer          itypfb(nfabor) , itrifb(nfabor)
integer          itepa(nbpmax,nivep)
integer          ia(*)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision vagaus(nbpmax,*) , croule(ncelet)
double precision auxl(nbpmax,3)
double precision distpa(ncelet) , distyp(ncelet)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision ra(*)

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
