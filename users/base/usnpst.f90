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

subroutine usnpst &
!================

 ( nvar   , nscal  , nvlsta ,                                     &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Force or inhibit post-processing output at the current time step.

! This subroutine is called at the end of each time step.

! We pass all the usual arguments to this routine to allow writing of complex
! tests if necessary (for example, output when a given variable reaches a
! given threshold).

! We may also use the following variables from optcal.h:

! ntpabs <-- Absolute number of the last time step of the previous calculation
!            in case of restart (0 otherwise)
! ntcabs <-- Absolute number of the current time step
! ntmabs <-- Absolute number of the last desired time step
! ttpabs <-- Absolute time of at the end of the previous calculation
!            in case of restart (0 otherwise)
! ttcabs <-- Absolute time at the current time step

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nvlsta           ! i  ! <-- ! number of Lagrangian statistical variables     !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
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
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal  , nvlsta

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision statis(ncelet,nvlsta)

! Local variables

integer          indwri, indact

!===============================================================================
! Activate or inhibit output
!===============================================================================

! For any given writer, default output is activated for time steps which are
! multiples of the writer's output frequency (likewise for the absolute
! time)

call pstntc(ntcabs, ttcabs)
!==========

! We may force the activation or deactivation of a given writer at a
! given time step:

! indwri = 0 to activate all writers, or a writer number for a specific writer
! indact = 1 to activate for the current time step, 0 to deactivate.

! By default, all writers are active at the last time step:

if (ntcabs .eq. ntmabs) then
  indwri = 0
  indact = 1
  call pstact(indwri, indact)
  !==========
endif

return
end subroutine
