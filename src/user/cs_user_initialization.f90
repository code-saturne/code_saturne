!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine cs_user_initialization &
!================================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfa , propfb , coefa  , coefb  )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Initialize variables

! This subroutine is called at beginning of the computation
! (restart or not) before the loop time step

! This subroutine enables to initialize or modify (for restart)
!     unkown variables and time step values

! rom and viscl values are equal to ro0 and viscl0 or initialize
! by reading the restart file
! viscls and cp variables (when there are defined) have no value
! excepted if they are read from a restart file

! Physical quantities are defined in the following arrays:
!  propce (physical quantities defined at cell center),
!  propfa (physical quantities defined at interior face center),
!  propfa (physical quantities defined at border face center).
!
! Examples:
!  propce(iel, ipproc(irom  )) means rom  (iel)
!  propce(iel, ipproc(iviscl)) means viscl(iel)
!  propce(iel, ipproc(icp   )) means cp   (iel)
!  propce(iel, ipproc(ivisls(iscal))) means visls(iel, iscal)
!  propfa(ifac, ipprof(ifluma(ivar))) means flumas(ifac, ivar)
!  propfb(ifac, ipprob(irom )) means romb  (ifac)
!  propfb(ifac, ipprob(ifluma(ivar))) means flumab(ifac, ivar)

! Modification of the behaviour law of physical quantities (rom, viscl,
! viscls, cp) is not done here. It is the purpose of the user subroutine
! usphyv

! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the
! 'cs_user_boundary_conditions' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp(ncelet, *)   ! ra ! <-- ! computed variables at cell centers at current  !
!                  !    !     ! time steps                                     !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
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
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use atincl
use ctincl
use elincl
use ppcpfu
use cs_coal_incl
use cs_fuel_incl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables

integer          iel, iutile

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) then
!       Tag to know if a call to this subroutine has already been done
  iusini = 0
  return
endif

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1.  Initialization of local variables
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))


!===============================================================================
! 2. Unknown variables initialization:
!      ONLY done if there is no restart computation
!===============================================================================

! --- Example:  isca(1) is the variable number in RTP related to the first
!               user-defined scalar variable
!               rtp(iel,isca(1)) is the value of this variable in cell number
!               iel.

if (isuite.eq.0) then

  do iel = 1, ncel
    rtp(iel,isca(1)) = 25.d0
  enddo

endif


!===============================================================================
! 3. Time step:
!===============================================================================

! --- Example: We do a computation restart with a variable in time and constant
!              in space time step or with a variable in time and space time
!              step. We want to modify the time step given by the reading of
!              the restart file (in order to overcome a too slow evolution
!              for instance)

! The test below allows checking that the following example compiles
! while disabling it by default.

iutile = 0

if (iutile.eq.0) return

! ----------------------------------------------


if (isuite.eq.1.and.(idtvar.eq.1.or.idtvar.eq.2)) then
  do iel = 1, ncel
    dt (iel) = 10.d0*dt(iel)
  enddo
endif

!--------
! Formats
!--------

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
