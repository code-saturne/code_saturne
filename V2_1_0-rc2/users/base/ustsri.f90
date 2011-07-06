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

subroutine ustsri &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smcelp , gamma  , grdvit , produc , &
   crvexp , crvimp )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for the equations of Rij and
!    epsilon in Rij-LRR model (ITURB=30) or Rij-SSG model (ITURB=31)

!
! Usage
! -----
! The routine is called for each variable Rij and epsilon. It is therefore
! necessary to test the value of the variable ivar to separate the treatments
! of the variables ivar=ir11, ir22, ir33, ir12,
! ir13, ir23, or iep.
!
! The additional source term is decomposed into an explicit part (crvexp) and
! an implicit part (crvimp) that must be provided here.
! The resulting equation solved by the code for a variable var is:
!
!  rho*volume*d(var)/dt + .... = crvimp*(var) + crvexp
!
! Note that crvexp and crvimp are defined after the Finite Volume integration
! over the cells, so they include the "volume" term. More precisely:
!   - crvexp is expressed in kg.m2/s3 for Rij
!   - crvexp is expressed in kg.m2/s4 for epsilon
!   - crvimp is expressed in kg/s for Rij and epsilon
!
! The crvexp and crvimp arrays are already initialized to 0 before entering the
! the routine. It is not needed to do it in the routine (waste of CPU time).
!
! For stability reasons, Code_Saturne will not add -crvimp directly to the
! diagonal of the matrix, but Max(-crvimp,0). This way, the crvimp term is
! treated implicitely only if it strengthens the diagonal of the matrix.
! However, when using the second-order in time scheme, this limitation cannot
! be done anymore and -crvimp is added directly. The user should therefore test
! the negativity of crvimp by himself.
!
! When using the second-order in time scheme, one should supply:
!   - crvexp at time n
!   - crvimp at time n+1/2
!
!
! When entering the routine, two additional work arrays are already set for
! potential user need:
!   produc is the turbulent production term
!          produc(6,ncelet) = (P11, P22, P33, P12, P13, P23)
!          with Pij=-Rik.dUj/dxk - Rjk.dUi/dxk
!
!   grdvit is the
!          grdvit(ncelet,i,j) = dUi/dxj


! WARNING
! =======
! produc is only allocated and defined for the Rij-LRR model (ITURB=30)
! grdvit is only allocated and defined for the Rij-SSG model (ITURB=31)
! DO NOT USE produc WITH THE SSG MODEL
! DO NOT USE grdvit WITH THE LRR MODEL


!
! The selection of cells where to apply the source terms is based on a getcel
! command. For more info on the syntax of the getcel command, refer to the
! user manual or to the comments on the similar command getfbr in the routine
! usclim.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! ncssmp           ! i  ! <-- ! number of cells with mass source terms         !
! ivar             ! i  ! <-- ! index number of the current variable           !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see ustsma)                                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time steps)                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smcelp(ncelet)   ! ra ! <-- ! value of variable ivar associated to mass      !
!                  ! ra !     !  source term (see ustsma)                      !
! gamma(ncelet)    ! ra ! <-- ! volumic rate of mass source term               !
! grdvit(ncelet,3,3! ra ! <-- ! velocity gradient (only for iturb=31)          !
! produc(6,ncelet) ! ra ! <-- ! turbulent production term (only for iturb=30)  !
! crvexp           ! ra ! --> ! explicit part of the source term               !
! crvimp           ! ra ! --> ! implicit part of the source term               !
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
use entsor
use optcal
use cstphy
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itpsmp(ncesmp)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6)
double precision smcelp(ncesmp), gamma(ncesmp)
double precision grdvit(ncelet,3,3), produc(6,ncelet)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables

integer          iel, ipcrom
double precision ff, tau, xx

integer, allocatable, dimension(:) :: lstelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate a temporary array for cells selection
allocate(lstelt(ncel))


! --- Index number of the density in the propce array
ipcrom = ipproc(irom)

if(iwarni(ir11).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Example of arbitrary additional source term for R11 and epsilon

!      Source term for R11 :
!         rho volume d(R11)/dt       = ...
!                        ... - rho*volume*ff*epsilon - rho*volume*R11/tau

!      Source term for epsilon :
!         rho volume d(epsilon)/dt = ...
!                        ... + rho*volume*xx

!      With xx = 2.d0, ff=3.d0 and tau = 4.d0

!===============================================================================


! ---  For R11
!      -------

if(ivar.eq.ir11) then

  ff  = 3.d0
  tau = 4.d0

!   -- Explicit source term

  do iel = 1, ncel
    crvexp(iel) = -propce(iel,ipcrom)*volume(iel)                 &
                                               *ff*rtpa(iel,iep)
  enddo

!    -- Implicit source term

  do iel = 1, ncel
    crvimp(iel) = -propce(iel,ipcrom)*volume(iel)/tau
  enddo


! ---  For epsilon
!      -----------

elseif(ivar.eq.iep) then

  xx  = 2.d0

!   -- Explicit source term

  do iel = 1, ncel
    crvexp(iel) =  propce(iel,ipcrom)*volume(iel)*xx
  enddo

!    -- Implicit source term
!        crvimp is already initialized to 0, no need to set it here


endif

!--------
! Formats
!--------

 1000 format(' User source terms for Rij and epsilon',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
