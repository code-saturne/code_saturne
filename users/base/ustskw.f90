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

subroutine ustskw &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , s2kw   , divukw ,          &
   gkgw   , ggrho  , xf1    ,                                     &
   crkexp , crwexp , crkimp , crwimp ,                            &
   viscf  , viscb  , xam    ,                                     &
   ra     )

!===============================================================================
!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for k and omega equations
!    when using k-omega SST (ITURB=60)
!
! Usage
! -----
!
! The additional source term is decomposed into an explicit part (crkexp,crwexp) and
! an implicit part (crkimp,crwimp) that must be provided here.
! The resulting equations solved by the code are:
!
!  rho*volume*d(k)/dt     + .... = crkimp*k     + crkexp

!  rho*volume*d(omega)/dt + .... = crwimp*omega + crwexp
!
! Note that crkexp, crkimp, crwexp and crwimp are defined after the Finite Volume
! integration over the cells, so they include the "volume" term. More precisely:
!   - crkexp is expressed in kg.m2/s3
!   - crwexp is expressed in kg/s2
!   - crkimp is expressed in kg/s
!   - crwimp is expressed in kg/s
!
! The crkexp, crkimp, crwexp and crwimp arrays are already initialized to 0 before
! entering the routine. It is not needed to do it in the routine (waste of CPU time).
!
! For stability reasons, Code_Saturne will not add -crkimp directly to the
! diagonal of the matrix, but Max(-crkimp,0). This way, the crkimp term is
! treated implicitely only if it strengthens the diagonal of the matrix.
! However, when using the second-order in time scheme, this limitation cannot
! be done anymore and -crkimp is added directly. The user should therefore test
! the negativity of crkimp by himself.
! The same mechanism applies to cvwimp.
!
! When using the second-order in time scheme, one should supply:
!   - crkexp and crwexp at time n
!   - crkimp and crwimp at time n+1/2
!
! When entering the routine, some additional work arrays are already set for
! potential user need:
!   s2kw   =  2 (S11)**2 + 2 (S22)**2 + 2 (S33)**2
!          +  (2 S12)**2 + (2 S13)**2 + (2 S23)**2
!
!            where Sij = (dUi/dxj+dUj/dxi)/2
!
!   divukw = du/dx + dv/dy + dw/dz

!   gkgw = grad(k).grad(omega)

!   ggrho = -g.grad(rho)/prdtur/rho if igrake>0
!         = 0                       if igrake=0
!            where prdtur is the turbulent Prandtl number

!   xf1   = k-eps/k-omega blending coefficient

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
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see ustsma)                                  !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time steps)                        !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc(ncepdp,6) ! ra ! <-- ! head loss coefficient                          !
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate (see ustsma)        !
! s2kw             ! ra ! <-- ! turbulent production term (see comment above)  !
! divukw           ! ra ! <-- ! velocity divergence (see comment above)        !
! gkgw             ! ra ! <-- ! grad(k).grad(omega)                            !
! ggrho            ! ra ! <-- ! gravity source term (see comment above)        !
! xf1              ! ra ! <-- ! k-eps/k-w blending function (see comment above)!
! crkexp           ! ra ! --> ! explicit part of the source term for k         !
! crwexp           ! ra ! --> ! explicit part of the source term for omega     !
! crkimp           ! ra ! --> ! implicit part of the source term for k         !
! crwimp           ! ra ! --> ! implicit part of the source term for omega     !
! viscf(nfac)      ! ra ! --- ! work array                                     !
! viscb(nfabor)    ! ra ! --- ! work array                                     !
! xam(nfac,2)      ! ra ! --- ! work array                                     !
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

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision s2kw(ncelet)  , divukw(ncelet)
double precision gkgw(ncelet)  , ggrho(ncelet), xf1(ncelet)
double precision crkexp(ncelet), crkimp(ncelet)
double precision crwexp(ncelet), crwimp(ncelet)
double precision viscf(nfac), viscb(nfabor), xam(nfac,2)
double precision ra(*)

! Local variables

integer          iel, ikiph, iomgip, ipcrom
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


! --- Index numbers of variables k and omega
ikiph  = ik
iomgip = iomg

! --- Index number of the density in the propce array
ipcrom = ipproc(irom)

if(iwarni(ikiph).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Example of arbitrary additional source term for k and omega

!      Source term for k :
!         rho volume d(k)/dt       = ...
!                        ... - rho*volume*ff*omega - rho*volume*k/tau

!      Source term for omega :
!         rho volume d(omega)/dt = ...
!                        ... + rho*volume*xx

!      With xx = 2.d0, ff=3.d0 and tau = 4.d0

!===============================================================================

! --- Explicit source terms

ff  = 3.d0
tau = 4.d0
xx  = 2.d0

do iel = 1, ncel
  crkexp(iel) = -propce(iel,ipcrom)*volume(iel)*ff*rtpa(iel,iomgip)
  crwexp(iel) =  propce(iel,ipcrom)*volume(iel)*xx
enddo

! --- Implicit source terms
!        crwimp is already initialized to 0, no need to set it here

do iel = 1, ncel
  crkimp(iel) = -propce(iel,ipcrom)*volume(iel)/tau
enddo

!--------
! Formats
!--------

 1000 format(' User source terms for k and omega',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
