!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
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

subroutine ustsv2 &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , produc , gphigk ,          &
   crvexp , crvimp )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for the equations of phi and f_bar
!    with the v2-f phi-model turbulence model (ITURB=50)
!
! Usage
! -----
! The routine is called for both phi and f_bar. It is therefore necessary
! to test the value of the variable ivar to separate the treatments of the
! ivar=iphi or ivar=ifb.
!
! The additional source term is decomposed into an explicit part (crvexp) and
! an implicit part (crvimp) that must be provided here.
! The resulting equation solved by the code are as follows:
!
! For f_bar (ivar=ifb):
!  volume*div(grad(f_bar))= ( volume*f_bar + ..... + crvimp*f_bar + crvexp )/L^2
!
! For phi (ivar=iphi)
!  rho*volume*d(phi)/dt + .... = crvimp*phi + crvexp

!
! Note that crvexp and crvimp are defined after the Finite Volume integration
! over the cells, so they include the "volume" term. More precisely:
!   - crvexp is expressed in m3/s for f_bar
!   - crvexp is expressed in kg/s for phi
!   - crvimp is expressed in m3 for f_bar
!   - crvimp is expressed in kg/s for phi
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

! When entering the routine, two additional work arrays are already set for
! potential user need:
!   produc = 2*mu_t*Sij*Sij -2/3*rho*k*div(u) -2/3*mu_t*div(u)**2
!          + gravity term (when activated)
!          (produc is the production term in the equation for k)
!
!   gphigk = grad(phi).grad(k)

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
! smacel           ! ra ! <-- ! value associated to each variable in the mass  !
!  (ncesmp,nvar)   !    !     !  source terms or mass rate (see ustsma)        !
! produc(ncelet)   ! ra ! <-- ! production term for k                          !
! gphigk(ncelet)   ! ra ! <-- ! grad(phi).grad(k)                              !
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
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
double precision produc(ncelet), gphigk(ncelet)

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

if(iwarni(ifb).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. Example of arbitrary additional source term for k and epsilon

!      Source term for f_bar :
!         volume div(grad f_barre) = (... + volume*xx )/L^2

!      Source term for phi :
!         rho volume d(phi)/dt       = ...
!                        ... + rho*volume*ff*f_bar - rho*volume*phi/tau

!      With xx = 2.d0, ff=3.d0 and tau = 4.d0

!===============================================================================

! ---  For f_bar
!      ---------

if(ivar.eq.ifb) then

  xx  = 2.d0

!    -- Explicit source term

  do iel = 1, ncel
    crvexp(iel) =  volume(iel)*xx
  enddo

!    -- Implicit source term
!
!       crvimp is already initialized to 0, no need to set it here


! ---  For phi
!      -------

elseif(ivar.eq.iphi) then

  ff  = 3.d0
  tau = 4.d0

!   -- Explicit source term

  do iel = 1, ncel
    crvexp(iel) = propce(iel,ipcrom)*volume(iel)*ff*rtpa(iel,ifb)
  enddo

!    -- Implicit source term

  do iel = 1, ncel
    crvimp(iel) = -propce(iel,ipcrom)*volume(iel)/tau
  enddo


endif

!--------
! Formats
!--------

 1000 format(' User source terms for phi and f_bar',/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
