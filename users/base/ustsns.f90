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

subroutine ustsns &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   crvexp , crvimp ,                                              &
   ra     )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for velocity components equation
!    (Navier-Stokes)

!
! Usage
! -----
! The routine is called for each velocity component. It is therefore necessary
! to test the value of the variable ivar to separate the treatments of the
! components iu, iv or iw.
!
! The additional source term is decomposed into an explicit part (crvexp) and
! an implicit part (crvimp) that must be provided here.
! The resulting equation solved by the code for a velocity component u is:
!
!  rho*volume*du/dt + .... = crvimp*u + crvexp
!
! Note that crvexp and crvimp are defined after the Finite Volume integration
! over the cells, so they include the "volume" term. More precisely:
!   - crvexp is expressed in kg.m/s2
!   - crvimp is expressed in kg/s
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
! crvexp           ! ra ! --> ! explicit part of the source term               !
! crvimp           ! ra ! --> ! implicit part of the source term               !
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
integer          ia(*)

double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
integer          iel, ipcrom, ipp, iutile
double precision ckp, qdm

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


ipp    = ipprtp(ivar)

if(iwarni(ivar).ge.1) then
  chaine = nomvar(ipp)
  write(nfecra,1000) chaine(1:8)
endif

ipcrom = ipproc(irom  )


!===============================================================================
! 2. Example of arbitrary source term for component u:

!                             S = A * u + B

!            appearing in the equation under the form:

!                       rho*du/dt = S (+ standard Navier-Stokes terms)


!In the following example:
!  A = -rho*CKP
!  B =  XMMT
!
!with:
!  CKP = 1.D0 [1/s       ] (return term on velocity)
!  MMT = 100.D0 [kg/m2/s2] (momentum production by volume and time unit)
!
!which yields:
!     crvimp(iel) = volume(iel)* A = - volume(iel)*(rho*CKP )
!     crvexp(iel) = volume(iel)* B =   volume(iel)*(XMMT    )

! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------

if (ivar.eq.iu) then

  ckp  = 10.d0
  qdm  = 100.d0

  do iel = 1, ncel
     crvimp(iel) = - volume(iel)*propce(iel,ipcrom)*ckp
  enddo

  do iel = 1, ncel
     crvexp(iel) =   volume(iel)*qdm
  enddo

endif


!--------
! Formats
!--------

 1000 format(' User source termes for variable ',A8,/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
