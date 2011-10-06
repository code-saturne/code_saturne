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

subroutine ustssc &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   crvexp , crvimp )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Additional right-hand side source terms for scalar equations (user
!     scalars and specific physics scalars).

!
! Usage
! -----
! The routine is called for each scalar, user or specific physisc. It is
! therefore necessary to test the value of the scalar number iscal to separate
! the treatments of the different scalars (if (iscal.eq.p) then ....).
!
! The additional source term is decomposed into an explicit part (crvexp) and
! an implicit part (crvimp) that must be provided here.
! The resulting equation solved by the code for a scalar f is:
!
!  rho*volume*df/dt + .... = crvimp*f + crvexp
!
!
! Note that crvexp and crvimp are defined after the Finite Volume integration
! over the cells, so they include the "volume" term. More precisely:
!   - crvexp is expressed in kg.[scal]/s, where [scal] is the unit of the scalar
!   - crvimp is expressed in kg/s
!
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

!
! STEEP SOURCE TERMS
!===================
! In case of a complex, non-linear source term, say F(f), for scalar f, the
! easiest method is to implement the source term explicitely.
!
!   df/dt = .... + F(f(n))
!   where f(n) is the value of f at time tn, the beginning of the time step.
!
! This yields :
!   crvexp = volume*F(f(n))
!   crvimp = 0
!
! However, if the source term is potentially steep, this fully explicit
! method will probably generate instabilities. It is therefore wiser to
! partially implicit the term by writing:
!
!   df/dt = .... + dF/df*f(n+1) - dF/df*f(n) + F(f(n))
!
! This yields:
!   crvexp = volume*( F(f(n)) - dF/df*f(n) )
!   crvimp = volume*dF/df

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss terms           !
! ncssmp           ! i  ! <-- ! number of cells with mass source terms         !
! iscal            ! i  ! <-- ! index number of the current scalar             !
! icepdc(ncepdp)   ! ia ! <-- ! index number of cells with head loss terms     !
! icetsm(ncesmp)   ! ia ! <-- ! index number of cells with mass source terms   !
! itypsm           ! ia ! <-- ! type of mass source term for each variable     !
!  (ncesmp,nvar)   !    !     !  (see ustsma)                                  !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (preceding time step)                         !
! rtp              ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (current time step)                           !
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
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)

! Local variables

character*80     chaine
integer          ivar, iiscvr, ipcrom, iel, iutile
integer          ilelt, nlelt

double precision tauf, prodf, volf, pwatt

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


! --- Index number of the variable associated to scalar iscal
ivar = isca(iscal)

! --- Name of the the variable associated to scalar iscal
chaine = nomvar(ipprtp(ivar))

! --- Indicateur of variance scalars
!         If iscavr(iscal) = 0:
!           the scalar iscal is not a variance
!         If iscavr(iscal) > 0 and iscavr(iscal) < nscal + 1 :
!           the scalar iscal is the variance of the scalar iscavr(iscal)
iiscvr = iscavr(iscal)

! --- Index number of the density in the propce array
ipcrom = ipproc(irom)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif


!===============================================================================
! 2. Example of arbitrary source term for the scalar f, 2nd scalar in the
!    calculation

!                             S = A * f + B

!            appearing in the equation under the form

!                       rho*df/dt = S (+ regular terms in the equation)


!In the following example:
!     A = - rho / tauf
!     B =   rho * prodf
!        AVEC
!     tauf   = 10.d0  [ s  ] (dissipation time for f)
!     prodf  = 100.d0 [ [f]/s ] (production of f by unit of time)

!which yields
!     crvimp(iel) = volume(iel)* A = - volume(iel)*rho/tauf
!     crvexp(iel) = volume(iel)* B =   volume(iel)*rho*prodf

!===============================================================================


! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------

!Source term applied to second scalar
if (iscal.eq.2) then

   tauf  = 10.d0
   prodf = 100.d0

   do iel = 1, ncel
      crvimp(iel) = - volume(iel)*propce(iel,ipcrom)/tauf
   enddo

   do iel = 1, ncel
      crvexp(iel) =   volume(iel)*propce(iel,ipcrom)*prodf
   enddo

endif

!===============================================================================
! 3. Example of arbitrary volumic heat term in the equation for enthalpy h

! In the considered example, a uniform volumic source of heating is imposed
! in the cells with coordinate X in [0;1.2] and Y in [3.1;4]

! The global heating power if Pwatt (in W) and the total volume of the concerned
! cells is volf (in m3)

! This yields
!     crvimp(iel) = 0
!     crvexp(iel) = volume(iel)* Pwatt/volf

!===============================================================================


! ----------------------------------------------

! It is quite frequent to forget to remove this example when it is
!  not needed. Therefore the following test is designed to prevent
!  any bad surprise.

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------

! WARNING :
! It is assumed here that the thermal scalar in an enthalpy.
!  If the scalar in a temperature, PWatt must be devided by Cp.

pwatt = 100.d0

! calculation of volf

volf  = 0.d0
CALL GETCEL('X > 0.0 and X < 1.2 and Y > 3.1 and'//               &
            'Y < 4.0',NLELT,LSTELT)

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
  volf = volf + volume(iel)
enddo

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
! No implicit source term
  crvimp(iel) = 0.d0
! Explicit source term
  crvexp(iel) = volume(iel)*pwatt/volf
enddo

!--------
! Formats
!--------

 1000 format(' User source terms for variable ',A8,/)

!----
! End
!----

! Deallocate the temporary array
deallocate(lstelt)

return
end subroutine
