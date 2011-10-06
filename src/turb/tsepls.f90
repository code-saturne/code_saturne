!-------------------------------------------------------------------------------

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

subroutine tsepls &
!================

 ( dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     )

!===============================================================================
! FONCTION :
! ----------

! CALCULATION OF THE E TERM OF THE EPSILON EQUATION (BL-V2/K MODEL)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1(ncelet)       ! ra ! --> ! work array to store the E-term                 !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use optcal
use cstphy
use entsor
use pointe, only: coefau, coefbu
use mesh

!===============================================================================

implicit none

! Arguments

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision w1(ncelet)

! Local variables

integer          iel, ifac, init, inc, iccocg, ivar, iphydp
integer          isou, ii, jj, nswrgp, imligp, iwarnp
double precision climgp, prdtur, extrap
double precision w1f, w2f, w3f, pfac
double precision pondi, flux, somsur, epsrgp

double precision, allocatable, dimension(:) :: w7
double precision, allocatable, dimension(:,:) :: grad

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate a temporary for the gradient calculation
allocate(grad(ncelet,3))

! Allocate work arrays
allocate(w7(ncelet))

!===============================================================================
! 2. CALCULATION OF THE TERM d2Ui/dxkdxj*d2Ui/dxkdxj
!===============================================================================

do iel = 1, ncel
  w1(iel) = 0.0d0
enddo

do isou = 1, 3

  if(isou.eq.1) ivar = iu
  if(isou.eq.2) ivar = iv
  if(isou.eq.3) ivar = iw

  do iel=1,ncel
    w7(iel) = 0.0d0
  enddo

  inc = 1
  iccocg = 1

  nswrgp = nswrgr(ivar)
  epsrgp = epsrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  iphydp = 0

  call grdcel &
  !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,ivar)   , coefa(1,ivar) , coefb(1,ivar) ,               &
   grad   )

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    pondi = pond(ifac)
    w1f = pondi * grad(ii,1) + (1.0d0 - pondi) * grad(jj,1)
    w2f = pondi * grad(ii,2) + (1.0d0 - pondi) * grad(jj,2)
    w3f = pondi * grad(ii,3) + (1.0d0 - pondi) * grad(jj,3)

    somsur = surfac(1,ifac) + surfac(2,ifac) + surfac(3,ifac)

    flux = (w1f + w2f + w3f)*somsur

    w7(ii) = w7(ii) + flux
    w7(jj) = w7(jj) - flux

  enddo

  do ifac = 1, nfabor

    ii = ifabor(ifac)
    w1f = grad(ii,1)
    w2f = grad(ii,2)
    w3f = grad(ii,3)
    somsur = surfbo(1,ifac) + surfbo(2,ifac) + surfbo(3,ifac)
    flux = (w1f + w2f + w3f)*somsur
    w7(ii) = w7(ii) + flux

  enddo

  do iel = 1, ncel
    w1(iel) = w1(iel) + (w7(iel)/volume(iel))**2
  enddo

enddo

! Free memory
deallocate(grad)
deallocate(w7)

!----
! FIN
!----

return

end subroutine
