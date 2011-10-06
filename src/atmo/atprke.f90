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

subroutine atprke &
!================

 ( nscal  ,                                                       &
   ipcvto,                                                        &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   tinstk , tinste )

!===============================================================================
! FONCTION :
! --------
!  SPECIFIQUE AU CAS ATMOSPHERIQUE :
!  CALCUL DU TERME DE PRODUCTION LIEE A LA FLOTTABILITE:
!  G = G*GRAD(THETA)/PRDTUR/THETA
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! irespr(ncelet    ! te ! --- ! tab entier multigrille                         !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !               !    faces de bord
! coefa, coefb     ! tr !  <- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !               !    faces de bord
!tinstk(ncelet)    ! tr ! <-- ! prod et terme de gravite pour eq k             !
!tinste(ncelet)    ! tr ! <-- ! prod et terme de gravite pour eq eps           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use cstnum
use cstphy
use optcal
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nscal
integer          ipcvto



double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rtp (ncelet,*), rtpa (ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision tinstk(ncelet), tinste(ncelet)

! Local variables
integer         iel
integer         itpp , icltpp
integer         iccocg, inc
integer         iivar
integer         nswrgp, imligp
integer         iwarnp

double precision gravke, prdtur
double precision epsrgp, climgp, extrap

double precision, allocatable, dimension(:,:) :: grad

!===============================================================================
!
!===============================================================================
! 1. Initialisation
!===============================================================================

! Allocate work arrays
allocate(grad(ncelet,3))


!===============================================================================
! 2. Calcul des derivees de la temperature potentielle
!===============================================================================

if (ippmod(iatmos).eq.1) then


  itpp = isca(iscalt)
  icltpp = iclrtp(itpp,icoef)

! ---- Options de calcul:

  iccocg = 1
  inc = 1

  nswrgp = nswrgr(itpp)
  epsrgp = epsrgr(itpp)
  imligp = imligr(itpp)
  iwarnp = iwarni(itpp)
  climgp = climgr(itpp)
  extrap = extrag(itpp)

  iivar = itpp

  call grdcel                                                     &
  !==========
 ( iivar  , imrgra , inc    , iccocg , nswrgp ,imligp,            &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,itpp), coefa(1,icltpp) , coefb(1,icltpp) ,              &
   grad   )


!      Production et terme de gravite
!         TINSTK=P+G et TINSTE=P+(1-CE3)*G

  if(iscalt.gt.0.and.nscal.ge.iscalt) then
    prdtur = sigmas(iscalt)
  else
    prdtur = 1.d0
  endif

!     En production lineaire, on multiplie tout de suite le terme
!     de gravite par VISCT, car le reste est deja multiplie.
!     Dans les autres cas, la multiplication est faite plus tard.
  if (iturb.eq.21) then
    do iel = 1, ncel
      gravke =   (grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz) &
               / (rtp(iel,itpp)*prdtur)
      tinste(iel) = tinstk(iel) + propce(iel,ipcvto)*max(gravke,zero)
      tinstk(iel) = tinstk(iel) + propce(iel,ipcvto)*gravke
    enddo
  else
    do iel = 1, ncel
      gravke =   (grad(iel,1)*gx + grad(iel,2)*gy + grad(iel,3)*gz) &
               / (rtp(iel,itpp)*prdtur)
      tinste(iel) = tinstk(iel) + max(gravke,zero)
      tinstk(iel) = tinstk(iel) + gravke
    enddo
  endif

endif

! Allocate work arrays
deallocate(grad)

!----
! FIN
!----
return
end subroutine
