!-------------------------------------------------------------------------------

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

subroutine atprke &
!================

 ( idbia0 , idbra0 ,                                              &
   nscal  , nphas  ,                                              &
   iphas  , ipcvto,                                               &
   ia     ,                                                       &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3    ,                                      &
   w4     , w5     , w6    ,                                      &
   tinstk , tinste ,                                              &
   ra     )


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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! iphas            ! i  ! <-- ! phase number                                   !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! irespr(ncelet    ! te ! --- ! tab entier multigrille                         !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !               !    faces de bord
! coefa, coefb     ! tr !  <- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !               !    faces de bord
! w1...6(ncelet    ! tr ! --- ! tableaux de travail                            !
!tinstk(ncelet)    ! tr ! <-- ! prod et terme de gravite pour eq k             !
!tinste(ncelet)    ! tr ! <-- ! prod et terme de gravite pour eq eps           !
! ra(*)            ! ra ! --- ! main real work array                           !
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

integer          idbia0 , idbra0
integer          nscal  , nphas
integer          iphas  , ipcvto


integer          ia(*)

double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rtp (ncelet,*), rtpa (ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision tinstk(ncelet), tinste(ncelet)
double precision ra(*)

! Local variables
integer         idebra, idebia
integer         iel
integer         itpp , icltpp
integer         iccocg, inc
integer         iivar, iphydp
integer         nswrgp, imligp
integer         iwarnp

double precision gravke, prdtur
double precision epsrgp, climgp, extrap


!
!===============================================================================
!
!===============================================================================
! 1. Initialisation
!===============================================================================

idebia = idbia0
idebra = idbra0

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

  iphydp = 0
  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   nphas  ,                                                       &
   iivar  , imrgra , inc    , iccocg , nswrgp ,imligp, iphydp,    &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w1     , w1     , w1     ,                                     &
   rtpa(1,itpp), coefa(1,icltpp) , coefb(1,icltpp) ,              &
   w4     , w5     , w6     ,                                     &
!        ------   ------   ------
   w1     , w2     , w3     ,                                     &
   ra     )


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
      gravke =   (w4(iel)*gx + w5(iel)*gy + w6(iel)*gz) &
               / (rtp(iel,itpp)*prdtur)
      tinste(iel) = tinstk(iel) + propce(iel,ipcvto)*max(gravke,zero)
      tinstk(iel) = tinstk(iel) + propce(iel,ipcvto)*gravke
    enddo
  else
    do iel = 1, ncel
      gravke =   (w4(iel)*gx + w5(iel)*gy + w6(iel)*gz) &
               / (rtp(iel,itpp)*prdtur)
      tinste(iel) = tinstk(iel) + max(gravke,zero)
      tinstk(iel) = tinstk(iel) + gravke
    enddo
  endif

endif
!----
! FIN
!----
return
end subroutine
