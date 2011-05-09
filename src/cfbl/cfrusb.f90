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

subroutine cfrusb &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   imodif , iphas  ,                                              &
   ia     ,                                                       &
   gammag ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   sorti1 , sorti2 , gamagr , masmor ,                            &
   ra     )

!===============================================================================
! FONCTION :
! ---------

! Flux de rusanov au bord pour euler + energie

! d rho   /dt + div rho u             = 0
! d rho u /dt + div rho u u + grad  P = 0
! d E     /dt + div rho u E + div u P = 0



!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! imodif           ! e  ! <-- ! modification directe de rtp (imodif=1          !
! iphas            ! i  ! <-- ! phase number                                   !
! ia(*)            ! ia ! --- ! main integer work array                        !
! gammag           ! r  ! <-- ! gamma du gaz                                   !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp,rtpa         ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! sorti1,2(*)      ! tr ! --> ! variables de sortie                            !
! gamagr(*)        ! tr ! --- ! constante gamma equivalent du gaz              !
! masmor(*)        ! tr ! --- ! masse molaire des constituants du gaz          !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use numvar
use optcal
use cstphy
use cstnum
use parall
use pointe
use entsor
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas
integer          imodif , iphas

integer          ia(*)

double precision gammag

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*),propfa(nfac,*),propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision sorti1(*), sorti2(*), gamagr(*), masmor(*)
double precision ra(*)


! Local variables

integer          idebia , idebra
integer          ifac0
integer          iel    , ifac
integer          ipriph , irhiph , ieniph
integer          iuiph  , iviph  , iwiph
integer          iclp   , iclr   , icle
integer          iclu   , iclv   , iclw
integer          iflmab
double precision und    , uni    , rund   , runi   , cd     , ci
double precision rrus   , runb


!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

ipriph = ipr
irhiph = isca(irho  )
ieniph = isca(ienerg)
iuiph  = iu
iviph  = iv
iwiph  = iw
iclp = iclrtp(ipriph,icoef)
iclr = iclrtp(irhiph,icoef)
icle = iclrtp(ieniph,icoef)
iclu = iclrtp(iuiph,icoef)
iclv = iclrtp(iviph,icoef)
iclw = iclrtp(iwiph,icoef)

iflmab = ipprob(ifluma(ieniph))

ifac0 = imodif
ifac  = ifac0
iel   = ifabor(ifac)

!===============================================================================
! 1. GRANDEURS LIEES A RUSANOV
!===============================================================================


und   = (coefa(ifac,iclu)*surfbo(1,ifac)                          &
       + coefa(ifac,iclv)*surfbo(2,ifac)                          &
       + coefa(ifac,iclw)*surfbo(3,ifac))/surfbn(ifac)
uni   = (rtp(iel,iuiph)*surfbo(1,ifac)                            &
       + rtp(iel,iviph)*surfbo(2,ifac)                            &
       + rtp(iel,iwiph)*surfbo(3,ifac))/surfbn(ifac)
rund  = coefa(ifac,iclr)*und
runi  = rtp(iel,irhiph)     *uni
cd    = sqrt(gammag*coefa(ifac,iclp)/coefa(ifac,iclr))
ci    = sqrt(gammag*rtp(iel,ipriph)/rtp(iel,irhiph))
rrus  = max(abs(und)+cd,abs(uni)+ci)

runb  = 0.5d0*(coefa(ifac,iclr)*und+rtp(iel,irhiph)*uni)          &
      - 0.5d0*rrus*(coefa(ifac,iclr)-rtp(iel,irhiph))

!===============================================================================
! 2. FLUX CONVECTIFS DE RUSANOV
!===============================================================================


!     Reperage de la face pour annuler les flux convectifs
!       calcules au bord par bilsc2 ou cfbsc2 pour la qdm (div(rho u u))
!       et l'energie (div(rho u E)) ainsi que les termes en
!       grad(P) et div(u P)

ia(iifbru+ifac-1+(iphas-1)*nfabor) = 1

!     Flux de masse
propfb(ifac,iflmab) = runb*surfbn(ifac)

!     Flux de Qdm (la partie centree en pression pourrait etre prise dans
!       la condition à la limite de pression, ce qui eviterait de retoucher
!       le gradient de pression de la qdm, mais qui donne moins de
!       flexibilité quant à la condition à la limite de pression utilisee
!       pour la reconstruction du gradient, si le maillage est non
!       orthogonal en entree)
propfb(ifac,ipprob(ifbrhu)) = surfbn(ifac)*                &
             0.5d0*(                                              &
             (rund*coefa(ifac,iclu)+runi*rtp(iel,iuiph))          &
             -rrus*( coefa(ifac,iclr)*coefa(ifac,iclu)            &
                    -rtp(iel,irhiph)     *rtp(iel,iuiph)     ))   &
                        + surfbo(1,ifac)*                         &
             0.5d0*(coefa(ifac,iclp)+rtp(iel,ipriph))
propfb(ifac,ipprob(ifbrhv)) = surfbn(ifac)*                &
             0.5d0*(                                              &
             (rund*coefa(ifac,iclv)+runi*rtp(iel,iviph))          &
             -rrus*( coefa(ifac,iclr)*coefa(ifac,iclv)            &
                    -rtp(iel,irhiph)     *rtp(iel,iviph)     ))   &
                        + surfbo(2,ifac)*                         &
             0.5d0*(coefa(ifac,iclp)+rtp(iel,ipriph))
propfb(ifac,ipprob(ifbrhw)) = surfbn(ifac)*                &
             0.5d0*(                                              &
             (rund*coefa(ifac,iclw)+runi*rtp(iel,iwiph))          &
             -rrus*( coefa(ifac,iclr)*coefa(ifac,iclw)            &
                    -rtp(iel,irhiph)     *rtp(iel,iwiph)     ))   &
                        + surfbo(3,ifac)*                         &
             0.5d0*(coefa(ifac,iclp)+rtp(iel,ipriph))
!     Flux de E
propfb(ifac,ipprob(ifbene)) = surfbn(ifac)*                &
             0.5d0*(                                              &
              rund*coefa(ifac,icle)+runi*rtp(iel,ieniph)          &
              +und*coefa(ifac,iclp)+ uni*rtp(iel,ipriph)          &
             -rrus*( coefa(ifac,iclr)*coefa(ifac,icle)            &
                    -rtp(iel,irhiph)     *rtp(iel,ieniph)     ))



!----
! FIN
!----

return

end subroutine
