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
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   imodif , iphas  ,                                              &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   gammag ,                                                       &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   sorti1 , sorti2 , gamagr , masmor ,                            &
   rdevel , rtuser , ra     )

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
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ndim             ! e  ! <-- ! dimension de l'espace                          !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nfml             ! e  ! <-- ! nombre de familles d entites                   !
! nprfml           ! e  ! <-- ! nombre de proprietese des familles             !
! nnod             ! e  ! <-- ! nombre de sommets                              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac (optionnel          !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr (optionnel          !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! imodif           ! e  ! <-- ! modification directe de rtp (imodif=1          !
! iphas            ! e  ! <-- ! numero de la phase courante                    !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! nfml  ,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! gammag           ! r  ! <-- ! gamma du gaz                                   !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! <-- ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! <-- ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! <-- ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! <-- ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds (optionnel)              !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp,rtpa         ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! sorti1,2(*)      ! tr ! --> ! variables de sortie                            !
! gamagr(*)        ! tr ! --- ! constante gamma equivalent du gaz              !
! masmor(*)        ! tr ! --- ! masse molaire des constituants du gaz          !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "pointe.h"
include "entsor.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "cfpoin.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          imodif , iphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision gammag
double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*),propfa(nfac,*),propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision sorti1(*), sorti2(*), gamagr(*), masmor(*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! VARIABLES LOCALES

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

ipriph = ipr(iphas)
irhiph = isca(irho  (iphas))
ieniph = isca(ienerg(iphas))
iuiph  = iu(iphas)
iviph  = iv(iphas)
iwiph  = iw(iphas)
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
       + coefa(ifac,iclw)*surfbo(3,ifac))/ra(isrfbn+ifac-1)
uni   = (rtp(iel,iuiph)*surfbo(1,ifac)                            &
       + rtp(iel,iviph)*surfbo(2,ifac)                            &
       + rtp(iel,iwiph)*surfbo(3,ifac))/ra(isrfbn+ifac-1)
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
propfb(ifac,iflmab) = runb*ra(isrfbn+ifac-1)

!     Flux de Qdm (la partie centree en pression pourrait etre prise dans
!       la condition à la limite de pression, ce qui eviterait de retoucher
!       le gradient de pression de la qdm, mais qui donne moins de
!       flexibilité quant à la condition à la limite de pression utilisee
!       pour la reconstruction du gradient, si le maillage est non
!       orthogonal en entree)
propfb(ifac,ipprob(ifbrhu(iphas))) = ra(isrfbn+ifac-1)*           &
             0.5d0*(                                              &
             (rund*coefa(ifac,iclu)+runi*rtp(iel,iuiph))          &
             -rrus*( coefa(ifac,iclr)*coefa(ifac,iclu)            &
                    -rtp(iel,irhiph)     *rtp(iel,iuiph)     ))   &
                        + surfbo(1,ifac)*                         &
             0.5d0*(coefa(ifac,iclp)+rtp(iel,ipriph))
propfb(ifac,ipprob(ifbrhv(iphas))) = ra(isrfbn+ifac-1)*           &
             0.5d0*(                                              &
             (rund*coefa(ifac,iclv)+runi*rtp(iel,iviph))          &
             -rrus*( coefa(ifac,iclr)*coefa(ifac,iclv)            &
                    -rtp(iel,irhiph)     *rtp(iel,iviph)     ))   &
                        + surfbo(2,ifac)*                         &
             0.5d0*(coefa(ifac,iclp)+rtp(iel,ipriph))
propfb(ifac,ipprob(ifbrhw(iphas))) = ra(isrfbn+ifac-1)*           &
             0.5d0*(                                              &
             (rund*coefa(ifac,iclw)+runi*rtp(iel,iwiph))          &
             -rrus*( coefa(ifac,iclr)*coefa(ifac,iclw)            &
                    -rtp(iel,irhiph)     *rtp(iel,iwiph)     ))   &
                        + surfbo(3,ifac)*                         &
             0.5d0*(coefa(ifac,iclp)+rtp(iel,ipriph))
!     Flux de E
propfb(ifac,ipprob(ifbene(iphas))) = ra(isrfbn+ifac-1)*           &
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
