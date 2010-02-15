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

subroutine laggra &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtp    , propce , coefa  , coefb  ,                            &
   gradpr , gradvf ,                                              &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    1)   Calcul de (- (GRADIENT DE PRESSION) / ROM )

!    2)   Calcul du gradient de Vitesse

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
!  (nfml,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac                      !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr                      !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !                                                !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! coefa,coefb      ! tr ! <-- ! tableaux des cond lim pour pvar                !
!   (nfabor)       !    !     !  sur la normale a la face de bord              !
! gradpr(ncel,3    ! tr ! --> ! gradient de pression                           !
! gradvf(ncel,9    ! tr ! --> ! gradient de vitesse fluide                     !
! w1...w3(ncel)    ! tr ! --- ! tableau de travail                             !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "dimfbr.h"
include "paramx.h"
include "numvar.h"
include "optcal.h"
include "entsor.h"
include "cstphy.h"
include "pointe.h"
include "period.h"
include "parall.h"
include "lagpar.h"
include "lagran.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndnod , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision coefa(ndimfb,*) , coefb(ndimfb,*)
double precision rtp(ncelet,*)
double precision propce(ncelet,*)
double precision gradpr(ncelet,3) , gradvf(ncelet,9)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra

integer          inc , iccocg , ipriph , iclipr
integer          iuiph , iviph , iwiph
integer          ipcliu , ipcliv , ipcliw
integer          iromf  , idimte , itenso
integer          iel    , iphydp , iphas
double precision unsrho

!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 0. PARAMETRES
!===============================================================================

!-->Rappel de differents paramtrage possibles pour le calcul du gradient

!     INC     = 1               ! 0 RESOL SUR INCREMENT 1 SINON
!     ICCOCG  = 1               ! 1 POUR RECALCUL DE COCG 0 SINON
!     NSWRGR  = 100             ! mettre 1 pour un maillage regulie (par defaut : 100)
!     IMLIGR  = -1              ! LIMITATION DU GRADIENT : < 0 PAS DE LIMITATION
!     IWARNI  = 1               ! NIVEAU D'IMPRESSION
!     EPSRGR  = 1.D-8           ! PRECISION RELATIVE POUR LA REC GRA 97
!     CLIMGR  = 1.5D0           ! COEF GRADIENT*DISTANCE/ECART
!     EXTRAG  = 0               ! COEF GRADIENT*DISTANCE/ECART
!     IPHYDP  = 0               ! =1 NE CONCERNE QUE LE GRADIENT DE PRESSION


!-->Parametrage des calculs des gradients

inc     = 1
iccocg  = 1
iphydp  = 0

!-->Numero de la phase continue qui est la phase porteuse

iphas = ilphas

!===============================================================================
! 1. CALCUL DE :  - (GRADIENT DE PRESSION)/ROM
!===============================================================================

ipriph = ipr(iphas)
iclipr = iclrtp(ipriph,icoef)

! Calcul du gradient de pression


! ---> TRAITEMENT DU PARALLELISME

if (irangp.ge.0) call parcom (rtp(1,ipriph))
                 !==========

! ---> TRAITEMENT DE LA PERIODICITE

if (iperio.eq.1) then
  idimte = 0
  itenso = 0
  call percom                                                     &
  !==========
  (idimte , itenso ,                                              &
   rtp(1,ipriph)  , rtp(1,ipriph) , rtp(1,ipriph) ,               &
   rtp(1,ipriph)  , rtp(1,ipriph) , rtp(1,ipriph) ,               &
   rtp(1,ipriph)  , rtp(1,ipriph) , rtp(1,ipriph)  )
endif

call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ipriph , imrgra , inc    , iccocg ,                            &
   nswrgr(ipriph)  , imligr(ipriph)  , iphydp ,                   &
   iwarni(ipriph)  , nfecra ,                                     &
   epsrgr(ipriph)  , climgr(ipriph)  , extrag(ipriph)  ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtp(1,ipriph)  , coefa(1,iclipr) , coefb(1,iclipr) ,           &
   gradpr(1,1)     , gradpr(1,2)     , gradpr(1,3)     ,          &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

! Pointeur sur la masse volumique en fonction de l'ecoulement

if ( ippmod(icp3pl).ge.0 .or. ippmod(icfuel).ge.0 ) then
  iromf = ipproc(irom1)
else
  iromf = ipproc(irom(iphas))
endif

! Calcul de -Grad P / Rom

do iel = 1,ncel
  unsrho = 1.d0 / propce(iel,iromf)
  gradpr(iel,1) = -gradpr(iel,1) * unsrho
  gradpr(iel,2) = -gradpr(iel,2) * unsrho
  gradpr(iel,3) = -gradpr(iel,3) * unsrho
enddo

!===============================================================================
! 2. CALCUL DU GRADIENT DE VITESSE :
!===============================================================================

if (modcpl.gt.0 .and. iplas.ge.modcpl) then

  iuiph = iu(iphas)
  iviph = iv(iphas)
  iwiph = iw(iphas)
  ipcliu = iclrtp(iuiph,icoef)
  ipcliv = iclrtp(iviph,icoef)
  ipcliw = iclrtp(iwiph,icoef)

! ---> TRAITEMENT DU PARALLELISME

  if(irangp.ge.0) then
    call parcom (rtp(1,iuiph))
    !==========
    call parcom (rtp(1,iviph))
    !==========
    call parcom (rtp(1,iwiph))
    !==========
  endif

! ---> TRAITEMENT DE LA PERIODICITE

  if(iperio.eq.1) then
    idimte = 1
    itenso = 0
    call percom                                                   &
    !==========
    (idimte , itenso ,                                            &
     rtp(1,iuiph),rtp(1,iuiph),rtp(1,iuiph),                      &
     rtp(1,iviph),rtp(1,iviph),rtp(1,iviph),                      &
     rtp(1,iwiph),rtp(1,iwiph),rtp(1,iwiph) )
  endif

!     COMPOSANTE X
!     ============

!    Sans prise en compte de la pression hydrostatique

  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iuiph  , imrgra , inc    , iccocg ,                            &
   nswrgr(iuiph)   , imligr(iuiph)  , iphydp ,                    &
   iwarni(iuiph)   , nfecra ,                                     &
   epsrgr(iuiph)   , climgr(iuiph)  , extrag(iuiph)  ,            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtp(1,iuiph)   , coefa(1,ipcliu) , coefb(1,ipcliu) ,           &
   gradvf(1,1)     , gradvf(1,2)     , gradvf(1,3)     ,          &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

!     COMPOSANTE Y
!     ============

!    Sans prise en compte de la pression hydrostatique

  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iviph  , imrgra , inc    , iccocg ,                            &
   nswrgr(iviph)   , imligr(iviph)  , iphydp ,                    &
   iwarni(iviph)   , nfecra ,                                     &
   epsrgr(iviph)   , climgr(iviph)  , extrag(iviph)  ,            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtp(1,iviph)   , coefa(1,ipcliv) , coefb(1,ipcliv) ,           &
   gradvf(1,4)     , gradvf(1,5)     , gradvf(1,6)     ,          &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

!     COMPOSANTE Z
!     ============

!    Sans prise en compte de la pression hydrostatique

  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iwiph  , imrgra , inc    , iccocg ,                            &
   nswrgr(iwiph)   , imligr(iwiph)  , iphydp ,                    &
   iwarni(iwiph)   , nfecra ,                                     &
   epsrgr(iwiph)   , climgr(iwiph)  , extrag(iwiph)  ,            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtp(1,iwiph)   , coefa(1,ipcliw) , coefb(1,ipcliw) ,           &
   gradvf(1,7)     , gradvf(1,8)     , gradvf(1,9)     ,          &
   w1     , w2     , w3     ,                                     &
   rdevel , rtuser , ra     )

endif

!----
! FIN
!----

return

end subroutine
