!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine cfdivs &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iphas  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtp    , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   diverg , ux     , uy     , uz     ,                            &
   vistot ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------
!                                   v
! CALCULE DIVERG = DIVERG + DIV(SIGMA .U)

!          v               t
! AVEC SIGMA = MU (GRAD(U) + GRAD(U)) + (KAPPA - 2/3 MU) DIV(U) Id

! ET MU = MU_LAMINAIRE + MU_TURBULENT

! DIV(U)            EST CALCULE A PARTIR DE FLUMAS/RHO
! GRAD_TRANSPOSE(U) EST UN GRADIENT CELLULE

! REMARQUES :
!  - Theoriquement le terme en div(u) devrait plutot etre calcule
!      par un gradient cellule, pour correspondre exactement au
!      terme en dUj/dxi. Mais comme la partie en dUi/dxj est
!      calculee completement autrement (gradient facette et implicitation)
!      de toute facon on n'aura jamais Trace(tau_ij)=0 exactement.
!  - Pour la meme raison, comme le terme en dUi/dxj est calcule sur les
!      elements de bord et pas celui en dUj/dxi, il est difficile de
!      traiter le terme en div(u) de maniere rigoureuse. Il est donc
!      conserve sur les elements de bord.
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
! ncepdp           ! e  ! <-- ! nombre de cellules avec pdc                    !
! ncesmp           ! e  ! <-- ! nombre de cellules a source de masse           !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! iphas            ! e  ! <-- ! phase courante                                 !
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
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
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
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant)                  !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! diverg(ncelet    ! tr ! --> ! div(sigma.u)                                   !
! ux,y,z(ncelet    ! tr ! <-- ! composantes du vecteur u                       !
! vistot(ncelet    ! tr ! --- ! tableau de travail pour mu                     !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
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
include "cstphy.h"
include "entsor.h"
include "numvar.h"
include "optcal.h"
include "vector.h"
include "period.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse , iphas

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision diverg(ncelet)
double precision ux(ncelet), uy(ncelet), uz(ncelet)
double precision vistot(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          iccocg, inc, iel, ifac, ivar, isou, ii, jj
integer          iuiph, iviph, iwiph
integer          iclvar
integer          nswrgp, imligp, iwarnp
integer          ipcvis, ipcvst, ipcvsv
integer          idimte, itenso, iphydp
double precision epsrgp, climgp, extrap
double precision vecfac, visttt

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

iuiph  = iu(iphas)
iviph  = iv(iphas)
iwiph  = iw(iphas)

ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))
if(iviscv(iphas).gt.0) then
  ipcvsv = ipproc(iviscv(iphas))
else
  ipcvsv = 0
endif


! --- Calcul de la viscosite totale

if (itytur(iphas).eq.3 ) then
  do iel = 1, ncel
    vistot(iel) = propce(iel,ipcvis)
  enddo
else
  do iel = 1, ncel
    vistot(iel) = propce(iel,ipcvis) + propce(iel,ipcvst)
  enddo
endif

!    Pour la periodicite de rotation, il faut avoir calcule
!      le gradient avec grdcel. La seule solution consiste donc a
!      echanger VISTOT puis a faire le produit, y compris sur les
!      cellules halo (calcul sur le halo, exceptionnellement).
!    Pour le parallelisme, on s'aligne sur la sequence ainsi definie.

! ---> TRAITEMENT DU PARALLELISME

if(irangp.ge.0) then
  call parcom (vistot)
  !==========
  if(ipcvsv.gt.0) then
    call parcom (propce(1,ipcvsv))
    !==========
  endif
endif

! ---> TRAITEMENT DE LA PERIODICITE

if(iperio.eq.1) then
  idimte = 0
  itenso = 0
  call percom                                                     &
  !==========
( idimte , itenso ,                                               &
  vistot , vistot , vistot ,                                      &
  vistot , vistot , vistot ,                                      &
  vistot , vistot , vistot )
  if(ipcvsv.gt.0) then
    call percom                                                   &
    !==========
( idimte , itenso ,                                               &
  propce(1,ipcvsv) , propce(1,ipcvsv) , propce(1,ipcvsv) ,        &
  propce(1,ipcvsv) , propce(1,ipcvsv) , propce(1,ipcvsv) ,        &
  propce(1,ipcvsv) , propce(1,ipcvsv) , propce(1,ipcvsv) )
  endif
endif


!===============================================================================
! 2.  CALCUL DES TERMES DE LA DIVERGENCE
!===============================================================================

! --- Boucle sur les composantes de vitesse Ui

do isou = 1, 3

  if (isou.eq.1) ivar = iuiph
  if (isou.eq.2) ivar = iviph
  if (isou.eq.3) ivar = iwiph

! Ceci pointe eventuellement sur ICLRTP(IVAR,ICOEF)
  iclvar = iclrtp(ivar,icoeff)

! --- Calcul du gradient de la vitesse

  iccocg = 1
  inc    = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w6     , w6     , w6     ,                                     &
   rtp(1,ivar)     , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )



! --- Assemblage sur les faces internes

! On a echange le gradient dans grdcel et vistot plus haut

  if(ipcvsv.gt.0) then
    if    (isou.eq.1) then
      do iel = 1, ncelet
        visttt = propce(iel,ipcvsv) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*( 2.d0*w1(iel)*ux(iel)              &
                                   + w2(iel)*uy(iel)              &
                                   + w3(iel)*uz(iel) )            &
                     + visttt*w1(iel)*ux(iel)
        w5(iel) = vistot(iel)*w2(iel)*ux(iel)                     &
                     + visttt*w1(iel)*uy(iel)
        w6(iel) = vistot(iel)*w3(iel)*ux(iel)                     &
                     + visttt*w1(iel)*uz(iel)
      enddo

    elseif(isou.eq.2) then
      do iel = 1, ncelet
        visttt = propce(iel,ipcvsv) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*w1(iel)*uy(iel)                     &
                     + visttt*w2(iel)*ux(iel)
        w5(iel) = vistot(iel)*(      w1(iel)*ux(iel)              &
                              + 2.d0*w2(iel)*uy(iel)              &
                                   + w3(iel)*uz(iel) )            &
                     + visttt*w2(iel)*uy(iel)
        w6(iel) = vistot(iel)*w3(iel)*uy(iel)                     &
                     + visttt*w2(iel)*uz(iel)
      enddo

    elseif(isou.eq.3) then
      do iel = 1, ncelet
        visttt = propce(iel,ipcvsv) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*w1(iel)*uz(iel)                     &
                     + visttt*w3(iel)*ux(iel)
        w5(iel) = vistot(iel)*w2(iel)*uz(iel)                     &
                     + visttt*w3(iel)*uy(iel)
        w6(iel) = vistot(iel)*(      w1(iel)*ux(iel)              &
                                   + w2(iel)*uy(iel)              &
                              + 2.d0*w3(iel)*uz(iel) )            &
                     + visttt*w3(iel)*uz(iel)
      enddo

    endif

  else

    if    (isou.eq.1) then
      do iel = 1, ncelet
        visttt = viscv0(iphas) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*( 2.d0*w1(iel)*ux(iel)              &
                                   + w2(iel)*uy(iel)              &
                                   + w3(iel)*uz(iel) )            &
                     + visttt*w1(iel)*ux(iel)
        w5(iel) = vistot(iel)*w2(iel)*ux(iel)                     &
                     + visttt*w1(iel)*uy(iel)
        w6(iel) = vistot(iel)*w3(iel)*ux(iel)                     &
                     + visttt*w1(iel)*uz(iel)
      enddo

    elseif(isou.eq.2) then
      do iel = 1, ncelet
        visttt = viscv0(iphas) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*w1(iel)*uy(iel)                     &
                     + visttt*w2(iel)*ux(iel)
        w5(iel) = vistot(iel)*(      w1(iel)*ux(iel)              &
                              + 2.d0*w2(iel)*uy(iel)              &
                                   + w3(iel)*uz(iel) )            &
                     + visttt*w2(iel)*uy(iel)
        w6(iel) = vistot(iel)*w3(iel)*uy(iel)                     &
                     + visttt*w2(iel)*uz(iel)
      enddo

    elseif(isou.eq.3) then
      do iel = 1, ncelet
        visttt = viscv0(iphas) - 2.d0/3.d0*vistot(iel)
        w4(iel) = vistot(iel)*w1(iel)*uz(iel)                     &
                     + visttt*w3(iel)*ux(iel)
        w5(iel) = vistot(iel)*w2(iel)*uz(iel)                     &
                     + visttt*w3(iel)*uy(iel)
        w6(iel) = vistot(iel)*(      w1(iel)*ux(iel)              &
                                   + w2(iel)*uy(iel)              &
                              + 2.d0*w3(iel)*uz(iel) )            &
                     + visttt*w3(iel)*uz(iel)
      enddo

    endif

  endif



! On initialise DIVERG(NCEL+1, NCELET)
!     (valeur bidon, mais pas NaN : les calculs sur le halo sont
!      par principe denue de sens, sauf exception)
  if(ncelet.gt.ncel) then
    do iel = ncel+1, ncelet
      diverg(iel) = 0.d0
    enddo
  endif



  if(ivecti.eq.1) then

!CDIR NODEP
    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
!MO             VECFAC = SURFAC(ISOU,IFAC)
!MO     &                  *(POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
      vecfac = surfac(1,ifac)*(w4(ii)+w4(jj))*0.5d0               &
             + surfac(2,ifac)*(w5(ii)+w5(jj))*0.5d0               &
             + surfac(3,ifac)*(w6(ii)+w6(jj))*0.5d0
      diverg(ii) = diverg(ii) + vecfac
      diverg(jj) = diverg(jj) - vecfac
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1, nfac
      ii = ifacel(1,ifac)
      jj = ifacel(2,ifac)
!MO             VECFAC = SURFAC(ISOU,IFAC)
!MO     &                  *(POND(IFAC)*W4(II)+(1.D0-POND(IFAC))*W4(JJ))
      vecfac = surfac(1,ifac)*(w4(ii)+w4(jj))*0.5d0               &
             + surfac(2,ifac)*(w5(ii)+w5(jj))*0.5d0               &
             + surfac(3,ifac)*(w6(ii)+w6(jj))*0.5d0
      diverg(ii) = diverg(ii) + vecfac
      diverg(jj) = diverg(jj) - vecfac
    enddo

  endif


! --- Assemblage sur les faces de bord

  if(ivectb.eq.1) then

!CDIR NODEP
    do ifac = 1, nfabor
      ii = ifabor(ifac)
      vecfac = surfbo(1,ifac)*w4(ii)                              &
             + surfbo(2,ifac)*w5(ii)                              &
             + surfbo(3,ifac)*w6(ii)
      diverg(ii) = diverg(ii) + vecfac
    enddo

  else

! VECTORISATION NON FORCEE
    do ifac = 1, nfabor
      ii = ifabor(ifac)
      vecfac = surfbo(1,ifac)*w4(ii)                              &
             + surfbo(2,ifac)*w5(ii)                              &
             + surfbo(3,ifac)*w6(ii)
      diverg(ii) = diverg(ii) + vecfac
    enddo

  endif

enddo


return

end
