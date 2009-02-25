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

subroutine mttycl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   itypfb , itrifb , icodcl , isostd ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl , frcxt  ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! TRAITEMENT DES CODES DE CONDITIONS AUX LIMITES DE MATISSE

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
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas)          !    !     !                                                !
! itrifb(nfabor    ! te ! --> ! tab d'indirection pour tri des faces           !
!  nphas)          !    !     !                                                !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! isostd           ! te ! --> ! indicateur de sortie standard                  !
!    (nfabor+1)    !    !     !  +numero de la face de reference               !
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
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! frcxt(ncelet,    ! tr ! <-- ! force exterieure generant la pression          !
!   3,nphas)       !    !     !  hydrostatique                                 !
! w1,2,3,4,5,6     ! tr ! --- ! tableaux de travail                            !
!  (ncelet         !    !     !  (calcul du gradient de pression)              !
! rijipb           ! tr ! --- ! tab de trav pour valeurs en iprime             !
! (nfabor,6   )    !    !     !  des rij au bord                               !
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

include "dimfbr.h"
include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstnum.h"
include "cstphy.h"
include "entsor.h"
include "pointe.h"
include "parall.h"
include "matiss.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          isostd(nfabor+1,nphas)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision rcodcl(nfabor,nvar,3)
double precision frcxt(ncelet,3,nphas)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,3)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          iphas , ifac  , ifml  , icoul
double precision dbm   , roe

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! --- Gestion memoire

idebia = idbia0
idebra = idbra0

! --- Une seule phase

iphas = 1

!===============================================================================
! 2.  BOUCLE SUR LES FACES DE BORD
!===============================================================================

do ifac = 1, nfabor

! --- Couleur de la face de bord

  ifml  = ifmfbr(ifac  )
  icoul = iprfml(ifml,1)


!===============================================================================
! 3.  ENTREE : classique en convection forcee, a pression imposee sinon
!===============================================================================

  if(icoul.eq.icmtfi) then


! --- Convection forcee
    if(icofor.eq.1)then

!     Type de base : IENTRE
      itypfb(ifac,iphas) = ientre

!     Temperature exterieure prise a TINIT
      rcodcl(ifac,isca(itaamt),1) = tinit
      rcodcl(ifac,isca(itpcmt),1) = tinit
      rcodcl(ifac,isca(itppmt),1) = tinit

!     Vitesses transverses nulles
      rcodcl(ifac,iu(iphas),1) = 0.d0
      rcodcl(ifac,iv(iphas),1) = 0.d0

!     Vitesse debitante W est calculee a partir du debit massique par
!       W = - debit massique / (masse volumique * surface)
!     . On suppose la vitesse verticale descendante (d'ou le signe).
!     . ROE est la masse volumique exterieure evaluee a partir de
!       la temperature TINIT (degres C) en utilisant la loi d'etat.
!     . La surface de l'entree du maillage est calculee par la formule
!       NPTRAN*PTRRES*EPCHEM/RCONVE. En effet, la surface au bas de la
!       cheminee est NPTRAN*PTRRES*EPCHEM et on divise par le rapport
!       du convergent represente sur le maillage (>=1) pour obtenir
!       la surface d'entree. (en 2D, il faut RCONVE=1, sinon, on impose
!       un debit plus grand que le debit souhaite).
!     . La correction du debit massique par FRDTRA correspond au
!       rapport d'echelle transverse eventuel entre le maillage et la
!       realite    .

      roe = rrfmat*(trfmat+tkelvi)/(tinit+tkelvi)
      dbm = debmas/frdtra
      rcodcl(ifac,iw(iphas),1) =                                  &
             -dbm*rconve/(roe*nptran*ptrres*epchem)


! --- Convection naturelle (ICOFOR=0)
    else

!     Type de base indefini
!       (dans une version ulterieure de Code_Saturne, on disposera de
!        sorties a pression imposee)
      itypfb(ifac,iphas) = iindef

!     Dirichlet sur les temperatures
      icodcl(ifac,isca(itaamt)  ) = 1
      rcodcl(ifac,isca(itaamt),1) = tinit
      icodcl(ifac,isca(itpcmt)  ) = 1
      rcodcl(ifac,isca(itpcmt),1) = tinit
      icodcl(ifac,isca(itppmt)  ) = 1
      rcodcl(ifac,isca(itppmt),1) = tinit

!     Dirichlet nul sur les vitesses transverses (nul par symetrie)
      icodcl(ifac,iu(iphas)  ) = 1
      rcodcl(ifac,iu(iphas),1) = 0.d0
      icodcl(ifac,iv(iphas)  ) = 1
      rcodcl(ifac,iv(iphas),1) = 0.d0

!     Neumann homogene sur la vitesse debitante (on impose la pression)
      icodcl(ifac,iw(iphas)  ) = 3

!     Pression nulle
      icodcl(ifac,ipr(iphas)  ) = 1
      rcodcl(ifac,ipr(iphas),1) = p0(iphas)

    endif


!===============================================================================
! 4.  SORTIE (pression imposee)
!===============================================================================

  elseif(icoul.eq.icmtfo) then

!     Type de base : sortie libre
    itypfb(ifac,iphas) = isolib

!     Temperature imposee en cas de reentree
    rcodcl(ifac,isca(itaamt),1) = tcrit
    rcodcl(ifac,isca(itpcmt),1) = tinit
    rcodcl(ifac,isca(itppmt),1) = tinit

!     Pression imposee (ecart hydrostatique par rapport a l'entree
!                       + decalage eventuel DPVENT)
    icodcl(ifac,ipr(iphas)  ) = 1
    rcodcl(ifac,ipr(iphas),1) = p0(iphas) + dpvent                &
         - dabs(gz) * ro0(iphas) * (hcheva - hchali)

!===============================================================================
! 5.  LE RESTE : glissement et Neumann homogene
!===============================================================================

! --- Sol
  elseif(icoul.eq.icmtfg) then
    itypfb(ifac,iphas)   = isymet

! --- Plafond
  elseif(icoul.eq.icmtfc) then
    itypfb(ifac,iphas)   = isymet

! --- Symetries
  elseif(icoul.eq.icmtfs) then
    itypfb(ifac,iphas)   = isymet

! --- Parois
  elseif(icoul.eq.icmtfw) then
    itypfb(ifac,iphas)   = isymet

  endif

enddo


!----
! FORMATS
!----

!----
! FIN
!----

return

end
