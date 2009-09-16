!-------------------------------------------------------------------------------

!VERS


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

subroutine usvort &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , iappel ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr , irepvo ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

! METHODE DES VORTEX POUR LES CONDITIONS AUX LIMITES D'ENTREE
!  EN L.E.S. :
!  DEFINITION DES ENTREES AVEC VORTEX
!  DEFINITION DES CARACTERISTIQUES DES VORTEX


! Boundary faces identification
! =============================

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


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
! iphas            ! e  ! <-- ! numero de la phase                             !
! iappel           ! e  ! <-- ! indique les donnes a renvoyer                  !
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
! maxelt           !  e ! <-- ! nb max d'elements (cell,fac,fbr)               !
! lstelt(maxelt) te ! --- ! tableau de travail                             !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! irepvo           ! te ! <-- ! numero de l'entree associe a chaque            !
!     (nfabor)     !    !     ! face de bord (=0 si pas de vortex)             !
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
include "optcal.h"
include "entsor.h"
include "vortex.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas  , iappel

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          irepvo(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          ifac, ient
integer          ilelt, nlelt

!===============================================================================
! 1. PARAMETRES GLOBAUX
!===============================================================================

! --- Nombre d'entrees avec la methode des vortex

nnent = 2

! --- Nombre de vortex a mettre dans chaque entree

!   NVORT min ~ Surface d'entree/(pi*SIGMA**2)

nvort(1) = 500
nvort(2) = 500

if (iappel.eq.1) then

!===============================================================================
! 2. DEFINITION DES ZONES D'ENTREE (AU PREMIER PASSAGE)
!===============================================================================

  do ifac = 1, nfabor
    irepvo(ifac) = 0
  enddo

! ------------------
!   ENTREE 1
! ------------------
  CALL GETFBR('3',NLELT,LSTELT)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    ient = 1
    irepvo(ifac) = ient

  enddo

! ------------------
!   ENTREE 2
! ------------------
  CALL GETFBR('1',NLELT,LSTELT)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    ient = 2
    irepvo(ifac) = ient

  enddo

elseif (iappel.eq.2) then

!===============================================================================
! 3. PARAMETRES GEOMETRIQUES ET CONDITIONS LIMITES
!===============================================================================

! --- Cas traité

! ICAS = 1...Conduite rectangulaire
!        2...Conduite circulaire
!        3...Geometrie quelconque sans traitement specifique des conditions aux limites
!        4...Geometrie quelconque sans traitement specifique des conditions aux limites
!            ni fichier de donnees (la vitesse moyenne, le niveau de k et de epsilon
!            sont fournis par l'utilisateur)

  ient = 1
  icas(ient) = 1

  ient = 2
  icas(ient) = 2


! --- Repere definissant le plan d'entree

!     Si ICAS = 4, le code se charge de ces donnees
!     Sinon il faut preciser les vecteurs DIR1 et DIR2 definissant
!     un repère directe tel que DIR3 soit un vecteur entrant normal
!     a la face d'entree.

  ient = 1
  if(icas(ient).eq.1.or.icas(ient).eq.2.or.icas(ient).eq.3) then
    dir1(1,ient) = 1.d0
    dir1(2,ient) = 0.d0
    dir1(3,ient) = 0.d0

    dir2(1,ient) = 0.d0
    dir2(2,ient) = 1.d0
    dir2(3,ient) = 0.d0
  endif

  ient = 2
  if(icas(ient).eq.1.or.icas(ient).eq.2.or.icas(ient).eq.3) then
    dir1(1,ient) = 0.d0
    dir1(2,ient) = 1.d0
    dir1(3,ient) = 0.d0

    dir2(1,ient) = 0.d0
    dir2(2,ient) = 0.d0
    dir2(3,ient) = 1.d0
  endif

! --- Centre du repere local dans le plan d'entree

!     Si ICAS = 1 ou ICAS = 2, le centre du repere doit correspondre
!                 au centre de gravite de la zone d'entree (rectangle ou cercle)

  ient = 1

  cen(1,ient) = 0.d0
  cen(2,ient) = 0.d0
  cen(3,ient) = -6.05d-1

  ient = 2

  cen(1,ient) = -3.664d-1
  cen(2,ient) = 0.d0
  cen(3,ient) = 0.d0

! --- Condition aux limites

! -> Si ICAS = 1...Il faut specifier le type de condition aux limite ICLVOR
!               dans les directions DIR1, DIR2, - DIR1, -DIR2

!               Ces conditions peuvent etre de 3 types :

! ICLVOR = 1...Condition de paroi
!          2...Condition de symetrie
!          3...Condition de periodicite

!                    y = LLY/2
!                    (ICLVOR 1)
!           +-----------------------+
!           |           ^ DIR1      |
!           |           |           |
!           |           |           |
! z=- LLZ/2 |           +----> DIR2 | z = LLZ/2
! (ICLVOR 4)|                       | (ICLVOR 2)
!           |                       |
!           |                       |
!           +-----------------------+
!                    y = -LLY/2
!                    (ICLVOR 3)


! -> Si ICAS = 2, les conditions sont necessairement de type paroi
! -> Si ICAS = 3 ou 4, pas de traitement particulier

  ient = 1

  if(icas(ient).eq.1) then
    iclvor(1,ient) = 1
    iclvor(2,ient) = 2
    iclvor(3,ient) = 1
    iclvor(4,ient) = 2
  endif

! LLY et LLZ sont les dimensions de l'entree dans les directions DIR1 et DIR2
! LDD est le diametre de la conduite


  ient = 1
  lly(ient) = 0.2d0
  llz(ient) = 0.1d0

  ient = 2
  lld(2) = 0.154d0

!===============================================================================
! 5. PARAMETRES PHYSIQUES ET MARCHE EN TEMPS
!===============================================================================

! --- " Temps de vie " limite du vortex

! ITLIVO = 1...Les vortex sont retire au bout du temps TLIMVO
!                donne par l'utilisateur
!                ( par exemple TLIMVO = 10*DTREF)

!          2...Chaque vortex a un temps d'exitence limite valant
!                5.Cmu.k^(3/2).U/epsilon
!               ( ou U est la vitesse principale suivant DIR3)

  ient = 1
  itlivo(ient) = 1

  if(itlivo(ient).eq.1) then
    tlimvo(ient) = 10.d0*dtref
  endif

  ient = 2
  itlivo(ient) = 2


! --- " Diametre " des vortex

! ISGMVO = 1...diametre constant XSGMVO donne par l'utilisateur
!          2...basee sur la formule sigma = Cmu^(3/4).k^(3/2)/epsilon
!          3...basee sur la formule sigma = max(Lt, Lk) avec
!                 Lt = (5 nu.k/epsilon)^(1/2)
!             et  Lk = 200.(nu^3/epsilon)^(1/4)

  ient = 1
  isgmvo(ient) = 1

  if(isgmvo(ient).eq.1) then
    xsgmvo(ient) = 0.01d0
  endif

  ient = 2
  isgmvo(ient) = 2


! --- Mode de deplacement des vortex

! IDEPVO = 1...Deplacement en r*UD (r aleatoire dans [0,1])
!              UD a fournir par l'utilisateur
!          2...Convection par les vortex

  ient = 1
  idepvo(ient) = 2

  ient = 2
  idepvo(ient) = 1

  if(idepvo(ient).eq.1) then
    ud(ient) = 0.7d0
  endif

!===============================================================================
! 6. PARAMETRES D'ENTREE / SORTIES ET DONNEES UTILISATEUR
!===============================================================================

! --- Fichier de donnees utilisateur

! NDAT ...Nombre de lignes du fichier de donnees contenant les donnees :
!          x | y | z | U | V | W | Grad[u.DIR3].n | k | epsilon

!         dans le plan d'entree du calcul

!         Grad[u.DIR3].n est le gradient dans la direction normale
!         a la paroi, de la vitesse principale dans le plan d'entree.
!         Cette données n'est utilisée qu'avec ICAS=2

! FICVOR...Nom du fichier de donnees utilisateur

  ient = 1
  ndat(ient) = 2080

  ient = 2
  ndat(ient) = 2080

! Par les defaut les fichiers sont nommes "vordat" affecté de l'indice
! d'entrée

  ient = 1
  FICVOR(IENT) = 'entree_1.dat'

  ient = 2
  FICVOR(IENT) = 'entree_2.dat'

! Pour ICAS = 4, on precise juste la valeur moyenne de U, k et de espilon
! a l'entree

  if(icas(ient).eq.4) then
    udebit(ient) = 10.d0
    kdebit(ient) = 1.d0
    edebit(ient) = 1.d0
  endif

! --- Relecture d'un fichier suite eventuel

! ISUIVO = 0...Pas de relecture (reinitialisation des vortex)
!          1...Relecture du fichier suite de methode des vortex

  isuivo = isuite


endif


return
end

!===============================================================================
! 7. DEFINTION DE LA FONCTION PERMETAT D'IMPOSER LES DONNEES D'ENTREE
!===============================================================================

                   function phidat                                &
!==============

 ( nfecra , icas   , ndat   ,                                     &
   yy     , zz     , ydat   , zdat   ,                            &
   vardat , iii    )

!===============================================================================
! FONCTION :
! --------

! FONCTION PERMETTANT D'INTERPOLER LES DONNEES D'ENTREE FOURNIES
! PAR L'UTILISATEUR AU CENTRE DES FACES D'ENTREE POUR LESQUELLES
! EST UTILISEE LA METHODE DES VORTEX

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nfecra           ! e  ! <-- ! unite                                          !
! icas             ! e  ! <-- ! type de geometrie du cas                       !
! ndat             ! e  ! <-- ! nbr de lignes du fichier de donnees            !
! yy               ! e  ! <-- ! coordoonnes dans le repere local du            !
! zz               ! e  ! <-- ! point ou l'on cherche a connaitre la           !
!                  !    !     ! variable vardat                                !
! ydat             ! e  ! <-- ! coordoonnes ou est connue la variable          !
! zdat             ! e  ! <-- ! vardat dans le fichier de donnees              !
! vardat           ! e  ! <-- ! valeur de la variable vardat                   !
! iii              ! e  ! --> ! ligne ou a ete trouvee la donnee la            !
!                  !    !     ! plus proche du point (yy,zz)                   !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

integer          nfecra, icas, ndat, iii
double precision zz, yy
double precision zdat(ndat), ydat(ndat)
double precision vardat(ndat)

integer          ii
double precision phidat, dist1


! Dans l'exemple suivant, on se contente de retourne la valeur situee
! dans le fichier de donnee a l'abscisse la plus proche du point de
! coordonnée (Y,Z) ou l'on cherche a connaitre la valeur de la
! variable numero VARDAT.


if(icas.eq.1.or.icas.eq.2.or.icas.eq.3) then

  if(iii.eq.0) then
    dist1 = 1.d20
    do ii = 1,ndat
      if(sqrt((yy-ydat(ii))**2+(zz-zdat(ii))**2).lt.dist1) then
        dist1 = sqrt((zz-zdat(ii))**2+(yy-ydat(ii))**2)
        iii   = ii
        phidat = vardat(ii)
      endif
    enddo
  elseif(iii.ne.0) then
    phidat =  vardat(iii)
  endif

elseif(icas.eq.4) then
  phidat = vardat(1)
endif

return
end
