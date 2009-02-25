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

subroutine usmpst &
!================

 ( idbia0 , idbra0 , ipart  ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   ncelps , nfacps , nfbrps ,                                     &
   nideve , nrdeve , nituse , nrtuse , imodif ,                   &
   itypps , ifacel , ifabor , ifmfbr , ifmcel , iprfml ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statis ,                                     &
   tracel , trafac , trafbr , rdevel , rtuser , ra     )
!===============================================================================
! FONCTION :
! --------

! ROUTINE UTILISATEUR POUR LA MODIFICATION DES LISTES DE CELLULES
! OU FACES INTERNES ET DE BORD DEFINISSANT UN MAILLAGE DE POST
! TRAITEMENT EXISTANT ; CETTE ROUTINE EST APPELEE AUX PAS DE
! TEMPS AUQUEL CE MAILLAGE EST ACTIF, ET UNIQUEMENT POUR LES
! MAILLAGES POST UTILISATEUR PRINCIPAUX (NON ALIAS), SI TOUS LES
! "WRITERS" ASSOCIES A CE MAILLAGE OU SES ALIAS PERMETTENT
! CETTE MODIFICATION
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ipart            ! e  ! <-- ! numero du maillage post                        !
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
! nvlsta           ! e  ! <-- ! nombre de variables stat. lagrangien           !
! ncelps           ! e  ! <-- ! nombre de cellules du maillage post            !
! nfacps           ! e  ! <-- ! nombre de faces interieur post                 !
! nfbrps           ! e  ! <-- ! nombre de faces de bord post                   !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! imodif           ! e  ! <-- ! 0 si maillage non modifie par cette            !
!                  !    !     ! fonction, 1 si modifie                         !
! itypps(3)        ! te ! <-- ! indicateur de presence (0 ou 1) de             !
!                  !    !     ! cellules (1), faces (2), ou faces de           !
!                  !    !     ! de bord (3) dans le maillage post              !
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
! lstcel(ncelps    ! te ! <-- ! liste des cellules du maillage post            !
! lstfac(nfacps    ! te ! <-- ! liste des faces interieures post               !
! lstfbr(nfbrps    ! te ! <-- ! liste des faces de bord post                   !
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! te ! --- ! macro tableau entier                           !
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
! (ncelet)         !    !     !                                                !
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
! statis           ! tr ! <-- ! statistiques (lagrangien)                      !
!ncelet,nvlsta)    !    !     !                                                !
! tracel(*)        ! tr ! <-- ! tab reel valeurs cellules post                 !
! trafac(*)        ! tr ! <-- ! tab reel valeurs faces int. post               !
! trafbr(*)        ! tr ! <-- ! tab reel valeurs faces bord post               !
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
include "pointe.h"
include "entsor.h"
include "optcal.h"
include "numvar.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ipart
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , nvlsta
integer          ncelps , nfacps , nfbrps
integer          nideve , nrdeve , nituse , nrtuse, imodif

integer          itypps(3)
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtpa(ncelet,*), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision statis(ncelet,nvlsta)
double precision tracel(ncelps*3)
double precision trafac(nfacps*3), trafbr(nfbrps*3)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! VARIABLES LOCALES

integer          ifac  , iphas
integer          ii   , jj
double precision vmin2, v2, w2


!===============================================================================

!     Remarque : le tableau ITYPPS permet de savoir si le maillage post
!                contient a l'origine des cellules, des faces internes,
!                ou des faces de bord (sur l'ensemble des processeurs).

!                Ceci permet d'avoir un traitement "generique" qui
!                peut fonctionner pour tous les numeros de maillage,
!                mais si le maillage post est vide a un instant de
!                post traitement donne, on ne saura plus s'il contenait
!                des cellules ou faces. Dans ce cas, il est preferable
!                d'utiliser explicitement le numero du maillage post
!                pour bien determiner s'il doit contenir des cellules
!                ou des faces.

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================
!     1. TRAITEMENT DES MAILLAGES POST A REDEFINIR
!         A RENSEIGNER PAR L'UTILISATEUR aux endroits indiques
!===============================================================================

!     Exemple :
!               pour les maillage post utilisateur, on ne conserve que
!               les mailles auxquelles la vitesse est superieure à
!               un seuil donne.


if (ipart.eq.3) then

  imodif = 1

  ncelps = 0
  nfacps = 0
  nfbrps = 0

  vmin2 = (0.5d0)**2

!       SI LE MAILLAGE POST CONTIENT DES CELLULES
!       -----------------------------------------

  if (itypps(1) .eq. 1) then

    do ii = 1, ncel

      iphas = 1

      v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2        &
           + rtp(ii, iw(iphas))**2
      if (v2 .ge. vmin2) then
        ncelps = ncelps + 1
        lstcel(ncelps) = ii
      endif

    enddo

!       SI LE MAILLAGE POST CONTIENT DES FACES INTERNES
!       -----------------------------------------------

  else if (itypps(2) .eq. 1) then

    do ifac = 1, nfac

      iphas = 1

      ii = ifacel(1, ifac)
      jj = ifacel(2, ifac)

      v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2        &
           + rtp(ii, iw(iphas))**2
      w2 =   rtp(jj, iu(iphas))**2 + rtp(jj, iv(iphas))**2        &
           + rtp(jj, iw(iphas))**2

      if (v2 .ge. vmin2 .or. w2 .ge. vmin2) then
        nfacps = nfacps + 1
        lstfac(nfacps) = ifac
      endif

    enddo

!       SI LE MAILLAGE POST CONTIENT DES FACES DE BORD
!       ----------------------------------------------

  else if (itypps(3) .eq. 1) then

    do ifac = 1, nfabor

      iphas = 1

      ii = ifabor(ifac)

      v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2        &
           + rtp(ii, iw(iphas))**2

      if (v2 .ge. vmin2) then
        nfbrps = nfbrps + 1
        lstfbr(nfbrps) = ifac
      endif

    enddo

  endif

!       Fin du test sur le type de mailles deja existantes

else if (ipart.eq.4) then

  imodif = 1

  ncelps = 0
  nfacps = 0
  nfbrps = 0

  vmin2 = (0.5d0)**2

!       SELECTION DES FACES INTERNES
!       ----------------------------

  do ifac = 1, nfac

    iphas = 1

    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)

    v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2          &
         + rtp(ii, iw(iphas))**2
    w2 =   rtp(jj, iu(iphas))**2 + rtp(jj, iv(iphas))**2          &
         + rtp(jj, iw(iphas))**2

    if (     (v2 .ge. vmin2 .and. w2 .lt. vmin2)                  &
        .or. (v2 .lt. vmin2 .and. w2 .ge. vmin2)) then
      nfacps = nfacps + 1
      lstfac(nfacps) = ifac
    endif

  enddo

!       SELECTION DES FACES DE BORD
!       ---------------------------

  do ifac = 1, nfabor

    iphas = 1

    ii = ifabor(ifac)

    v2 =   rtp(ii, iu(iphas))**2 + rtp(ii, iv(iphas))**2          &
         + rtp(ii, iw(iphas))**2

    if (v2 .ge. vmin2) then
      nfbrps = nfbrps + 1
      lstfbr(nfbrps) = ifac
    endif

  enddo

endif
!     Fin du test sur le numero de maillage post.


return

end
