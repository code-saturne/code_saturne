!-------------------------------------------------------------------------------

!                      Code_Saturne version 2.0.0-beta1
!                      --------------------------

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

subroutine usvpst &
!================

 ( idbia0 , idbra0 , ipart  ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   ncelps , nfacps , nfbrps ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
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

! ROUTINE UTILISATEUR POUR LA SORTIE DE VARIABLES SUR UN MAILLAGE
!   DE POST TRAITEMENT DEJA DEFINI

! PAR DEFAUT, DEUX MAILLAGES SONT DEFINIS AUTOMATIQUEMENT :
!      - LE MAILLAGE VOLUMIQUE (IPART=-1) SELON LA VALEUR DE ICHRVL
!      - LE MAILLAGE DE BORD   (IPART=-2) SELON LA VALEUR DE ICHRBO
! DES MAILLAGES SUPPLEMENTAIRES (CELLULES OU FACES INTERNES ET
!   DE BORD) PEUVENT ETRE DEFINIS ET PARAMETRES A TRAVERS
!   usdpst    .F ET usmpst.F.

!ETTE  CETTE ROUTINE EST APPELEE UNE FOIS PAR MAILLAGE POST
!   ET PAR PAS DE TEMPS AUQUEL CE MAILLAGE EST ACTIF

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
include "cstnum.h"
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
integer          nideve , nrdeve , nituse , nrtuse

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

character*32     namevr

integer          ntindp
integer          iel   , ifac  , iloc  , iphas, ivar , iclt
integer          idimt , ii   , jj
integer          idimte, itenso, ientla, ivarpr
integer          imom1, imom2, ipcmo1, ipcmo2, idtcm
double precision pond
double precision rbid(1)

integer          ipass
data             ipass /0/
save             ipass


!===============================================================================



!===============================================================================
!     1. TRAITEMENT DES VARIABLES A SORTIR
!         A RENSEIGNER PAR L'UTILISATEUR aux endroits indiques
!===============================================================================


!     Un maillage de posttraitement est une "part"
!       (au sens EnSight ; les équivalents MED et CGNS sont le maillage
!       et la base respectivement)
!     L'utilisateur aura defini ses maillages de posttraitement dans
!       usdpst (NBPART maillages de posttraitement)


!     La routine est appelee une fois pour chaque maillage IPART.
!     Pour chaque maillage et pour chacune des variables que l'on
!       souhaite posttraiter, on doit definir certains parametres et les
!       passer a la routine PSTEVA qui se charge de l'ecriture effective.
!     Ces parametres sont :
!       NAMEVR : nom de la variable
!       IDIMT  : dimension de la variable
!       IENTLA : dans le cas ou IDIMT est >1, IENTLA permet de specifier
!                si le tableau contenant la variable est range de maniere
!                entrelacee X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,... (IENTLA=1)
!                ou non entrelancee X1,X2,X3,...,Y1,Y2,Y3,...,Z1,Z2,Z3,...
!                                                               (IENTLA=0)
!       IVARPR : specifie si le tableau contenant la variable traitee est
!                defini sur le maillage "parent"/
!                En effet, meme si le maillage IPART considere contient
!                le meme nombre d'elements que le maillage complet "parent"
!                (NCELPS=NCEL), l'ordre de numerotation des elements n'est pas
!                forcement le meme. Le tableau TRACEL passe en argument de
!                PSTEVA est construit selon la numerotation du maillage IPART.
!                Pour posttraiter une variable contenue dans un tableau RTUSER
!                par exemple, il faut donc d'abord le reordonner dans le
!                tableau TRACEL :
!                  DO ILOC = 1, NCELPS
!                    IEL = LSTCEL(ILOC)
!                    TRACEL(ILOC) = RTUSER(IEL)
!                  ENDDO
!                Une alternative est cependant offerte, pour eviter des copies
!                inutiles. Si NCELPS=NCEL, on peut directement passer RTUSER
!                en argument de PSTEVA, en specifiant IVARPR=1, pour avertir
!                le code que la numerotation est celle du maillage "parent".
!                L'exemple ci-dessus concerne les cellules, mais la
!                problematique est la meme avec les faces internes ou faces de
!                bord.


!     Remarque : attention aux longueurs des noms de variables.

!                On autorise ici jusqu'à 32 caracteres , mais selon le
!                format utilise, les noms peuvent etre tronques :

!                  - a 19 caracteres pour une sortie EnSight
!                  - a 32 caracteres pour ue sortie MED 2.2

!                La longueur du nom n'est pas limitee en interne, et en
!                cas de deux variables aux noms tronques ne differant
!                qu'apres le 19ieme caractere, les lignes correspondantes
!                apparaitront normalement dans le fichier texte ".case"
!                EnSight, avec un meme champ de description ; il suffit
!                alors de renommer un de ces champs dans ce fichier
!                texte pour corriger le probleme.

!                Les caracteres blancs en debut ou fin de chaine sont
!                supprimes automatiquement. Selon le format utilise,
!                les caracteres interdits (sous EnSight, les caracteres
!                (  ) ] [ + - @           ! # * ^ $ / ainsi que les blancs et les
!                tabulations seront remplaces par le caractere ___________.
!                Ceci ne necessite aucune intervention utilisateur
!                (i.e. inutile d'utiliser ici les anciennes fonctions
!                Fortran VERLON et UNDSCR).


!     Exemples :
!               pour le maillage post 1, on sort
!                 la temperature interpolee a la face
!                   (et 0 sur les eventuelles faces interieures)

!               pour le maillage post 2, on sort
!                 la temperature


if (ipart.eq.1 ) then


!       Initialisation
!         pas d'intervention utilisateur requise
  do ii = 1, 32
    NAMEVR (II:II) = ' '
  enddo

!       Nom de la variable
!         a renseigner par l'utilisateur
  NAMEVR = 'Temperature interpolee'

!       Dimension de la variable (3 = vecteur, 1=scalaire)
!         a renseigner par l'utilisateur
  idimt = 1

!       Valeurs entrelacées
  ientla = 0


!       Calcul des valeurs de la variable sur les faces internes
!         Pour simplifier l'exemple, on se contente ici d'une simple
!           interpolation lineaire.
!         Dans les calculs paralleles, si l'on utilise les voisins,
!           il faut faire un echange au prealable, comme d'habitude.
!         Dans les calculs avec periodicite, il faut egalement le faire
!         Pour les calculs ou periodicite et parallelisme coexistent,
!           l'appel a ces routines doit etre fait dans l'ordre
!           PARCOM puis PERCOM

!         a renseigner par l'utilisateur

  if(irangp.ge.0) then
    call parcom (rtp(1,isca(1)))
  endif

  if(iperio.eq.1) then
    ivar = isca(1)
    idimte = 0
    itenso = 0
    call percom                                                   &
    !==========
         ( idimte , itenso ,                                      &
           rtp(1,ivar), rtp(1,ivar), rtp(1,ivar),                 &
           rtp(1,ivar), rtp(1,ivar), rtp(1,ivar),                 &
           rtp(1,ivar), rtp(1,ivar), rtp(1,ivar))
  endif

  do iloc = 1, nfacps

    ifac = lstfac(iloc)
    ii = ifacel(1, ifac)
    jj = ifacel(2, ifac)
    pond = ra(ipond-1+ifac)

    trafac(iloc)                                                  &
         =         pond  * rtp(ii, isca(1))                       &
         + (1.d0 - pond) * rtp(jj, isca(1))

  enddo


!       Ecriture effective des valeurs calculees

  ivarpr = 0

  call psteva(ipart , namevr, idimt, ientla, ivarpr,              &
  !==========
              ntcabs, ttcabs, tracel, trafac, trafbr)




else if  (ipart.eq.2) then


!       1.4.1 TRAITEMENT DE LA VITESSE
!       ------------------------------

!       Initialisation
!         pas d'intervention utilisateur requise
  do ii = 1, 32
    NAMEVR (II:II) = ' '
  enddo

!       Nom de la variable
!         a renseigner par l'utilisateur
  NAMEVR = 'Temperature'

!       Dimension de la variable (3 = vecteur, 1=scalaire)
!         a renseigner par l'utilisateur
  idimt = 1

!       Valeurs non entrelacées
  ientla = 0

  do iloc = 1, ncelps
    iel = lstcel(iloc)
    tracel(iloc) = rtp(iel,isca(1))
  enddo

  ivarpr = 0

  call psteva(ipart , namevr, idimt, ientla, ivarpr,              &
  !==========
              ntcabs, ttcabs, tracel, trafac, trafbr)

endif
!     Fin du test sur le numero de maillage post.


return

end
