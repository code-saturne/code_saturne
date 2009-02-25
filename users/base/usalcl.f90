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

subroutine usalcl &
!================

 ( idbia0 , idbra0 , itrale ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itypfb , ialtyb , impale ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl , xyzno0 , depale ,                   &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!    ROUTINE UTILISATEUR
!    REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
!      (IALTYB,ICODCL,RCODCL) POUR LA VITESSE DE MAILLAGE

!    POSSIBILITE DE FIXER DIRECTEMENT LE DEPLACEMENT DES
!      NOEUDS



!    CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE EN ALE
!    ====================================================



! INTRODUCTION
! ============

! On donne ici les conditions aux limites par face de bord, pour
!    pour la vitesse de maillage uniquement (pour les autres
!    variables, se reporter a l'Interface Graphique ou a usclim)


! On balaye les faces de bord et selon un critere on affecte tel
!    ou tel type de conditions aux limites. Dans l'exemple donne
!    ci dessous, c'est la couleur (propriete 1 de la famille)
!    qui permet de distinguer les differents types de bord. On
!    aurait pu aussi travailler avec les coordonnees du centre
!    des faces, mais "c'eut ete moins pratique".


! TYPE DE CONDITIONS AUX LIMITES
! ==============================

! On peut affecter les conditions aux limites de deux manieres.


!    Pour les conditions "standards" :
!    --------------------------------

! Trois types de conditions "standards" sont disponibles. On
!     specifie pour chaque face le type de condition choisi.
!     Le choix se fait en remplissant le tableau IALTYB.


! * IALTYB(IFAC) = IBFIXE : la face IFAC correspond a un bord
!     de maillage fixe. Une condition de Dirichlet nulle sera
!     affectee automatiquement a la vitesse de maillage. En
!     outre, le deplacement des noeuds sera automatiquement
!     impose a 0 (cf. plus base sur IMPALE), sauf si
!     l'utilisateur a modifie la condition pour au moins une
!     des conposantes de la vitesse de maillage (ICODCL modifie,
!     cf. conditions non standards)


! * IALTYB(IFAC) = IGLISS : la face IFAC est une face a maillage
!     glissant. La composante normale a la face de la vitesse
!     de maillage sera forcee a 0, les autres composantes seront
!     traitees en Neumann homogene (pour les faces non alignees
!     avec les axes, la condition de Neumann homogene est
!     seulement partiellement implicitee), comme pour les
!     conditions de symetrie de la vitesse fluide.


! * IALTYB(IFAC) = IVIMPO : la face IFAC a une vitesse de maillage
!     imposee. Cette vitesse est specifiee en remplissant le
!     tableau RCODCL correspondant :
!     RCODCL(IFAC,IUMA,1) = vitesse de maillage selon X
!     RCODCL(IFAC,IVMA,1) = vitesse de maillage selon Y
!     RCODCL(IFAC,IWMA,1) = vitesse de maillage selon Z
!     Les composantes non specifiees de RCODCL(.,I.MA,1) seront mises
!     a 0 (seules les composantes non nulles sont donc a specifier)



!    Pour les conditions "non standards" :
!    ------------------------------------

!     Autres que fixe, glissement ou vitesse imposee, on donne

!      - pour chaque face et chaque composante IVAR=IUMA, IVMA, IWMA :
!        -> un code     ICODCL(IFAC,IVAR)
!        -> trois reels RCODCL(IFAC,IVAR,1)
!                       RCODCL(IFAC,IVAR,2)
!                       RCODCL(IFAC,IVAR,3)
!     La valeur de ICODCL est prise parmi les suivantes :
!       1 : Dirichlet
!       3 : Neumann
!       4 : Symetrie
!     Les valeurs des 3 reels RCODCL sont les suivantes

!      RCODCL(IFAC,IVAR,1) :
!         Dirichlet sur la variable si ICODCL(IFAC,IVAR)=  1
!         La dimension de RCODCL(IFAC,IVAR,1) est en m/s

!      RCODCL(IFAC,IVAR,2) :
!         coefficient d'echange "exterieur" en kg/(m2 s) (entre la
!                    valeur imposee et la valeur au bord du domaine)
!                    RINFIN = infini par defaut
!           RCODCL(IFAC,IVAR,2) =            (VISCMA) / D
!              (D a la dimension d'une distance en m, VISCMA est
!               la viscosite de maillage)
!         NB : la definition de RCODCL(.,.,2) est calquee sur le cas
!              des autres variables standards de l'ecoulement. Il n'a pas
!              grand sens physique dans le cas de la vitesse de maillage.

!      RCODCL(IFAC,IVAR,3) :
!        Densite de flux en kg/(m s2) = J si ICODCL(IFAC,IVAR)= 3
!                       (< 0 si gain, n normale orientee vers l'exterieur)
!           RCODCL(IFAC,IVAR,3) =           -(VISCMA) * (GRAD Um).n
!            (ou Um represente la vitesse de maillage)
!         NB : la definition de RCODCL(.,.,3) est calquee sur le cas
!              des autres variables standards de l'ecoulement.
!              RCODCL(.,.,3) = 0.D0 permet de specifier un Neumann
!              homogene sur la vitesse de maillage. Toute autre valeur
!              aura une signification physique plus discutable.



!      Noter bien que si l'utilisateur affecte une valeur a IALTYB
!       parmi IBFIXE, IGLISS, IVIMPO,
!       et qu'il ne modifie pas ICODCL (valeur nulle par defaut),
!       c'est ITYPFB qui imposera la condition limite.

!      Par contre, si l'utilisateur impose
!        ICODCL(IFAC,IVAR) (non nul),
!        ce sont alors les valeurs de ICODCL et RCODCL qu'il aura
!        fournies qui sont retenues pour la face et la composante
!        consideree (s'il ne precise pas RCODCL, ce sont les valeurs
!        par defaut qui sont retenues pour la face et
!        la variable consideree soit :
!                                 RCODCL(IFAC,IVAR,1) = 0.D0
!                                 RCODCL(IFAC,IVAR,2) = RINFIN
!                                 RCODCL(IFAC,IVAR,3) = 0.D0)


!      Si l'utilisateur specifie de lui-meme la valeur de ICODCL
!         et de RCODCL pour TOUTES les composantes de la vitesse
!         de maillage, il n'est pas necessaire de specifier IALTYB
!         pour la face concernee (sa valeur ne sera pas utilisee)


! REGLE DE COHERENCE
! ==================

!       Si une condition de symetrie (ICODCL=4) a ete imposee sur
!          une des composantes, elle doit l'etre sur toutes les
!          composantes.


! DEPLACEMENT FORCE DES NOEUDS
! ============================

! Pour plus de precision dans le deplacement du maillage, on peut
!   aussi forcer directement le deplacement de certains noeuds,
!   internes ou de bord. Pour cela on remplit les tableaux DEPALE
!   et IMPALE :
!   DEPALE(INOD,1) = deplacement du noeud INOD dans la direction X
!   DEPALE(INOD,2) = deplacement du noeud INOD dans la direction Y
!   DEPALE(INOD,3) = deplacement du noeud INOD dans la direction Z
!   Ce deplacement s'entend comme le deplacement absolu du noeud
!   a partir de sa position dans le maillage d'origine.
!   IMPALE(INOD) = 1 indique que le noeud INOD est a deplacement
!   impose (IMPALE est initialise a 0 ; si sa valeur est laissee
!   nulle, DEPALE ne sera pas pris en compte).

! Lors de la mise a jour du maillage, les noeuds tels que IMPALE=1
!   (internes ou de bord) ne seront pas deplaces a partir de la
!   vitesse de maillage mais directement a partir de DEPALE.
! Dans le cas ou tous les noeuds d'une face sont a deplacement
!   impose, il n'est pas necessaire de remplir les tableaux de
!   conditions aux limites de vitesse de maillage pour cette faces,
!   ils seront ecrases :
!    -> ICODCL sera mis a 1 (Dirichlet)
!    -> RCODCL sera mis a la valeur moyenne des vitesses des noeuds
!       de la face (vitesses calculees a partir de DEPALE)

! Dans le cas de faces specifies a maillage fixes (IALTYB(IFAC)=IBFIXE),
!   tous les noeuds de la face sont automatiquement mis en
!   deplacement impose, avec une valeur nulle de DEPALE.


! CARACTERISTIQUES DES NOEUDS
! ===========================
! Le nombre total de noeuds est stocke dans la variable NNOD.
! Les coordonnees des noeuds sont accessibles par le tableau
!   XYZNOD(3,NNOD).
! Les coordonnees des noeuds dans le maillage initiale sont
!   accessibles par le tableau XYZNO0(3,NNOD).

! Le reperage des noeuds est possible a partir des faces internes et
!   de bord, grace aux tableaux IPNFAC, NODFAC, IPNFBR, NODFBR.
!   NODFAC (NODFBR) contient sequentiellement la liste des noeuds de
!     toutes les faces internes (de bord).
!   IPNFAC (IPNFBR) contient pour chaque face interne (de bord) le
!     numero de la premiere case de NODFAC (NODFBR) lui correspondant.

! Par exemple, pour recuperer sequentiellement tous les noeuds de la
!   face interne IFAC, il suffit d'ecrire la boucle suivante :
!   DO II = IPNFAC(IFAC), IPNFAC(IFAC+1)-1 <- indices des elements
!                                             de NODFAC correspondant
!                                             a IFAC
!     INOD = NODFAC(II)                    <- recuperation du numero
!                                             du IIeme noeud de la face
!                                             IFAC
!     ...
!   ENDDO


! INFLUENCE SUR LES CONDITIONS AUX LIMITES DE VITESSE FLUIDE EN PAROI
! ===================================================================
! Dans le cas de faces de paroi pour le fluide (ITYPFB=IPAROI ou IPARUG),
!  l'influence de la vitesse de maillage depend de sa signification
!  physique.
!  En effet, dans le cas d'une structure immergee par exemple, le
!  mouvement des faces de paroi de la structure correspond à un mouvement
!  physique et doit donc entrainer le fluide.
!  Dans le cas d'un piston au contraire, les bords latéraux sont des
!  parois ou le mouvement des noeuds n'a aucune signification physique et
!  ne doit pas entrainer le fluide.
!  Dans tous les cas, la vitesse de maillage normale a la face est prise
!  en compte (u.n=w.n a la face), c'est le traitement de la vitesse
!  tangentielle qui differe suivant les cas.

! Par defaut, Code_Saturne gere la relation entre vitesse fluide et vitesse
!  de maillage pour les parois de la maniere suivante :
!  - Si IALTYB(IFAC) = IBFIXE, la vitesse de maillage est nulle, il n'y a
!    pas de probleme (et si la paroi est defilante, l'utilisateur le
!    specifiera dans l'interface ou dans usclim).
!  - Si IALTYB(IFAC) = IGLISS, la vitesse tangentielle de maillage n'est
!    pas prise en compte dans les conditions aux limites de paroi pour le
!    fluide (et si la paroi est defilante, l'utilisateur le specifiera
!    dans l'interface ou dans usclim).
!  - Si IALTYB(IFAC) = IVIMPO, la vitesse tangentielle de maillage est
!    prise en compte comme une vitesse de defilement dans les conditions
!    aux limites de paroi pour le fluide, sauf si une vitesse de
!    defilement de paroi a ete specifiee par l'utilisateur dans l'interface
!    ou dans usclim (auquel cas c'est cette vitesse qui est consideree).
!  - Si IMPALE(INOD) = 1 pour tous les noeuds d'une face, la vitesse
!    tangentielle de maillage deduite de ce deplacement sera prise en
!    compte comme une vitesse de defilement dans les conditions
!    aux limites de paroi pour le fluide, sauf si une vitesse de
!    defilement de paroi a ete specifiee par l'utilisateur dans l'interface
!    ou dans usclim (auquel cas c'est cette vitesse qui est consideree).

! Pour les autres types de conditions aux limites pour le fluide (IENTRE et
!  ISOLIB), la vitesse de maillage n'a pas d'influence.

! Dans le cas de conditions aux limites non standards, c'est a l'utilisateur
!  de gerer directement la relation entre les conditions sur la vitesse de
!  maillage et celles sur la vitesse fluide (les conditions aux limites du
!  fluide pouvant etre modifiees dans cette routine).



! IDENTIFICATION DES FACES DE BORD
! ================================
! L'identification des faces de bord concernees se fait grace
! a la commande GETFBR.

!  GETFBR(CHAINE,NLELT,LSTELT) :
!  - CHAINE est une chaine de caractere fournie par l'utilisateur
!    qui donne les criteres de selection
!  - NLTELT est renvoye par la commande. C'est un entier qui
!    correspond au nombre de faces de bord trouveees repondant au
!    critere
!  - LSTELT est renvoye par la commande. C'est un tableau d'entiers
!    de taille NLTELT donnant la liste des faces de bord trouvees
!    repondant au critere.

!  CHAINE peut etre constitue de :
!  - references de couleurs (ex. : 1, 8, 26, ...
!  - references de groupes (ex. : entrees, groupe1, ...)
!  - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
!  Ces criteres peuvent etre combines par des operateurs logiques
!  (AND et OR) et des parentheses
!  ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
!  les faces de bord de couleur 1, appartenant aux groupes 'groupe2'
!  ou 'groupe3' et de coordonnee Y inferieure a 1.


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! itrale           ! e  ! <-- ! numero d'iteration pour l'ale                  !
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
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord pour le fluide          !
!  nphas      )    !    !     !                                                !
! ialtyb(nfabor    ! te ! --> ! type des faces de bord pour la                 !
!                  ! te !     !                  vitesse de maillage           !
! impale(nnod)     ! te ! <-- ! indicateur de delacement impose                !
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
!                  !    !     !  flux (negatif si gain) w/m2                   !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! depale(nnod,3    ! tr ! <-- ! deplacement aux noeuds                         !
! xyzno0(3,nnod    ! tr ! <-- ! coordonnees noeuds maillage initial            !
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
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "period.h"
include "ihmpre.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0 , itrale
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itypfb(nfabor,nphas), ialtyb(nfabor)
integer          impale(nnod)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision depale(nnod,3), xyzno0(3,nnod)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ifac, iel, ii
integer          inod
integer          ilelt, nlelt

double precision delta, deltaa

!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
!       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
!       LA BIBLIOTHEQUE)
!===============================================================================

write(nfecra,9000)
call csexit (1)

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@     LA METHODE ALE A ETE ENCLENCHEE                        ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usalcl DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================
! 1.  INITIALISATIONS

!===============================================================================

idebia = idbia0
idebra = idbra0


!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR LES FACES DE BORD
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LA CONDITION LIMITE

!          IMPOSER ICI LES CONDITIONS LIMITES SUR LES FACES DE BORD

!===============================================================================

!     Calcul du deplacement au temps courant
deltaa = sin(3.141596d0*(ntcabs-1)/50.d0)
delta  = sin(3.141596d0*ntcabs/50.d0)

! --- On impose en couleur 4 une vitesse de deplacement des faces

CALL GETFBR('4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!     ELEMENT ADJACENT A LA FACE DE BORD
  iel = ifabor(ifac)

  ialtyb(ifac) = ivimpo
  rcodcl(ifac,iuma,1) = 0.d0
  rcodcl(ifac,ivma,1) = 0.d0
  rcodcl(ifac,iwma,1) = (delta-deltaa)/dt(iel)

enddo

! --- On impose en couleur 5 un deplacement des noeuds

CALL GETFBR('5',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
    inod = nodfbr(ii)
    if (impale(inod).eq.0) then
      depale(inod,1) = 0.d0
      depale(inod,2) = 0.d0
      depale(inod,3) = delta
      impale(inod) = 1
    endif
  enddo

enddo

! --- On impose en couleur 6 un glissement

CALL GETFBR('6',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ialtyb(ifac) = igliss

enddo

! --- On impose ailleurs une condition de bord fixe

CALL GETFBR( 'not (4 or 5 or 6)',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  ialtyb(ifac) = ibfixe

enddo

!----
! FORMATS
!----

!----
! FIN
!----

return
end
