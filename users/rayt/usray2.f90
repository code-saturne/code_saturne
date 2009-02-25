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

                  subroutine usray2                               &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , iphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , itypfb ,          &
   maxelt , lstelt ,                                              &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , izfrdp , isothp ,                                     &
   tmin   , tmax   , tx     ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb , rcodcl , &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser ,                                              &
   tparop , qincid , hfcnvp , flcnvp ,                            &
   xlamp  , epap   , epsp   , textp  , tintp  ,                   &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
!   -----------------------------------------

!   CARACTERISTIQUES ET COMPORTEMENT RADIATIF DES FACES DE PAROI

!   Ce sous-programme est appele autant de fois qu'il y a de phases
!   pour lesquelles il faut faire un calcul de rayonnement
!   semi-transparent

!   IPHAS contient le numero absolu de la phase concernee



!  EXTERIEUR       PAROI       DOMAINE FLUIDE
!            |               |
!            |               |
!            |               |
!            | EPSP Qincid<--|<----- Qincid             [flux incident]
!            |               |
!            | [absorption]  |
!            |               |
!            |               |-----> (1-EPSP) Qincid    [reflexion]
!            |               |
!            |               |
!            |               |                       4
!            |               |-----> EPSP sigma Tparop  [emission]
!            |               |
!            |               |
!            |               |
!            | [conduction]  |     .....o Tfluide       [convection]
!            |               |    .
!            |     XLAMP     |   .    Hfluide
!            |               |  .
!            |            ...o..
!            |          ..   | Tparop (temperature en paroi)
!            |         .     |
!            |       ..      |
!            |      .        |
!            |    ..         |
!            |   .           |
!            o...            |
!       Textp|               |
!            |               |
!            |<------------->|
!            |     EPAP      |
!            |               |


!  LE CALCUL DES TEMPERATURES EST REALISEE PAR UN BILAN DE FLUX
!  COMME SUIT :

!  Q           =  Q           + (Qrayt           - Qrayt        )
!   conduction     convection         absorption        emission

!  SOIT (flux positif si sortant du domaine de calcul) :

!   XLAMP
!   -----(Tparop-Textp) =
!   EPAP
!                                                             4
!         Hfluide (Tfluide-Tparop) + EPSP (QINCID - SIGMA Tparop )


!      CORPS                        EPSP
!      ------------------------------------
!      acier poli                   0,06
!      acier oxydé                  0,80
!      aceir rugueux                0,94
!      aluminium poli               0,04
!      aluminium oxydé (intérieur)  0,09
!      aluminium oxydé (air humide) 0,90
!      brique                       0,93
!      béton                        0,93
!      papier                     0,8 à 0,9
!      eau                          0,96


! IDENTIFICATION DES FACES DE PAROI
! =================================
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
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
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
! icodcl           ! te ! <-- ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! izfrdp(nfabor    ! te ! --> ! numero de zone pour les faces de bord          !
! isothp(nfabor    ! te ! --> ! liste des frontieres isothermes                !
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
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! w1...6(ncelet    ! tr ! --- ! tableau de travail                             !
! tparop(nfabor    ! tr ! <-- ! temperature de paroi en kelvin                 !
! qincid(nfabor    ! tr ! <-- ! densite de flux radiatif aux bords             !
! hfcnvp(nfabor    ! tr ! <-- ! coefficient d'echange fluide aux               !
!                  !    !     ! faces de bord                                  !
! flcnvp(nfabor    ! tr ! <-- ! densite de flux convectif aux faces            !
!                  !    !     ! de bord                                        !
! xlamp(nfabor)    ! tr ! --> ! coefficient de conductivite thermique          !
!                  !    !     ! des facettes de paroi (w/m/k)                  !
! epap(nfabor)     ! tr ! --> ! epaisseur des facettes de paroi (m)            !
! epsp (nfabor)    ! tr ! --> ! emissivite des facettes de bord                !
! textp(nfabor)    ! tr ! --> ! temperature de bord externe                    !
!                  !    !     ! en kelvin                                      !
! tintp(nfabor)    ! tr ! --> ! temperature de bord interne                    !
!                  !    !     ! en kelvin                                      !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "pointe.h"
include "parall.h"
include "period.h"
include "radiat.h"
include "ihmpre.h"


!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , iphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml) , itypfb(nfabor,*)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse), ia(*)

integer          icodcl(nfabor,nvar)
integer          izfrdp(nfabor), isothp(nfabor)

double precision tmin , tmax , tx

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision coefa(nfabor,*), coefb(nfabor,*)

double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)

double precision tparop(nfabor), qincid(nfabor)
double precision hfcnvp(nfabor),flcnvp(nfabor)
double precision xlamp(nfabor), epap(nfabor)
double precision epsp(nfabor)
double precision textp(nfabor), tintp(nfabor)

double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! VARIABLES LOCALES

integer          idebia , idebra
integer          ifac , ivar, iok
integer          ilelt, nlelt

!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if (iihmpr.eq.1) then
  return
else
  write(nfecra,9000)
  call csexit (1)
  !==========
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usray2 DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
!===============================================================================
! 0. GESTION MEMOIRE
!       Aucune modification requise
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 1. NUMERO DE LA VARIABLE THERMIQUE

!       Noter que TKELVI, disponible ici, contient 273.16D0
!       Noter que SIG   , disponible ici, contient 5.6703D-8
!===============================================================================

!     IVAR est le numero de la variable thermique.

ivar = isca(iscalt(iphas))


!===============================================================================
!  2. VALEURS MIN ET MAX ADMISSIBLES POUR LA TEMPERATURE DE PAROI
!                            (clipping)



!                            EN KELVIN


!        (quel que soit le choix fait pour la variable thermique)




!   TMIN et TMAX sont respectivement les valeurs minimale et maximale
!     autorisees pour les temperatures de paroi. A chaque nouveau
!     calcul de temperature de paroi, on verifie qu'il n'y a pas
!     depassement, sinon il y a clipping.

!   L'exemple ci-dessous correspond a un clipping tres peu contraignant.
!     (valeurs par defaut)

!===============================================================================

tmin = 0.d0
tmax = grand + tkelvi

!===============================================================================
! 3. ZONES FRONTIERE ET PROPRIETES DES PAROIS
!===============================================================================


!     DEFINITION DES ZONES FRONTIERES
!     ===============================

!     On definit des zones de faces de paroi, et on leur affecte un type.
!       Ceci permet d'appliquer les conditions aux limites et de realiser
!       des bilans en traitant separement les differents sous ensembles
!       ainsi constitues (on peut ainsi connaitre par exemple le flux au
!       travers des differentes zones definies par l'utilisateur).

!     Pour chaque face de bord IFAC (pas uniquement les faces de paroi),
!       l'utilisateur definit selon son propre choix un numero de zone
!       frontiere IZFRDP(IFAC) a partir
!         de la couleur des faces de bord,
!         ou, plus generalement, de leurs proprietes (couleurs, groupes...),
!         ou des conditions aux limites fixees dans usclim,
!         ou meme de leur coordonnees.
!     Attention : il est indispensable que TOUTES les faces de bord
!       aient ete affectees a une zone.
!     Le numero des zones (la valeur de IZFRDP(IFAC)) est
!       arbitrairement choisi par l'utilisateur, mais doit etre un
!       entier strictement positif et inferieur ou egal a NBZRDM
!       (valeur fixee en parametre dans radiat.h, a ne pas modifier).



!     PROPRIETES DES PAROIS
!     =====================


!      ATTENTION :
!      ---------

!             L'unite des temperatures est le KELVIN

!        quel que soit le choix de resolution fait par ailleurs



!      DONNEES OBLIGATOIRES :
!      --------------------

!      ISOTHP = ITPIMP
!                repere les facettes de paroi a temperature imposee
!                  ("Indicateur Temperature Paroi IMPosee" ou "isotherme")
!                  la temperature a imposer est TINTP(Kelvin)
!             = IPGRNO ou IPREFL
!                repere une facette de paroi dont il faut calculer
!                  la temperature (un bilan de flux sera fait).
!                  IPGRNO : Indicateur de Paroi GRise ou NOire (EPSP > 0)
!                  IPREFL : Indicateur de Paroi REFLechissante (EPSP = 0)
!             = IFGRNO ou IFREFL
!                repere une facette de paroi pour laquelle on impose un
!                  flux (conductif en paroi ou total sur le fluide)
!                  IFGRNO : Indicateur de Flux en paroi GRise ou NOire
!                  IFREFL : Indicateur de Flux en paroi REFLechissante

!      TINTP = Temperature de peau interne (KELVIN)
!               sert a l'initialisation de TPAROP en debut de calcul
!               (lorsqu'il n'y a pas de fichier suite rayonnement)
!                 Si ISOTHP = ITPIMP, la valeur de TPAROP est egale a
!                   celle imposee dans TINTP a chaque pas de temps
!                 Dans les autres cas, TINTP ne sert qu'a l'initialisation.


!      RENSEIGNEMENTS COMLEMENTAIRES A FOURNIR SUIVANT ISOTHP (VOIR EXEMPLES)
!      --------------------------------------------------------------------

!      RCODCL = valeur de flux (voir USCLIM pour la definition)

!      EPSP   = Emissivite aux facettes de paroi (dans l'intervalle [0;1])

!      XLAMP  = Conductivite thermique des parois (W/m/K)
!               > 0 quand il est renseigne

!      EPAP   = Epaisseur des parois (m)
!               > 0 quand il est renseigne

!      TEXTP  = Temperature de peau externe en KELVIN
!               > 0 quand il est renseigne



!     EXEMPLE
!     =======

!      Attention a la coherence avec USCLIM

!      Dans cet EXEMPLE,
!               -------

!        Les faces de paroi (IPAROI et IPARUG), sont reparties en
!          5 sous ensembles (zones) reperes par IFRFAC(IFAC) variant de
!          51 a 55 (arbitraire) a chacun desquels on applique une
!          condition (ISOTHP) differente.
!        Les autres faces de bord (entree, sortie, sous-ensembles
!          arbitraires sur lesquels on souhaite voir imprimer des
!          informations) sont repartis en zones reperees
!          par IFRFAC(IFAC), dont la valeur peut etre arbitrairement
!          choisie (>0 et <NBZRDM+1).

!        Remarque : une maniere simple (mais non generale) d'affecter
!          rapidement des valeurs a IFRFAC est IFRFAC(IFAC) = ICOUL,
!          ou ICOUL est la couleur de la face.










!             DANS TOUS LES CAS, IL EST INTERDIT
!                                       ========

!                       DE MODIFIER TPAROP ou QINCID ICI.








!    Indicateur : nombre de faces oubliees.
iok = 0

!   -------------------------------------------------------------------
!-->  Exemple 1 :
!      Pour les faces PAROI de couleur 1 :
!           Profil de temperature imposee
!       ------------------------------------

CALL GETFBR('1',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac,iphas).eq.iparoi ) then

!      Numero de zone
    izfrdp(ifac) = 51

!      Type de condition : Paroi grise ou noire et
!                          Profil de temperature imposee TINTP
    isothp(ifac) = itpimp

!      Donnees complementaires necessaires et suffisantes
!        Emissivite
    epsp  (ifac) = 0.1d0
!        Temperature de paroi imposee ("Profil" simple ici = 473.16 K)
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo

!   -------------------------------------------------------------------
!-->  Exemple 2 :
!      Pour les faces PAROI RUGUEUSE de couleur 2 :
!           Paroi grise ou noire et
!           Profil de temperature imposee TEXTP
!       ------------------------------------

CALL GETFBR('2',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac,iphas).eq.iparug ) then

!      Numero de zone
    izfrdp(ifac) = 52

!      Type de condition : Paroi grise ou noire et
!                          Profil de temperature imposee TEXTP
    isothp(ifac) = ipgrno

!      Donnees complementaires necessaires et suffisantes
!        Emissivite (>0)
    epsp  (ifac) = 0.9d0
!        Conductivite (W/m/K)
    xlamp (ifac) = 3.0d0
!        Epaisseur    (m)
    epap  (ifac) = 0.1d0
!        Temperature externe imposee : 473.16 K
    textp (ifac) = 200.d0 + tkelvi
!        Temperature de paroi initiale : 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo

!   -------------------------------------------------------------------
!-->  Exemple 3 :
!      Pour les faces PAROI de couleur 3 :
!           Paroi reflechissante (EPSP = 0) et
!           Profil de temperature imposee TEXTP
!       ------------------------------------

CALL GETFBR('3',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac,iphas).eq.iparoi ) then

!      Numero de zone
    izfrdp(ifac) = 53

!      Type de condition : Paroi reflechissante (EPSP = 0) et
!                          Profil de temperature imposee TEXTP
    isothp(ifac) = iprefl

!      Donnees complementaires necessaires et suffisantes
!        Conductivite (W/m/K)
    xlamp (ifac) = 3.0d0
!        Epaisseur    (m)
    epap  (ifac) = 0.1d0
!        Temperature externe imposee : 473.16 K
    textp (ifac) = 200.d0 + tkelvi
!        Temperature de paroi initiale : 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo

!   -------------------------------------------------------------------
!-->  Exemple 4 :
!      Pour les faces PAROI de couleur 4 :
!           Paroi grise ou noire et
!           Flux de conduction impose dans la paroi

!        XLAMP
!        -----(Tparop-Textp) = Flux de conduction impose (W/m2)
!        EPAP
!                         = RODCL(IFAC,IVAR,3)

!       Si le flux de conduction est nul la paroi est adiabatique
!       Le tableau RCODCL(IFAC,IVAR,3) recoit la valeur du flux.
!       Densite de flux (< 0 si gain pour le fluide)
!       Reprise de USCLIM :
!         Pour les temperatures T,         en Watt/m2       :
!            RCODCL(IFAC,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
!         Pour les enthalpies H,           en Watt/m2       :
!            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
!       ------------------------------------

CALL GETFBR('4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac,iphas).eq.iparoi ) then

!      Numero de zone
    izfrdp(ifac) = 54

!      Type de condition : Paroi grise ou noire et
!                          Flux de conduction impose dans la paroi
    isothp(ifac) = ifgrno

!      Donnees complementaires necessaires et suffisantes
!        Emissivite (>0)
    epsp  (ifac) = 0.9d0
!        Flux de conduction (W/m2)
    rcodcl(ifac,ivar,3) = 0.d0
!        Temperature de paroi initiale : 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo

!   -------------------------------------------------------------------
!-->  Exemple 5 :
!      Pour les faces PAROI de couleur 5 :
!           Paroi reflechissante (EPSP = 0) et
!           Flux de conduction impose dans la paroi

!      Equivalent a imposer une condition de flux au fluide

!        XLAMP
!        -----(Tparop-Textp) = Flux de conduction impose ET EPSP = 0
!        EPAP
!                         = RODCL(IFAC,IVAR,3)

!       Si le flux de conduction est nul la paroi est adiabatique
!       Le tableau RCODCL(IFAC,IVAR,3) contient la valeur du flux (W/m2)
!       imposee comme condition au limite pour le fluide
!       Densite de flux (< 0 si gain pour le fluide)
!       Reprise de USCLIM :
!         Pour les temperatures T,         en Watt/m2       :
!            RCODCL(IFAC,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
!         Pour les enthalpies H,           en Watt/m2       :
!            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
!       ------------------------------------

CALL GETFBR('5',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  if ( itypfb(ifac,iphas).eq.iparoi ) then

!      Numero de zone
    izfrdp(ifac) = 55

!      Type de condition : Paroi reflechissante (EPSP = 0) et
!                          Flux de conduction impose dans la paroi
    isothp(ifac) = ifrefl

!      Donnees complementaires necessaires et suffisantes
!        Flux de conduction (W/m2)
    rcodcl(ifac,ivar,3) = 0.d0
!        Temperature de paroi initiale : 473.16 K
    tintp (ifac) = 200.d0 + tkelvi

  endif

enddo







!                           ATTENTION











!   -------------------------------------------------------------------
!-->   Pour les faces non paroi : reperage obligatoire de la zone
!        Les informations (flux par exemple) seront donnees zone par zone
!        et il est donc conseille de regrouper entre elles les faces
!        de meme nature.

!                             PAR SECURITE
!                        DANS CE SOUS-PROGRAMME
!          ON DEMANDE A L'UTILISATEUR DE REPERER TOUTES LES FACES
!                             PAROI OU NON

!        (il est indispensable d'affecter une valeur a tous les
!            IZFRDP(IFAC), IFAC = 1, NFABOR)




!       ------------------------------------

do ifac = 1, nfabor

  if     ( itypfb(ifac,iphas).eq.isolib                  ) then
    izfrdp(ifac) = 61
  elseif ( itypfb(ifac,iphas).eq.ientre.and.                      &
           cdgfbo(2,ifac)    .gt.0.d0                    ) then
    izfrdp(ifac) = 62
  elseif ( itypfb(ifac,iphas).eq.ientre.and.                      &
           cdgfbo(2,ifac)    .le.0.d0                    ) then
    izfrdp(ifac) = 63
  elseif ( itypfb(ifac,iphas).eq.isymet                  ) then
    izfrdp(ifac) = 64



!   -------------------------------------------------------------------
!-->  Exemple 7 :
!      Verification que toutes les faces de paroi ont ete vues
!       Cette precaution est souhaitable dans tous les cas.
!       ------------------------------------

  elseif ( itypfb(ifac,iphas).eq.iparoi .or.                      &
           itypfb(ifac,iphas).eq.iparug     ) then
    if (izfrdp(ifac) .eq. -1) then
      write(nfecra,1000)ifac
      iok = iok + 1
    endif
  endif




!     Fin de la boucle sur les faces de bord
!     --------------------------------------

enddo

! Stop si des faces ont ete oubliees
if(iok.ne.0) then
  call csexit (1)
  !==========
endif

! -------
! FORMAT
! -------


 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     LES DONNEES RAYONNEMENT NE SONT PAS RENSEIGNEES POUR   ',/,&
'@       LA FACE DE PAROI ',I10                                ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usray2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
! ---
! FIN
! ---


return

end
