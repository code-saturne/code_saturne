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

subroutine uslag2                                                 &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   itypfb , itrifb , itepa  , ifrlag ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   ,                                              &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    SOUS-PROGRAMME UTILISATEUR (INTERVENTION OBLIGATOIRE)

!    ROUTINE UTILISATEUR POUR LES CONDITIONS AUX LIMITES RELATIVES
!      AUX PARTICULES (ENTREE ET TRAITEMENT AUX AUTRES BORDS)


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
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! e  ! <-- ! nombre d'elements ayant au moins une           !
!                  !    !     ! face de bord                                   !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
!                  !    !     ! le module lagrangien                           !
! ifacel           ! te ! <-- ! elements voisins d'une face interne            !
! (2, nfac)        !    !     !                                                !
! ifabor           ! te ! <-- ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! <-- ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! <-- ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
!  (nfml,nprfml    !    !     !                                                !
! maxelt           !  e ! <-- ! nb max d'elements (cell,fac,fbr)               !
! lstelt(maxelt) te ! --- ! tableau de travail                             !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac                      !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr                      !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !                                                !
! itrifb(nfabor    ! te ! <-- ! indirection pour tri des faces de brd          !
!  nphas      )    !    !     !                                                !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
! ifrlag(nfabor    ! te ! --> ! type des faces de bord lagrangien              !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
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
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet          !    !     !                                                !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec)                     !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
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
include "cstnum.h"
include "cstphy.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"
include "ppppar.h"
include "ppthch.h"
include "cpincl.h"
include "ihmpre.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ifac , izone, nbclas, iclas
integer          icha
integer          ilelt, nlelt

double precision pis6 , mp0 , temp

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
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
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@     MODULE LAGRANGIEN :                                    ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uslag2 DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2. INITIALISATION
!===============================================================================

pis6 = pi / 6.d0

!===============================================================================
! 3. CONSTRUCTION DES ZONES FRONTIERE
!===============================================================================


!     DEFINITION DES ZONES FRONTIERES

!     Pour le module Lagrangien, l'utilisateur definit NFRLAG zones
!       frontieres a partir de la couleur des faces de bord, plus
!       generalement de leurs proprietes (couleurs, groupes...),
!       ou des conditions aux limites fixees dans usclim
!       ou meme de leur coordonnees. Pour ce faire, on renseigne
!       le tableau IFRLAG(NFABOR) qui donne pour
!       chaque face de bord IFAC le numero de la zone a laquelle elle
!       appartient IFRLAG(IFAC).
!     Attention : il est indispensable que TOUTES les faces aient ete
!       affectees a une zone.
!     Le numero des zones (donc les valeurs de IFRLAG(IFAC)) est
!       arbitrairement choisi par l'utilisateur, mais doit etre un
!       entier strictement positif et inferieur ou egal a NFLAGM
!       (valeur fixee en parametre dans lagpar.h).
!     On affecte ensuite a chaque zone un type ITYLAG qui sera utilise
!       pour imposer des conditions aux limites globales.


izone = -1

! ---> Premiere zone, numerotee IZONE = 1 ( = couleur 10)
CALL GETFBR('10',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 1
  ifrlag(ifac) = izone

enddo

! ---> Deuxieme zone, numerotee IZONE = 2 ( = partie des couleur 4)
CALL GETFBR('4 and Y < 1.0',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 2
  ifrlag(ifac) = izone

enddo

! ---> Troisieme zone, numerotee IZONE = 4 ( = entree phase 1)
do ifac = 1, nfabor
  if(itypfb(ifac,1).eq.ientre) then
    izone        = 4
    ifrlag(ifac) = izone
  endif
enddo

! ---> Nieme zone, numerotee IZONE = 5 ( = couleur 3)
CALL GETFBR('3',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  izone        = 5
  ifrlag(ifac) = izone

enddo


!===============================================================================
! 4. INJECTION PAR CLASSES DE PARTICULES DANS LE DOMAINE DE CALCUL
!===============================================================================

!   POUR DONNER LES INFORMATIONS SUR LES CLASSES DE PARTICULES,
!   ON PROCEDE EN DEUX ETAPES :

!   1) TOUT D'ABORD, ON RENSEIGNE LE NOMBRE DE CLASSES DE PARTICULES
!      PAR ZONE FRONTIERE : IUSNCL (ce nombre est nul par defaut)

!   2) ENSUITE PAR ZONE DE COULEUR ET PAR CLASSE,
!      ON DONNE LES CARACTERISTIQUES D'ENTREE ET
!      LES CARACTERISTIQUE PHYSIQUES DES PARTICULES




! --> Nombre de classes de particules entrantes
!   On renseigne ici le nombre de classes pour toutes les zones
!     identifiees precedemment.
!   Par defaut le nombre de classes de particules est nul.
!   Le nombre max de classes est NCLAGM donne dans lagpar

! ---> Premiere zone, numerotee IZONE = 1 : 1 classe injectee
izone     = 1
nbclas    = 1
iusncl(izone) = nbclas

! ---> Deuxieme zone, numerotee IZONE = 2 : 0 classe injectee
izone     = 2
nbclas    = 0
iusncl(izone) = nbclas

! ---> Troisieme zone, numerotee IZONE = 4 : 0 classe injectee
izone     = 4
nbclas    = 0
iusncl(izone) = nbclas

! ---> Nieme zone,     numerotee IZONE = 5 : 0 classe injectee
izone     = 5
nbclas    = 0
iusncl(izone) = nbclas


! --> Pour chaque classe de particules associee a une zone,
!     il faut fournir des informations.


!     IUSNCL nbr de classes par zones
!     IUSCLB conditions au bord pour les particules
!     = IENTRL -> zone d'injection de particules
!     = ISORTL -> sortie du domaine
!     = IREBOL -> rebond des particules
!     = IDEPO1 -> deposition definitive
!     = IDEPO2 -> deposition definitive mais la particule reste en
!                 memoire (utile si IENSI2 = 1 uniquement)
!     = IDEPO3 -> deposition et remise en suspension possible
!                 suivant les condition de l'ecoulement
!     = IDEPFA -> deposition de la particule avec force d'attachement,
!                 vitesse conservee et re-entrainement possible
!                 (Possible si LADLVO = 1 )
!     = IENCRL -> encrassement (Charbon uniquement IPHYLA = 2)
!     = JBORD1 -> interaction part/frontiere utilisateur (cf. USLABO)
!     = JBORD2 -> interaction part/frontiere utilisateur (cf. USLABO)
!     = JBORD3 -> interaction part/frontiere utilisateur (cf. USLABO)
!     = JBORD4 -> interaction part/frontiere utilisateur (cf. USLABO)
!     = JBORD5 -> interaction part/frontiere utilisateur (cf. USLABO)



!     Tableau IUSLAG :
!     ================
!        IJNBP : nbr de part par classe et zones
!        IJFRE : frequence d'injection, si la frequence est nulle alors
!                  il y a injection uniquement a la 1ere iteration absolue
!        ICLST : numero de groupe auquel appartient la particule
!                 (uniquement si on souhaite des statistiques par groupe)
!        IJUVW : type de condition vitesse
!                  = -1 vitesse fluide imposee
!                  =  0 vitesse imposee selon la direction normale a la face
!                        de bord et de norme RUSLAG(ICLAS,IZONE,IUNO)
!                  =  1 vitesse imposee : on donne RUSLAG(ICLAS,IZONE,IUPT)
!                                                  RUSLAG(ICLAS,IZONE,IVPT)
!                                                  RUSLAG(ICLAS,IZONE,IWPT)
!                  =  2 profil utilisateur
!        IJPRTP : type de condition temperature
!                  =  1 temperature imposee : on donne RUSLAG(ICLAS,IZONE,ITPT)
!                  =  2 profil utilisateur
!        IJPRDP : type de condition diametre
!                  =  1 vitesse imposee : on donne RUSLAG(ICLAS,IZONE,IDPT)
!                                                  RUSLAG(ICLAS,IZONE,IVDPT)
!                  =  2 profil utilisateur
!        INUCHL : numero du charbon de la particule (uniquement si IPHYLA=2)

!     Tableau RUSLAG :
!     ================
!        IUNO  : norme de la vitesse (m/s)
!        IUPT  : U par classe et zones (m/s)
!        IVPT  : V par classe et zones (m/s)
!        IWPT  : W par classe et zones (m/s)
!        IDEBT : debit massique (kg/s)
!        IPOIT : poids statistique (nombre d'echantillons) associe
!                  a la particule (il est calcule automatiquement pour
!                  respecter un debit massique ce dernier est fournis)

!        En fonction de la physique
!          IDPT   : diametre (m)
!          IVDPT  : ecart-type du diametre (m)
!          ITPT   : temperature en degres Celsius (Pas d'enthalpie)
!          ICPT   : chaleur specifique (J/kg/K)
!          IEPSI  : emissivite (si =0 alors aucun effet radiatif n'est pris en compte)
!          IROPT  : masse volumique (kg/m3)

!          Si Charbon (IPHYLA=2)
!            IHPT  : temperature en degres Celsius (Pas d'enthalpie)
!            IMCHT : masse de charbon reactif (kg)
!            IMCKT : masse de coke (kg)


! ---> Premiere zone, numerotee IZONE = 1 (NBCLAS classes)
!        IUSCLB : adherence de la particule a une face de paroi
!        IJNBP  : 10 particules par classe,
!        IJFRE  : injection tous les 2 pas de temps
!        IJUVW, IUPT, IVPT, IWPT : vitesse imposee a 1.1D0, 0.0D0, 0.0D0
!        ICPT   : cp de 10000
!        ITPT   : temperature de 25 degres Celsius
!        IDPT   : diametre de 50.E-6 m
!        IEPSI  : emissivite a 0.7
!        IVDPT  : diametre constant ==> Ecart-type nul
!        IROPT  : masse volumique
!        IPOIT  : poids statistique (nombre de particules physiques
!                 represente par une particule-echantillon statistique)
!        IDEBT  : debit massique


izone     = 1
nbclas    = iusncl(izone)
iusclb (izone)         =  ientrl
do iclas  = 1, nbclas

  iuslag (iclas,izone,ijnbp) = 10
  iuslag (iclas,izone,ijfre) = 2

  if (nbclst.gt.0) then
    iuslag(iclas,izone,iclst) = 1
  endif

  iuslag (iclas,izone,ijuvw) = -1
  ruslag (iclas,izone,iupt)  = 1.1d0
  ruslag (iclas,izone,ivpt)  = 0.0d0
  ruslag (iclas,izone,iwpt)  = 0.0d0

  iuslag (iclas,izone,ijprpd)= 1
  ruslag (iclas,izone,ipoit) = 1.d0
  ruslag (iclas,izone,idebt) = 0.d0

!    Si physique simple

  if ( iphyla.eq.0 .or. iphyla.eq.1 ) then

!        Diametre et ecart-type du diametre

    iuslag (iclas,izone,ijprdp)= 1
    ruslag (iclas,izone,idpt)  = 50.d-6
    ruslag (iclas,izone,ivdpt) = 0.d0

!        Masse Volumique

    ruslag(iclas,izone,iropt) = 2500.d0

    if ( iphyla.eq.1 ) then

!        Temperature et Cp

      if ( itpvar.eq.1 ) then
        iuslag (iclas,izone,ijprtp) = 1
        ruslag(iclas,izone,itpt)    = 20.d0

        ruslag(iclas,izone,icpt)    = 1400.d0
        ruslag(iclas,izone,iepsi)   = 0.7d0
      endif

    endif

!    Charbon

  else if ( iphyla.eq.2 ) then

!    ATTENTION : 1) Pour transporter et bruler des particules de
!                   charbon en Lagrangien, une physique particuliere
!                   liee au charbon pulverise doit etre active pour
!                   la phase porteuse.

!                2) Les proprietes physiques des grains de charbon
!                   sont connus a partir du fichier de thermochimie :
!                                      dp_FCP.

!                3) Pour la classe Lagrangienne courante ICLAS, et pour
!                   la zone frontiere NB courante, on donne aux grains
!                   de charbon les proprietes du charbon ICHA de classe
!                   ICLAS lu dans le fichier dp_FCP.

!                4) ICHA : numero du Charbon compris entre 1 et NCHARB
!                   defini dans le fichier dp_FCP par l'utilisateur.
!                   (NCHARB = NCHARM = 3 au maximum)


    icha = ichcor(iclas)
    temp = 800.d0

!        Numero du charbon

    iuslag(iclas,izone,inuchl) = icha

!        Temperature et Cp

    ruslag(iclas,izone,ihpt) = temp
    ruslag(iclas,izone,icpt) = cp2ch(icha)

!        Diametre et son ecart-type (nul)

    ruslag (iclas,izone,idpt)  = diam20(iclas)
    ruslag (iclas,izone,ivdpt) = 0.d0

!        Masse Volumique

    ruslag(iclas,izone,iropt) =  rho0ch(icha)

!        Masse de charbon actif et
!        Masse de Coke (nulle si le charbon n'a jamais brule)

    mp0 = pis6 * ( ruslag(iclas,izone,idpt)**3 )                  &
               * ruslag(iclas,izone,iropt)
    ruslag(iclas,izone,imcht) = mp0 * (1.d0-xashch(icha))
    ruslag(iclas,izone,imckt) = 0.d0

  endif

enddo

! ---> Deuxieme zone, numerotee IZONE = 2 (NBCLAS classes)
!        IUSCLB : rebond de la particule

izone     = 2
iusclb (izone)         =  irebol


izone     = 4
iusclb (izone)         =  irebol


! de meme pour les autres zones ...

!===============================================================================

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end
