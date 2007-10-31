c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS                  Code_Saturne version 1.3
C                      ------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2007 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
                        SUBROUTINE USLAG2
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBPMAX , NVP    , NVP1   , NVEP   , NIVEP  ,
     &   NTERSL , NVLSTA , NVISBR ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ITYPFB , ITRIFB , ITEPA  , IFRLAG ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   ETTP   , TEPA   ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC     SOUS-PROGRAMME UTILISATEUR (INTERVENTION OBLIGATOIRE)
CFONC
CFONC     ROUTINE UTILISATEUR POUR LES CONDITIONS AUX LIMITES RELATIVES
CFONC       AUX PARTICULES (ENTREE ET TRAITEMENT AUX AUTRES BORDS)
CFONC
CFONC
CFONC  IDENTIFICATION DES FACES DE BORD
CFONC  ================================
CFONC L'identification des faces de bord concernees se fait grace
CFONC a la commande GETFBR.
CFONC
CFONC   GETFBR(CHAINE,NLELT,LSTELT) :
CFONC   - CHAINE est une chaine de caractere fournie par l'utilisateur
CFONC     qui donne les criteres de selection
CFONC   - NLTELT est renvoye par la commande. C'est un entier qui
CFONC     correspond au nombre de faces de bord trouveees repondant au
CFONC     critere
CFONC   - LSTELT est renvoye par la commande. C'est un tableau d'entiers
CFONC     de taille NLTELT donnant la liste des faces de bord trouvees
CFONC     repondant au critere.
CFONC
CFONC   CHAINE peut etre constitue de :
CFONC   - references de couleurs (ex. : 1, 8, 26, ...
CFONC   - references de groupes (ex. : entrees, groupe1, ...)
CFONC   - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
CFONC   Ces criteres peuvent etre combines par des operateurs logiques
CFONC   (AND et OR) et des parentheses
CFONC   ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
CFONC   les faces de bord de couleur 1, appartenant aux groupes 'groupe2'
CFONC   ou 'groupe3' et de coordonnee Y inferieure a 1.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC           !
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR           !
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NBPMAX       ! E  !  -> ! NOMBRE MAX DE PARTICULIES AUTORISE   !
CARGU ! NVP          ! E  !  -> ! NOMBRE DE VARIABLES PARTICULAIRES    !
CARGU ! NVP1         ! E  !  -> ! NVP SANS POSITION, VFLUIDE, VPART    !
CARGU ! NVEP         ! E  !  -> ! NOMBRE INFO PARTICULAIRES (REELS)    !
CARGU ! NIVEP        ! E  !  -> ! NOMBRE INFO PARTICULAIRES (ENTIERS)  !
CARGU ! NTERSL       ! E  !  -> ! NBR TERMES SOURCES DE COUPLAGE RETOUR!
CARGU ! NVLSTA       ! E  !  -> ! NOMBRE DE VAR STATISTIQUES LAGRANGIEN!
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU !              !    !     ! LE MODULE LAGRANGIEN                 !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU !  (NFML,NPRFML!    !     !                                      !
CARGU ! MAXELT       !  E !  -> ! NB MAX D'ELEMENTS (CELL,FAC,FBR)     !
CARGU ! LSTELT(MAXELT) TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC            !
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !                                      !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR            !
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !                                      !
CARGU ! ITRIFB(NFABOR! TE !  -> ! INDIRECTION POUR TRI DES FACES DE BRD!
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! IFRLAG(NFABOR! TE ! <-  ! TYPE DES FACES DE BORD LAGRANGIEN    !
CARGU ! ITEPA        ! TE !  -> ! INFO PARTICULAIRES (ENTIERS)         !
CARGU ! (NBPMAX,NIVEP!    !     !   (CELLULE DE LA PARTICULE,...)      !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! VECTEUR SURFACE DES FACES INTERNES   !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! SURFBO       ! TR !  -> ! VECTEUR SURFACE DES FACES DE BORD    !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS                !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTPA         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT PREC)           !
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! ETTP         ! TR !  -> ! TABLEAUX DES VARIABLES LIEES         !
CARGU !  (NBPMAX,NVP)!    !     !   AUX PARTICULES ETAPE COURANTE      !
CARGU ! TEPA         ! TR !  -> ! INFO PARTICULAIRES (REELS)           !
CARGU ! (NBPMAX,NVEP)!    !     !   (POIDS STATISTIQUES,...)           !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "cpincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NBPMAX , NVP    , NVP1   , NVEP  , NIVEP
      INTEGER          NTERSL , NVLSTA , NVISBR
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1) , NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1) , NODFBR(LNDFBR)
      INTEGER          ITYPFB(NFABOR,NPHAS) , ITRIFB(NFABOR,NPHAS)
      INTEGER          ITEPA(NBPMAX,NIVEP) , IFRLAG(NFABOR)
      INTEGER          IDEVEL(NIDEVE) , ITUSER(NITUSE)
      INTEGER          IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC) , SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC) , CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD) , VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET) , RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*) , PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*) , COEFB(NFABOR,*)
      DOUBLE PRECISION ETTP(NBPMAX,NVP) , TEPA(NBPMAX,NVEP)
      DOUBLE PRECISION RDEVEL(NRDEVE) , RTUSER(NRTUSE)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC , IFML , ICOUL , IZONE, NBCLAS, ICLAS
      INTEGER          IIII , ICHA
      INTEGER          ILELT, NLELT
C
      DOUBLE PRECISION PIS6 , MP0 , TEMP
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C=======================================================================
C
      IF(1.EQ.1) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,
     &'@    *********                                               ',/,
     &'@     MODULE LAGRANGIEN :                                    ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR uslag2 DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C=======================================================================
C 1.  GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 2. INITIALISATION
C=======================================================================
C
      PIS6 = PI / 6.D0
C
C=======================================================================
C 3. CONSTRUCTION DES ZONES FRONTIERE
C=======================================================================
C
C
C     DEFINITION DES ZONES FRONTIERES
C
C     Pour le module Lagrangien, l'utilisateur definit NFRLAG zones
C       frontieres a partir de la couleur des faces de bord, plus
C       generalement de leurs proprietes (couleurs, groupes...),
C       ou des conditions aux limites fixees dans usclim
C       ou meme de leur coordonnees. Pour ce faire, on renseigne
C       le tableau IFRLAG(NFABOR) qui donne pour
C       chaque face de bord IFAC le numero de la zone a laquelle elle
C       appartient IFRLAG(IFAC).
C     Attention : il est indispensable que TOUTES les faces aient ete
C       affectees a une zone.
C     Le numero des zones (donc les valeurs de IFRLAG(IFAC)) est
C       arbitrairement choisi par l'utilisateur, mais doit etre un
C       entier strictement positif et inferieur ou egal a NFLAGM
C       (valeur fixee en parametre dans lagpar.h).
C     On affecte ensuite a chaque zone un type ITYLAG qui sera utilise
C       pour imposer des conditions aux limites globales.
C
C
      IZONE = -1
C
C ---> Premiere zone, numerotee IZONE = 1 ( = couleur 10)
      CALL GETFBR('10',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IZONE        = 1
        IFRLAG(IFAC) = IZONE
C
      ENDDO
C
C ---> Deuxieme zone, numerotee IZONE = 2 ( = partie des couleur 4)
      CALL GETFBR('4 and Y < 1.0',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IZONE        = 2
        IFRLAG(IFAC) = IZONE
C
      ENDDO
C
C ---> Troisieme zone, numerotee IZONE = 4 ( = entree phase 1)
      DO IFAC = 1, NFABOR
        IF(ITYPFB(IFAC,1).EQ.IENTRE) THEN
          IZONE        = 4
          IFRLAG(IFAC) = IZONE
        ENDIF
      ENDDO
C
C ---> Nieme zone, numerotee IZONE = 5 ( = couleur 3)
      CALL GETFBR('3',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IZONE        = 5
        IFRLAG(IFAC) = IZONE
C
      ENDDO
C
C
C=======================================================================
C 4. INJECTION PAR CLASSES DE PARTICULES DANS LE DOMAINE DE CALCUL
C=======================================================================
C
C   POUR DONNER LES INFORMATIONS SUR LES CLASSES DE PARTICULES,
C   ON PROCEDE EN DEUX ETAPES :
C
C   1) TOUT D'ABORD, ON RENSEIGNE LE NOMBRE DE CLASSES DE PARTICULES
C      PAR ZONE FRONTIERE : IUSNCL (ce nombre est nul par defaut)
C
C   2) ENSUITE PAR ZONE DE COULEUR ET PAR CLASSE,
C      ON DONNE LES CARACTERISTIQUES D'ENTREE ET
C      LES CARACTERISTIQUE PHYSIQUES DES PARTICULES
C
C
C
C
C --> Nombre de classes de particules entrantes
C   On renseigne ici le nombre de classes pour toutes les zones
C     identifiees precedemment.
C   Par defaut le nombre de classes de particules est nul.
C   Le nombre max de classes est NCLAGM donne dans lagpar
C
C ---> Premiere zone, numerotee IZONE = 1 : 1 classe injectee
      IZONE     = 1
      NBCLAS    = 1
      IUSNCL(IZONE) = NBCLAS
C
C ---> Deuxieme zone, numerotee IZONE = 2 : 0 classe injectee
      IZONE     = 2
      NBCLAS    = 0
      IUSNCL(IZONE) = NBCLAS
C
C ---> Troisieme zone, numerotee IZONE = 4 : 0 classe injectee
      IZONE     = 4
      NBCLAS    = 0
      IUSNCL(IZONE) = NBCLAS
C
C ---> Nieme zone,     numerotee IZONE = 5 : 0 classe injectee
      IZONE     = 5
      NBCLAS    = 0
      IUSNCL(IZONE) = NBCLAS
C
C
C --> Pour chaque classe de particules associee a une zone,
C     il faut fournir des informations.
C
C
C     IUSNCL nbr de classes par zones
C     IUSCLB conditions au bord pour les particules
C     = IENTRL -> zone d'injection de particules
C     = ISORTL -> sortie du domaine
C     = IREBOL -> rebond des particules
C     = IDEPO1 -> deposition definitive
C     = IDEPO2 -> deposition definitive mais la particule reste en
C                 memoire (utile si IENSI2 = 1 uniquement)
C     = IDEPO3 -> deposition et remise en suspension possible
C                 suivant les condition de l'ecoulement
C     = IDEPFA -> deposition de la particule avec force d'attachement,
C                 vitesse conservee et re-entrainement possible
C                 (Possible si LADLVO = 1 )
C     = IENCRL -> encrassement (Charbon uniquement IPHYLA = 2)
C     = JBORD1 -> interaction part/frontiere utilisateur (cf. USLABO)
C     = JBORD2 -> interaction part/frontiere utilisateur (cf. USLABO)
C     = JBORD3 -> interaction part/frontiere utilisateur (cf. USLABO)
C     = JBORD4 -> interaction part/frontiere utilisateur (cf. USLABO)
C     = JBORD5 -> interaction part/frontiere utilisateur (cf. USLABO)
C
C
C
C     Tableau IUSLAG :
C     ================
C        IJNBP : nbr de part par classe et zones
C        IJFRE : frequence d'injection, si la frequence est nulle alors
C                  il y a injection uniquement a la 1ere iteration absolue
C        ICLST : numero de groupe auquel appartient la particule
C                 (uniquement si on souhaite des statistiques par groupe)
C        IJUVW : type de condition vitesse
C                  = -1 vitesse fluide imposee
C                  =  0 vitesse imposee selon la direction normale a la face
C                        de bord et de norme RUSLAG(ICLAS,IZONE,IUNO)
C                  =  1 vitesse imposee : on donne RUSLAG(ICLAS,IZONE,IUPT)
C                                                  RUSLAG(ICLAS,IZONE,IVPT)
C                                                  RUSLAG(ICLAS,IZONE,IWPT)
C                  =  2 profil utilisateur
C        IJPRTP : type de condition temperature
C                  =  1 temperature imposee : on donne RUSLAG(ICLAS,IZONE,ITPT)
C                  =  2 profil utilisateur
C        IJPRDP : type de condition diametre
C                  =  1 vitesse imposee : on donne RUSLAG(ICLAS,IZONE,IDPT)
C                                                  RUSLAG(ICLAS,IZONE,IVDPT)
C                  =  2 profil utilisateur
C        INUCHL : numero du charbon de la particule (uniquement si IPHYLA=2)
C
C     Tableau RUSLAG :
C     ================
C        IUNO  : norme de la vitesse (m/s)
C        IUPT  : U par classe et zones (m/s)
C        IVPT  : V par classe et zones (m/s)
C        IWPT  : W par classe et zones (m/s)
C        IDEBT : debit massique (kg/s)
C        IPOIT : poids statistique (nombre d'echantillons) associe
C                  a la particule (il est calcule automatiquement pour
C                  respecter un debit massique ce dernier est fournis)
C
C        En fonction de la physique
C          IDPT   : diametre (m)
C          IVDPT  : ecart-type du diametre (m)
C          ITPT   : temperature en degres Celsius (Pas d'enthalpie)
C          ICPT   : chaleur specifique (J/kg/K)
C          IEPSI  : emissivite (si =0 alors aucun effet radiatif n'est pris en compte)
C          IROPT  : masse volumique (kg/m3)
C
C          Si Charbon (IPHYLA=2)
C            IHPT  : temperature en degres Celsius (Pas d'enthalpie)
C            IMCHT : masse de charbon reactif (kg)
C            IMCKT : masse de coke (kg)
C
C
C ---> Premiere zone, numerotee IZONE = 1 (NBCLAS classes)
C        IUSCLB : adherence de la particule a une face de paroi
C        IJNBP  : 10 particules par classe,
C        IJFRE  : injection tous les 2 pas de temps
C        IJUVW, IUPT, IVPT, IWPT : vitesse imposee a 1.1D0, 0.0D0, 0.0D0
C        ICPT   : cp de 10000
C        ITPT   : temperature de 25 degres Celsius
C        IDPT   : diametre de 50.E-6 m
C        IEPSI  : emissivite a 0.7
C        IVDPT  : diametre constant ==> Ecart-type nul
C        IROPT  : masse volumique
C        IPOIT  : poids statistique (nombre de particules physiques
C                 represente par une particule-echantillon statistique)
C        IDEBT  : debit massique
C
C
      IZONE     = 1
      NBCLAS    = IUSNCL(IZONE)
      IUSCLB (IZONE)         =  IENTRL
      DO ICLAS  = 1, NBCLAS
C
        IUSLAG (ICLAS,IZONE,IJNBP) = 10
        IUSLAG (ICLAS,IZONE,IJFRE) = 2
C
        IF (NBCLST.GT.0) THEN
          IUSLAG(ICLAS,IZONE,ICLST) = 1
        ENDIF
C
        IUSLAG (ICLAS,IZONE,IJUVW) = -1
        RUSLAG (ICLAS,IZONE,IUPT)  = 1.1D0
        RUSLAG (ICLAS,IZONE,IVPT)  = 0.0D0
        RUSLAG (ICLAS,IZONE,IWPT)  = 0.0D0
C
        IUSLAG (ICLAS,IZONE,IJPRPD)= 1
        RUSLAG (ICLAS,IZONE,IPOIT) = 1.D0
        RUSLAG (ICLAS,IZONE,IDEBT) = 0.D0
C
C    Si physique simple
C
        IF ( IPHYLA.EQ.0 .OR. IPHYLA.EQ.1 ) THEN
C
C        Diametre et ecart-type du diametre
C
          IUSLAG (ICLAS,IZONE,IJPRDP)= 1
          RUSLAG (ICLAS,IZONE,IDPT)  = 50.D-6
          RUSLAG (ICLAS,IZONE,IVDPT) = 0.D0
C
C        Masse Volumique
C
          RUSLAG(ICLAS,IZONE,IROPT) = 2500.D0
C
          IF ( IPHYLA.EQ.1 ) THEN
C
C        Temperature et Cp
C
            IF ( ITPVAR.EQ.1 ) THEN
              IUSLAG (ICLAS,IZONE,IJPRTP) = 1
              RUSLAG(ICLAS,IZONE,ITPT)    = 20.D0
C
              RUSLAG(ICLAS,IZONE,ICPT)    = 1400.D0
              RUSLAG(ICLAS,IZONE,IEPSI)   = 0.7D0
            ENDIF
C
          ENDIF
C
C    Charbon
C
        ELSE IF ( IPHYLA.EQ.2 ) THEN
C
C    ATTENTION : 1) Pour transporter et bruler des particules de
C                   charbon en Lagrangien, une physique particuliere
C                   liee au charbon pulverise doit etre active pour
C                   la phase porteuse.
C
C                2) Les proprietes physiques des grains de charbon
C                   sont connus a partir du fichier de thermochimie :
C                                      dp_FCP.
C
C                3) Pour la classe Lagrangienne courante ICLAS, et pour
C                   la zone frontiere NB courante, on donne aux grains
C                   de charbon les proprietes du charbon ICHA de classe
C                   ICLAS lu dans le fichier dp_FCP.
C
C                4) ICHA : numero du Charbon compris entre 1 et NCHARB
C                   defini dans le fichier dp_FCP par l'utilisateur.
C                   (NCHARB = NCHARM = 3 au maximum)
C
C
          ICHA = ICHCOR(ICLAS)
          TEMP = 800.D0
C
C        Numero du charbon
C
          IUSLAG(ICLAS,IZONE,INUCHL) = ICHA
C
C        Temperature et Cp
C
          RUSLAG(ICLAS,IZONE,IHPT) = TEMP
          RUSLAG(ICLAS,IZONE,ICPT) = CP2CH(ICHA)
C
C        Diametre et son ecart-type (nul)
C
          RUSLAG (ICLAS,IZONE,IDPT)  = DIAM20(ICLAS)
          RUSLAG (ICLAS,IZONE,IVDPT) = 0.D0
C
C        Masse Volumique
C
          RUSLAG(ICLAS,IZONE,IROPT) =  RHO0CH(ICHA)
C
C        Masse de charbon actif et
C        Masse de Coke (nulle si le charbon n'a jamais brule)
C
          MP0 = PIS6 * ( RUSLAG(ICLAS,IZONE,IDPT)**3 )
     &               * RUSLAG(ICLAS,IZONE,IROPT)
          RUSLAG(ICLAS,IZONE,IMCHT) = MP0 * (1.D0-XASHCH(ICHA))
          RUSLAG(ICLAS,IZONE,IMCKT) = 0.D0
C
        ENDIF
C
      ENDDO
C
C ---> Deuxieme zone, numerotee IZONE = 2 (NBCLAS classes)
C        IUSCLB : rebond de la particule
C
      IZONE     = 2
      IUSCLB (IZONE)         =  IREBOL
C
C
      IZONE     = 4
      IUSCLB (IZONE)         =  IREBOL
C

C de meme pour les autres zones ...
C
C=======================================================================
C
C--------
C FORMATS
C--------
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
