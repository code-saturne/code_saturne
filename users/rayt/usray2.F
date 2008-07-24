c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2008 EDF S.A., France
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
                        SUBROUTINE USRAY2
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , IZFRDP , ISOTHP ,
     &   TMIN   , TMAX   , TX     ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
     &   TPAROP , QINCID , HFCNVP , FLCNVP ,
     &   XLAMP  , EPAP   , EPSP   , TEXTP  , TINTP  ,
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE DE RAYONNEMENT :
CFONC    -----------------------------------------
CFONC
CFONC    CARACTERISTIQUES ET COMPORTEMENT RADIATIF DES FACES DE PAROI
CFONC
CFONC    Ce sous-programme est appele autant de fois qu'il y a de phases
CFONC    pour lesquelles il faut faire un calcul de rayonnement
CFONC    semi-transparent
CFONC
CFONC    IPHAS contient le numero absolu de la phase concernee
CFONC
CFONC
CFONC
CFONC   EXTERIEUR       PAROI       DOMAINE FLUIDE
CFONC             |               |
CFONC             |               |
CFONC             |               |
CFONC             | EPSP Qincid<--|<----- Qincid             [flux incident]
CFONC             |               |
CFONC             | [absorption]  |
CFONC             |               |
CFONC             |               |-----> (1-EPSP) Qincid    [reflexion]
CFONC             |               |
CFONC             |               |
CFONC             |               |                       4
CFONC             |               |-----> EPSP sigma Tparop  [emission]
CFONC             |               |
CFONC             |               |
CFONC             |               |
CFONC             | [conduction]  |     .....o Tfluide       [convection]
CFONC             |               |    .
CFONC             |     XLAMP     |   .    Hfluide
CFONC             |               |  .
CFONC             |            ...o..
CFONC             |          ..   | Tparop (temperature en paroi)
CFONC             |         .     |
CFONC             |       ..      |
CFONC             |      .        |
CFONC             |    ..         |
CFONC             |   .           |
CFONC             o...            |
CFONC        Textp|               |
CFONC             |               |
CFONC             |<------------->|
CFONC             |     EPAP      |
CFONC             |               |
CFONC
CFONC
CFONC   LE CALCUL DES TEMPERATURES EST REALISEE PAR UN BILAN DE FLUX
CFONC   COMME SUIT :
CFONC
CFONC   Q           =  Q           + (Qrayt           - Qrayt        )
CFONC    conduction     convection         absorption        emission
CFONC
CFONC   SOIT (flux positif si sortant du domaine de calcul) :
CFONC
CFONC    XLAMP
CFONC    -----(Tparop-Textp) =
CFONC    EPAP
CFONC                                                              4
CFONC          Hfluide (Tfluide-Tparop) + EPSP (QINCID - SIGMA Tparop )
CFONC
CFONC
CFONC       CORPS                        EPSP
CFONC       ------------------------------------
CFONC       acier poli                   0,06
CFONC       acier oxyde                  0,80
CFONC       aceir rugueux                0,94
CFONC       aluminium poli               0,04
CFONC       aluminium oxyde (interieur)  0,09
CFONC       aluminium oxyde (air humide) 0,90
CFONC       brique                       0,93
CFONC       beton                        0,93
CFONC       papier                     0,8 a 0,9
CFONC       eau                          0,96
CFONC
CFONC
CFONC  IDENTIFICATION DES FACES DE PAROI
CFONC  =================================
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
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE COURANTE          !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! MAXELT       !  E !  -> ! NB MAX D'ELEMENTS (CELL,FAC,FBR)     !
CARGU ! LSTELT(MAXELT) TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ICODCL       ! TE !  -> ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
CARGU ! IZFRDP(NFABOR! TE ! <-  ! NUMERO DE ZONE POUR LES FACES DE BORD!
CARGU ! ISOTHP(NFABOR! TE ! <-  ! LISTE DES FRONTIERES ISOTHERMES      !
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
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS (OPTIONNEL)    !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! RCODCL       ! TR ! <-  ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2 OU      !
CARGU !              !    !     !  HAUTEUR DE RUGOSITE (m) si ICODCL=6 !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! W1...6(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! TPAROP(NFABOR! TR !  -> ! TEMPERATURE DE PAROI EN KELVIN       !
CARGU ! QINCID(NFABOR! TR !  -> ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU ! HFCNVP(NFABOR! TR !  -> ! COEFFICIENT D'ECHANGE FLUIDE AUX     !
CARGU !              !    !     ! FACES DE BORD                        !
CARGU ! FLCNVP(NFABOR! TR !  -> ! DENSITE DE FLUX CONVECTIF AUX FACES  !
CARGU !              !    !     ! DE BORD                              !
CARGU ! XLAMP(NFABOR)! TR ! <-  ! COEFFICIENT DE CONDUCTIVITE THERMIQUE!
CARGU !              !    !     ! DES FACETTES DE PAROI (W/m/K)        !
CARGU ! EPAP(NFABOR) ! TR ! <-  ! EPAISSEUR DES FACETTES DE PAROI (m)  !
CARGU ! EPSP (NFABOR)! TR ! <-  ! EMISSIVITE DES FACETTES DE BORD      !
CARGU ! TEXTP(NFABOR)! TR ! <-  ! TEMPERATURE DE BORD EXTERNE          !
CARGU !              !    !     ! EN KELVIN                            !
CARGU ! TINTP(NFABOR)! TR ! <-  ! TEMPERATURE DE BORD INTERNE          !
CARGU !              !    !     ! EN KELVIN                            !
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
C-----------------------------------------------------------------------
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
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "radiat.h"
      INCLUDE "ihmpre.h"
C
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , IPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML) , ITYPFB(NFABOR,*)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          IZFRDP(NFABOR), ISOTHP(NFABOR)
C
      DOUBLE PRECISION TMIN , TMAX , TX
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
C
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
C
      DOUBLE PRECISION TPAROP(NFABOR), QINCID(NFABOR)
      DOUBLE PRECISION HFCNVP(NFABOR),FLCNVP(NFABOR)
      DOUBLE PRECISION XLAMP(NFABOR), EPAP(NFABOR)
      DOUBLE PRECISION EPSP(NFABOR)
      DOUBLE PRECISION TEXTP(NFABOR), TINTP(NFABOR)
C
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IFAC , IVAR, IOK
      INTEGER          ILELT, NLELT
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
      IF (IIHMPR.EQ.1) THEN
        RETURN
      ELSE
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usray2 DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C=======================================================================
C 0. GESTION MEMOIRE
C       Aucune modification requise
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. NUMERO DE LA VARIABLE THERMIQUE
C
C       Noter que TKELVI, disponible ici, contient 273.16D0
C       Noter que SIG   , disponible ici, contient 5.6703D-8
C=======================================================================
C
C     IVAR est le numero de la variable thermique.
C
      IVAR = ISCA(ISCALT(IPHAS))
C
C
C=======================================================================
C  2. VALEURS MIN ET MAX ADMISSIBLES POUR LA TEMPERATURE DE PAROI
C                            (clipping)
C
C
C
C                            EN KELVIN
C
C
C        (quel que soit le choix fait pour la variable thermique)
C
C
C
C
C   TMIN et TMAX sont respectivement les valeurs minimale et maximale
C     autorisees pour les temperatures de paroi. A chaque nouveau
C     calcul de temperature de paroi, on verifie qu'il n'y a pas
C     depassement, sinon il y a clipping.
C
C   L'exemple ci-dessous correspond a un clipping tres peu contraignant.
C     (valeurs par defaut)
C
C=======================================================================
C
      TMIN = 0.D0
      TMAX = GRAND + TKELVI
C
C=======================================================================
C 3. ZONES FRONTIERE ET PROPRIETES DES PAROIS
C=======================================================================
C
C
C     DEFINITION DES ZONES FRONTIERES
C     ===============================
C
C     On definit des zones de faces de paroi, et on leur affecte un type.
C       Ceci permet d'appliquer les conditions aux limites et de realiser
C       des bilans en traitant separement les differents sous ensembles
C       ainsi constitues (on peut ainsi connaitre par exemple le flux au
C       travers des differentes zones definies par l'utilisateur).
C
C     Pour chaque face de bord IFAC (pas uniquement les faces de paroi),
C       l'utilisateur definit selon son propre choix un numero de zone
C       frontiere IZFRDP(IFAC) a partir
C         de la couleur des faces de bord,
C         ou, plus generalement, de leurs proprietes (couleurs, groupes...),
C         ou des conditions aux limites fixees dans usclim,
C         ou meme de leur coordonnees.
C     Attention : il est indispensable que TOUTES les faces de bord
C       aient ete affectees a une zone.
C     Le numero des zones (la valeur de IZFRDP(IFAC)) est
C       arbitrairement choisi par l'utilisateur, mais doit etre un
C       entier strictement positif et inferieur ou egal a NBZRDM
C       (valeur fixee en parametre dans radiat.h, a ne pas modifier).
C
C
C
C     PROPRIETES DES PAROIS
C     =====================
C
C
C      ATTENTION :
C      ---------
C
C             L'unite des temperatures est le KELVIN
C
C        quel que soit le choix de resolution fait par ailleurs
C
C
C
C      DONNEES OBLIGATOIRES :
C      --------------------
C
C      ISOTHP = ITPIMP
C                repere les facettes de paroi a temperature imposee
C                  ("Indicateur Temperature Paroi IMPosee" ou "isotherme")
C                  la temperature a imposer est TINTP(Kelvin)
C             = IPGRNO ou IPREFL
C                repere une facette de paroi dont il faut calculer
C                  la temperature (un bilan de flux sera fait).
C                  IPGRNO : Indicateur de Paroi GRise ou NOire (EPSP > 0)
C                  IPREFL : Indicateur de Paroi REFLechissante (EPSP = 0)
C             = IFGRNO ou IFREFL
C                repere une facette de paroi pour laquelle on impose un
C                  flux (conductif en paroi ou total sur le fluide)
C                  IFGRNO : Indicateur de Flux en paroi GRise ou NOire
C                  IFREFL : Indicateur de Flux en paroi REFLechissante
C
C      TINTP = Temperature de peau interne (KELVIN)
C               sert a l'initialisation de TPAROP en debut de calcul
C               (lorsqu'il n'y a pas de fichier suite rayonnement)
C                 Si ISOTHP = ITPIMP, la valeur de TPAROP est egale a
C                   celle imposee dans TINTP a chaque pas de temps
C                 Dans les autres cas, TINTP ne sert qu'a l'initialisation.
C
C
C      RENSEIGNEMENTS COMLEMENTAIRES A FOURNIR SUIVANT ISOTHP (VOIR EXEMPLES)
C      --------------------------------------------------------------------
C
C      RCODCL = valeur de flux (voir USCLIM pour la definition)
C
C      EPSP   = Emissivite aux facettes de paroi (dans l'intervalle [0;1])
C
C      XLAMP  = Conductivite thermique des parois (W/m/K)
C               > 0 quand il est renseigne
C
C      EPAP   = Epaisseur des parois (m)
C               > 0 quand il est renseigne
C
C      TEXTP  = Temperature de peau externe en KELVIN
C               > 0 quand il est renseigne
C
C
C
C     EXEMPLE
C     =======
C
C      Attention a la coherence avec USCLIM
C
C      Dans cet EXEMPLE,
C               -------
C
C        Les faces de paroi (IPAROI et IPARUG), sont reparties en
C          5 sous ensembles (zones) reperes par IFRFAC(IFAC) variant de
C          51 a 55 (arbitraire) a chacun desquels on applique une
C          condition (ISOTHP) differente.
C        Les autres faces de bord (entree, sortie, sous-ensembles
C          arbitraires sur lesquels on souhaite voir imprimer des
C          informations) sont repartis en zones reperees
C          par IFRFAC(IFAC), dont la valeur peut etre arbitrairement
C          choisie (>0 et <NBZRDM+1).
C
C        Remarque : une maniere simple (mais non generale) d'affecter
C          rapidement des valeurs a IFRFAC est IFRFAC(IFAC) = ICOUL,
C          ou ICOUL est la couleur de la face.
C
C
C
C
C
C
C
C
C
C
C             DANS TOUS LES CAS, IL EST INTERDIT
C                                       ========
C
C                       DE MODIFIER TPAROP ou QINCID ICI.
C
C
C
C
C
C
C
C
C    Indicateur : nombre de faces oubliees.
      IOK = 0
C
C   -------------------------------------------------------------------
C-->  Exemple 1 :
C      Pour les faces PAROI de couleur 1 :
C           Profil de temperature imposee
C       ------------------------------------
C
      CALL GETFBR('1',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IF ( ITYPFB(IFAC,IPHAS).EQ.IPAROI ) THEN
C
C      Numero de zone
          IZFRDP(IFAC) = 51
C
C      Type de condition : Paroi grise ou noire et
C                          Profil de temperature imposee TINTP
          ISOTHP(IFAC) = ITPIMP
C
C      Donnees complementaires necessaires et suffisantes
C        Emissivite
          EPSP  (IFAC) = 0.1D0
C        Temperature de paroi imposee ("Profil" simple ici = 473.16 K)
          TINTP (IFAC) = 200.D0 + TKELVI
C
        ENDIF
C
      ENDDO
C
C   -------------------------------------------------------------------
C-->  Exemple 2 :
C      Pour les faces PAROI RUGUEUSE de couleur 2 :
C           Paroi grise ou noire et
C           Profil de temperature imposee TEXTP
C       ------------------------------------
C
      CALL GETFBR('2',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IF ( ITYPFB(IFAC,IPHAS).EQ.IPARUG ) THEN
C
C      Numero de zone
          IZFRDP(IFAC) = 52
C
C      Type de condition : Paroi grise ou noire et
C                          Profil de temperature imposee TEXTP
          ISOTHP(IFAC) = IPGRNO
C
C      Donnees complementaires necessaires et suffisantes
C        Emissivite (>0)
          EPSP  (IFAC) = 0.9D0
C        Conductivite (W/m/K)
          XLAMP (IFAC) = 3.0D0
C        Epaisseur    (m)
          EPAP  (IFAC) = 0.1D0
C        Temperature externe imposee : 473.16 K
          TEXTP (IFAC) = 200.D0 + TKELVI
C        Temperature de paroi initiale : 473.16 K
          TINTP (IFAC) = 200.D0 + TKELVI
C
        ENDIF
C
      ENDDO
C
C   -------------------------------------------------------------------
C-->  Exemple 3 :
C      Pour les faces PAROI de couleur 3 :
C           Paroi reflechissante (EPSP = 0) et
C           Profil de temperature imposee TEXTP
C       ------------------------------------
C
      CALL GETFBR('3',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IF ( ITYPFB(IFAC,IPHAS).EQ.IPAROI ) THEN
C
C      Numero de zone
          IZFRDP(IFAC) = 53
C
C      Type de condition : Paroi reflechissante (EPSP = 0) et
C                          Profil de temperature imposee TEXTP
          ISOTHP(IFAC) = IPREFL
C
C      Donnees complementaires necessaires et suffisantes
C        Conductivite (W/m/K)
          XLAMP (IFAC) = 3.0D0
C        Epaisseur    (m)
          EPAP  (IFAC) = 0.1D0
C        Temperature externe imposee : 473.16 K
          TEXTP (IFAC) = 200.D0 + TKELVI
C        Temperature de paroi initiale : 473.16 K
          TINTP (IFAC) = 200.D0 + TKELVI
C
        ENDIF
C
      ENDDO
C
C   -------------------------------------------------------------------
C-->  Exemple 4 :
C      Pour les faces PAROI de couleur 4 :
C           Paroi grise ou noire et
C           Flux de conduction impose dans la paroi
C
C        XLAMP
C        -----(Tparop-Textp) = Flux de conduction impose (W/m2)
C        EPAP
C                         = RODCL(IFAC,IVAR,3)
C
C       Si le flux de conduction est nul la paroi est adiabatique
C       Le tableau RCODCL(IFAC,IVAR,3) recoit la valeur du flux.
C       Densite de flux (< 0 si gain pour le fluide)
C       Reprise de USCLIM :
C         Pour les temperatures T,         en Watt/m2       :
C            RCODCL(IFAC,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
C         Pour les enthalpies H,           en Watt/m2       :
C            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
C       ------------------------------------
C
      CALL GETFBR('4',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IF ( ITYPFB(IFAC,IPHAS).EQ.IPAROI ) THEN
C
C      Numero de zone
          IZFRDP(IFAC) = 54
C
C      Type de condition : Paroi grise ou noire et
C                          Flux de conduction impose dans la paroi
          ISOTHP(IFAC) = IFGRNO
C
C      Donnees complementaires necessaires et suffisantes
C        Emissivite (>0)
          EPSP  (IFAC) = 0.9D0
C        Flux de conduction (W/m2)
          RCODCL(IFAC,IVAR,3) = 0.D0
C        Temperature de paroi initiale : 473.16 K
          TINTP (IFAC) = 200.D0 + TKELVI
C
        ENDIF
C
      ENDDO
C
C   -------------------------------------------------------------------
C-->  Exemple 5 :
C      Pour les faces PAROI de couleur 5 :
C           Paroi reflechissante (EPSP = 0) et
C           Flux de conduction impose dans la paroi
C
C      Equivalent a imposer une condition de flux au fluide
C
C        XLAMP
C        -----(Tparop-Textp) = Flux de conduction impose ET EPSP = 0
C        EPAP
C                         = RODCL(IFAC,IVAR,3)
C
C       Si le flux de conduction est nul la paroi est adiabatique
C       Le tableau RCODCL(IFAC,IVAR,3) contient la valeur du flux (W/m2)
C       imposee comme condition au limite pour le fluide
C       Densite de flux (< 0 si gain pour le fluide)
C       Reprise de USCLIM :
C         Pour les temperatures T,         en Watt/m2       :
C            RCODCL(IFAC,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
C         Pour les enthalpies H,           en Watt/m2       :
C            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
C       ------------------------------------
C
      CALL GETFBR('5',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IF ( ITYPFB(IFAC,IPHAS).EQ.IPAROI ) THEN
C
C      Numero de zone
          IZFRDP(IFAC) = 55
C
C      Type de condition : Paroi reflechissante (EPSP = 0) et
C                          Flux de conduction impose dans la paroi
          ISOTHP(IFAC) = IFREFL
C
C      Donnees complementaires necessaires et suffisantes
C        Flux de conduction (W/m2)
          RCODCL(IFAC,IVAR,3) = 0.D0
C        Temperature de paroi initiale : 473.16 K
          TINTP (IFAC) = 200.D0 + TKELVI
C
        ENDIF
C
      ENDDO
C
C
C
C
C
C
C
C                           ATTENTION
C
C
C
C
C
C
C
C
C
C
C
C   -------------------------------------------------------------------
C-->   Pour les faces non paroi : reperage obligatoire de la zone
C        Les informations (flux par exemple) seront donnees zone par zone
C        et il est donc conseille de regrouper entre elles les faces
C        de meme nature.
C
C                             PAR SECURITE
C                        DANS CE SOUS-PROGRAMME
C          ON DEMANDE A L'UTILISATEUR DE REPERER TOUTES LES FACES
C                             PAROI OU NON
C
C        (il est indispensable d'affecter une valeur a tous les
C            IZFRDP(IFAC), IFAC = 1, NFABOR)
C
C
C
C
C       ------------------------------------
C
      DO IFAC = 1, NFABOR
C
        IF     ( ITYPFB(IFAC,IPHAS).EQ.ISOLIB                  ) THEN
          IZFRDP(IFAC) = 61
        ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.IENTRE.AND.
     &           CDGFBO(2,IFAC)    .GT.0.D0                    ) THEN
          IZFRDP(IFAC) = 62
        ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.IENTRE.AND.
     &           CDGFBO(2,IFAC)    .LE.0.D0                    ) THEN
          IZFRDP(IFAC) = 63
        ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.ISYMET                  ) THEN
          IZFRDP(IFAC) = 64
C
C
C
C   -------------------------------------------------------------------
C-->  Exemple 7 :
C      Verification que toutes les faces de paroi ont ete vues
C       Cette precaution est souhaitable dans tous les cas.
C       ------------------------------------
C
        ELSEIF ( ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &           ITYPFB(IFAC,IPHAS).EQ.IPARUG     ) THEN
          IF (IZFRDP(IFAC) .EQ. -1) THEN
            WRITE(NFECRA,1000)IFAC
            IOK = IOK + 1
          ENDIF
        ENDIF
C
C
C
C
C     Fin de la boucle sur les faces de bord
C     --------------------------------------
C
      ENDDO
C
C Stop si des faces ont ete oubliees
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C -------
C FORMAT
C -------
C
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     LES DONNEES RAYONNEMENT NE SONT PAS RENSEIGNEES POUR   ',/,
     &'@       LA FACE DE PAROI ',I10                                ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C ---
C FIN
C ---
C
C
      RETURN
C
      END
