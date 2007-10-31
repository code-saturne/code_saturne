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
                        SUBROUTINE USCFCL
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ITRIFB , ITYPFB , IZFPPP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     ROUTINE UTILISATEUR
CFONC     REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
CFONC     (ICODCL,RCODCL) POUR LES VARIABLES INCONNUES
CFONC
CFONC
CFONC
CFONC     CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE
CFONC     =============================================
CFONC
CFONC
CFONC
CFONC INTRODUCTION
CFONC ============
CFONC
CFONC On donne ici les conditions aux limites par face de bord.
CFONC
CFONC L'identification des faces de bord concernees se fait grace
CONC a la commande GETFBR.
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
CFONC
CFONC TYPE DE CONDITIONS AUX LIMITES
CFONC ==============================
CFONC
CFONC En compressible, on ne peut affecter que les conditions aux
CFONC   limites predefinies
CFONC
CFONC     IPAROI, ISYMET, IESICF, ISSPCF, ISOPCF, IERUCF, IEQHCF
CFONC
CFONC     IPAROI : paroi standard
CFONC     ISYMET : symetrie standard
CFONC
CFONC     IESICF, ISSPCF, ISOPCF, IERUCF, IEQHCF : entree/sortie
CFONC
CFONC Pour les entrees/sorties, on peut
CFONC   imposer une valeur pour la turbulence et les scalaires
CFONC   passifs dans RCODCL(.,.,1) pour le cas ou le flux de masse
CFONC   serait entrant. Si on ne le fait pas, une condition de
CFONC   nul est appliquee.
CFONC
CFONC IESICF : entree sortie imposee (par exemple entree supersonique)
CFONC          l'utilisateur impose la vitesse et toutes les
CFONC            variables thermodynamiques
CFONC ISSPCF : sortie supersonique
CFONC          l'utilisateur n'impose rien
CFONC ISOPCF : sortie subsonique a pression imposee
CFONC          l'utilisateur impose la pression
CFONC IERUCF : entree subsonique a vitesse et rho imposes
CFONC          l'utilisateur impose la vitesse et la masse volumique
CFONC IEQHCF : entree subsonique a debit et debit enthalpique imposes
CFONC          a implementer
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
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
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
CARGU ! ICODCL       ! TE ! <-  ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
CARGU ! ITRIFB(NFABOR! TE !  -> ! INDIRECTION POUR TRI DES FACES DE BRD!
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! ITYPFB(NFABOR! TE ! <-  ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! IZFPPP       ! TE ! <-  ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU ! (NFABOR)     !    !     !  POUR LE MODULE PHYS. PART.          !
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
CARGU ! COEFA, COEFB ! TR ! <-> ! CONDITIONS AUX LIMITES AUX           !
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
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! COEFU        ! TR !  -  ! TAB DE TRAV                          !
CARGU !  (NFABOR,3)  !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
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
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
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
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ITRIFB(NFABOR,NPHAS), ITYPFB(NFABOR,NPHAS)
      INTEGER          IZFPPP(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION COEFU(NFABOR,NDIM)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC  , IEL   , II    , IPHAS
      INTEGER          ICOUL , IZONE , IFML  , IUTILE
      INTEGER          ILELT, NLELT
C
      DOUBLE PRECISION UREF2 , DHYD  , RHOMOY
      DOUBLE PRECISION USTAR2, XKENT , XEENT , D2S3
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
        WRITE(NFECRA,9001)
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
 9001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,
     &'@    *********                                               ',/,
     &'@     COMPRESSIBLE                                           ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR uscfcl DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 1.  INITIALISATIONS
C
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      D2S3 = 2.D0/3.D0
C
C=======================================================================
C 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
C       ON BOUCLE SUR LES FACES DE BORD
C         ON DETERMINE LA FAMILLE ET SES PROPRIETES
C           ON IMPOSE LA CONDITION LIMITE
C
C          IMPOSER ICI LES CONDITIONS LIMITES SUR LES FACES DE BORD
C
C=======================================================================
C
C      FAMILLE ET PREMIERE PROPRIETE
C
C        En general, la couleur des faces est la premiere propriete :
C          ICOUL = IPRFML(IFML,1).
C
C        Dans le cas plus specifique de maillages contenant des groupes
C          (ensembles dont l'intersection peut ne pas etre vide),
C          il est possible que l'on obtienne plusieurs proprietes par
C          famille. Ainsi, une face appartenant aux groupes 3 et 5 et
C          portant la couleur 8 sera dans une famille dont les
C          proprietes seront -3, -5 et 8 (ordre arbitraire).
C          Il faut donc, si l'on veut savoir si la face IFAC porte
C          la couleur 8, balayer toutes les proprietes de la
C          famille :
C            IFML  = IFMFBR(IFAC   )
C            DO IPROP = 1, NPRFML
C              IF(IPRFML(IFML,IPROP).EQ.8) THEN
C                La face porte la couleur 8 ...
C              ENDIF
C            ENDDO
C
C
      IPHAS = 1
C
C
C --- Exemple d'entree/sortie pour laquelle tout est connu
C
C       sans presumer du caractere subsonique ou supersonique,
C         l'utilisateur souhaite imposer toutes les caracteristiques
C         de l'ecoulement
C       une entree supersonique est un cas particulier
C
C       La turbulence et les scalaires utilisateur prennent un flux nul
C         si la vitesse est sortante.
C
      CALL GETFBR('1 and X <= 1.0 ',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C      ELEMENT ADJACENT A LA FACE DE BORD
C
        IEL = IFABOR(IFAC)
C
C       On numerote les zones de 1 a n...
        IZONE = 1
        IZFPPP(IFAC) = IZONE
C
        ITYPFB(IFAC,IPHAS) = IESICF
C
C   - Vitesse
        RCODCL(IFAC,IU(IPHAS),1) = 5.0D0
        RCODCL(IFAC,IV(IPHAS),1) = 0.0D0
        RCODCL(IFAC,IW(IPHAS),1) = 0.0D0
C
C   - Pression, Masse Volumique, Temperature, Energie Totale Specifique
C
C       Seules 2 variables sur les 4 sont independantes
C       On peut donc fixer le couple de variables que l'on veut
C       (sauf Temperature-Energie) et les 2 autres variables seront
C       calculees automatiquement
C
C    ** Choisir les 2 variables a imposer et effacer les autres
C       (elles sont calculées par les lois thermodynamiques dans uscfth)
C
C     Pression (en Pa)
        RCODCL(IFAC,IPR(IPHAS),1) = 5.D5
C
C     Masse Volumique (en kg/m3)
c        RCODCL(IFAC,ISCA(IRHO  (IPHAS)),1) = 1.D0
C
C     Temperature (en K)
        RCODCL(IFAC,ISCA(ITEMPK(IPHAS)),1) = 300.D0
C
C     Energie Totale Specifique (en J/kg)
c        RCODCL(IFAC,ISCA(IENERG(IPHAS)),1) = 355.D3
C
C
C   - Turbulence
C
        UREF2 = RCODCL(IFAC,IU(IPHAS),1)**2
     &         +RCODCL(IFAC,IV(IPHAS),1)**2
     &         +RCODCL(IFAC,IW(IPHAS),1)**2
        UREF2 = MAX(UREF2,1.D-12)
C
C
C       Exemple de turbulence calculee a partir
C         de formules valables pour une conduite
C
C       On veillera a specifier le diametre hydraulique
C         adapte a l'entree courante.
C
C       On s'attachera egalement a utiliser si besoin une formule
C         plus precise pour la viscosite dynamique utilisee dans le
C         calcul du nombre de Reynolds (en particulier, lorsqu'elle
C         est variable, il peut etre utile de reprendre ici la loi
C         imposee dans USCFPV. On utilise ici par defaut la valeur
C         VISCL0 donnee dans USINI1
C       En ce qui concerne la masse volumique, on peut utiliser directement
C         sa valeur aux faces de bord si elle est connue (et imposée
C         ci-dessus). Dans le cas général, on propose d'utiliser la valeur
C         de la cellule adjacente.
C
C         Diametre hydraulique
        DHYD   = 0.075D0
C
C         Vitesse de frottement (au carre = USTAR2)
        RHOMOY = RTP(IEL,ISCA(IRHO(IPHAS)))
        USTAR2 = 0.D0
        XKENT  = EPZERO
        XEENT  = EPZERO
C
        CALL KEENDB
C       ===========
     &    ( UREF2, DHYD, RHOMOY, VISCL0(IPHAS), CMU, XKAPPA,
     &      USTAR2, XKENT, XEENT )
C
        IF    (ITYTUR(IPHAS).EQ.2) THEN
C
          RCODCL(IFAC,IK(IPHAS),1)  = XKENT
          RCODCL(IFAC,IEP(IPHAS),1) = XEENT
C
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
          RCODCL(IFAC,IR11(IPHAS),1) = 2.D0/3.D0*XKENT
          RCODCL(IFAC,IR22(IPHAS),1) = 2.D0/3.D0*XKENT
          RCODCL(IFAC,IR33(IPHAS),1) = 2.D0/3.D0*XKENT
          RCODCL(IFAC,IR12(IPHAS),1) = 0.D0
          RCODCL(IFAC,IR13(IPHAS),1) = 0.D0
          RCODCL(IFAC,IR23(IPHAS),1) = 0.D0
          RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
C
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
C
          RCODCL(IFAC,IK(IPHAS),1)   = XKENT
          RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
          RCODCL(IFAC,IPHI(IPHAS),1) = D2S3
          RCODCL(IFAC,IFB(IPHAS),1)  = 0.D0
C
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
C
          RCODCL(IFAC,IK(IPHAS),1)   = XKENT
          RCODCL(IFAC,IOMG(IPHAS),1) = XEENT/CMU/XKENT
C
        ENDIF
C
C   - On traite les scalaires rattaches a la phase courante
C       (ne pas boucler sur NSCAL sous peine de modifier rho et energie)
        IF(NSCAUS.GT.0) THEN
          DO II = 1, NSCAUS
            IF(IPHSCA(II).EQ.IPHAS) THEN
              RCODCL(IFAC,ISCA(II),1) = 1.D0
            ENDIF
          ENDDO
        ENDIF
C
      ENDDO
C
C --- Exemple de sortie supersonique
C
C       toutes les caractéristiques sortent
C       on ne doit rien imposer (ce sont les valeurs internes qui sont
C         utilisees pour le calcul des flux de bord)
C
C       pour la turbulence et les scalaires, si on fournit ici des
C         valeurs de RCODCL, on les impose en Dirichlet si le flux
C         de masse est entrant ; sinon, on impose un flux nul (flux de
C         masse sortant ou RCODCL renseigné ici).
C         Noter que pour la turbulence, il faut renseigner RCODCL pour
C         toutes les variables turbulentes (sinon, on applique un flux
C         nul).
C
C
      CALL GETFBR('2',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C       On numerote les zones de 1 a n...
        IZONE = 2
        IZFPPP(IFAC) = IZONE
C
        ITYPFB(IFAC,IPHAS) = ISSPCF
C
      ENDDO
C
C --- Exemple d'entree subsonique (debit, debit enthalpique)
C
C       2 caracteristiques sur 3 entrent : il faut donner 2 informations
C         la troisieme est deduite par un scenario de 2-contact et
C         3-detente dans le domaine
C       ici on choisit de donner (rho*(U.n), rho*(U.n)*H)
C         avec H = 1/2 U*U + P/rho + e
C              n la normale unitaire entrante
C
C       ATTENTION, on donne des DENSITES de debit (par unite de surface)
C
      CALL GETFBR('3',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C       On numerote les zones de 1 a n...
        IZONE = 3
        IZFPPP(IFAC) = IZONE
C
        ITYPFB(IFAC,IPHAS) = IEQHCF
C
C   - Densite de debit massique (en kg/(m2 s))
        RCODCL(IFAC,IRUN(IPHAS),1) = 5.D5
C
C   - Densite de debit enthalpique (en J/(m2 s))
        RCODCL(IFAC,IRUNH(IPHAS),1) = 5.D5
C
C
C     Condition non disponible dans la version presente
        CALL CSEXIT (1)
C       ===========
C
      ENDDO
C
C --- Exemple d'entree subsonique (masse volumique, vitesse)
C
C       2 caracteristiques sur 3 entrent : il faut donner 2 informations
C         la troisieme est deduite par un scenario de 2-contact et
C         3-detente dans le domaine
C       ici on choisit de donner (rho, U)
C
      CALL GETFBR('4',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C       On numerote les zones de 1 a n...
        IZONE = 4
        IZFPPP(IFAC) = IZONE
C
        ITYPFB(IFAC,IPHAS) = IERUCF
C
C   - Vitesse d'entree
        RCODCL(IFAC,IU(IPHAS),1) = 5.0D0
        RCODCL(IFAC,IV(IPHAS),1) = 0.0D0
        RCODCL(IFAC,IW(IPHAS),1) = 0.0D0
C
C   - Masse Volumique (en kg/m3)
        RCODCL(IFAC,ISCA(IRHO  (IPHAS)),1) = 1.D0
C
C
C   - Turbulence
C
        UREF2 = RCODCL(IFAC,IU(IPHAS),1)**2
     &         +RCODCL(IFAC,IV(IPHAS),1)**2
     &         +RCODCL(IFAC,IW(IPHAS),1)**2
        UREF2 = MAX(UREF2,1.D-12)
C
C
C       Exemple de turbulence calculee a partir
C         de formules valables pour une conduite
C
C       On veillera a specifier le diametre hydraulique
C         adapte a l'entree courante.
C
C       On s'attachera egalement a utiliser si besoin une formule
C         plus precise pour la viscosite dynamique utilisee dans le
C         calcul du nombre de Reynolds (en particulier, lorsqu'elle
C         est variable, il peut etre utile de reprendre ici la loi
C         imposee dans USCFPV. On utilise ici par defaut la valeur
C         VISCL0 donnee dans USINI1
C       En ce qui concerne la masse volumique, on peut utiliser directement
C         sa valeur aux faces de bord si elle est connue (et imposée
C         ci-dessus). Dans le cas général, on propose d'utiliser la valeur
C         de la cellule adjacente.
C
C         Diametre hydraulique
        DHYD   = 0.075D0
C
C         Calcul de la vitesse de frottement au carre (USTAR2)
C           et de k et epsilon en entree (XKENT et XEENT) a partir
C           de lois standards en conduite circulaire
C           (leur initialisation est inutile mais plus propre)
        RHOMOY = PROPFB(IFAC,IPPROB(IROM(IPHAS)))
        USTAR2 = 0.D0
        XKENT  = EPZERO
        XEENT  = EPZERO
C
        CALL KEENDB
C       ===========
     &    ( UREF2, DHYD, RHOMOY, VISCL0(IPHAS), CMU, XKAPPA,
     &      USTAR2, XKENT, XEENT )
C
        IF    (ITYTUR(IPHAS).EQ.2) THEN
C
          RCODCL(IFAC,IK(IPHAS),1)  = XKENT
          RCODCL(IFAC,IEP(IPHAS),1) = XEENT
C
        ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
          RCODCL(IFAC,IR11(IPHAS),1) = 2.D0/3.D0*XKENT
          RCODCL(IFAC,IR22(IPHAS),1) = 2.D0/3.D0*XKENT
          RCODCL(IFAC,IR33(IPHAS),1) = 2.D0/3.D0*XKENT
          RCODCL(IFAC,IR12(IPHAS),1) = 0.D0
          RCODCL(IFAC,IR13(IPHAS),1) = 0.D0
          RCODCL(IFAC,IR23(IPHAS),1) = 0.D0
          RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
C
        ELSEIF(ITURB(IPHAS).EQ.50) THEN
C
          RCODCL(IFAC,IK(IPHAS),1)   = XKENT
          RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
          RCODCL(IFAC,IPHI(IPHAS),1) = D2S3
          RCODCL(IFAC,IFB(IPHAS),1)  = 0.D0
C
        ELSEIF(ITURB(IPHAS).EQ.60) THEN
C
          RCODCL(IFAC,IK(IPHAS),1)   = XKENT
          RCODCL(IFAC,IOMG(IPHAS),1) = XEENT/CMU/XKENT
C
        ENDIF
C
C   - On traite les scalaires rattaches a la phase courante
C       (ne pas boucler sur NSCAL sous peine de modifier rho et energie)
        IF(NSCAUS.GT.0) THEN
          DO II = 1, NSCAUS
            IF(IPHSCA(II).EQ.IPHAS) THEN
              RCODCL(IFAC,ISCA(II),1) = 1.D0
            ENDIF
          ENDDO
        ENDIF
C
      ENDDO
C
C --- Exemple de sortie subsonique
C
C       1 caracteristique sur 3 sort : il faut donner 1 information
C         les deux autres sont deduites par un scenario de 2-contact et
C         3-detente dans le domaine
C       ici on choisit de donner P
C
C       La turbulence et les scalaires utilisateur prennent un flux nul.
C
      CALL GETFBR('5',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C       On numerote les zones de 1 a n...
        IZONE = 5
        IZFPPP(IFAC) = IZONE
C
        ITYPFB(IFAC,IPHAS) = ISOPCF
C
C     Pression (en Pa)
        RCODCL(IFAC,IPR(IPHAS),1) = 5.D5
C
      ENDDO
C
C --- Exemple de paroi
C
      CALL GETFBR('7',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C       On numerote les zones de 1 a n...
        IZONE = 7
        IZFPPP(IFAC) = IZONE
C
        ITYPFB(IFAC,IPHAS) = IPAROI
C
C
C     Paroi défilante
C       Par défaut, la paroi n'est pas défilante
C       Si la paroi est défilante, donner les composantes non nulles
C         de la vitesse. La vitesse sera projetée dans le plan
C         tangeant à la paroi. Dans l'exemple suivant, on impose
C         Ux = 1. (l'exemple est activé si IUTILE=1)
C
        IUTILE = 0
        IF(IUTILE.EQ.1) THEN
          RCODCL(IFAC,IU(IPHAS),1) = 1.D0
        ENDIF
C
C     Température imposée
C       Par défaut, la paroi est adiabatique
C       Si la paroi est à température imposée, l'indiquer par
C         ICODCL = 5 et donner la valeur en Kelvin dans RCODCL(.,.,1)
C         Dans l'exemple suivant, on impose T = 293.15 K (l'exemple
C         est activé si IUTILE=1)
C
        IUTILE = 0
        IF(IUTILE.EQ.1) THEN
          ICODCL(IFAC,ISCA(ITEMPK(IPHAS)))   = 5
          RCODCL(IFAC,ISCA(ITEMPK(IPHAS)),1) = 20.D0 + 273.15D0
        ENDIF
C
C     Flux imposé
C       Par défaut, la paroi est adiabatique
C       Si la paroi est à flux imposé, l'indiquer par
C         ICODCL = 3 et donner la valeur en Watt/m2 dans RCODCL(.,.,3)
C         Dans l'exemple suivant, on impose un flux de 1000 W/m2
C         - la plage en été - (l'exemple est activé si IUTILE=1)
C
        IUTILE = 0
        IF(IUTILE.EQ.1) THEN
          ICODCL(IFAC,ISCA(ITEMPK(IPHAS)))   = 3
          RCODCL(IFAC,ISCA(ITEMPK(IPHAS)),3) = 1000.D0
        ENDIF
C
      ENDDO
C
C --- Exemple de symetrie
C
      CALL GETFBR('8',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C       On numerote les zones de 1 a n...
        IZONE = 8
        IZFPPP(IFAC) = IZONE
C
        ITYPFB(IFAC,IPHAS) = ISYMET
C
      ENDDO
C
C     Il est deconseille d'utiliser d'autres types de conditions
C     aux limites que ceux proposes ci-dessus.
C
C----
C FORMATS
C----
C
C----
C FIN
C----
C
      RETURN
      END
c@z
