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
                        SUBROUTINE USVORT
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IREPVO ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC METHODE DES VORTEX POUR LES CONDITIONS AUX LIMITES D'ENTREE
CFONC   EN L.E.S. :
CFONC   DEFINITION DES ENTREES AVEC VORTEX
CFONC   DEFINITION DES CARACTERISTIQUES DES VORTEX
CFONC
CFONC
CFONC IDENTIFICATION DES FACES DE BORD
CFONC ================================
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
C                                  ARGUMENTS
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
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE                   !
CARGU ! IAPPEL       ! E  !  -> ! INDIQUE LES DONNES A RENVOYER        !
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
CARGU ! IREPVO       ! TE ! <-> ! NUMERO DE L'ENTREE ASSOCIE A CHAQUE  !
CARGU !     (NFABOR) !    !     ! FACE DE BORD (=0 SI PAS DE VORTEX)   !
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
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
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
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
      INCLUDE "vortex.h"
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
      INTEGER          IPHAS  , IAPPEL
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IREPVO(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IFAC, IENT, IEL, ICOUL, IFML
      INTEGER          ILELT, NLELT
C
      DOUBLE PRECISION XX
C
C=======================================================================
C 1. PARAMETRES GLOBAUX
C=======================================================================
C
C --- Nombre d'entrees avec la methode des vortex
C
      NNENT = 2
C
C --- Nombre de vortex a mettre dans chaque entree
C
C   NVORT min ~ Surface d'entree/(pi*SIGMA**2)
C
      NVORT(1) = 500
      NVORT(2) = 500
C
      IF (IAPPEL.EQ.1) THEN
C
C=======================================================================
C 2. DEFINITION DES ZONES D'ENTREE (AU PREMIER PASSAGE)
C=======================================================================
C
        DO IFAC = 1, NFABOR
          IREPVO(IFAC) = 0
        ENDDO
C
C ------------------
C   ENTREE 1
C ------------------
        CALL GETFBR('3',NLELT,LSTELT)
C       ===========
C
        DO ILELT = 1, NLELT
C
          IFAC = LSTELT(ILELT)
C
          IENT = 1
          IREPVO(IFAC) = IENT
C
        ENDDO
C
C ------------------
C   ENTREE 2
C ------------------
        CALL GETFBR('1',NLELT,LSTELT)
C       ===========
C
        DO ILELT = 1, NLELT
C
          IFAC = LSTELT(ILELT)
C
          IENT = 2
          IREPVO(IFAC) = IENT
C
        ENDDO
C
      ELSEIF (IAPPEL.EQ.2) THEN
C
C=======================================================================
C 3. PARAMETRES GEOMETRIQUES ET CONDITIONS LIMITES
C=======================================================================
C
C --- Cas traité
C
C ICAS = 1...Conduite rectangulaire
C        2...Conduite circulaire
C        3...Geometrie quelconque sans traitement specifique des conditions aux limites
C        4...Geometrie quelconque sans traitement specifique des conditions aux limites
C            ni fichier de donnees (la vitesse moyenne, le niveau de k et de epsilon
C            sont fournis par l'utilisateur)
C
        IENT = 1
        ICAS(IENT) = 1
C
        IENT = 2
        ICAS(IENT) = 2
C
C
C --- Repere definissant le plan d'entree
C
C     Si ICAS = 4, le code se charge de ces donnees
C     Sinon il faut preciser les vecteurs DIR1 et DIR2 definissant
C     un repère directe tel que DIR3 soit un vecteur entrant normal
C     a la face d'entree.
C
        IENT = 1
        IF(ICAS(IENT).EQ.1.OR.ICAS(IENT).EQ.2.OR.ICAS(IENT).EQ.3) THEN
          DIR1(1,IENT) = 1.D0
          DIR1(2,IENT) = 0.D0
          DIR1(3,IENT) = 0.D0
C
          DIR2(1,IENT) = 0.D0
          DIR2(2,IENT) = 1.D0
          DIR2(3,IENT) = 0.D0
        ENDIF
C
        IENT = 2
        IF(ICAS(IENT).EQ.1.OR.ICAS(IENT).EQ.2.OR.ICAS(IENT).EQ.3) THEN
          DIR1(1,IENT) = 0.D0
          DIR1(2,IENT) = 1.D0
          DIR1(3,IENT) = 0.D0
C
          DIR2(1,IENT) = 0.D0
          DIR2(2,IENT) = 0.D0
          DIR2(3,IENT) = 1.D0
        ENDIF
C
C --- Centre du repere local dans le plan d'entree
C
C     Si ICAS = 1 ou ICAS = 2, le centre du repere doit correspondre
C                 au centre de gravite de la zone d'entree (rectangle ou cercle)
C
        IENT = 1
C
        CEN(1,IENT) = 0.D0
        CEN(2,IENT) = 0.D0
        CEN(3,IENT) = -6.05D-1
C
        IENT = 2
C
        CEN(1,IENT) = -3.664D-1
        CEN(2,IENT) = 0.D0
        CEN(3,IENT) = 0.D0
C
C --- Condition aux limites
C
C -> Si ICAS = 1...Il faut specifier le type de condition aux limite ICLVOR
C               dans les directions DIR1, DIR2, - DIR1, -DIR2
C
C               Ces conditions peuvent etre de 3 types :
C
C ICLVOR = 1...Condition de paroi
C          2...Condition de symetrie
C          3...Condition de periodicite
C
C                    y = LLY/2
C                    (ICLVOR 1)
C           +-----------------------+
C           |           ^ DIR1      |
C           |           |           |
C           |           |           |
C z=- LLZ/2 |           +----> DIR2 | z = LLZ/2
C (ICLVOR 4)|                       | (ICLVOR 2)
C           |                       |
C           |                       |
C           +-----------------------+
C                    y = -LLY/2
C                    (ICLVOR 3)
C
C
C -> Si ICAS = 2, les conditions sont necessairement de type paroi
C -> Si ICAS = 3 ou 4, pas de traitement particulier
C
        IENT = 1
C
        IF(ICAS(IENT).EQ.1) THEN
          ICLVOR(1,IENT) = 1
          ICLVOR(2,IENT) = 2
          ICLVOR(3,IENT) = 1
          ICLVOR(4,IENT) = 2
        ENDIF
C
C LLY et LLZ sont les dimensions de l'entree dans les directions DIR1 et DIR2
C LDD est le diametre de la conduite
C
C
        IENT = 1
        LLY(IENT) = 0.2D0
        LLZ(IENT) = 0.1D0
C
        IENT = 2
        LLD(2) = 0.154D0
C
C=======================================================================
C 5. PARAMETRES PHYSIQUES ET MARCHE EN TEMPS
C=======================================================================
C
C --- " Temps de vie " limite du vortex
C
C ITLIVO = 1...Les vortex sont retire au bout du temps TLIMVO
C                donne par l'utilisateur
C                ( par exemple TLIMVO = 10*DTREF)
C
C          2...Chaque vortex a un temps d'exitence limite valant
C                5.Cmu.k^(3/2).U/epsilon
C               ( ou U est la vitesse principale suivant DIR3)
C
        IENT = 1
        ITLIVO(IENT) = 1
C
        IF(ITLIVO(IENT).EQ.1) THEN
          TLIMVO(IENT) = 10.D0*DTREF
        ENDIF
C
        IENT = 2
        ITLIVO(IENT) = 2
C
C
C --- " Diametre " des vortex
C
C ISGMVO = 1...diametre constant XSGMVO donne par l'utilisateur
C          2...basee sur la formule sigma = Cmu^(3/4).k^(3/2)/epsilon
C          3...basee sur la formule sigma = max(Lt, Lk) avec
C                 Lt = (5 nu.k/epsilon)^(1/2)
C             et  Lk = 200.(nu^3/epsilon)^(1/4)
C
        IENT = 1
        ISGMVO(IENT) = 1
C
        IF(ISGMVO(IENT).EQ.1) THEN
          XSGMVO(IENT) = 0.01D0
        ENDIF
C
        IENT = 2
        ISGMVO(IENT) = 2
C
C
C --- Mode de deplacement des vortex
C
C IDEPVO = 1...Deplacement en r*UD (r aleatoire dans [0,1])
C              UD a fournir par l'utilisateur
C          2...Convection par les vortex
C
        IENT = 1
        IDEPVO(IENT) = 2
C
        IENT = 2
        IDEPVO(IENT) = 1
C
        IF(IDEPVO(IENT).EQ.1) THEN
          UD(IENT) = 0.7D0
        ENDIF
C
C=======================================================================
C 6. PARAMETRES D'ENTREE / SORTIES ET DONNEES UTILISATEUR
C=======================================================================
C
C --- Fichier de donnees utilisateur
C
C NDAT ...Nombre de lignes du fichier de donnees contenant les donnees :
C          x | y | z | U | V | W | Grad[u.DIR3].n | k | epsilon
C
C         dans le plan d'entree du calcul
C
C         Grad[u.DIR3].n est le gradient dans la direction normale
C         a la paroi, de la vitesse principale dans le plan d'entree.
C         Cette données n'est utilisée qu'avec ICAS=2
C
C FICVOR...Nom du fichier de donnees utilisateur
C
        IENT = 1
        NDAT(IENT) = 2080
C
        IENT = 2
        NDAT(IENT) = 2080
C
C Par les defaut les fichiers sont nommes "vordat" affecté de l'indice
C d'entrée
C
        IENT = 1
        FICVOR(IENT) = 'entree_1.dat'
C
        IENT = 2
        FICVOR(IENT) = 'entree_2.dat'
C
C Pour ICAS = 4, on precise juste la valeur moyenne de U, k et de espilon
C a l'entree
C
        IF(ICAS(IENT).EQ.4) THEN
          UDEBIT(IENT) = 10.D0
          KDEBIT(IENT) = 1.D0
          EDEBIT(IENT) = 1.D0
        ENDIF
C
C --- Relecture d'un fichier suite eventuel
C
C ISUIVO = 0...Pas de relecture (reinitialisation des vortex)
C          1...Relecture du fichier suite de methode des vortex
C
        ISUIVO = ISUITE
C
C
      ENDIF
C
C
      RETURN
      END
C
c@z
C=======================================================================
C 7. DEFINTION DE LA FONCTION PERMETAT D'IMPOSER LES DONNEES D'ENTREE
C=======================================================================
C
c@a
                         FUNCTION PHIDAT
C                        ***************
C     -------------------------------------------------------------
     & ( NFECRA , ICAS   , NDAT   ,
     &   YY     , ZZ     , YDAT   , ZDAT   ,
     &   VARDAT , III    )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC FONCTION PERMETTANT D'INTERPOLER LES DONNEES D'ENTREE FOURNIES
CFONC PAR L'UTILISATEUR AU CENTRE DES FACES D'ENTREE POUR LESQUELLES
CFONC EST UTILISEE LA METHODE DES VORTEX
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NFECRA       ! E  !  -> ! UNITE                                !
CARGU ! ICAS         ! E  !  -> ! TYPE DE GEOMETRIE DU CAS             !
CARGU ! NDAT         ! E  !  -> ! NBR DE LIGNES DU FICHIER DE DONNEES  !
CARGU ! YY           ! E  !  -> ! COORDOONNES DANS LE REPERE LOCAL DU  !
CARGU ! ZZ           ! E  !  -> ! POINT OU L'ON CHERCHE A CONNAITRE LA !
CARGU !              !    !     ! VARIABLE VARDAT                      !
CARGU ! YDAT         ! E  !  -> ! COORDOONNES OU EST CONNUE LA VARIABLE!
CARGU ! ZDAT         ! E  !  -> ! VARDAT DANS LE FICHIER DE DONNEES    !
CARGU ! VARDAT       ! E  !  -> ! VALEUR DE LA VARIABLE VARDAT         !
CARGU ! III          ! E  ! <-  ! LIGNE OU A ETE TROUVEE LA DONNEE LA  !
CARGU !              !    !     ! PLUS PROCHE DU POINT (YY,ZZ)         !
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
      INTEGER          NFECRA, ICAS, NDAT, III
      DOUBLE PRECISION ZZ, YY
      DOUBLE PRECISION ZDAT(NDAT), YDAT(NDAT)
      DOUBLE PRECISION VARDAT(NDAT)
C
      INTEGER          II
      DOUBLE PRECISION PHIDAT, DIST1

C
C Dans l'exemple suivant, on se contente de retourne la valeur situee
C dans le fichier de donnee a l'abscisse la plus proche du point de
C coordonnée (Y,Z) ou l'on cherche a connaitre la valeur de la
C variable numero VARDAT.
C
C
      IF(ICAS.EQ.1.OR.ICAS.EQ.2.OR.ICAS.EQ.3) THEN
C
        IF(III.EQ.0) THEN
          DIST1 = 1.D20
          DO II = 1,NDAT
            IF(SQRT((YY-YDAT(II))**2+(ZZ-ZDAT(II))**2).LT.DIST1) THEN
              DIST1 = SQRT((ZZ-ZDAT(II))**2+(YY-YDAT(II))**2)
              III   = II
              PHIDAT = VARDAT(II)
            ENDIF
          ENDDO
        ELSEIF(III.NE.0) THEN
          PHIDAT =  VARDAT(III)
        ENDIF
C
      ELSEIF(ICAS.EQ.4) THEN
        PHIDAT = VARDAT(1)
      ENDIF
C
      RETURN
      END
c@z
