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
                        SUBROUTINE USKPDC
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   NCEPDP , NCKPDP , IPHAS  , IAPPEL ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , ICEPDC ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , RTP    , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC                     PERTES DE CHARGE (PDC)
CFONC                        POUR LA PHASE IPHAS
CFONC
CFONC IAPPEL = 1 :
CFONC              CALCUL DU NOMBRE DE CELLULES OU L'ON IMPOSE UNE PDC
CFONC              DETERMINATION DU NOMBRE DE COEF DU TENSEUR DE PDC
CFONC IAPPEL = 2 :
CFONC              REPERAGE DES CELLULES OU L'ON IMPOSE UNE PDC
CFONC IAPPEL = 3 :
CFONC              CALCUL DES VALEURS DES COEFS DE PDC
CFONC
CFONC
CFONC CKUPDC EST LE COEFF DE PDC CALCULE.
CFONC
CFONC   IL INTERVIENT DANS LA QDM COMME SUIT :
CFONC     RHO DU/DT = - GRAD P + TSPDC        (+ AUTRES TERMES)
CFONC                       AVEC TSPDC = - RHO CKUPDC U ( en kg/(m2 s))
CFONC
CFONC
CFONC   POUR UNE PDC REPARTIE,
CFONC
CFONC     SOIT KSIL = DHL/(0.5 RHO U**2) DONNE DANS LA LITTERATURE
CFONC     (DHL EST LA PERTE DE CHARGE PAR UNITE DE LONGUEUR)
CFONC
CFONC     LE TERME SOURCE TSPDC VAUT DHL = - KSIL *(0.5 RHO U**2)
CFONC
CFONC     ON A CKUPDC = 0.5 KSIL ABS(U)
CFONC
CFONC
CFONC   POUR UNE PDC SINGULIERE,
CFONC
CFONC     SOIT KSIS = DHS/(0.5 RHO U**2) DONNE DANS LA LITTERATURE
CFONC     (DHS EST LA PERTE DE CHARGE SINGULIERE)
CFONC
CFONC     LE TERME SOURCE TSPDC VAUT DHS/L = - KSIS/L *(0.5 RHO U**2)
CFONC
CFONC     ON A CKUPDC = 0.5 KSIS/L ABS(U)
CFONC
CFONC     OU L DESIGNE LA LONGUEUR SUR LAQUELLE
CFONC            ON A CHOISI DE REPRESENTER LA ZONE DE PDC SINGULIERE
CFONC
CFONC
CFONC IDENTIFICATION DES CELLULES/FACES DE BORD/FACES INTERNES
CFONC ========================================================
CFONC
CFONC   Les commandes GETCEL, GETFBR et GETFAC permettent d'identifier
CFONC   respectivement les cellules, faces ou faces de bord en fonction
CFONC   de differents criteres.
CFONC
CFONC   GETCEL(CHAINE,NLELT,LSTELT) :
CFONC   - CHAINE est une chaine de caractere fournie par l'utilisateur
CFONC     qui donne les criteres de selection
CFONC   - NLTELT est renvoye par la commande. C'est un entier qui
CFONC     correspond au nombre de cellules trouveees repondant au
CFONC     critere
CFONC   - LSTELT est renvoye par la commande. C'est un tableau d'entiers
CFONC     de taille NLTELT donnant la liste des cellules trouvees
CFONC     repondant au critere.
CFONC
CFONC   CHAINE peut etre constitue de :
CFONC   - references de couleurs (ex. : 1, 8, 26, ...
CFONC   - references de groupes (ex. : entrees, groupe1, ...)
CFONC   - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
CFONC   Ces criteres peuvent etre combines par des operateurs logiques
CFONC   (AND et OR) et des parentheses
CFONC   ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
CFONC   les cellules de couleur 1, appartenant aux groupes 'groupe2'
CFONC   ou 'groupe3' et de coordonnee Y inferieure a 1.
CFONC
CFONC   La syntaxe des commandes GETFBR et GETFAC est identique.
CFONC
CFONC
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
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NCEPDP       ! E  ! <-> ! NOMBRE DE CELLULES AVEC PDC          !
CARGU ! NCKPDP       ! E  ! <-> ! NBR DE COEF DU TENSEUR DE PDC (3 OU 6!
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
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
CARGU ! ICEPDC(NCEPDP! TE ! <-> ! NUMERO DES NCEPDP CELLULES AVEC PDC  !
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
CARGU ! RTP, RTPA    ! TR ! <-> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! CKUPDC(NCEPDP! TR ! <-> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !     , NCKPDP)!    !     !                                      !
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
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
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
      INTEGER          NCEPDP , NCKPDP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
      INTEGER          IAPPEL
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICEPDC(NCEPDP)
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
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL, IELPDC, IPHAS, IKPDC
      INTEGER          ILELT, NLELT
C
      DOUBLE PRECISION ALPHA, COSALP, SINALP, VIT, CK1, CK2
C
C***********************************************************************
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IF(IAPPEL.EQ.1.OR.IAPPEL.EQ.2) THEN
C
C=======================================================================
C
C 1. POUR CHAQUE PHASE : UN OU DEUX APPELS
C
C      PREMIER APPEL :
C
C        IAPPEL = 1 : NCEPDP : CALCUL DU NOMBRE DE CELLULES
C                                AVEC PERTES DE CHARGE
C                     NCKPDP : NOMBRE DE COEFFS DU TENSEUR
C
C
C      DEUXIEME APPEL (POUR LES PHASES AVEC NCEPDP > 0) :
C
C        IAPPEL = 2 : ICEPDC : REPERAGE DU NUMERO DES CELLULES
C                                AVEC PERTES DE CHARGE
C
C REMARQUES :
C
C        Ne pas utiliser CKUPDC dans cette section
C          (il est rempli au troisieme appel, IAPPEL = 3)
C
C        Ne pas utiliser ICEPDC dans cette section
C           au premier appel (IAPPEL = 1)
C
C        On passe ici a chaque pas de temps
C           (ATTENTION au cout calcul de vos developpements)
C
C=======================================================================
C
C ---------------------------------------------------------------------
C  1.1 A completer par l'utilisateur : selection des cellules
C  -----------------------------------------------------------
C
C --- Exemple 1 : Aucune pdc (defaut)
C
        IELPDC = 0
C
C
C --- Exemple 2 : Pdc definies par coordonnees pour la phase 1
C                 Pas de pertes de charge pour la phase 2
C                 Le traitement etant different pour les phases
C                   un test est necessaire.
C
C       Ce test permet de desactiver l'exemple
        IF(1.EQ.0) THEN
          IF(IPHAS.EQ.1) THEN
            IELPDC = 0
C
            CALL GETCEL('X <= 6.0 and X >= 4.0 and Y >= 2.0 and'//
     &                  'Y <= 8.0',NLELT,LSTELT)
C
            DO ILELT = 1, NLELT
              IEL = LSTELT(ILELT)
              IELPDC = IELPDC + 1
              IF (IAPPEL.EQ.2) ICEPDC(IELPDC) = IEL
            ENDDO
C
          ELSE
            IELPDC = 0
          ENDIF
        ENDIF
C
C ---------------------------------------------------------------------
C  1.2 Sous section generique a ne pas modifier
C  ---------------------------------------------
C
C --- Pour IAPPEL = 1,
C      Renseigner NCEPDP, nombre de cellules avec pdc
C      Le bloc ci dessous est valable pourles 2 exemples ci dessus
C
        IF (IAPPEL.EQ.1) THEN
          NCEPDP = IELPDC
        ENDIF
C
C ---------------------------------------------------------------------
C  1.3 A completer par l'utilisateur : nombre de coeffs
C  -----------------------------------------------------
C
C --- Pour IAPPEL = 1,
C      Renseigner l'entier NCKPDP qui indique si le tenseur de pdc est
C        diagonal (NCKPDP=3) ou non (NCKPDP=6)
C      Utile uniquement si NCEPDP > 0. Si NCEPDP = 0, la valeur de
C        NCKPDP est ignoree.
C      Le bloc ci dessous est valable pourles 2 exemples ci dessus
C
        IF (IAPPEL.EQ.1) THEN
C
          NCKPDP = 0
          IF(IPHAS.EQ.1) THEN
            NCKPDP = 6
          ELSE
            NCKPDP = 3
          ENDIF
C
        ENDIF
C
C-----------------------------------------------------------------------
C
      ELSEIF(IAPPEL.EQ.3) THEN
C
C=======================================================================
C
C 2. POUR CHAQUE PHASE AVEC NCEPDP > 0 , TROISIEME APPEL
C
C      TROISIEME APPEL (POUR LES PHASES AVEC NCEPDP > 0) :
C
C       IAPPEL = 3 : CKUPDC : CALCUL DES COEFFICIENTS DE PERTE DE CHARGE
C                             DANS LE REPERE DE CALCUL
C                             STOCKES DANS L'ORDRE
C                             K11, K22, K33, K12, K13, K23
C
C
C    REMARQUE :
C
C        Veillez a ce que les coefs diagonaux soient positifs.
C
C        Vous risquez un PLANTAGE si ce n'est pas le cas.
C
C        AUCUN controle ulterieur ne sera effectue.
C
C      ===========================================================
C
C ---------------------------------------------------------------------
C  2.1 A completer par l'utilisateur : valeur des coefs
C  -----------------------------------------------------
C
C --- Attention
C   Il est important que les CKUPDC soient completes (par des valeurs
C     nulles eventuellement) dans la mesure ou ils seront utilises pour
C     calculer un terme source dans les cellules identifiees precedemment.
C
C   On les initialise tous par des valeurs nulles.
C       Et on demande a l'utilisateur de conserver cette initialisation.
C                                        =========
C
        IF(NCKPDP.EQ.3.OR.NCKPDP.EQ.6) THEN
          DO IKPDC = 1, NCKPDP
            DO IELPDC = 1, NCEPDP
              CKUPDC(IELPDC,IKPDC) = 0.D0
            ENDDO
          ENDDO
        ENDIF
C
        IF(IPHAS.EQ.1) THEN
C
C --- Tenseur diagonal
C   Exemple de pertes de charges dans la direction x
C
C
          IF(NCKPDP.EQ.3) THEN
            DO IELPDC = 1, NCEPDP
              IEL=ICEPDC(IELPDC)
              VIT = SQRT(  RTPA(IEL,IU(IPHAS))**2
     &                   + RTPA(IEL,IV(IPHAS))**2
     &                   + RTPA(IEL,IW(IPHAS))**2)
              CKUPDC(IELPDC,1) = 10.D0*VIT
              CKUPDC(IELPDC,2) =  0.D0*VIT
              CKUPDC(IELPDC,3) =  0.D0*VIT
            ENDDO
C
C --- Tenseur 3x3
C   Exemple de pertes de charges a ALPHA = 45 degres x,y
C      la direction x resiste par ck1 et y par ck2
C      ck2 nul represente des ailettes comme ceci :  ///////
C      dans le repere de calcul X Y
C
C                 Y|    /y
C                  |  /
C                  |/
C                  \--------------- X
C                   \ / ALPHA
C                    \
C                     \ x
C
          ELSEIF(NCKPDP.EQ.6) THEN
C
            ALPHA  = PI/4.D0
            COSALP = COS(ALPHA)
            SINALP = SIN(ALPHA)
            CK1 = 10.D0
            CK2 =  0.D0
C
            DO IELPDC = 1, NCEPDP
              IEL=ICEPDC(IELPDC)
              VIT = SQRT(  RTPA(IEL,IU(IPHAS))**2
     &                   + RTPA(IEL,IV(IPHAS))**2
     &                   + RTPA(IEL,IW(IPHAS))**2)
              CKUPDC(IELPDC,1) = (COSALP**2*CK1 + SINALP**2*CK2)*VIT
              CKUPDC(IELPDC,2) = (SINALP**2*CK1 + COSALP**2*CK2)*VIT
              CKUPDC(IELPDC,3) =  0.D0
              CKUPDC(IELPDC,4) = COSALP*SINALP*(-CK1+CK2)*VIT
              CKUPDC(IELPDC,5) =  0.D0
              CKUPDC(IELPDC,6) =  0.D0
            ENDDO
C
          ENDIF
C
        ENDIF
C
C
C-----------------------------------------------------------------------
C
      ENDIF
C
      RETURN
C
      END
c@z
