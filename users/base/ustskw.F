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
                        SUBROUTINE USTSKW
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , NCEPDP , NCKPDP , NCESMP ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICEPDC , ICETSM , ITYPSM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , CKUPDC , SMACEL , S2KW   , DIVUKW ,
     &   GKGW   , GGRHO  , XF1    ,
     &   CRKEXP , CRWEXP , CRKIMP , CRWIMP ,
     &   VISCF  , VISCB  , XAM    , W1     , W2     ,
     &   W3     , W4     , W5     , W6     , W7     ,
     &   RDEVEL , RTUSER , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC ROUTINE UTILISATEUR ON PRECISE LES TERMES SOURCES UTILISATEURS
CFONC    POUR K-OMEGA SUR UN PAS DE TEMPS (PHASE IPHAS)
CFONC
CFONC ON RESOUT RHO*VOLUME*D(VAR)/DT = CRVIMP*VAR + CRVEXP
CFONC
CFONC ON FOURNIT ICI CRVIMP ET CRVEXP (ILS CONTIENNENT RHO*VOLUME)
CFONC     CRVEXP en kg variable/s :
CFONC      ex : pour k                     kg.m2/s3
CFONC           pour omega                 kg/s2
CFONC     CRVIMP en kg/s :
CFONC
CFONC VEILLER A UTILISER UN CRVIMP NEGATIF
CFONC  (ON IMPLICITERA CRVIMP
CFONC   IE SUR LA DIAGONALE DE LA MATRICE, LE CODE AJOUTERA :
CFONC    MAX(-CRVIMP,0) EN SCHEMA STANDARD EN TEMPS
CFONC        -CRVIMP    SI LES TERMES SOURCES SONT A L'ORDRE 2
CFONC
CFONC CES TABLEAUX SONT INITIALISES A ZERO AVANT APPEL A CE SOUS
CFONC    PROGRAMME ET AJOUTES ENSUITE AUX TABLEAUX PRIS EN COMPTE
CFONC    POUR LA RESOLUTION
CFONC
CFONC  EN CAS D'ORDRE 2 DEMANDE SUR LES TERMES SOURCES, ON DOIT
CFONC    FOURNIR CR.EXP A L'INSTANT N     (IL SERA EXTRAPOLE) ET
CFONC            CR.IMP A L'INSTANT N+1/2 (IL EST  DANS LA MATRICE,
CFONC                                      ON LE SUPPOSE NEGATIF)
CFONC
CFONC S2KW   = 2 (S11)**2 + 2 (S22)**2 + 2 (S33)**2
CFONC        + (2 S12)**2 + (2 S13)**2 + (2 S23)**2 )
CFONC
CFONC      Sij = (dUi/dxj+dUj/dxi)/2
CFONC
CFONC DIVUKW = DUDX + DVDY + DWDZ
CFONC
CFONC GKGW   = GRAD(K).GRAD(OMEGA)
CFONC
CFONC GGRHO  = -G.GRAD(RHO)/PRDTUR/RHO SI IGRAKE>0
CFONC GGRHO  = 0                       SI IGRAKE=0
CFONC
CFONC XF1    = COEFFICIENT F1 DE PONDERATION K-EPS/K-OMEGA
CFONC
CFONC L'identification des cellules concernees peut s'appuyer sur la
CFONC commande GETCEL.
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
CARGU ! NCEPDP       ! E  !  -> ! NOMBRE DE CELLULES AVEC PDC          !
CARGU ! NCKPDP       ! E  !  -> ! NBR DE COEF DU TENSEUR DE PDC (3 OU 6!
CARGU ! NCESMP       ! E  !  -> ! NOMBRE DE CELLULES A SOURCE DE MASSE !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE PHASE                      !
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
CARGU ! ICEPDC(NCELET! TE !  -> ! NUMERO DES NCEPDP CELLULES AVEC PDC  !
CARGU ! ICETSM(NCESMP! TE !  -> ! NUMERO DES CELLULES A SOURCE DE MASSE!
CARGU ! ITYPSM       ! TE !  -> ! TYPE DE SOURCE DE MASSE POUR LES     !
CARGU ! (NCESMP,NVAR)!    !     !  VARIABLES (cf. USTSMA)              !
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
CARGU ! RTPA         ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT            PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! CKUPDC(NCEPDP! TR !  -> ! TABLEAU DE TRAVAIL POUR PDC          !
CARGU !     , NCKPDP)!    !     !                                      !
CARGU ! SMACEL       ! TR !  -> ! VALEUR DES VARIABLES ASSOCIEE A LA   !
CARGU ! (NCESMP,*   )!    !     !  SOURCE DE MASSE                     !
CARGU !              !    !     !  POUR IVAR=IPR, SMACEL=FLUX DE MASSE !
CARGU ! S2KW(NCELET) ! TR !  -> ! 2.SijSij  (EN 1/S2)                  !
CARGU ! DIVUKW(NCELET! TR !  -> ! DIV(U)     (EN 1/S )                 !
CARGU ! GKGW(NCELET) ! TR !  -> ! GRAD(K).GRAD(OMEGA)                  !
CARGU ! GGRHO(NCELET)! TR !  -> ! -G*GRAD(RHO)/PRDTUR/RHO SI IGRAKE>0  !
CARGU ! XF1(NCELET)  ! TR !  -> ! COEFFICIENT DE PONDERATION K-EPS/K-W !
CARGU ! CRKEXP(NCELET! TR ! <-  ! TABLEAU POUR K     POUR PART EXPLICIT!
CARGU ! CRKIMP(NCELET! TR ! <-  ! TABLEAU POUR K     POUR TERME INSTAT !
CARGU ! CRWEXP(NCELET! TR ! <-  ! TABLEAU POUR OMEGA POUR PART EXPLICIT!
CARGU ! CRWIMP(NCELET! TR ! <-  ! TABLEAU POUR OMEGA POUR TERME INSTAT !
CARGU ! VISCF(NFAC)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! XAM(NFAC,2)  ! TR !  -  ! TABLEAU DE TRAVAIL    FACES DE BORD  !
CARGU ! W1..10(NCELET! TR !  -  ! TABLEAU DE TRAVAIL    CELLULES       !
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
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
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
      INTEGER          NCEPDP , NCKPDP , NCESMP
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE , IPHAS
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML), MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICEPDC(NCEPDP)
      INTEGER          ICETSM(NCESMP), ITYPSM(NCESMP,NVAR)
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
      DOUBLE PRECISION CKUPDC(NCEPDP,NCKPDP), SMACEL(NCESMP,NVAR)
      DOUBLE PRECISION S2KW(NCELET)  , DIVUKW(NCELET)
      DOUBLE PRECISION GKGW(NCELET)  , GGRHO(NCELET), XF1(NCELET)
      DOUBLE PRECISION CRKEXP(NCELET), CRKIMP(NCELET)
      DOUBLE PRECISION CRWEXP(NCELET), CRWIMP(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR), XAM(NFAC,2)
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET), W7(NCELET)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IEL, IKIPH, IOMGIP, IPCROM
      DOUBLE PRECISION FF, TAU, XX
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Numero des variables k et epsilon de la phase IPHAS courante
      IKIPH  = IK (IPHAS)
      IOMGIP = IOMG(IPHAS)
C
C --- Numero des grandeurs physiques (voir usclim) : masse volumique
      IPCROM = IPPROC(IROM(IPHAS))
C
      IF(IWARNI(IKIPH).GE.1) THEN
        WRITE(NFECRA,1000) IPHAS
      ENDIF
C
C=======================================================================
C 2. EXEMPLE FICTIF :
C
C    Pour la phase 2
C
C      Terme source k :
C         rho volume d(k)/dt       = ...
C                        ... - rho*volume*ff*epsilon - rho*volume*k/tau
C
C      Terme source omega :
C         rho volume d(omega)/dt = ...
C                        ... + rho*volume*xx
C
C      Avec, pour l'exemple,
C                 xx = 2.d0, ff=3.d0, tau = 4.d0
C
C=======================================================================
C
C
C --- Si la phase est bien la phase selectionnee
C
      IF(IPHAS.NE.2) THEN
C
C
C --- Termes sources explicites
C
        FF  = 3.D0
        TAU = 4.D0
        XX  = 2.D0
C
        DO IEL = 1, NCEL
          CRKEXP(IEL) = -PROPCE(IEL,IPCROM)*VOLUME(IEL)
     &                                          *FF*RTPA(IEL,IOMGIP)
          CRWEXP(IEL) =  PROPCE(IEL,IPCROM)*VOLUME(IEL)*XX
        ENDDO
C
C --- Termes sources implicites (diagonale)
C        CRKIMP et CREIMP sont initialises a zero avant l'entree dans ce
C          sous-programme : il est donc inutile de completer CRWIMP
C          s'il est nul.
C
        DO IEL = 1, NCEL
          CRKIMP(IEL) = -PROPCE(IEL,IPCROM)*VOLUME(IEL)/TAU
        ENDDO
C
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(' TERMES SOURCES UTILISATEURS K-EPSILON PHASE ',I4,/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
