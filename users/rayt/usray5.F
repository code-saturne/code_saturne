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
                        SUBROUTINE USRAY5
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IZFRDP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   COFRUA , COFRUB ,
     &   W1     , W2     , W3     , W4     , W5     ,  W6     ,
     &   RDEVEL , RTUSER ,
     &   TPAROI , QINCID , FLUNET , XLAM   , EPA    , EPS     ,  CK   ,
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
CFONC    Ce sous-programme est appele 2 fois pour chaque phase IPHAS
CFONC    pour laquelle il faut faire un calcul de rayonnement
CFONC    semi-transparent.
CFONC
CFONC
CFONC
CFONC   1. PREMIER APPEL (IAPPEL = 1)
CFONC   =============================
CFONC
CFONC
CFONC
CFONC     1.1 Conditions aux limites pour la luminance
CFONC     --------------------------------------------
CFONC
CFONC        Il faut completer COFRUA qui fournit la luminance au bord
CFONC          selon le type de frontiere (condition de Dirichlet).
CFONC        La luminance est la puissance surfacique par unite
CFONC          d'angle solide
CFONC
CFONC
CFONC        Par exemple, on a
CFONC
CFONC
CFONC  1/ Paroi grise : rayonnement isotrope.
CFONC                                     4
CFONC                       eps.sig.tparoi         (1-eps).qincid
CFONC         cofrua   =    --------------    +    --------------
CFONC                             pi                     pi
CFONC   luminance de bord   emission propre         flux reflechi.
CFONC
CFONC      (eps=1 : paroi noire ; eps=0 : paroi reflechissante )
CFONC
CFONC  Pour une paroi , la luminance (i.e l'energie rayonnee par la
CFONC  paroi ) comprend son emission propre  et l'energie reflechie.
CFONC   (CF LE SOUS-PROGRAMME usray2)
CFONC
CFONC
CFONC  2/ Milieu libre : luminance rentrante nulle
CFONC
CFONC         cofrua   =   0.D0
CFONC
CFONC     (si l'utilisateur a plus d'informations, il peut ameliorer
CFONC      la situation)
CFONC
CFONC
CFONC
CFONC     L'exemple fourni ci-apres represente les conditions "par defaut"
CFONC       et suffit generalement (si l'utilisateur a utilise
CFONC       les types standard de conditions aux limites dans usclim)
CFONC
CFONC
CFONC
CFONC     1.2 Conditions aux limites pour le modele P-1
CFONC     ---------------------------------------------
CFONC
CFONC
CFONC
CFONC
CFONC   2. DEUXIEME APPEL (IAPPEL = 2)
CFONC   =============================
CFONC
CFONC       La densite de flux net radiatif doit etre calculee
CFONC         de maniere coherente avec les conditions aux limites
CFONC         de la luminance. La densite de flux net radiatif est
CDONF         le bilan entre le rayonnement qu'une face de bord emet
CFONC         (et non pas le rayonnement qu'elle reflechit) et celui
CFONC         qu'elle absorbe.
CFONC
CFONC       L'exemple fourni est coherent avec l'exemple fourni pour
CFONC         conditions aux limites sur la luminance au premier appel
CFONC         et suffi donc en general.
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
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! IZFRDP(NFABOR! TE !  -> ! NUMERO DE ZONE POUR LES FACES DE BORD!
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
CARGU ! COFRUA,COFRUB! TR ! <-  ! CONDITIONS AUX LIMITES AUX           !
CARGU !(NFABOR)      !    !     !    FACES DE BORD POUR LA LUMINANCES  !
CARGU ! W1...6(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! TPAROI(NFABOR! TR !  -> ! TEMPERATURE DE PAROI EN KELVIN       !
CARGU ! QINCID(NFABOR! TR !  -> ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU ! FLUNET(NFABOR! TR ! <-  ! DENSITE DE FLUX NET RADIATIF         !
CARGU ! CK (NCELET)  ! TR ! <-  ! COEFFICIENT D'ABSORPTION DU MILIEU   !
CARGU !              !    !     ! (NUL SI TRANSPARENT)                 !
CARGU ! XLAM(NFABOR) ! TR !  -> ! COEFFICIENT DE CONDUCTIVITE THERMIQUE!
CARGU !              !    !     ! DES FACETTES DE PAROI (W/m/K)        !
CARGU ! EPA (NFABOR) ! TR !  -> ! EPAISSEUR DES FACETTES DE PAROI (m)  !
CARGU ! EPS (NFABOR) ! TR !  -> ! EMISSIVITE DES FACETTES DE BORD      !
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
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
      INCLUDE "ihmpre.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , IPHAS  , IAPPEL
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML) , ITYPFB(NFABOR)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR),IZFRDP(NFABOR)
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
C
      DOUBLE PRECISION COFRUA(NFABOR), COFRUB(NFABOR)
C
      DOUBLE PRECISION W1(NCELET), W2(NCELET), W3(NCELET)
      DOUBLE PRECISION W4(NCELET), W5(NCELET), W6(NCELET)
C
      DOUBLE PRECISION TPAROI(NFABOR), QINCID(NFABOR)
      DOUBLE PRECISION XLAM(NFABOR), EPA(NFABOR)
      DOUBLE PRECISION EPS(NFABOR), FLUNET(NFABOR)
      DOUBLE PRECISION CK(NCELET)
C
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA , IEL, IFAC, IOK
      DOUBLE PRECISION UNSPI, VV, SF, XLC, XIT, XKMIN, DISTBF, PP
C
C***********************************************************************
C
C=======================================================================
C 0 - GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C Indicateur d'arret (pour savoir si des faces ont ete oubliees)
      IOK = 0
C
      UNSPI = 1.D0/PI
C
C=======================================================================
C  1. PREMIER APPEL
C  ================
C=======================================================================
C
      IF (IAPPEL.EQ.1) THEN
C
C
C=======================================================================
C  1.1 - CONDITIONS AUX LIMITES :
C        MODELE DOM : COFRUA CONTIENT LA LUMINANCE
C        MODELE P-1 : COFRUA ET COFRUB SONT A REMPLIR
C      LES EXEMPLES DONNES ICI SONT LES INITIALISATIONS FAITES
C      PAR DEFAUT ET SONT SUFFISANTES EN GENERAL
C=======================================================================
C
C
C
C
C      A - MODELE DOM
C      ^^^^^^^^^^^^^^
C
C
C
C
        IF (IRAYON(IPHAS).EQ.1) THEN
C
          DO IFAC = 1,NFABOR
C
C      1.1.1 - SYMETRIE :
C              ----------
C          REFLEXION TOTALE DU RAYONNEMENT ( EPS=0 )
C          -----------------------------------------
C
            IF (ITYPFB(IFAC).EQ.ISYMET) THEN
C
              COFRUA(IFAC) = QINCID(IFAC) * UNSPI
C
C
C      1.1.2 - PAROIS 'FLUIDES' : LUMINANCES RENTRANTES "NULLES"
C              (ATTENTION LOGIQUE DIFFERENTE DU MODELE P-1)
C          -------------------------------------------------
C
            ELSE IF (ITYPFB(IFAC).EQ.IENTRE
     &          .OR. ITYPFB(IFAC).EQ.ISOLIB) THEN
C
              COFRUA(IFAC) = EPZERO
C
C
C      1.1.3. - PAROIS 'SOLIDES' DE TEMPERATURE TPAROI ET D'EMISSIVITE EPS
C               ----------------------------------------------------------
C
            ELSE IF (ITYPFB(IFAC).EQ.IPAROI
     &          .OR. ITYPFB(IFAC).EQ.IPARUG) THEN
C
              COFRUA(IFAC)  = EPS(IFAC)*STEPHN*(TPAROI(IFAC)**4)*UNSPI
     &                          + (1.D0-EPS(IFAC))* QINCID(IFAC)*UNSPI
C
            ELSE
C
C      1.1.4 - SI DES FACES N'ONT PAS ETE TRAITEES, IL FAUT S'ARRETER
C              ------------------------------------------------------
C
C           ==============================================
C
C             CONSERVER IMPERATIVEMENT LE TEST D'ARRET
C
C           ==============================================
C
              WRITE (NFECRA,1000) IFAC,IZFRDP(IFAC),ITYPFB(IFAC)
              IOK = IOK + 1
            ENDIF
C
          ENDDO
C
C
C
C
C   B - MODELE P-1
C   ^^^^^^^^^^^^^^
C
C
C
C
C
        ELSE IF (IRAYON(IPHAS).EQ.2) THEN
C
          DO IFAC = 1,NFABOR
C
C      1.1.1 - SYMETRIE ET PAROI REFLECHISSANTE (EPS = 0) :
C              CONDITION DE FLUX NUL
C              ------------------------------------------
C
            IF (ITYPFB(IFAC).EQ.ISYMET     .OR.
     &         ((ITYPFB(IFAC).EQ.IPAROI.OR.
     &           ITYPFB(IFAC).EQ.IPARUG).AND.EPS(IFAC).EQ.0D0)) THEN
C
              COFRUA(IFAC) = 0.D0
              COFRUB(IFAC) = 1.D0
C
C
C      1.1.2 - PAROIS 'FLUIDES' : CONDITION DE FLUX NUL
C              (ATTENTION LOGIQUE DIFFERENTE DU MODELE DOM)
C              --------------------------------------------
C
            ELSE IF (ITYPFB(IFAC).EQ.IENTRE
     &          .OR. ITYPFB(IFAC).EQ.ISOLIB) THEN
C
              COFRUA(IFAC) = 0.D0
              COFRUB(IFAC) = 1.D0
C
C
C      1.1.3 - PAROIS 'SOLIDES' DE TEMPERATURE TPAROI ET D'EMISSIVITE EPS
C              (EPS NON NUL)
C              ----------------------------------------------------------
C
            ELSE IF (ITYPFB(IFAC).EQ.IPAROI .OR.
     &               ITYPFB(IFAC).EQ.IPARUG ) THEN
C
              DISTBF = RA(IDISTB-1+IFAC)
C
              XIT = 1.5D0 *DISTBF *CK(IFABOR(IFAC))
     &            * (2.D0 /(2.D0-EPS(IFAC)) -1.D0)
C
              COFRUB(IFAC) = 1.D0 / (1.D0 + XIT)
              COFRUA(IFAC) = XIT * TPAROI(IFAC)**4 * COFRUB(IFAC)
C
            ELSE
C
C      1.1.4 - SI DES FACES N'ONT PAS ETE TRAITEES, IL FAUT S'ARRETER
C              ------------------------------------------------------
C
C           ==============================================
C
C             CONSERVER IMPERATIVEMENT LE TEST D'ARRET
C
C           ==============================================
C
              WRITE (NFECRA,1000) IFAC,IZFRDP(IFAC),ITYPFB(IFAC)
            IOK = IOK + 1
          ENDIF
C
        ENDDO
C
        ENDIF
C
        IF (IOK.NE.0) THEN
          WRITE (NFECRA,1100) IPHAS
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C=======================================================================
C  2 - DEUXIEME APPEL
C  ===================
C=======================================================================
C
      ELSE IF (IAPPEL.EQ.2) THEN
C
C=======================================================================
C  2.1 - DENSITE DE FLUNET RADIATIF AUX DIFFERENTES FRONTIERES
C      L'EXEMPLE DONNE ICI EST L'INITIALISATION FAITE PAR DEFAUT
C=======================================================================
C
C    DANS LE CAS OU LES CONDITIONS A LA LIMITES CI-DESSUS
C      AURAIENT ETE MODIFIEES, IL EST NECESSAIRE DE MODIFIER
C      LA MANIERE DONT EST CALCULE LA DENSITE DE FLUX NET RADIATIF,
C      DE MANIERE COHERENTE.
C    LA REGLE EST LA SUIVANTE :
C      LA DENSITE DE FLUX NET EST UN BILAN ENTRE CE QU'UNE FACE
C      DE BORD EMET COMME RAYONNEMENT (ET NON CE QU'ELLE REFLECHIT)
C      ET CE QU'ELLE ABSORBE  (ORIENTATION DE LA NORMALE SORTANTE)
C      AINSI, SI UNE PAROI CHAUFFE LE FLUIDE, FLUNET < 0
C
C
C
C
        DO IFAC = 1,NFABOR
C
          IF (ITYPFB(IFAC).EQ.IPAROI .OR.
     &        ITYPFB(IFAC).EQ.IPARUG) THEN
C
C      2.1.1 - PAROIS 'SOLIDES' DE TEMPERATURE TPAROI ET D'EMISSIVITE EPS
C              ----------------------------------------------------------
C
            FLUNET(IFAC) =
     &      EPS(IFAC) *(QINCID(IFAC) - STEPHN*TPAROI(IFAC)**4)
C
C
C      2.1.2 - SYMETRIE :
C              ----------
C          REFLEXION TOTALE DU RAYONNEMENT ( FLUNET = 0 )
C          ----------------------------------------------
          ELSE IF (ITYPFB(IFAC).EQ.ISYMET) THEN
C
            FLUNET(IFAC)= ZERO
C
C
C      2.1.3 - PAROIS 'FLUIDES'
C              ----------------
C
          ELSE IF (ITYPFB(IFAC).EQ.IENTRE
     &        .OR. ITYPFB(IFAC).EQ.ISOLIB) THEN
C
            IF (IRAYON(IPHAS).EQ.1) THEN
C
            FLUNET(IFAC)= QINCID(IFAC) -PI*COFRUA(IFAC)
C
            ELSE IF (IRAYON(IPHAS).EQ.2) THEN
C
              FLUNET(IFAC)= 0.D0
C
            ENDIF
C
C
C      2.1.4 - SI DES FACES N'ONT PAS ETE TRAITEES, IL FAUT S'ARRETER
C              ------------------------------------------------------
          ELSE
C
C           ==============================================
C
C             CONSERVER IMPERATIVEMENT LE TEST D'ARRET
C
C           ==============================================
C
            WRITE (NFECRA,2000) IFAC,IZFRDP(IFAC),ITYPFB(IFAC)
            IOK = IOK + 1
C
          ENDIF
C
        ENDDO
C
C
        IF (IOK.NE.0) THEN
          WRITE (NFECRA,2100) IPHAS
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C
      ENDIF
C
C -------
C FORMATS
C -------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@                CONDITIONS AUX LIMITES NON RENSEIGNEES      ',/,
     &'@                                                            ',/,
     &'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )
C
 1100 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@    LES CONDITIONS AUX LIMITES NE SONT PAS RENSEIGNEES POUR ',/,
     &'@     CERTAINES FACES DE BORD (Phase ',I10   ,')             ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier le codage de usray3.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT (FLUNET    NON RENSEIGNE)       ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )
C
 2100 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@    LE FLUNET    N''EST PAS RENSEIGNEE POUR CERTAINES       ',/,
     &'@        FACES DE BORD (Phase ',I10   ,')                    ',/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier le codage de usray3.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C
      END
