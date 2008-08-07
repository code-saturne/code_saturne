c@a
c@versb
C-----------------------------------------------------------------------
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
                        SUBROUTINE RAYDOM
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IZFRAD ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
C
     &   COFRUA , COFRUB ,
     &   FLURDS , FLURDB ,
C
     &   DTR    , VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , W10    ,
     &   RDEVEL , RTUSER ,
C
     &   RAYEXP , RAYIMP , QX     , QY     , QZ     ,
     &   RAYABS , RAYEMI , CK     , TEMPK  ,
     &   TPAROI , QINCID , XLAM   , EPA    , EPS    ,
     &   FLUNET , FLCONV , HFCONV ,
C
     &   RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE RAYONNEMENT :
CFONC    --------------------------------------
CFONC
CFONC   Enveloppe principale du module de résolution de l'équation
CFONC   des transferts radiatifs
CFONC
CFONC   Deux methodes sont disponibles :
CFONC
CFONC     1) La methode : "Discretes Ordinates Methods" (DOM)
CFONC     2) L'approximation P-1 (recommandé uniquement pour le CP)
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
CARGU ! IZFRAD(NFABOR! TE !  -> ! NUMERO DE ZONE DES FACES DE BORD     !
CARGU !   ,NPHAST)   !    !     !                                      !
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
CARGU ! COFRUA,COFRUB! TR !  -  ! CONDITIONS AUX LIMITES AUX           !
CARGU !(NFABOR)      !    !     !    FACES DE BORD POUR LA LUMINANCES  !
CARGU ! FLURDS,FLURDB! TR !  -  ! PSEUDO FLUX DE MASSE (FACES INTERNES !
CARGU !(NFAC)(NFABOR)!    !     !    ET FACES DE BORD )                !
CARGU ! DTR(NCELET)  ! TR !  -  ! DT*CDTVAR                            !
CARGU ! VISCF(NFAC)  ! TR !  -  ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -  ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! DAM(NCELET   ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! XAM(NFAC,*)  ! TR !  -  ! TABLEAU DE TRAVAIL POUR MATRICE      !
CARGU ! DRTP(NCELET  ! TR !  -  ! TABLEAU DE TRAVAIL POUR INCREMENT    !
CARGU ! SMBRS(NCELET ! TR !  -  ! TABLEAU DE TRAVAIL POUR SEC MEM      !
CARGU ! ROVSDT(NCELET! TR !  -  ! TABLEAU DE TRAVAIL POUR TERME INSTAT !
CARGU ! W1...9(NCELET! TR !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! RAYEXP(NCELET! TR ! <-  ! TERME SOURCE RADIATIF EXPLICITE      !
CARGU !   ,NPHASC)   !    !     !                                      !
CARGU ! RAYIMP(NCELET! TR ! <-  ! TERME SOURCE RADIATIF IMPLICITE      !
CARGU !   ,NPHASC)   !    !     !                                      !
CARGU ! QXQYQZ(NCELET! TR ! <-  ! COMPOSANTE DU VECTEUR DENSITE DE FLUX!
CARGU !   ,NPHAST)   !    !     ! RADIATIF EXPLICITE                   !
CARGU ! RAYABS(NCELET! TR ! <-  ! PART D'ABSORPTION DU TERME SOURCE    !
CARGU !   ,NPHASC)   !    !     ! RADIATIF                             !
CARGU ! RAYEMI(NCELET! TR ! <-  ! PART D'EMISSION DU TERME SOURCE      !
CARGU !   ,NPHASC)   !    !     ! RADIATIF EXPLICITE                   !
CARGU ! CK (NCELET   ! TR ! <-  ! COEFFICIENT D'ABSORPTION DU MILIEU   !
CARGU !   ,NPHASC)   !    !     ! (NUL SI TRANSPARENT)                 !
CARGU ! TEMPK(NCELET)! TR ! <-  ! TEMPERATURE EN KELVIN                !
CARGU !   ,NPHASC)   !    !     !                                      !
CARGU ! TPAROI(NFABOR! TR !  -  ! TEMPERATURE DE PAROI EN KELVIN       !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! QINCID(NFABOR! TR ! <-  ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! XLAM(NFABOR  ! TR ! <-  ! COEFFICIENT DE CONDUCTIVITE THERMIQUE!
CARGU !   ,NPHAST)   !    !     ! DES FACETTES DE PAROI (W/m/K)        !
CARGU ! EPA (NFABOR  ! TR ! <-  ! EPAISSEUR DES FACETTES DE PAROI (m)  !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! EPS (NFABOR  ! TR ! <-  ! EMISSIVITE DES FACETTES DE BORD      !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! FLUNET(NFABOR! TR ! <-  ! DENSITE DE FLUX NET RADIATIF AUX     !
CARGU !   ,NPHAST)   !    !     ! FACES DE BORD                        !
CARGU ! FLCONV(NFABOR! TR ! <-  ! DENSITE DE FLUX CONVECTIF AUX FACES  !
CARGU !   ,NPHAST)   !    !     ! DE BORD                              !
CARGU ! HFCONV(NFABOR! TR ! <-  ! COEFFICIENT D'ECHANGE FLUIDE AUX     !
CARGU !   ,NPHAST)   !    !     ! FACES DE BORD                        !
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
      INCLUDE "lagpar.h"
      INCLUDE "lagdim.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
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
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML) , ITYPFB(NFABOR,NPHAS)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          IZFRAD(NFABOR,NPHAST)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE)
      INTEGER          IA(*)
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
      DOUBLE PRECISION FLURDS(NFAC), FLURDB(NFABOR)
C
      DOUBLE PRECISION DTR(NCELET)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION DAM(NCELET), XAM(NFAC,2)
      DOUBLE PRECISION DRTP(NCELET), SMBRS(NCELET)
      DOUBLE PRECISION ROVSDT(NCELET)
      DOUBLE PRECISION W1(NCELET) , W2(NCELET) , W3(NCELET)
      DOUBLE PRECISION W4(NCELET) , W5(NCELET) , W6(NCELET)
      DOUBLE PRECISION W7(NCELET) , W8(NCELET) , W9(NCELET)
      DOUBLE PRECISION W10(NCELET)
C
      DOUBLE PRECISION RAYEXP(NCELET,NPHASC)
      DOUBLE PRECISION RAYIMP(NCELET,NPHASC)
      DOUBLE PRECISION CK(NCELET,NPHASC)
      DOUBLE PRECISION TEMPK(NCELET,NPHASC)
      DOUBLE PRECISION RAYABS(NCELET,NPHASC)
      DOUBLE PRECISION RAYEMI(NCELET,NPHASC)
C
      DOUBLE PRECISION QX(NCELET,NPHAST), QY(NCELET,NPHAST)
      DOUBLE PRECISION QZ(NCELET,NPHAST)
      DOUBLE PRECISION TPAROI(NFABOR,NPHAST), QINCID(NFABOR,NPHAST)
      DOUBLE PRECISION XLAM(NFABOR,NPHAST), EPA(NFABOR,NPHAST)
      DOUBLE PRECISION EPS(NFABOR,NPHAST), FLUNET(NFABOR,NPHAST)
      DOUBLE PRECISION FLCONV(NFABOR,NPHAST), HFCONV(NFABOR,NPHAST)
C
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IPH    , IPHAS  , IAPPEL
      INTEGER          IFAC   , IEL    , IOK    , IZONE
      INTEGER          INC    , ICCOCG , IWARNP , IMLIGP , NSWRGP
      INTEGER          MODE   , ICLA   , IPCLA  , IVAR0
      INTEGER          ISCAT  , IVART  , IPHYDP
      INTEGER          IDIMTE , ITENSO
      INTEGER          IFLUX(NOZRDM)
      DOUBLE PRECISION EPSRGP, CLIMGP, EXTRAP
      DOUBLE PRECISION SURFBN
      DOUBLE PRECISION AA, BB, CKMIN, UNSPI, XLIMIT, COFRMN, FLUNMN
      DOUBLE PRECISION FLUX(NOZRDM)
      DOUBLE PRECISION VV, SF, XLC, XKMIN, PP
C
      INTEGER    IPADOM
      DATA       IPADOM /0/
      SAVE       IPADOM
C
C**********************************************************************
C=======================================================================
C 0. GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 1. INITIALISATIONS GENERALES
C=======================================================================
C
C---> NUMERO DE PASSAGE RELATIF
C
      IPADOM = IPADOM + 1
      IF (IPADOM.GT.1 .AND. MOD(NTCABS,NFREQR).NE.0) RETURN
C
      WRITE(NFECRA,1000)
C
C---> INITIALISATION DES CONSTANTES
C
      UNSPI = 1.D0/PI
C
C=======================================================================
C 2. BOUCLE SUR LES PHASES...
C=======================================================================
C
      DO IPH = 1, NPHAST
C
        IPHAS = IRAPHA(IPH)
C
C---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE
C
        ISCAT = ISCALT(IPHAS)
        IVART = ISCA(ISCALT(IPHAS))
C
C=======================================================================
C 3.1 COEFFICIENT D'ABSORPTION DU MILIEU SEMI-TRANSPARENT
C=======================================================================
C
C--> INITIALISATION NON ADMISSIBLE POUR TEST APRES USRAY3
C
        DO IEL = 1,NCEL
          CK(IEL,IPH) = -GRAND
        ENDDO
C
C--> COEFFICIENT D'ABSORPTION POUR LA PHYSIQUE PARTICULIERE
C
C  ATTENTION  : DANS LE CAS DE L'APPROXIMATION P-1, LE COEFFICIENT
C    D'ABSORPTION EST UTILISE POUR CALCULER LES CONDITIONS AUX
C    LIMITES DE L'EQUATION A RESOUDRE.
C
        IF (NSCAPP.GT.0) THEN
C
          CALL PPCABS
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     ,
     &   RDEVEL , RTUSER ,
     &   CK     ,  RA    )
C
C-----> W10 SERT A STOCKER TEMPORAIREMENT
C       LE COEFFICIENT D'ABSORPTION DU MELANGE GAZ-PARTICULE
C
          IF (IPPMOD(ICP3PL).GE.0 .OR. IPPMOD(ICFUEL).GE.0) THEN
C
            DO IEL = 1,NCEL
              W10(IEL) = CK(IEL,IPH)
            ENDDO
C
            IF (IPPMOD(ICP3PL).GE.0 ) THEN
              DO ICLA = 1,NCLACP
                IPCLA = 1+ICLA
                DO IEL = 1,NCEL
                  W10(IEL) = W10(IEL)
     &                   + ( PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                   *   CK(IEL,IPCLA)                 )
                ENDDO
              ENDDO
            ELSE IF ( IPPMOD(ICFUEL) .GE.0 ) THEN
              DO ICLA = 1,NCLAFU
                IPCLA = 1+ICLA
                DO IEL = 1,NCEL
                  W10(IEL) = W10(IEL)
     &                     + RTPA(IEL,ISCA(IYFOL(ICLA)))*CK(IEL,IPCLA)
                ENDDO
              ENDDO
            ENDIF
C
            DO IEL = 1,NCEL
              CK(IEL,IPH) = W10(IEL)
            ENDDO
          ENDIF
C
        ELSE
C
C
C---> LECTURES DES DONNEES UTILISATEURS
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
C
          IF (IIHMPR.EQ.1) THEN
C
            CALL UIRAY3 (CK, IPH, NCELET, NCEL)
C           ===========
C
            IF (IRAYON(IPHAS).EQ.2 .AND. NSCAPP.LE.0
     &          .AND. IPADOM.LE.3) THEN
              SF = 0.D0
              VV = 0.D0
C
C         Calcul de la longueur caractéristique du domaine de calcul
              DO IFAC = 1,NFABOR
                SF = SF + SQRT(SURFBO(1,IFAC)**2 +
     &                         SURFBO(2,IFAC)**2 +
     &                         SURFBO(3,IFAC)**2 )
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARSOM(SF)
C               ===========
              ENDIF
C
              DO IEL = 1,NCEL
                VV = VV + VOLUME(IEL)
              ENDDO
              IF (IRANGP.GE.0) THEN
                CALL PARSOM(VV)
C               ===========
              ENDIF
C
              XLC = 3.6D0 * VV / SF
C
C             Clipping pour la variable CK
C
              XKMIN = 1.D0 / XLC
C
              IOK = 0
              DO IEL = 1,NCEL
                IF (CK(IEL,IPH).LT.XKMIN) THEN
                  IOK = IOK +1
                ENDIF
              ENDDO
C
C     Alerte si epaisseur optique trop grande

              PP = XNP1MX/100.0D0
              IF (DBLE(IOK).GT.PP*DBLE(NCEL)) THEN
                WRITE(NFECRA,6000) XKMIN, DBLE(IOK)/DBLE(NCEL)*100.D0,
     &                             XNP1MX
              ENDIF
            ENDIF

          ENDIF
#endif
C
          CALL USRAY3
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPH    , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IZFRAD(1,IPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   CK(1,IPH), W1   , W2     , W3     , W4     , W5     ,  W6    ,
     &   RDEVEL , RTUSER ,
     &   RA     )
C
        ENDIF
C
C--> VERIFICATIONS GENERALES :
C
C---> P-1 : VERIFICATION QUE CK DU MILIEU EST STRICTEMENT
C       SUPERIEUR A ZERO POUR TOUTES CELLULES
C
        IF (IRAYON(IPHAS).EQ.2) THEN
C
          CKMIN = CK(1,IPH)
          DO IEL = 1, NCEL
            CKMIN = MIN(CKMIN,CK(IEL,IPH))
          ENDDO
          IF (CKMIN.LT.0.D0) THEN
            WRITE(NFECRA,2020)
            CALL CSEXIT (1)
C           ===========
          ENDIF
C
        ELSE IF (IRAYON(IPHAS).EQ.1) THEN
C
C---> DOM : VERIFICATION QUE CK DU MILIEU EST SUPERIEUR OU EGAL A -1D-12
C
          CKMIN = CK(1,IPH)
          DO IEL = 1, NCEL
            CKMIN = MIN(CKMIN,CK(IEL,IPH))
          ENDDO
          IF (CKMIN.LE.-EPZERO) THEN
            WRITE(NFECRA,2010) CKMIN
            CALL CSEXIT (1)
C           ===========
          ENDIF
C
        ENDIF
C
C---> VERIFICATION D'UN CAS TRANSPARENT
C
        AA = ZERO
        DO IEL = 1,NCEL
          AA = AA + CK(IEL,IPH)
        ENDDO
        IF (IRANGP.GE.0) THEN
          CALL PARMAX(AA)
C         ===========
        ENDIF
        IF (AA.LE.EPZERO) THEN
          WRITE(NFECRA,1100)
          IDIVER = -1
        ENDIF
C
C=======================================================================
C 3.2 CONDITIONS AUX LIMITES POUR LES EQNS DE LA DOM ET DE L'APPROX P-1
C
C     REMPLISSAGE DES CONDITIONS AUX LIMITES POUR LA LUMINANCE
C     (TABLEAUX  COFRUA & COFRUB)
C=======================================================================
C
C-----> INITIALISATIONS NON ADMISSIBLES POUR TEST APRES USRAY5
C
        DO IFAC = 1,NFABOR
          COFRUA(IFAC) = -GRAND
          COFRUB(IFAC) = -GRAND
        ENDDO
C
        IAPPEL = 1
C
        CALL USRAY5
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IZFRAD(1,IPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   COFRUA , COFRUB ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
     &   TPAROI(1,IPH) , QINCID(1,IPH) , FLUNET(1,IPH) ,
     &   XLAM(1,IPH)   , EPA(1,IPH)    , EPS(1,IPH)    , CK(1,IPH)    ,
     &   RA     )
C
C
C---> BLINDAGE POUR UN DIRICHLET SUR LA LUMINANCE
C      (UNIQUEMENT POUR LA METHODE DOM, SANS OBJET POUR L'APPROX P-1)
C
        IF (IRAYON(IPHAS).EQ.1) THEN
          DO IFAC = 1,NFABOR
            COFRUB(IFAC) = ZERO
          ENDDO
        ENDIF
C
C---> VERIFICATION DU REMPLISSAGE DE COFRUA & COFRUB
C
C     Attention : dans le cas de l'approx P-1 la valeur de COFRUA peut
C     etre grande (de l'ordre de TPAROI**4), d'ou la valeur de COFRMN
C
        IOK = 0
        XLIMIT = -GRAND*0.1D0
C GRAND n'est pas assez grand...
        COFRMN = RINFIN
C
        DO IFAC = 1,NFABOR
          IF (COFRUA(IFAC).LE.XLIMIT) THEN
            IOK = IOK + 1
            COFRMN = MIN(COFRMN,COFRUA(IFAC))
            WRITE(NFECRA,3000)IFAC,IZFRAD(IFAC,IPH),ITYPFB(IFAC,IPHAS)
          ENDIF
        ENDDO
C
        IF (IOK.NE.0) THEN
          WRITE(NFECRA,3100) IPHAS, COFRMN
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
        COFRMN = RINFIN
C
        IF (IRAYON(IPHAS).EQ.2) THEN
C
          DO IFAC = 1,NFABOR
            IF (COFRUB(IFAC).LE.XLIMIT) THEN
              IOK = IOK + 1
              COFRMN = MIN(COFRMN,COFRUB(IFAC))
              WRITE(NFECRA,3000)IFAC,IZFRAD(IFAC,IPH),ITYPFB(IFAC,IPHAS)
            ENDIF
          ENDDO
C
          IF (IOK.NE.0) THEN
            WRITE(NFECRA,3200) IPHAS,COFRMN
            CALL CSEXIT (1)
C           ===========
          ENDIF
C
        ENDIF
C
C=======================================================================
C 4. STOCKAGE DE LA TEMPERATURE (en Kelvin) dans TEMPK(IEL,IPH)
C=======================================================================
C
        IF (IDIVER.GE.0) THEN
C
          IF(ABS(ISCSTH(ISCAT)).EQ.1) THEN
C
C---> TRANSPORT DE LA TEMPERATURE
C
            IF (ISCSTH(ISCAT).EQ.-1) THEN
              DO IEL = 1, NCEL
                TEMPK(IEL,IPH) = RTPA(IEL,IVART) + TKELVI
              ENDDO
            ELSE
              DO IEL = 1, NCEL
                TEMPK(IEL,IPH) = RTPA(IEL,IVART)
              ENDDO
            ENDIF
C
C---> TRANSPORT DE L'ENTHALPIE (FLURDB est un auxiliaire)
C
          ELSE IF (ISCSTH(ISCAT).EQ.2) THEN
C
            MODE = 1
C
            IF (NSCAPP.EQ.0) THEN
C
              CALL USRAY4
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   MODE   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
C
     &   TPAROI(1,IPH)   , FLURDB , TEMPK(1,IPH)   ,
C
     &   RA     )
C
            ELSE
C
              CALL PPRAY4
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   MODE   ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
     &   TPAROI(1,IPH)   , FLURDB , TEMPK(1,IPH)   ,
     &   RA     )
C
            ENDIF
C
            IF ( IPPMOD(ICP3PL).GE.0 ) THEN
C
C Temperature des particules
C
              DO ICLA = 1, NCLACP
                IPCLA = 1+ICLA
                DO IEL = 1,NCEL
                  TEMPK(IEL,IPCLA) = PROPCE(IEL,IPPROC(ITEMP2(ICLA)))
                ENDDO
              ENDDO
C
C         Fuel
C
            ELSE IF ( IPPMOD(ICFUEL).GE.0 ) THEN
C
              DO ICLA = 1, NCLAFU
                IPCLA = 1+ICLA
                DO IEL = 1,NCEL
                  TEMPK(IEL,IPCLA) = PROPCE(IEL,IPPROC(ITEMP3(ICLA)))
                ENDDO
              ENDDO
C
            ENDIF
C
          ELSE
            WRITE(NFECRA,3500)ISCAT,ISCSTH(ISCAT)
            CALL CSEXIT (1)
C           ===========
          ENDIF
C
C---> ON SE SERT DE RAYIMP(IEL,IPH) COMME UN AUXILIAIRE POUR
C       STOCKER STEPHN*CK*TEMPK**4 ICI
C       PLUS BAS ON JUSTIFIERA LE NOM.
C
          IF ( IPPMOD(ICOD3P).EQ.-1 .AND.
     &         IPPMOD(ICOEBU).EQ.-1       ) THEN
C
C           Rayonnement standard, flamme CP ou fuel
C
            DO IEL = 1,NCEL
              RAYIMP(IEL,IPH) = STEPHN*CK(IEL,IPH)*(TEMPK(IEL,IPH)**4)
            ENDDO
C
          ELSE
C
C           Flamme de diffusion ou flamme de premelange
C
            DO IEL = 1,NCEL
              RAYIMP(IEL,IPH) = STEPHN*CK(IEL,IPH)
     &                             *PROPCE(IEL,IPPROC(IT4M))
            ENDDO
C
          ENDIF
C
C     Charbon
C
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
            DO ICLA = 1,NCLACP
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                RAYIMP(IEL,IPCLA) =
     &                  STEPHN*CK(IEL,IPCLA)*(TEMPK(IEL,IPCLA)**4)
              ENDDO
            ENDDO
C
C       Fuel
C
          ELSE IF ( IPPMOD(ICFUEL).GE.0 ) THEN
            DO ICLA = 1,NCLAFU
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                RAYIMP(IEL,IPCLA) =
     &                  STEPHN*CK(IEL,IPCLA)*(TEMPK(IEL,IPCLA)**4)
              ENDDO
            ENDDO
          ENDIF
C
        ELSE
          DO IEL = 1,NCEL
            RAYIMP(IEL,IPH) = ZERO
          ENDDO
C fin de IF (IDIVER.GE.0) THEN
        ENDIF
C
C=======================================================================
C 5.1 MODELE DE RAYONNEMENT P-1
C=======================================================================
C
        IF (IRAYON(IPHAS).EQ.2) THEN
C
C--> Terme source explicite de l'equation sur Theta4
C
          DO IEL = 1, NCEL
            SMBRS(IEL) = 3.D0 * CK(IEL,IPH) * ( TEMPK(IEL,IPH) ** 4)
     &                   * VOLUME(IEL)
          ENDDO
C
C Tenir compte de l'absorption des particules
C
C       Charbon
C
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
            DO ICLA = 1,NCLACP
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                SMBRS(IEL) = SMBRS(IEL)
     &                         + ( 3.D0 * PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                         *   CK(IEL,IPCLA)*(TEMPK(IEL,IPCLA)**4)
     &                         *   VOLUME(IEL) )
              ENDDO
            ENDDO
C
C       FUEL
C
          ELSE IF ( IPPMOD(ICFUEL).GE.0 ) THEN
            DO ICLA = 1,NCLAFU
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                SMBRS(IEL) = SMBRS(IEL)
     &                      +( 3.D0 * RTPA(IEL,ISCA(IYFOL(ICLA)))
     &                        *CK(IEL,IPCLA)*(TEMPK(IEL,IPCLA)**4)
     &                        *VOLUME(IEL) )
              ENDDO
            ENDDO
          ENDIF
C
C--> Terme source implicite de l'equation sur Theta4
C
          DO IEL = 1, NCEL
            ROVSDT(IEL) =  3.D0 * CK(IEL,IPH) * VOLUME(IEL)
          ENDDO
C
C Tenir compte de l'absorption des particules
C
C        Charbon
C
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
            DO ICLA = 1,NCLACP
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                ROVSDT(IEL) = ROVSDT(IEL)
     &                          + (3.D0 * PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                          *  CK(IEL,IPCLA) * VOLUME(IEL) )
              ENDDO
            ENDDO
C
C        Fuel
C
          ELSE IF ( IPPMOD(ICFUEL).GE.0 ) THEN
            DO ICLA = 1,NCLAFU
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                ROVSDT(IEL) = ROVSDT(IEL)
     &                       + (3.D0*RTPA(IEL,ISCA(IYFOL(ICLA)))
     &                              *CK(IEL,IPCLA)*VOLUME(IEL) )
              ENDDO
            ENDDO
C
          ENDIF
C
C--> Inverse du coefficient de diffusion de l'equation sur Theta4
C       A priori W10 contient deja la bonne info, mais pour plus de
C       securite  on le re-remplit
C
          DO IEL = 1, NCEL
            W10(IEL) =  CK(IEL,IPH)
          ENDDO
C
C Tenir compte de l'absorption des particules
C
C        Charbon
C
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
            DO ICLA = 1,NCLACP
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                W10(IEL) = W10(IEL)
     &                    + ( PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                       *CK(IEL,IPCLA)                )
              ENDDO
            ENDDO

C        Fuel
C
          ELSE IF ( IPPMOD(ICFUEL).GE.0 ) THEN
            DO ICLA = 1,NCLAFU
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                W10(IEL) = W10(IEL)
     &                    +RTPA(IEL,ISCA(IYFOL(ICLA)))*CK(IEL,IPCLA)
              ENDDO
            ENDDO
C
          ENDIF
C
          CALL RAYPUN
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   COFRUA , COFRUB ,
     &   FLURDS , FLURDB ,
     &   DTR    , VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
C
     &   RAYABS(1,IPH) , RAYEMI(1,IPH) , RAYEXP(1,IPH) ,
     &   QX(1,IPH)     , QY(1,IPH)     , QZ(1,IPH)     ,
     &   QINCID(1,IPH) , EPS(1,IPH)    , TPAROI(1,IPH) ,
C
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , W10    ,
     &   RDEVEL , RTUSER , RA     )
C
C=======================================================================
C 5.2 RESOLUTION DE L'EQUATION DES TRANSFERTS RADIATIFS
C=======================================================================
C
        ELSE IF (IRAYON(IPHAS).EQ.1) THEN
C
C--> Terme source explicite de l'equation sur la luminance
C
C
          DO IEL = 1, NCEL
            SMBRS(IEL) = RAYIMP(IEL,IPH) *VOLUME(IEL) *UNSPI
          ENDDO

C       Charbon
C
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
            DO ICLA = 1,NCLACP
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                SMBRS(IEL) = SMBRS(IEL)
     &                      + PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                       *RAYIMP(IEL,IPCLA)*VOLUME(IEL)*UNSPI
              ENDDO
            ENDDO
C
C       Fuel
C
          ELSEIF ( IPPMOD(ICFUEL).GE.0 ) THEN
            DO ICLA = 1,NCLAFU
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                SMBRS(IEL) = SMBRS(IEL)
     &                      + RTPA(IEL,ISCA(IYFOL(ICLA)))
     &                       *RAYIMP(IEL,IPCLA)*VOLUME(IEL)*UNSPI
              ENDDO
            ENDDO
C
          ENDIF
C
C--> Terme source implicite de l'equation sur la luminance
C      KL + div(LS) = KL0 integre sur le volume de controle
C
          DO IEL = 1, NCEL
            ROVSDT(IEL) =  CK(IEL,IPH) * VOLUME(IEL)
          ENDDO
C
C        Charbon
C
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
            DO ICLA = 1,NCLACP
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                ROVSDT(IEL) = ROVSDT(IEL) +
     &                          PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                          * CK(IEL,IPCLA) * VOLUME(IEL)
              ENDDO
            ENDDO
C
C        Fuel
C
          ELSEIF ( IPPMOD(ICFUEL).GE.0 ) THEN
            DO ICLA = 1,NCLAFU
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                ROVSDT(IEL) = ROVSDT(IEL)
     &                       + RTPA(IEL,ISCA(IYFOL(ICLA)))
     &                        *CK(IEL,IPCLA)*VOLUME(IEL)
              ENDDO
            ENDDO
C
          ENDIF
C
          CALL RAYSOL
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   COFRUA , COFRUB ,
     &   FLURDS , FLURDB ,
     &   DTR    , VISCF  , VISCB  ,
     &   DAM    , XAM    ,
     &   DRTP   , SMBRS  , ROVSDT ,
     &   W1     , W2     , W3     , W4     , W5     ,
     &   W6     , W7     , W8     , W9     , W10    ,
     &   RDEVEL , RTUSER ,
     &   RAYABS(1,IPH) , RAYEMI(1,IPH) , RAYEXP(1,IPH) ,
     &   QX(1,IPH)     , QY(1,IPH)     , QZ(1,IPH)     ,
     &   QINCID(1,IPH) , FLUNET(1,IPH) ,
C
     &   RA     )
C
C
        ENDIF
C
C=======================================================================
C 5.3 STOCKAGE DE L'INTEGRALE DE LA LUMINANCE POUR LE LANGRANGIEN
C=======================================================================
C
C  Si dans le module lagrangien on resout une equation de la temperature
C    sur les particules (IPHYLA=1 et ITPVAR=1) ou si les particules
C    sont des grains de charbon (IPHYLA=2) alors on stocke
C                                     /    ->  ->
C    l'integrale de la luminance SA= /  L( X , S ). DOMEGA
C                                   /4.PI
C
        IF ( IPPMOD(ICP3PL).GE.0 .OR. IPPMOD(ICFUEL).GE.0 ) THEN
          DO IEL = 1,NCEL
            PROPCE(IEL,IPPROC(ILUMI)) = RAYEXP(IEL,IPH)
          ENDDO
        ELSE IF (IILAGR.GT.0 .AND. IPHYLA.EQ.2 .AND. ITPVAR.EQ.1) THEN
          DO IEL = 1,NCEL
            PROPCE(IEL,IPPROC(ILUMN)) = RAYEXP(IEL,IPH)
          ENDDO
        ENDIF
C
C=======================================================================
C 6. FLUX NET RADIATIF AUX PAROIS : CALCUL ET INTEGRATION
C=======================================================================
C
C
C---> INITIALISATION NON ADMISSIBLE POUR TEST APRES USRAY5
C
        DO IFAC = 1,NFABOR
          FLUNET(IFAC,IPH) = -GRAND
        ENDDO
C
C---> LECTURES DES DONNEES UTILISATEURS
C
        IAPPEL = 2
C
        CALL USRAY5
C       ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  , IAPPEL ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR , IZFRAD(1,IPH) ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   COFRUA , COFRUB ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
     &   TPAROI(1,IPH) , QINCID(1,IPH) , FLUNET (1,IPH)  ,
     &   XLAM(1,IPH)   , EPA(1,IPH)    , EPS(1,IPH)    ,  CK(1,IPH)  ,
     &   RA     )
C
C---> VERIFICATION DE FLUNET
C
        IOK = 0
        XLIMIT = -GRAND*0.1D0
        FLUNMN = GRAND
C
        DO IFAC = 1,NFABOR
          IF (FLUNET(IFAC,IPH).LE.XLIMIT) THEN
            IOK = IOK + 1
            FLUNMN = MIN(FLUNMN,FLUNET(IFAC,IPH))
            WRITE(NFECRA,4000)IFAC,IZFRAD(IFAC,IPH), ITYPFB(IFAC,IPHAS)
          ENDIF
        ENDDO
C
        IF (IOK.NE.0) THEN
          WRITE(NFECRA,4100) IPHAS,FLUNMN
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C--> Intégration du flux net sur les différentes zones de frontieres
C     IFLUX sert en parallele pour reperer les zones existantes
C
        DO IZONE = 1, NOZRDM
          FLUX(IZONE) = 0.D0
          IFLUX(IZONE) = 0
        ENDDO
        DO IFAC = 1,NFABOR
          SURFBN = RA(ISRFBN-1+IFAC)
          IZONE = IZFRAD(IFAC,IPH)
          FLUX(IZONE) = FLUX(IZONE)
     &         + FLUNET(IFAC,IPH)*SURFBN
          IFLUX(IZONE) = 1
        ENDDO
        IF(IRANGP.GE.0) THEN
          CALL PARRSM(NOZARM,FLUX )
          CALL PARIMX(NOZARM,IFLUX)
        ENDIF
C
C
        WRITE(NFECRA,5000)
        WRITE(NFECRA,5010)
        DO IZONE = 1, NOZARM(IPHAS)
          IF(IFLUX(IZONE).EQ.1) THEN
            WRITE(NFECRA,5020) IZONE,FLUX(IZONE)
          ENDIF
        ENDDO
        WRITE(NFECRA,5000)
C
C
C--> Intégration de la densité de flux net aux frontieres
C
        AA = ZERO
        DO IFAC = 1,NFABOR
          SURFBN = RA(ISRFBN-1+IFAC)
          AA =  AA + FLUNET(IFAC,IPH) * SURFBN
        ENDDO
        IF(IRANGP.GE.0) THEN
          CALL PARSOM(AA)
        ENDIF
        WRITE(NFECRA,5030) AA
C
C=======================================================================
C 7. TERMES SOURCES RADIATIFS IMPLICITE ET EXPLICITE
C=======================================================================
C
C
C=======================================================================
C 7.1 TERMES SOURCES RADIATIFS SEMI-ANALYTIQUES
C=======================================================================
C
        IF (IDIVER.GE.0) THEN
C
C--> On stocke dans le tableau de travail W9 le CP
C    Attention : il faut conserver W9 dans la suite de la routine,
C    car son contenu est utilisé plus loin
C
          IF (ICP(IPHAS).GT.0) THEN
            DO IEL = 1,NCEL
              W9(IEL) = 1.D0/PROPCE(IEL,IPPROC(ICP(IPHAS)))
            ENDDO
          ELSE
            DO IEL = 1,NCEL
              W9(IEL) = 1.D0/CP0(IPHAS)
            ENDDO
          ENDIF
C
          DO IEL = 1,NCEL
C
C--> part d'absorption du terme source explicite
C
            RAYABS(IEL,IPH) = CK(IEL,IPH) * RAYEXP(IEL,IPH)
C
C--> part d'émission du terme source explicite
C
            RAYEMI(IEL,IPH) = -4.D0 * RAYIMP(IEL,IPH)
C
          ENDDO
C
C Combustion CP : On rajoute la contribution des particules
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
            DO ICLA = 1,NCLACP
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
C         Fluide
                RAYABS(IEL,IPH) = RAYABS(IEL,IPH)
     &                           + PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                            *CK(IEL,IPCLA)*RAYEXP(IEL,IPH)
                RAYEMI(IEL,IPH) = RAYEMI(IEL,IPH)
     &                           - 4.0D0*PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                                  * RAYIMP(IEL,IPCLA)
C         Particule
                RAYABS(IEL,IPCLA) = CK(IEL,IPCLA)*RAYEXP(IEL,IPH)
                RAYEMI(IEL,IPCLA) = - 4.0D0 * RAYIMP(IEL,IPCLA)
                RAYEXP(IEL,IPCLA) = RAYABS(IEL,IPCLA)
     &                             +RAYEMI(IEL,IPCLA)
              ENDDO
            ENDDO
C
C Combustion Fuel : On rajoute la contribution des particules
C
          ELSEIF ( IPPMOD(ICFUEL).GE.0 ) THEN
            DO ICLA = 1,NCLAFU
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
C         Fluide
                RAYABS(IEL,IPH) = RAYABS(IEL,IPH)
     &                           + RTPA(IEL,ISCA(IYFOL(ICLA)))
     &                            *CK(IEL,IPCLA)*RAYEXP(IEL,IPH)
                RAYEMI(IEL,IPH) = RAYEMI(IEL,IPH)
     &                         - 4.0D0*RTPA(IEL,ISCA(IYFOL(ICLA)))
     &                                * RAYIMP(IEL,IPCLA)
C         Particule
                RAYABS(IEL,IPCLA) = CK(IEL,IPCLA)*RAYEXP(IEL,IPH)
                RAYEMI(IEL,IPCLA) = -4.D0*RAYIMP(IEL,IPCLA)
                RAYEXP(IEL,IPCLA) = RAYABS(IEL,IPCLA)
     &                             +RAYEMI(IEL,IPCLA)
              ENDDO
            ENDDO
C
          ENDIF
C
C--> Première méthode pour le calcul du terme source explicite :
C    il est calculé comme la somme des termes d'absorption et d'émission
C    (il faudra multiplier ce terme par VOLUME(IEL) dans COVOFI->RAYSCA)
C
          DO IEL = 1,NCEL
            RAYEXP(IEL,IPH) = RAYABS(IEL,IPH) + RAYEMI(IEL,IPH)
          ENDDO
C
C--> Terme source implicite,
C    (il faudra multiplier ce terme par VOLUME(IEL) dans COVOFI->RAYSCA)
C
        IF ( IPPMOD(ICOD3P).EQ.-1 .AND.
     &       IPPMOD(ICOEBU).EQ.-1       ) THEN
C
C         Rayonnement standard, flamme CP ou fuel
C
          DO IEL = 1,NCEL
            RAYIMP(IEL,IPH) =
     &       -16.D0*CK(IEL,IPH) *STEPHN *(TEMPK(IEL,IPH)**3) * W9(IEL)
          ENDDO
C
        ELSE
C
C         Flamme de diffusion ou flamme de premelange
C
          DO IEL = 1,NCEL
            RAYIMP(IEL,IPH) =
     &       -16.D0*STEPHN*CK(IEL,IPH)*PROPCE(IEL,IPPROC(IT3M))
     &                              *W9(IEL)
          ENDDO
C
        ENDIF
C
C
C Combustion CP : On rajoute la contribution des particules
C
          IF ( IPPMOD(ICP3PL).GE.0 ) THEN
            DO ICLA = 1,NCLACP
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                RAYIMP(IEL,IPH) = RAYIMP(IEL,IPH)
     &       -16.D0*CK(IEL,IPCLA)*PROPCE(IEL,IPPROC(IX2(ICLA)))
     &                          *STEPHN*(TEMPK(IEL,IPCLA)**3)
     &                              / CP2CH(ICHCOR(ICLA))
                RAYIMP(IEL,IPCLA) =
     &       -16.D0*CK(IEL,IPCLA) *STEPHN *(TEMPK(IEL,IPCLA)**3)
     &                                / CP2CH(ICHCOR(ICLA))
              ENDDO
            ENDDO
C
C Combustion FUEL : On rajoute la contribution des particules
C
          ELSEIF ( IPPMOD(ICFUEL).GE.0 ) THEN
            DO ICLA = 1,NCLAFU
              IPCLA = 1+ICLA
              DO IEL = 1,NCEL
                RAYIMP(IEL,IPH) = RAYIMP(IEL,IPH)
     &           -16.D0*CK(IEL,IPCLA)*RTPA(IEL,ISCA(IYFOL(ICLA)))
     &                 *STEPHN*(TEMPK(IEL,IPCLA)**3)
     &                 /CP2FOL
                RAYIMP(IEL,IPCLA) =
     &           -16.D0*CK(IEL,IPCLA)*STEPHN*(TEMPK(IEL,IPCLA)**3)
     &                           /CP2FOL
              ENDDO
            ENDDO
          ENDIF
C
        ELSE
          DO IEL = 1,NCEL
            RAYABS(IEL,IPH) = ZERO
            RAYEMI(IEL,IPH) = ZERO
            RAYEXP(IEL,IPH) = ZERO
            RAYIMP(IEL,IPH) = ZERO
          ENDDO
        ENDIF
C
C=======================================================================
C 7.2 TERME SOURCE RADIATIF EXPLICITE CONSERVATIF
C=======================================================================
C
C A partir d'ici COFRUA et COFRUB deviennent les CL pour la divergence
C
        IF (IDIVER.EQ.1 .OR. IDIVER.EQ.2) THEN
C
          DO IFAC = 1,NFABOR
            COFRUB(IFAC) = ZERO
          ENDDO
C
C--> calcul de la divergence
C
C    En periodique et parallele, echange avant calcul du gradient
C
C    Parallele
          IF(IRANGP.GE.0) THEN
            CALL PARCOM (QX(1,IPH))
C           ===========
            CALL PARCOM (QY(1,IPH))
C           ===========
            CALL PARCOM (QZ(1,IPH))
C           ===========
          ENDIF
C
C    Periodique
          IF(IPERIO.EQ.1) THEN
            IDIMTE = 1
            ITENSO = 0
            CALL PERCOM
C           ===========
     &    ( IDIMTE , ITENSO ,
     &      QX(1,IPH) , QX(1,IPH) , QX(1,IPH) ,
     &      QY(1,IPH) , QY(1,IPH) , QY(1,IPH) ,
     &      QZ(1,IPH) , QZ(1,IPH) , QZ(1,IPH) )
          ENDIF
C
C
C    Donnees pour le calcul de la divergence
C
          INC     = 1
          ICCOCG  = 1
          IMLIGP  = -1
          IWARNP  = IIMLUM
          EPSRGP  = 1.D-8
          CLIMGP  = 1.5D0
          EXTRAP  = 0.D0
          NSWRGP  = 100
C
C------->>direction X
C
          DO IFAC = 1,NFABOR
            SURFBN = RA(ISRFBN-1+IFAC)
            COFRUA(IFAC) =
     &      FLUNET(IFAC,IPH) *SURFBO(1,IFAC) /SURFBN
          ENDDO
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
C    sera a revoir pour la periodicite de rotation
          IVAR0 = 0
          IPHYDP = 0
          CALL GRDCEL
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   QX(1,IPH)       , COFRUA , COFRUB ,
     &   W1     , W2     , W3     ,
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          DO IEL = 1,NCEL
            RAYEXP(IEL,IPH) = - W1(IEL)
          ENDDO
C
C------->>direction Y
C
          DO IFAC = 1,NFABOR
            SURFBN = RA(ISRFBN-1+IFAC)
            COFRUA(IFAC) =
     &      FLUNET(IFAC,IPH) *SURFBO(2,IFAC) /SURFBN
          ENDDO
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
C    sera a revoir pour la periodicite de rotation
          IVAR0 = 0
          IPHYDP = 0
          CALL GRDCEL
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   QY(1,IPH)       , COFRUA , COFRUB ,
     &   W1     , W2     , W3     ,
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          DO IEL = 1,NCEL
            RAYEXP(IEL,IPH) = RAYEXP(IEL,IPH) - W2(IEL)
          ENDDO
C
C------->>direction Z
C
          DO IFAC = 1,NFABOR
            SURFBN = RA(ISRFBN-1+IFAC)
            COFRUA(IFAC) =
     &      FLUNET(IFAC,IPH) *SURFBO(3,IFAC) /SURFBN
          ENDDO
C
C  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
C     n'est pas la vitesse ni Rij)
C    sera a revoir pour la periodicite de rotation
          IVAR0 = 0
          IPHYDP = 0
          CALL GRDCEL
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IVAR0  , IMRGRA , INC    , ICCOCG , NSWRGP , IMLIGP , IPHYDP ,
     &   IWARNP , NFECRA , EPSRGP , CLIMGP , EXTRAP ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   W1     , W1     , W1     ,
     &   QZ(1,IPH)       , COFRUA , COFRUB ,
     &   W1     , W2     , W3     ,
     &   W4     , W5     , W6     ,
     &   RDEVEL , RTUSER , RA     )
C
          DO IEL = 1,NCEL
            RAYEXP(IEL,IPH) = RAYEXP(IEL,IPH) - W3(IEL)
          ENDDO
C
C Fin du calcul de la divergence
        ENDIF
C
C
C=======================================================================
C 7.3 TERME SOURCE RADIATIF EXPLICITE SEMI-ANALYTIQUE CORRIGE
C=======================================================================
C
C
        IF (IDIVER.EQ.2) THEN
C
C---> comparaison des termes sources semi-analytique et conservatif
C
          AA = ZERO
          DO IEL = 1,NCEL
            AA = AA + RAYEXP(IEL,IPH) * VOLUME(IEL)
          ENDDO
C
          BB = ZERO
          DO IEL = 1,NCEL
            BB = BB + (RAYABS(IEL,IPH) + RAYEMI(IEL,IPH)) * VOLUME(IEL)
          ENDDO
C
          IF(IRANGP.GE.0) THEN
            CALL PARSOM(AA)
            CALL PARSOM(BB)
          ENDIF
C
          AA = AA/BB
C
C---> correction du terme source semi-analytique par le conservatif
C
          DO IEL = 1,NCEL
            RAYEXP(IEL,IPH) = (RAYABS(IEL,IPH) + RAYEMI(IEL,IPH)) * AA
          ENDDO
C
        ENDIF
C
C=======================================================================
C 7.4 FINALISATION DU TERME SOURCE EXPLICITE
C=======================================================================
C
        IF (IDIVER.GE.0) THEN
C
C--> Integration volumique du terme source explicite
C    Le resultat de cette integration DOIT etre le meme que l'integration
C    surfacique de la densite de flux net radiatif faite plus haut
C    si  IDIVER = 1 ou 2
C
          AA = ZERO
          DO IEL = 1,NCEL
            AA = AA + RAYEXP(IEL,IPH) * VOLUME(IEL)
          ENDDO
          IF(IRANGP.GE.0) THEN
            CALL PARSOM(AA)
          ENDIF
          WRITE(NFECRA,5040) AA
          WRITE(NFECRA,5050)
          WRITE(NFECRA,5000)
C
C--> Correction du terme source explicite si
C    la variable transportee est la temperature
C    (il faudra multiplier ce terme par VOLUME(IEL) dans COVOFI->RAYSCA)
C
          IF (ABS(ISCSTH(ISCALT(IPHAS))).EQ.1) THEN
            DO IEL = 1,NCEL
              RAYEXP(IEL,IPH) = RAYEXP(IEL,IPH) * W9(IEL)
            ENDDO
          ENDIF
C
        ELSE
          WRITE(NFECRA,5000)
        ENDIF
C
C=======================================================================
C 8. FIN DE LA BOUCLE SUR LES PHASES
C=======================================================================
C
      ENDDO
C
C--------
C FORMATS
C--------
C
 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LE TERME SOURCE RADIATIF',/,
     &           3X,'   -----------------------------------------' )
 1100 FORMAT (/, 3X,'   Calcul effectue en rayonnement transparent'  ,/)
C
 2010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LE RAYONNEMENT EST ACTIVE.                              ',/,
     &'@      LA VALEUR MINIMALE DU COEFFICIENT D ABSORPTION A EST  ',/,
     &'@      EGALE A ', E14.5                                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LE RAYONNEMENT EST ACTIVE AVEC LE MODELE P-1.           ',/,
     &'@      LE COEFFICIENT D''ABSORBTION DOIT ETRE STRICTEMENT    ',/,
     &'@      SUPERIEUR A ZERO.                                     ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3000 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@                CONDITIONS AUX LIMITES MAL RENSEIGNEES      ',/,
     &'@                                                            ',/,
     &'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )
 3100 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@    LES COEFFICIENTS DE CONDITIONS AUX LIMITES (COFRUA)     ',/,
     &'@    NE SONT PAS RENSEIGNES POUR CERTAINES                   ',/,
     &'@        FACES DE BORD (Phase ',I10   ,')                    ',/,
     &'@                                                            ',/,
     &'@        Valeur minimale COFRUA ',E14.5                       ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier le codage de usray5.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3200 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@    LES COEFFICIENTS DE CONDITIONS AUX LIMITES (COFRUB)     ',/,
     &'@    NE SONT PAS RENSEIGNES POUR CERTAINES                   ',/,
     &'@        FACES DE BORD (Phase ',I10   ,')                    ',/,
     &'@                                                            ',/,
     &'@        Valeur minimale COFRUB ',E14.5                       ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier le codage de usray5.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 3500 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LE RAYONNEMENT EST ACTIVE.                              ',/,
     &'@                                                            ',/,
     &'@    Le scalaire ',I10   ,' devrait etre la temperature ou   ',/,
     &'@      l''enthalpie. On a ISCSTH = ',I10                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4000 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT (FLUNET    NON RENSEIGNE)       ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    Face = ',I10   ,' Zone = ',I10   ,' Type = ',I10           )
 4100 FORMAT(
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
     &'@        Valeur minimale ',E14.5                              ,/,
     &'@                                                            ',/,
     &'@    Le calcul ne sera pas execute.                          ',/,
     &'@                                                            ',/,
     &'@    Verifier le codage de usray5.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 5000 FORMAT('-----------------------------------------------------',
     &          '--------------')
C
 5010 FORMAT('Zone         Flux net radiatif (Watt) (normale',
     &          ' unitaire sortante)')
C
 5020 FORMAT(I6,13X,E10.4)
C
 5030 FORMAT('Flux net radiatif sur toutes les frontieres  Fnet = ',
     &           E10.4,' Watt')
C
 5040 FORMAT('Intégrale volumique du terme source radiatif Srad = ',
     &           E10.4,' Watt')
C
 5050 FORMAT('(Si IDIVER = 1 ou 2 alors on doit avoir Srad = -Fnet)')
C
 6000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT APPROXIMATION P-1  (RAYDOM)     ',/,
     &'@    *********                                               ',/,
     &'@                                                            ',/,
     &'@    LA LONGUEUR OPTIQUE DU MILIEU SEMI-TRANSPARENT          ',/,
     &'@      DOIT AU MOINS ETRE DE L''ORDRE DE L''UNITE POUR ETRE  ',/,
     &'@      DANS LE DOMAINE D''APPLICATION DE L''APPROXIMATION P-1',/,
     &'@    CELA NE SEMBLE PAS ETRE LE CAS ICI.                     ',/,
     &'@                                                            ',/,
     &'@    LE COEFFICIENT D''ABSORPTION MINIMUM POUR ASSURER CETTE ',/,
     &'@      LONGUEUR OPTIQUE EST XKMIN = ',E10.4                   ,/,
     &'@    CETTE VALEUR N''EST PAS ATTEINTE POUR ', E10.4,'%       ',/,
     &'@      DES CELLULES DU MAILLAGE.                             ',/,
     &'@    LE POURCENTAGE DE CELLULES DU MAILLAGE POUR LESQUELLES  ',/,
     &'@      ON ADMET QUE CETTE CONDITION SOIT VIOLEE EST IMPOSE   ',/,
     &'@      PAR DEFAUT OU DANS USINI1 A XNP1MX = ', E10.4,'%      ',/,
     &'@                                                            ',/,
     &'@    Verifier les valeurs du coefficient d''absorption CK    ',/,
     &'@      dans l''interface ou le modifier dans USRAY3.         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      END
c@z

