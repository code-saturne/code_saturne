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
                        SUBROUTINE RAYCLI
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , ISVHB  , ISVTB  ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ITRIFB , ITYPFB ,
     &   IZFRAD , ISOTHM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFA  , COEFB  , HBORD  , TBORD  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RAYEXP , RAYIMP ,
     &   TPAROI , QINCID , XLAM   , EPA    , EPS    ,
     &   FLUNET , FLCONV , HFCONV , TEXT   , TINT   , TEMPK  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE RAYONNEMENT :
CFONC    --------------------------------------
CFONC
CFONC   1) Calcul des temperatures de paroi
CFONC   2) Mise a jours des conditions aux limites de la variable
CFONC      energetique
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
CARGU ! ISVHB        ! E  !  -> ! INDICATEUR DE SAUVEGARDE DES         !
CARGU !              !    !     !  COEFFICIENTS D'ECHANGE AUX BORDS    !
CARGU ! ISVTB        ! E  !  -> ! INDICATEUR DE SAUVEGARDE DES         !
CARGU !              !    !     !  TEMPERATURES AUX BORDS              !
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
CARGU ! IZFRAD(NFABOR! TE ! <-> ! NUMERO DE ZONE DES FACES DE BORD     !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! ISOTHM(NFABOR! TE ! <-> ! TYPE DE CONDITION DE PAROI           !
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
CARGU ! COEFA, COEFB ! TR ! <-> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! HBORD        ! TR ! <-  ! COEFFICIENTS D'ECHANGE AUX BORDS     !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! TBORD        ! TR ! <-  ! TEMPERATURE AUX BORDS           I    !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! RAYEXP(NCELET! TR ! <-  ! TERME SOURCE RADIATIF EXPLICITE      !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! RAYIMP(NCELET! TR ! <-  ! TERME SOURCE RADIATIF IMPLICITE      !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! TEMPK(NCELET)! TR ! <-  ! TEMPERATURE EN KELVIN                !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! TPAROI(NFABOR! TR !  -  ! TEMPERATURE DE PAROI EN KELVIN       !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! QINCID(NFABOR! TR ! <-  ! DENSITE DE FLUX RADIATIF AUX BORDS   !
CARGU !   ,NPHAST)   !    !     !                                      !
CARGU ! XLAM (NFABOR ! TR ! <-  ! COEFFICIENT DE CONDUCTIVITE THERMIQUE!
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
CARGU ! TEXT (NFABOR ! TR ! <-  ! TEMPERATURE DE BORD EXTERNE          !
CARGU !   ,NPHAST)   !    !     ! EN DEGRES CELCIUS                    !
CARGU ! TINT (NFABOR ! TR ! <-  ! TEMPERATURE DE BORD INTERNE          !
CARGU !   ,NPHAST)   !    !     ! EN DEGRES CELCIUS                    !
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "radiat.h"
      INCLUDE "ihmpre.h"
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
      INTEGER          ISVHB  , ISVTB
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ITRIFB(NFABOR,NPHAS), ITYPFB(NFABOR,NPHAS)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
      INTEGER          IZFRAD(NFABOR,NPHAST),ISOTHM(NFABOR,NPHAST)
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
      DOUBLE PRECISION HBORD(NFABOR),TBORD(NFABOR)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
C
      DOUBLE PRECISION RAYEXP(NCELET,NPHAST), RAYIMP(NCELET,NPHAST)
      DOUBLE PRECISION TEMPK(NCELET,NPHAST)
      DOUBLE PRECISION TPAROI(NFABOR,NPHAST), QINCID(NFABOR,NPHAST)
      DOUBLE PRECISION XLAM(NFABOR,NPHAST), EPA(NFABOR,NPHAST)
      DOUBLE PRECISION EPS(NFABOR,NPHAST), FLUNET(NFABOR,NPHAST)
      DOUBLE PRECISION FLCONV(NFABOR,NPHAST), HFCONV(NFABOR,NPHAST)
      DOUBLE PRECISION TEXT(NFABOR,NPHAST), TINT(NFABOR,NPHAST)
C
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IEL, IPHAS, IPH, IDEB, IVART, ISCAT
      INTEGER          MODE, IOK, IFVU, II, IZONEM, IZONE
      INTEGER          MAXELT, IDBIA1, ILS
C
      DOUBLE PRECISION TMIN , TMAX   , TX
      DOUBLE PRECISION CPP, XMTK
C
      INTEGER    IPACLI
      DATA       IPACLI /0/
      SAVE       IPACLI
C
C***********************************************************************
C=======================================================================
C 0 - GESTION MEMOIRE
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      MAXELT = MAX(NCELET,NFAC,NFABOR)
C
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
C---> NUMERO DE PASSAGE RELATIF
C
      IPACLI = IPACLI + 1
      IDEB = 0
C
C---> VALEURS MIN ET MAX ADMISSIBLES POUR LA TEMPERATURE DE PAROI
C         EN KELVIN
C
      TMIN = 0.D0
      TMAX = GRAND + TKELVI
C
C---> COEFF DE RELAX
C
C      TX est strictement superieur a 0 et inferieur ou egal a 1
C
C      Pour calculer la temperature de paroi, on calcule un increment
C      de temperature DeltaT entre l'etape courante n et l'etape
C      precedente n-1, puis on calcule :
C           n    n-1                                 n-1
C          T  = T    + DeltaT si le rapport DeltaT/T    =< TX, sinon
C
C           n    n-1                      n-1             n-1
C          T  = T    * (1 + TX *((DeltaT/T   ) / |DeltaT/T   |))
C
      TX = 0.1D0
C
C
C---> INITIALISATIONS PAR DEFAUT BIDON
C
      DO IPH = 1,NPHAST
        DO IFAC = 1,NFABOR
          IZFRAD(IFAC,IPH) = -1
          ISOTHM(IFAC,IPH) = -1
          XLAM  (IFAC,IPH) = -GRAND
          EPA   (IFAC,IPH) = -GRAND
          EPS   (IFAC,IPH) = -GRAND
          TEXT  (IFAC,IPH) = -GRAND
          TINT  (IFAC,IPH) = -GRAND
        ENDDO
      ENDDO
C
C
C=======================================================================
C 2. SI PAS DE FICHIER SUITE ALORS INITIALISATION AU PREMIER PASSAGE
C    DE TPAROI ET QINCID :
C      LECTURE DE L'INITIALISATION DE TPAROI A TINT
C      QINCID EST INITIALISE A STEPHN*TINT**4 (SI ON INITIALISE QINCID
C      A ZERO, ON AURA UN DEFICIT SUR LA CONDITION LIMITE DE LUMINANCE
C      AUX PAROIS AU 1er PAS DE TEMPS EN DOM)
C=======================================================================
C
      IF (IPACLI.EQ.1 .AND. ISUIRD.EQ.0) THEN
C
C Indicateur : si non suite et premier pas de temps.
          IDEB = 1
C
          DO IPH = 1, NPHAST
C
            IPHAS = IRAPHA(IPH)
C
            DO IEL = 1,NCELET
              RAYIMP(IEL,IPH) = ZERO
              RAYEXP(IEL,IPH) = ZERO
            ENDDO
C
            DO IFAC = 1,NFABOR
              HFCONV(IFAC,IPH) = ZERO
              FLCONV(IFAC,IPH) = ZERO
            ENDDO
C
C     On utilise TBORD comme auxiliaire pour l'appel a USRAY2
C       pour être sur que TPAROI ne sera pas modifié
C       (puisqu'on a TBORD libre)
C     On utilise FLUNET comme auxiliaire pour l'appel a USRAY2
C       pour être sur que QINCID ne sera pas modifié
C       (puisqu'on a FLUNET libre)
C
            DO IFAC = 1,NFABOR
              TBORD(IFAC)      = ZERO
              FLUNET(IFAC,IPH) = ZERO
            ENDDO
C
C         - Interface Code_Saturne
C           ======================
C
            IF (IIHMPR.EQ.1) THEN
C
C---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE
              IVART = ISCA(ISCALT(IPHAS))
C
              CALL UIRAY2
C             ===========
     &       ( ITYPFB, IPAROI, IPARUG, IVART, IPH, NPHAST,  IZFRAD,
     &         ISOTHM, ITPIMP, IPGRNO, IPREFL, IFGRNO, IFREFL,
     &         NFABOR, NFML, IFMFBR , IPRFML , NVAR,
     &         EPS, EPA, TINT, TEXT,
     &         XLAM, RCODCL)
C
            ENDIF
C
            ILS    = IDEBIA
            IDBIA1 = ILS + MAXELT
            CALL IASIZE('RAYCLI',IDBIA1)
C
            CALL USRAY2
C           ===========
     & ( IDBIA1 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , IZFRAD(1,IPH) , ISOTHM(1,IPH) ,
     &   TMIN   , TMAX   , TX     ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
     &   TBORD  , FLUNET(1,IPH) , HFCONV(1,IPH) ,  FLCONV(1,IPH) ,
     &   XLAM(1,IPH)   , EPA(1,IPH)    , EPS(1,IPH)    ,
     &   TEXT(1,IPH)   , TINT(1,IPH)   ,
     &   RA     )
C
            WRITE(NFECRA,1000)
C
C Tparoi en Kelvin et QINCID en W/m2
            DO IFAC = 1,NFABOR
              IF ( ITYPFB(IFAC,IPHAS).EQ.IPAROI .OR.
     &             ITYPFB(IFAC,IPHAS).EQ.IPARUG ) THEN
                TPAROI(IFAC,IPH) = TINT(IFAC,IPH)
                QINCID(IFAC,IPH) = STEPHN*TINT(IFAC,IPH)**4
              ELSE
                TPAROI(IFAC,IPH) = 0.D0
                QINCID(IFAC,IPH) = 0.D0
              ENDIF
            ENDDO
C
C Fin de la boucle sur NPHAST
          ENDDO
C
C Fin détection premier passage
      ENDIF
C
C=======================================================================
C 3. BOUCLE SUR LES PHASES...
C=======================================================================
C
      DO IPH = 1,NPHAST
C
        IPHAS = IRAPHA(IPH)
C
C---> NUMERO DU SCALAIRE ET DE LA VARIABLE THERMIQUE
        ISCAT = ISCALT(IPHAS)
        IVART = ISCA(ISCALT(IPHAS))
C
C=======================================================================
C 3.1 DONNEES SUR LES FACES FRONTIERES
C=======================================================================
C
C     On utilise TBORD comme auxiliaire pour l'appel a USRAY2
C       pour être sur que TPAROI ne sera pas modifié
C       (puisqu'on a TBORD libre)
C     On utilise FLUNET comme auxiliaire pour l'appel a USRAY2
C       pour être sur que QINCID ne sera pas modifié
C       (puisqu'on a FLUNET libre)
C
        DO IFAC = 1,NFABOR
          TBORD (IFAC)     = TPAROI(IFAC,IPH)
          FLUNET(IFAC,IPH) = QINCID(IFAC,IPH)
        ENDDO
C
C     - Interface Code_Saturne
C       ======================
C
        IF (IIHMPR.EQ.1) THEN
C
          CALL UIRAY2
C         ===========
     &  ( ITYPFB, IPAROI, IPARUG, IVART, IPH, NPHAST, IZFRAD,
     &    ISOTHM, ITPIMP, IPGRNO, IPREFL, IFGRNO, IFREFL,
     &    NFABOR, NFML, IFMFBR , IPRFML , NVAR,
     &    EPS, EPA, TINT, TEXT,
     &    XLAM, RCODCL)
C
        ENDIF
C
        ILS    = IDEBIA
        IDBIA1 = ILS + MAXELT
        CALL IASIZE('RAYCLI',IDBIA1)
C
        CALL USRAY2
C       ===========
     & ( IDBIA1 , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB ,
     &   MAXELT , IA(ILS),
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , IZFRAD(1,IPH) , ISOTHM(1,IPH) ,
     &   TMIN   , TMAX   , TX     ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
     &   TBORD  , FLUNET(1,IPH) ,  HFCONV(1,IPH) ,  FLCONV(1,IPH) ,
     &   XLAM(1,IPH)   , EPA(1,IPH)    , EPS(1,IPH)    ,
     &   TEXT(1,IPH)   , TINT(1,IPH)   ,
     &   RA     )
C
C=======================================================================
C 3.2 CONTROLE DES DONNEES UTILISATEUR
C=======================================================================
C
C--> Arret si le numero de zone est non renseigne ou mal renseigne
C
        IOK = 0
C
        DO IFAC = 1, NFABOR
          IF (IZFRAD(IFAC,IPH).LE.0.OR.IZFRAD(IFAC,IPH).GT.NOZRDM) THEN
            IOK = IOK + 1
            WRITE(NFECRA,2000)IFAC,NOZRDM,IZFRAD(IFAC,IPH)
          ENDIF
        ENDDO
C
        IF(IOK.NE.0) THEN
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
C --> On construit une liste des numeros des zones frontieres.
C           (liste locale au processeur, en parallele)
C     Stop si depassement.
C
        NZFRAD(IPHAS) = 0
        DO IFAC = 1, NFABOR
          IFVU = 0
          DO II = 1, NZFRAD(IPHAS)
            IF (ILZRAD(II,IPHAS).EQ.IZFRAD(IFAC,IPH)) THEN
              IFVU = 1
            ENDIF
          ENDDO
          IF(IFVU.EQ.0) THEN
            NZFRAD(IPHAS) = NZFRAD(IPHAS) + 1
            IF(NZFRAD(IPHAS).LE.NBZRDM) THEN
              ILZRAD(NZFRAD(IPHAS),IPHAS) = IZFRAD(IFAC,IPH)
            ELSE
              WRITE(NFECRA,2001) NBZRDM
              WRITE(NFECRA,2002)(ILZRAD(II,IPHAS),II=1,NBZRDM)
              CALL CSEXIT (1)
C             ===========
            ENDIF
          ENDIF
        ENDDO
C
C ---> Plus grand numero de zone atteint
C
        IZONEM = 0
        DO II = 1, NZFRAD(IPHAS)
          IZONE = ILZRAD(II,IPHAS)
          IZONEM = MAX(IZONEM,IZONE)
        ENDDO
        IF(IRANGP.GE.0) THEN
          CALL PARCMX(IZONEM)
C         ===========
        ENDIF
        NOZARM(IPHAS) = IZONEM
C
C
C
C
C On verra si ca coute cher ou non.
C   Pour le moment on le fait tout le temps.
c        IF(IWARNI(IVART).GE.-1.OR.IPACLI.LE.3) THEN
        IF(1.EQ.1) THEN
C
          IOK = 0
C
C--> Si en paroi ISOTHM non renseignee : stop
          DO IFAC = 1, NFABOR
            IF( (ITYPFB(IFAC,IPHAS).EQ.IPAROI  .OR.
     &           ITYPFB(IFAC,IPHAS).EQ.IPARUG) .AND.
     &           ISOTHM(IFAC,IPH)  .EQ.-1    ) THEN
              IOK = IOK + 1
              WRITE(NFECRA,2110) IPHAS,IFAC,IZFRAD(IFAC,IPH)
            ENDIF
          ENDDO
C
C--> Si ISOTHM renseignee en non paroi : stop
          DO IFAC = 1, NFABOR
            IF( ITYPFB(IFAC,IPHAS).NE.IPAROI .AND.
     &          ITYPFB(IFAC,IPHAS).NE.IPARUG .AND.
     &          ISOTHM(IFAC,IPH)  .NE.-1         ) THEN
              IOK = IOK + 1
              WRITE(NFECRA,2111)
     &             IPHAS,IFAC,IZFRAD(IFAC,IPH),ISOTHM(IFAC,IPH)
            ENDIF
          ENDDO
C
C--> Si valeur physique erronee : stop
          DO IFAC = 1, NFABOR
            IF(ISOTHM(IFAC,IPH).EQ.ITPIMP ) THEN
              IF(EPS(IFAC,IPH) .LT.0.D0.OR.EPS(IFAC,IPH).GT.1.D0.OR.
     &           TINT(IFAC,IPH).LE.0.D0                          ) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2120) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &               EPS(IFAC,IPH),
     &                              TINT(IFAC,IPH)
              ENDIF
            ELSEIF(ISOTHM(IFAC,IPH).EQ.IPGRNO ) THEN
              IF(EPS(IFAC,IPH) .LT.0.D0.OR.EPS(IFAC,IPH).GT.1.D0.OR.
     &           XLAM(IFAC,IPH).LE.0.D0.OR.
     &           EPA (IFAC,IPH).LE.0.D0.OR.
     &           TEXT(IFAC,IPH).LE.0.D0.OR.
     &           TINT(IFAC,IPH).LE.0.D0                          ) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2130) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &               EPS(IFAC,IPH),XLAM(IFAC,IPH),EPA(IFAC,IPH),
     &               TEXT(IFAC,IPH),TINT(IFAC,IPH)
              ENDIF
            ELSEIF(ISOTHM(IFAC,IPH).EQ.IPREFL ) THEN
              IF(XLAM(IFAC,IPH).LE.0.D0.OR.
     &           EPA (IFAC,IPH).LE.0.D0.OR.
     &           TEXT(IFAC,IPH).LE.0.D0.OR.
     &           TINT(IFAC,IPH).LE.0.D0                          ) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2140) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &                             XLAM(IFAC,IPH),EPA(IFAC,IPH),
     &               TEXT(IFAC,IPH),TINT(IFAC,IPH)
              ENDIF
            ELSEIF(ISOTHM(IFAC,IPH).EQ.IFGRNO ) THEN
              IF(EPS(IFAC,IPH) .LT.0.D0.OR.EPS(IFAC,IPH).GT.1.D0.OR.
     &           TINT(IFAC,IPH).LE.0.D0                          ) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2150) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &               EPS(IFAC,IPH),
     &                              TINT(IFAC,IPH)
              ENDIF
            ELSEIF(ISOTHM(IFAC,IPH).EQ.IFREFL ) THEN
              IF(TINT(IFAC,IPH).LE.0.D0                          ) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2160) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &                              TINT(IFAC,IPH)
              ENDIF
            ELSEIF(ISOTHM(IFAC,IPH).NE.-1) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2170) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &                             ISOTHM(IFAC,IPH)
            ENDIF
          ENDDO
C
C--> Si valeur renseignee sans raison : stop
          DO IFAC = 1, NFABOR
           IF(ISOTHM(IFAC,IPH).EQ.ITPIMP ) THEN
              IF(XLAM(IFAC,IPH).GT.0.D0.OR.EPA(IFAC,IPH).GT.0.D0.OR.
     &           TEXT(IFAC,IPH).GT.0.D0                          ) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2220) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &               XLAM(IFAC,IPH),EPA(IFAC,IPH),TEXT(IFAC,IPH)
              ENDIF
            ELSEIF(ISOTHM(IFAC,IPH).EQ.IPREFL ) THEN
              IF(EPS(IFAC,IPH).GE.0.D0                           ) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2240) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &               EPS(IFAC,IPH)
              ENDIF
            ELSEIF(ISOTHM(IFAC,IPH).EQ.IFGRNO ) THEN
              IF(XLAM(IFAC,IPH).GT.0.D0.OR.EPA(IFAC,IPH).GT.0.D0.OR.
     &           TEXT(IFAC,IPH).GT.0.D0                          ) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2250) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &               XLAM(IFAC,IPH),EPA(IFAC,IPH),TEXT(IFAC,IPH)
              ENDIF
            ELSEIF(ISOTHM(IFAC,IPH).EQ.IFREFL ) THEN
              IF(EPS(IFAC,IPH).GE.0.D0.OR.
     &           XLAM(IFAC,IPH).GT.0.D0.OR.EPA(IFAC,IPH).GT.0.D0.OR.
     &           TEXT(IFAC,IPH).GT.0.D0                          ) THEN
                IOK = IOK + 1
                WRITE(NFECRA,2260) IPHAS,IFAC,IZFRAD(IFAC,IPH),
     &               EPS(IFAC,IPH),
     &               XLAM(IFAC,IPH),EPA(IFAC,IPH),TEXT(IFAC,IPH)
              ENDIF
            ENDIF
          ENDDO
C
C--> Stop si erreur
          IF(IOK.NE.0) THEN
            CALL CSEXIT (1)
C           ===========
          ENDIF
C
        ENDIF
C
C
C=======================================================================
C 3.2 COMPLETION DES DONNEES UTILISATEUR
C=======================================================================
C
C ICODCL et EPS (quand il est nul)
C
        DO IFAC = 1, NFABOR
          IF(    ISOTHM(IFAC,IPH).EQ.ITPIMP ) THEN
            ICODCL(IFAC,IVART) = 5
          ELSEIF(ISOTHM(IFAC,IPH).EQ.IPGRNO ) THEN
            ICODCL(IFAC,IVART) = 5
          ELSEIF(ISOTHM(IFAC,IPH).EQ.IPREFL ) THEN
            ICODCL(IFAC,IVART) = 5
            EPS(IFAC,IPH) = 0.D0
          ELSEIF(ISOTHM(IFAC,IPH).EQ.IFGRNO ) THEN
            ICODCL(IFAC,IVART) = 5
          ELSEIF(ISOTHM(IFAC,IPH).EQ.IFREFL ) THEN
            ICODCL(IFAC,IVART) = 3
            EPS(IFAC,IPH) = 0.D0
          ENDIF
        ENDDO
C
C
C=======================================================================
C 4. STOCKAGE DE LA TEMPERATURE (en Kelvin) dans TEMPK(IEL,IPH)
C=======================================================================
C
        IF (ABS(ISCSTH(ISCAT)).EQ.1) THEN
C
C---> ON REMPLIT TEMPK
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
        ELSE IF (ISCSTH(ISCAT).EQ.2) THEN
C
C---> LECTURES DES DONNEES UTILISATEURS (TBORD est un auxiliaire)
C
          MODE = 1
C
          IF (IPPMOD(IPHPAR).LE.1) THEN
C
            CALL USRAY4
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
C
     &   MODE   ,
C
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
C
     &   TPAROI(1,IPH)   , TBORD  , TEMPK(1,IPH)   ,
C                                   Resultat : T en K
C
     &   RA     )
C
          ELSE
C
            CALL PPRAY4
C           ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
C
     &   MODE   ,
C
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,
C
     &   TPAROI(1,IPH)   , TBORD  , TEMPK(1,IPH)   ,
C                                   Resultat : T en K
C
     &   RA     )
C
          ENDIF
C
        ENDIF
C
C=======================================================================
C 5.  CALCUL DES TEMPERATURES DE PAROIS
C=======================================================================
C
C
C DANS TOUS LES CAS HFCONV CONTIENT Lambda * Hturb / distance
C   (HFCONV : W/(m2 K) ; Hturb est sans dimension)
C  (au premier passage, il est nul)
C
C--> CALCUL DU FLUX CONVECTIF
C      Par flux convectif, on entend bien sur
C        flux convectif parallele a la paroi,
C        on suppose que la paroi est etanche...
C      Le flux est calcule dans condli clptur, sauf au premier
C        passage sans suite de calcul, puisque raycli est appele avant.
C
C
      IF (IDEB.EQ.1) THEN
C
        DO IFAC = 1,NFABOR
          IF (ISOTHM(IFAC,IPH).NE.-1) THEN
            FLCONV(IFAC,IPH) =
     &      HFCONV(IFAC,IPH)*(TEMPK(IFABOR(IFAC),IPH)-TPAROI(IFAC,IPH))
          ENDIF
        ENDDO
C
      ENDIF
C
C
C--> Les cas ou il faut calculer TPAROI sont, au premier passage sans suite
C      des cas a temperature imposee TPAROI = TINT
C
        IF (IDEB.EQ.1) THEN
C
          DO IFAC = 1,NFABOR
            IF (ISOTHM(IFAC,IPH).EQ.IPGRNO .OR.
     &          ISOTHM(IFAC,IPH).EQ.IPREFL .OR.
     &          ISOTHM(IFAC,IPH).EQ.IFGRNO    ) THEN
              ISOTHM(IFAC,IPH) = ITPIMP
            ENDIF
          ENDDO
C
        ENDIF
C
        IF(IDEB.EQ.0) THEN
C
          CALL RAYPAR
C         ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
C
     &   ICODCL , ISOTHM(1,IPH)  , IZFRAD(1,IPH) ,
C
     &   IDEVEL , ITUSER , IA     ,
C
     &   TMIN   , TMAX   , TX     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB , RCODCL ,
     &   COEFA  , COEFB  ,
     &   RDEVEL , RTUSER ,
C
     &   TPAROI(1,IPH) , QINCID(1,IPH) , TEXT(1,IPH)   , TINT(1,IPH)   ,
     &   XLAM(1,IPH)   , EPA(1,IPH)    , EPS(1,IPH)    ,
     &   HFCONV(1,IPH) , FLCONV(1,IPH) , TEMPK(1,IPH)  ,
C
     &   RA     )
C
        ENDIF
C
C
C=======================================================================
C 6.  CHANGEMENT DES CONDITIONS LIMITES UTILISATEUR
C=======================================================================
C
C=======================================================================
C 6.1  LA VARIABLE TRANSPORTEE EST LA TEMPERATURE
C=======================================================================
C
        IF (ABS(ISCSTH(ISCAT)).EQ.1) THEN
C
          IF(ISCSTH(ISCAT).EQ.-1) THEN
            XMTK = -TKELVI
          ELSE
            XMTK = 0.D0
          ENDIF
C
          DO IFAC = 1,NFABOR
C
            IF (ISOTHM(IFAC,IPH).EQ.ITPIMP .OR.
     &          ISOTHM(IFAC,IPH).EQ.IPGRNO .OR.
     &          ISOTHM(IFAC,IPH).EQ.IFGRNO    ) THEN
              RCODCL(IFAC,IVART,1) = TPAROI(IFAC,IPH)+XMTK
              RCODCL(IFAC,IVART,2) = RINFIN
              RCODCL(IFAC,IVART,3) = 0.D0
C
            ELSE IF (ISOTHM(IFAC,IPH).EQ.IPREFL) THEN
              RCODCL(IFAC,IVART,1) = TEXT(IFAC,IPH)+XMTK
              RCODCL(IFAC,IVART,2) = XLAM(IFAC,IPH)/EPA(IFAC,IPH)
              RCODCL(IFAC,IVART,3) = 0.D0
C
            ELSE IF (ISOTHM(IFAC,IPH).EQ.IFREFL) THEN
              ICODCL(IFAC,IVART) = 3
              RCODCL(IFAC,IVART,1) = 0.D0
              RCODCL(IFAC,IVART,2) = RINFIN
            ENDIF
C
          ENDDO
C
C=======================================================================
C 6.2  LA VARIABLE TRANSPORTEE EST L'ENTHALPIE
C=======================================================================
C
        ELSEIF (ISCSTH(ISCAT).EQ.2) THEN
C
C---> LECTURES DES DONNEES UTILISATEURS
C     ON CONVERTIT TPAROI EN ENTHALPIE DE BORD, STOCKEE DANS FLUNET,
C     QUI EST UTILISE COMME AUXILIAIRE
C
          MODE = 0
C
          DO IFAC = 1,NFABOR
            IF (ISOTHM(IFAC,IPH).EQ.ITPIMP.OR.
     &          ISOTHM(IFAC,IPH).EQ.IPGRNO.OR.
     &          ISOTHM(IFAC,IPH).EQ.IFGRNO    ) THEN
              MODE = -1
            ENDIF
          ENDDO
C
          IF (MODE.EQ.-1) THEN
C
            IF (IPPMOD(IPHPAR).LE.1) THEN
C
              CALL USRAY4
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
C
     &   MODE   ,
C
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,

     &   TPAROI(1,IPH)   , FLUNET(1,IPH) , TEMPK(1,IPH)  ,
C                          HPAROI
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
C
     &   MODE   ,
C
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,

     &   TPAROI(1,IPH)   , FLUNET(1,IPH) , TEMPK(1,IPH)  ,
C                          HPAROI
     &   RA     )
C
            ENDIF
C
          ENDIF
C
          MODE = 0
C
          DO IFAC = 1,NFABOR
            IF (ISOTHM(IFAC,IPH).EQ.IPREFL) THEN
              MODE = -1
            ENDIF
          ENDDO
C
          IF (MODE.EQ.-1) THEN
C
            IF (IPPMOD(IPHPAR).LE.1) THEN
C
              CALL USRAY4
C             ===========
     & ( IDEBIA , IDEBRA ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , IPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
C
     &   MODE   ,
C
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,

     &   TEXT(1,IPH)  , TBORD  , TEMPK(1,IPH)  ,
C                       HEXT
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
C
     &   MODE   ,
C
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , ITYPFB(1,IPHAS) ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  ,
     &   W1     , W2     , W3     , W4     , W5     , W6     ,
     &   RDEVEL , RTUSER ,

     &   TEXT(1,IPH)  , TBORD  , TEMPK(1,IPH)  ,
C                       HEXT
     &   RA     )
C
            ENDIF
C
          ENDIF
C
          DO IFAC = 1,NFABOR
C
            IF (ISOTHM(IFAC,IPH).EQ.ITPIMP.OR.
     &          ISOTHM(IFAC,IPH).EQ.IPGRNO.OR.
     &          ISOTHM(IFAC,IPH).EQ.IFGRNO    ) THEN
              RCODCL(IFAC,IVART,1) = FLUNET(IFAC,IPH)
              RCODCL(IFAC,IVART,2) = RINFIN
              RCODCL(IFAC,IVART,3) = 0.D0
C
            ELSE IF (ISOTHM(IFAC,IPH).EQ.IPREFL) THEN
C
              IF (ICP(IPHAS).GT.0) THEN
                IEL = IFABOR(IFAC)
                CPP = PROPCE(IEL,IPPROC(ICP(IPHAS)))
              ELSE
                CPP = CP0(IPHAS)
              ENDIF
C
              RCODCL(IFAC,IVART,1) = TBORD(IFAC)
              RCODCL(IFAC,IVART,2) =
     &             XLAM(IFAC,IPH)/(EPA(IFAC,IPH)*CPP)
              RCODCL(IFAC,IVART,3) = 0.D0
C
            ELSE IF (ISOTHM(IFAC,IPH).EQ.IFREFL) THEN
              ICODCL(IFAC,IVART) = 3
              RCODCL(IFAC,IVART,1) = 0.D0
              RCODCL(IFAC,IVART,2) = RINFIN
            ENDIF
C
          ENDDO
C
        ENDIF
C
C=======================================================================
C 7.  FIN DE LA BOUCLE SUR LES PHASES
C=======================================================================
C
      ENDDO
C
C--------
C FORMATS
C--------
C
 1000 FORMAT (/, 3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT ',/,
     &           3X,'   ------------------------------------------',/,
     &           3X,' Initialisation de la temperature de paroi   ',/,
     &           3X,' (TPAROI) avec le profil utilisateur (TINTP) ',/,
     &           3X,' et du flux incident aux parois (QINCID).    ',/)
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@    LES CONDITIONS AUX LIMITES SONT INCOMPLETES OU ERRONEES ',/,
     &'@                                                            ',/,
     &'@  Le numero de zone associee a la face ',I10   ,' doit etre ',/,
     &'@    un entier strictement positif et inferieur ou egal a    ',/,
     &'@    NOZRDM = ',I10                                           ,/,
     &'@  Ce numero (IZFRDP(IFAC)) vaut ici ',I10                    ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les conditions aux limites dans usray2.          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : RAYONNEMENT                                 ',/,
     &'@    *********                                               ',/,
     &'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,
     &'@                                                            ',/,
     &'@  Le nombre maximal de zones frontieres qui peuvent etre    ',/,
     &'@    definies par l''utilisateur est NBZRDM = ',I10           ,/,
     &'@    Il a ete depasse.                                       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier les conditions aux limites dans usray2.          ',/,
     &'@                                                            ',/,
     &'@  Les NBZRDM premieres zones frontieres                     ',/,
     &'@    portent ici les numeros suivants :                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2002 FORMAT(I10)
C
 2110 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    ISOTHP DOIT ETRE RENSEIGNE SUR TOUTES LES FACES DE PAROI',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Il ne l''a pas ete pour la face ',I10                      ,/,
     &'@                    zone         ',I10                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2111 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    ISOTHP A ETE RENSEIGNE SUR UNE FACE NON PAROI           ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Sur la face ',I10   ,', zone  ',I10   ,', ISOTHP a ete    ',/,
     &'@    renseigne dans usray2 (ISOTHP = ',I10   ,') alors que   ',/,
     &'@    la face n''a pas ete declaree de type IPAROI ou IPARUG  ',/,
     &'@    dans usclim.                                            ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2 et usclim.                                ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2120 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    EPSP  DOIT ETRE UN REEL INCLUS DANS [0.; 1.]            ',/,
     &'@    TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF             ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ITPIMP       ',/,
     &'@  EPSP  = ',E14.5                                            ,/,
     &'@  TINTP = ',E14.5                                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2130 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    EPSP DOIT ETRE UN REEL INCLUS DANS [0.; 1.]             ',/,
     &'@    XLAMP, EPAP, TINTP, TEXTP DOIVENT ETRE DES REELS        ',/,
     &'@                                      STRICTEMENT POSITIFS  ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPGRNO       ',/,
     &'@  EPSP  = ',E14.5                                            ,/,
     &'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,
     &'@  TEXTP = ',E14.5    ,' TINTP = ',E14.5                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2140 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    XLAMP, EPAP, TINTP, TEXTP DOIVENT ETRE DES REELS        ',/,
     &'@                                      STRICTEMENT POSITIFS  ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPREFL       ',/,
     &'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,
     &'@  TEXTP = ',E14.5    ,' TINTP = ',E14.5                      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2150 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    EPSP  DOIT ETRE UN REEL INCLUS DANS [0.; 1.]            ',/,
     &'@    TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF             ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFGRNO       ',/,
     &'@  EPSP  = ',E14.5                                            ,/,
     &'@  TINTP = ',E14.5                                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2160 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@  TINTP DOIT ETRE UN REEL STRICTEMENT POSITIF               ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFREFL       ',/,
     &'@  TINTP = ',E14.5                                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2170 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@   VALEUR NON ADMISSIBLE DE ISOTHP                          ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ',I10         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2220 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    XLAMP, EPAP ET TEXTP NE DOIVENT PAS ETRE RENSEIGNES     ',/,
     &'@                                     AVEC ISOTHP = ITPIMP   ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = ITPIMP       ',/,
     &'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,
     &'@  TEXTP = ',E14.5                                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2240 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    EPSP NE DOIT PAS ETRE RENSEIGNE AVEC ISOTHP = IPREFL    ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IPREFL       ',/,
     &'@  EPSP  = ',E14.5                                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2250 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    XLAMP, EPAP, TEXTP NE DOIVENT PAS ETRE RENSEIGNES       ',/,
     &'@                                       AVEC ISOTHP = IFGRNO ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFGRNO       ',/,
     &'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,
     &'@  TEXTP = ',E14.5                                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 2260 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES RAYONNEMENT   ',/,
     &'@    *********                                               ',/,
     &'@    XLAMP, EPAP, TEXTP NE DOIVENT PAS ETRE RENSEIGNES       ',/,
     &'@                                       AVEC ISOTHP = IFREFL ',/,
     &'@                                                            ',/,
     &'@  Phase ',I10                                                ,/,
     &'@  Face = ',I10   ,' Zone = ',I10   ,' ISOTHP = IFREFL       ',/,
     &'@  EPSP  = ',E14.5                                            ,/,
     &'@  XLAMP = ',E14.5    ,' EPAP  = ',E14.5                      ,/,
     &'@  TEXTP = ',E14.5                                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usray2.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
