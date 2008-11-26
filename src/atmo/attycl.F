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
                        SUBROUTINE ATTYCL
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NBMETD , NBMETT , NBMETM ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ITRIFB , ITYPFB , IZFPPP , IPROFM ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL ,
     &   TMPROM , ZTPROM , ZDPROM , XMET   , YMET   , PMER   ,
     &   TTPROM , QVPROM , UPROM  , VPROM  , EKPROM , EPPROM ,
     &   RPROM  , TPPROM , PHPROM ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     CONDITIONS AUX LIMITES AUTOMATIQUES
CFONC
CFONC            ECOULEMENTS ATMOSPHERIQUES
CFONC
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
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! IZFPPP       ! TE !  -> ! NUMERO DE ZONE DE LA FACE DE BORD    !
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
CARGU ! COEFA, COEFB ! TR !  -> ! CONDITIONS AUX LIMITES AUX           !
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
C ARGUMENTS
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
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
      INTEGER          NBMETD , NBMETT , NBMETM
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ITRIFB(NFABOR,NPHAS), ITYPFB(NFABOR,NPHAS)
      INTEGER          IZFPPP(NFABOR), IPROFM(NOZPPM)
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
      DOUBLE PRECISION TMPROM(NBMETM)
      DOUBLE PRECISION ZTPROM(NBMETT) , ZDPROM(NBMETD)
      DOUBLE PRECISION XMET(NBMETM)   , YMET(NBMETM)  , PMER(NBMETM)
      DOUBLE PRECISION TTPROM(NBMETT,NBMETM) , QVPROM(NBMETT,NBMETM)
      DOUBLE PRECISION UPROM(NBMETD,NBMETM)  , VPROM(NBMETD,NBMETM)
      DOUBLE PRECISION EKPROM(NBMETD,NBMETM) , EPPROM(NBMETD,NBMETM)
      DOUBLE PRECISION RPROM(NBMETT,NBMETM)  , TPPROM(NBMETT,NBMETM)
      DOUBLE PRECISION PHPROM(NBMETT,NBMETM)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION COEFU(NFABOR,NDIM)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IZONE, IPHAS
      DOUBLE PRECISION D2S3, ZENT, VS, XUENT, XVENT
      DOUBLE PRECISION XKENT, XEENT, TPENT
C
C***********************************************************************
C=======================================================================
C 1.  INITIALISATIONS
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      D2S3 = 2.D0/3.D0
C
      XUENT = 0.D0
      XVENT = 0.D0
      XKENT = 0.D0
      XEENT = 0.D0
      TPENT = 0.D0
C
C=======================================================================
C 2.  SI IPROFM = 1 : CHOIX ENTREE/SORTIE SUIVANT LE PROFIL METEO SI
C                       ITYPFB N'A PAS ETE MODIFIE
C                     VARIABLES TIREES DU PROFIL METEO SI
C                       RCODCL(IFAC,IVAR,1) N'A PAS ETE MODIFIE
C
C       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
C                     =========================
C=======================================================================
C
C
      DO IFAC = 1, NFABOR
C
        IZONE = IZFPPP(IFAC)
C
        IF (IPROFM(IZONE).EQ.1) THEN
C
C     On recupere les valeurs du profil et on met a jour RCODCL s'il n'a pas
C       ete modifie. Il servira si la face est une face d'entree ou si c'est une
C       face de sortie (si le flux est rentrant).
          ZENT=CDGFBO(3,IFAC)
C
          CALL INTPRF
C         ===========
     &   (NBMETD, NBMETM,
     &    ZDPROM, TMPROM, UPROM , ZENT  , TTCABS, XUENT )
C
          CALL INTPRF
C         ===========
     &   (NBMETD, NBMETM,
     &    ZDPROM, TMPROM, VPROM , ZENT  , TTCABS, XVENT )
C
          CALL INTPRF
C         ===========
     &   (NBMETD, NBMETM,
     &    ZDPROM, TMPROM, EKPROM, ZENT  , TTCABS, XKENT )
C
          CALL INTPRF
C         ===========
     &   (NBMETD, NBMETM,
     &    ZDPROM, TMPROM, EPPROM, ZENT  , TTCABS, XEENT )

          CALL INTPRF
C         ===========
     &   (NBMETT, NBMETM,
     &    ZTPROM, TMPROM, TPPROM, ZENT  , TTCABS, TPENT )
C
          DO IPHAS = 1, NPHAS
C
            IF (RCODCL(IFAC,IU(IPHAS),1).GT.RINFIN*0.5D0)
     &           RCODCL(IFAC,IU(IPHAS),1) = XUENT
            IF (RCODCL(IFAC,IV(IPHAS),1).GT.RINFIN*0.5D0)
     &           RCODCL(IFAC,IV(IPHAS),1) = XVENT
            IF (RCODCL(IFAC,IW(IPHAS),1).GT.RINFIN*0.5D0)
     &           RCODCL(IFAC,IW(IPHAS),1) = 0.D0
C
            IF    (ITYTUR(IPHAS).EQ.2) THEN
C
              IF (RCODCL(IFAC,IK(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IK(IPHAS),1) = XKENT
              IF (RCODCL(IFAC,IEP(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IEP(IPHAS),1) = XEENT
C
            ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
              IF (RCODCL(IFAC,IR11(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IR11(IPHAS),1) = D2S3*XKENT
              IF (RCODCL(IFAC,IR22(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IR22(IPHAS),1) = D2S3*XKENT
              IF (RCODCL(IFAC,IR33(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IR33(IPHAS),1) = D2S3*XKENT
              IF (RCODCL(IFAC,IR12(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IR12(IPHAS),1) = 0.D0
              IF (RCODCL(IFAC,IR13(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IR13(IPHAS),1) = 0.D0
              IF (RCODCL(IFAC,IR23(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IR23(IPHAS),1) = 0.D0
              IF (RCODCL(IFAC,IEP(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IEP(IPHAS),1) = XEENT
C
            ELSEIF(ITURB(IPHAS).EQ.50) THEN
C
              IF (RCODCL(IFAC,IK(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IK(IPHAS),1) = XKENT
              IF (RCODCL(IFAC,IEP(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IEP(IPHAS),1) = XEENT
              IF (RCODCL(IFAC,IPHI(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IPHI(IPHAS),1) = D2S3
              IF (RCODCL(IFAC,IFB(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IFB(IPHAS),1) = 0.D0
C
            ELSEIF(ITURB(IPHAS).EQ.60) THEN
C
              IF (RCODCL(IFAC,IK(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IK(IPHAS),1) = XKENT
              IF (RCODCL(IFAC,IOMG(IPHAS),1).GT.RINFIN*0.5D0)
     &             RCODCL(IFAC,IOMG(IPHAS),1) = XEENT/CMU/XKENT
C
            ENDIF
C
          ENDDO
C
          IF (NSCAL.GE.1) THEN
C
            IF (RCODCL(IFAC,ISCA(1),1).GT.RINFIN*0.5D0)
     &           RCODCL(IFAC,ISCA(1),1) = TPENT
C
          ENDIF
C
C     On met a jour le type de face de bord s'il n'a pas ete specifie
C       par l'utilisateur
C
          VS = XUENT*SURFBO(1,IFAC) + XVENT*SURFBO(2,IFAC)
C
          IF (VS.GT.0) THEN
            DO IPHAS = 1, NPHAS
              IF (ITYPFB(IFAC,IPHAS).EQ.0) ITYPFB(IFAC,IPHAS) = ISOLIB
            ENDDO
          ELSE
            DO IPHAS = 1, NPHAS
              IF (ITYPFB(IFAC,IPHAS).EQ.0) ITYPFB(IFAC,IPHAS) = IENTRE
            ENDDO
          ENDIF
C
        ENDIF
C
      ENDDO
C
C----
C FORMATS
C----
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
