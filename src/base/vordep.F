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
                        SUBROUTINE VORDEP
C                       *****************
C     -------------------------------------------------------------
     & ( NCEVOR , NVOR   , IENT   , DTREF  ,
     &   IVORCE , YZCEL  , XU     , XV     , XW     ,
     &   YZVOR  , YZVORA , SIGNV  , TEMPS  , TPSLIM )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ----------
c@foncb
CFONC
CFONC  METHODE DES VORTEX POUR LES ENTREES EN L.E.S. :
CFONC    - DEPLACEMENT DES VORTEX
CFONC    - TRAITEMENT DES VORTEX MORTS ET SORTANTS
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCEVOR       ! E  !  -> ! NOMBRE DE FACE A L'ENTREE OU EST     !
CARGU !              !    !     ! UTILISE LA METHODE                   !
CARGU ! NVOR         ! E  !  -> ! NOMBRE DE VORTEX A L'ENTREE          !
CARGU ! IENT         ! E  !  -> ! NUMERO DE L'ENTREE                   !
CARGU ! DTREF        ! R  !  -> ! PAS DE TEMPS                         !
CARGU ! IVORCE       ! TE ! <-> ! NUMERO DU VORTEX LE PLUS PROCHE D'UNE!
CARGU !     (NVOMAX) !    !     ! FACE DONNEE                          !
CARGU ! YZCEL        ! TR !  -> ! COORDONNEES DES FACES D'ENTREE DANS  !
CARGU !   (ICVMAX ,2)!    !     ! LE REFERENTIEL LOCAL                 !
CARGU ! XU(ICVMAX)   ! TR !  -  ! COMPOSANTE DE VITESSE PRINCIPALE     !
CARGU ! XV(ICVMAX)   ! TR ! <-> ! COMPOSANTES DE VITESSE TRANSVERSES   !
CARGU ! XW(ICVMAX)   ! TR ! <-> !                                      !
CARGU ! YZVOR        ! TR ! <-  ! NOUVELLES COORDONNEES DU CENTRE      !
CARGU !   (NVOMAX,2) !    !     ! DES VORTEX                           !
CARGU ! YZVORA       ! TR !  -> ! ANCIENNES COORDONNEES DU CENTRE      !
CARGU !   (NVOMAX,2) !    !     ! DES VORTEX                           !
CARGU ! SIGNV(NVOMAX)! TR ! <-> ! SENS DE ROTATION DES VORTEX          !
CARGU ! TEMPS        ! TR ! <-> ! TEMPS ECOULE DEPUIS LA CREATION      !
CARGU !     (NVOMAX) !    !     ! DU VORTEX                            !
CARGU ! TPSLIM       ! TR ! <-> ! DUREE DE VIE DU VORTEX               !
CARGU !     (NVOMAX) !    !     !                                      !
CARGU .______________.____._____.______________________________________.
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
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "vortex.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCEVOR , NVOR   , IENT
      INTEGER          IVORCE(NVOMAX)
C
      DOUBLE PRECISION DTREF
      DOUBLE PRECISION YZCEL(ICVMAX ,2)
      DOUBLE PRECISION XU(ICVMAX)      , XV(ICVMAX)      , XW(ICVMAX)
      DOUBLE PRECISION YZVOR(NVOMAX,2) , YZVORA(NVOMAX,2)
      DOUBLE PRECISION SIGNV(NVOMAX)
      DOUBLE PRECISION TEMPS(NVOMAX)   , TPSLIM(NVOMAX)
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ, KK, III, IUN
C
      DOUBLE PRECISION SENS, DY, DZ
      DOUBLE PRECISION XXV, XXW
      DOUBLE PRECISION DRAND(1), U_VOR, EK_VOR, EE_VOR
      DOUBLE PRECISION DD, YY, ZZ
      DOUBLE PRECISION PHIDAT
C
C=======================================================================
C 1. DEPLACEMENT DES VORTEX
C=======================================================================
      IUN = 1
      DO II = 1, NVOR
        TEMPS(II) = TEMPS(II) + DTREF
      ENDDO
C
C - Deplacement aleatoire
C
      IF(IDEPVO(IENT).EQ.1)THEN
C
        DO II = 1, NVOR
          YZVORA(II,1) = YZVOR(II,1)
          YZVORA(II,2) = YZVOR(II,2)
        ENDDO
        DO II = 1, NVOR
          SENS = 1.D0
          CALL ZUFALL(IUN,DRAND(1))
          IF (DRAND(1).LT.0.5D0) SENS = -1.D0
          CALL ZUFALL(IUN,DRAND(1))
          DY = DRAND(1) * UD(IENT) * DTREF
          YZVOR(II,1) = YZVOR(II,1) + SENS * DY
          SENS = 1.D0
          CALL ZUFALL(IUN,DRAND(1))
          IF (DRAND(1).LT.0.5D0) SENS = -1.D0
          CALL ZUFALL(IUN,DRAND(1))
          DZ = DRAND(1) * UD(IENT) * DTREF
          YZVOR(II,2) = YZVOR(II,2) + SENS * DZ
        ENDDO
C
C - Convection des vortex
C
      ELSEIF(IDEPVO(IENT).EQ.2) THEN
C
        DO II = 1, NVOR
          YZVORA(II,1) = YZVOR(II,1)
          YZVORA(II,2) = YZVOR(II,2)
        ENDDO
C
        DO II = 1, NVOR
          KK = IVORCE(II)
          XXV = XV(KK)
          XXW = XW(KK)
          YZVOR(II,1) = YZVOR(II,1) + DTREF * XXV
          YZVOR(II,2) = YZVOR(II,2) + DTREF * XXW
        ENDDO
C
      ENDIF
C
C=======================================================================
C 2. GESTION DES VORTEX SORTANT DU DOMAINE
C=======================================================================
C
      IF(ICAS(IENT).EQ.1) THEN
C
        IF(ICLVOR(1,IENT).EQ.1.OR.ICLVOR(1,IENT).EQ.2) THEN
          DO II = 1, NVOR
            IF(YZVOR(II,1).GT.(LLY(IENT)/2.D0)) THEN
              YZVOR(II,1) = YZVORA(II,1)
            ENDIF
          ENDDO
        ELSEIF(ICLVOR(1,IENT).EQ.3) THEN
          DO II = 1, NVOR
            IF(YZVOR(II,1).GT.(LLY(IENT)/2.D0).AND.
     &        YZVOR(II,1).LT.(3.D0*LLY(IENT)/2.D0)) THEN
              YZVOR(II,1) = YZVOR(II,1) - LLY(IENT)
            ELSEIF(YZVOR(II,1).LT.-(LLY(IENT)/2.D0).AND.
     &        YZVOR(II,1).GT.-(3.D0*LLY(IENT)/2.D0)) THEN
              YZVOR(II,1) = YZVOR(II,1) + LLY(IENT)
            ELSEIF(YZVOR(II,1).GT.(3.D0*LLY(IENT)/2.D0).OR.
     &        YZVOR(II,1).LT.-(3.D0*LLY(IENT)/2.D0)) THEN
              YZVOR(II,1) = YZVORA(II,1)
            ENDIF
          ENDDO
        ENDIF
C
        IF(ICLVOR(2,IENT).EQ.1.OR.ICLVOR(2,IENT).EQ.2) THEN
          DO II = 1, NVOR
            IF(YZVOR(II,2).GT.(LLZ(IENT)/2.D0)) THEN
              YZVOR(II,2) = YZVORA(II,2)
            ENDIF
          ENDDO
        ELSEIF(ICLVOR(2,IENT).EQ.3) THEN
          DO II = 1, NVOR
            IF(YZVOR(II,2).GT.(LLZ(IENT)/2.D0).AND.
     &         YZVOR(II,2).LT.(3.D0*LLZ(IENT)/2.D0)) THEN
              YZVOR(II,2) = YZVOR(II,2) - LLZ(IENT)
            ELSEIF(YZVOR(II,2).LT.-(LLZ(IENT)/2.D0).AND.
     &        YZVOR(II,2).GT.-(3.D0*LLZ(IENT)/2.D0)) THEN
              YZVOR(II,2) = YZVOR(II,2) + LLZ(IENT)
            ELSEIF(YZVOR(II,2).GT.(3.D0*LLZ(IENT)/2.D0).OR.
     &        YZVOR(II,2).LT.-(3.D0*LLZ(IENT)/2.D0)) THEN
              YZVOR(II,2) = YZVORA(II,2)
            ENDIF
          ENDDO
        ENDIF
C
        IF(ICLVOR(3,IENT).EQ.1.OR.ICLVOR(3,IENT).EQ.2) THEN
          DO II = 1, NVOR
            IF(YZVOR(II,1).LT.-(LLY(IENT)/2.D0)) THEN
              YZVOR(II,1) = YZVORA(II,1)
            ENDIF
          ENDDO
        ENDIF
C
        IF(ICLVOR(4,IENT).EQ.1.OR.ICLVOR(4,IENT).EQ.2) THEN
          DO II = 1, NVOR
            IF(YZVOR(II,2).LT.-(LLZ(IENT)/2.D0)) THEN
              YZVOR(II,2) = YZVORA(II,2)
            ENDIF
          ENDDO
        ENDIF
C
      ELSEIF(ICAS(IENT).EQ.2) THEN
        DO II = 1, NVOR
          IF((YZVOR(II,1)**2+YZVOR(II,2)**2).GT.
     &      (LLD(IENT)/2.0D0)**2)THEN
            YZVOR(II,1) = YZVORA(II,1)
            YZVOR(II,2) = YZVORA(II,2)
          ENDIF
        ENDDO
C
      ELSEIF(ICAS(IENT).EQ.3.OR.ICAS(IENT).EQ.4) THEN
        DO II = 1, NVOR
          IF(YZVOR(II,1).LT.YMIN(IENT).OR.
     &      YZVOR(II,1).GT.YMAX(IENT)) THEN
            YZVOR(II,1) = YZVORA(II,1)
          ENDIF
          IF(YZVOR(II,2).LT.ZMIN(IENT).OR.
     &      YZVOR(II,2).GT.ZMAX(IENT)) THEN
            YZVOR(II,2) = YZVORA(II,2)
          ENDIF
        ENDDO
      ENDIF
C
C=======================================================================
C 3. REGENERATION DES VORTEX MORTS
C=======================================================================
C
      DO II = 1, NVOR
        IF(TEMPS(II).GT.TPSLIM(II))THEN
          TEMPS(II) = 0.D0
C
C - Position
C
          IF(ICAS(IENT).EQ.1) THEN
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,1) = LLY(IENT) * DRAND(1) - LLY(IENT) / 2.D0
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,2) = LLZ(IENT) * DRAND(1) - LLZ(IENT) / 2.D0
          ELSEIF(ICAS(IENT).EQ.2) THEN
 15   CONTINUE
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,1) = LLD(IENT) * DRAND(1) - LLD(IENT) / 2.0D0
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,2) = LLD(IENT) * DRAND(1) - LLD(IENT) / 2.0D0
            IF((YZVOR(II,1)**2+YZVOR(II,2)**2).GT.
     &         (LLD(IENT)/2.0D0)**2) THEN
              GOTO 15
            ENDIF
          ELSEIF(ICAS(IENT).EQ.3.OR.ICAS(IENT).EQ.4) THEN
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,1) = YMIN(IENT) + LLY(IENT) * DRAND(1)
            CALL ZUFALL(IUN,DRAND(1))
            YZVOR(II,2) = ZMIN(IENT) + LLZ(IENT) * DRAND(1)
          ENDIF
C
C - Duree de vie
C
          IF(ITLIVO(IENT).EQ.1) THEN
            TPSLIM(II) = TLIMVO(IENT)
          ELSEIF(ITLIVO(IENT).EQ.2) THEN
            YY = YZVOR(II,1)
            ZZ = YZVOR(II,2)
            III = 0
            U_VOR  = PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &               YDAT(1,IENT),ZDAT(1,IENT),UDAT(1,IENT),III)
            EK_VOR = PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &               YDAT(1,IENT),ZDAT(1,IENT),KDAT(1,IENT),III)
            EE_VOR = PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &               YDAT(1,IENT),ZDAT(1,IENT),EPSDAT(1,IENT),III)
            TPSLIM(II) = 5.D0*CMU*EK_VOR**(3.D0/2.D0)/EE_VOR
            TPSLIM(II) = TPSLIM(II)/U_VOR
          ENDIF
C
C - Sens de rotation
C
          CALL ZUFALL(IUN,DRAND(1))
          IF(DRAND(1).GT.0.5D0) SIGNV(II) = -1.0D0*SIGNV(II)
        ENDIF
      ENDDO
C
C=======================================================================
C 4. RECHERCHE DE LA FACE LA PLUS PROCHE DE CHAQUE VORTEX
C=======================================================================
C
      DO II = 1, NVOR
         KK = 0
         DD = GRAND
         DO JJ = 1, NCEVOR
           IF(((YZCEL(JJ,1)-YZVOR(II,1))**2+
     &         (YZCEL(JJ,2)-YZVOR(II,2))**2).LT.DD)THEN
             DD = (YZCEL(JJ,1)-YZVOR(II,1))**2
     &           +(YZCEL(JJ,2)-YZVOR(II,2))**2
             KK = JJ
           ENDIF
         ENDDO
         IVORCE(II) = KK
      ENDDO
C
C ---
C FIN
C ---
C
      RETURN
C
      END
c@z
