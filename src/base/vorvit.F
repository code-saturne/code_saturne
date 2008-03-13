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
                        SUBROUTINE VORVIT
C                       *****************
C     -------------------------------------------------------------
     & ( NCEVOR , NVOR   , IENT   , IVORCE , VISCO  ,
     &   YZCEL  , XU     , XV     , XW     ,
     &   YZVOR  , SIGNV  , SIGMA  , GAMMA  , TEMPS )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  METHODE DES VORTEX POUR LES CONDITIONS AUX LIMITES D'ENTREE EN
CFONC    L.E.S. : CALCUL DES LA VITESSE DES VORTEX
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCEVOR       ! E  !  -> ! NOMBRE DE FACE A L'ENTREE OU EST     !
CARGU !              !    !     ! UTILISE LA METHODE                   !
CARGU ! NVOR         ! E  !     ! NOMBRE DE VORTEX A L'ENTREE OU EST   !
CARGU !              !    !     ! UTILISEE LA METHODE                  !
CARGU ! IENT         ! E  !  -> ! NUMERO DE L'ENTREE OU EST UTILISEE   !
CARGU !              !    !     ! LA METHODE                           !
CARGU ! IVORCE       ! TE ! <-> ! NUMERO DU VORTEX LE PLUS PROCHE D'UNE!
CARGU !     (NVOMAX) !    !     ! FACE DONNEE                          !
CARGU ! VISCO        ! TR !  -> ! VISCOSITE CINEMATIQUE SUR LES FACES  !
CARGU !(ICVMAX,NNENT)!    !     ! D'ENTREE                             !
CARGU ! YZCEL        ! TR !  -> ! COORDONNEES DES FACES D'ENTREE DANS  !
CARGU !   (ICVMAX ,2)!    !     ! LE REFERENTIEL LOCAL                 !
CARGU ! XU(ICVMAX)   ! TR ! <-> ! COMPOSANTE DE VITESSE PRINCIPALE     !
CARGU ! XV(ICVMAX)   ! TR ! <-> ! COMPOSANTES DE VITESSE TRANSVERSES   !
CARGU ! XW(ICVMAX)   ! TR ! <-> !                                      !
CARGU ! YZVOR        ! TR ! <-> ! COORDONNEES DU CENTRE DES VORTEX     !
CARGU !   (NVOMAX,2) !    !     !                                      !
CARGU ! SIGNV(NVOMAX)! TR ! <-> ! SENS DE ROTATION DES VORTEX          !
CARGU ! SIGMA        ! TR ! <-> ! TAILLE DES VORTEX                    !
CARGU !(NVOMAX,NNENT)!    !     !                                      !
CARGU ! GAMMA        ! TR ! <-> ! INTENSITE (CIRCULATION) DES VORTEX   !
CARGU !(NVOMAX,2,NNEN!    !     ! DANS LES DEUX DIRECTIONS DU PLAN     !
CARGU ! TEMPS        ! TR ! <-> ! TEMPS ECOULE DEPUIS LA CREATION      !
CARGU !     (NVOMAX) !    !     ! DU VORTEX                            !
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
      DOUBLE PRECISION YZCEL(ICVMAX,2) , VISCO(ICVMAX)
      DOUBLE PRECISION XU(ICVMAX )     , XV(ICVMAX )     , XW(ICVMAX )
      DOUBLE PRECISION YZVOR(NVOMAX,2) , SIGNV(NVOMAX)   , SIGMA(NVOMAX)
      DOUBLE PRECISION GAMMA(NVOMAX,2) , TEMPS(NVOMAX)
C
C VARIABLES LOCALES
C
C
      INTEGER          II, JJ, III
C
      DOUBLE PRECISION YY, ZZ, XVISC
      DOUBLE PRECISION NORME, VV, WW, THETA, DV, DW
      DOUBLE PRECISION YPAROI, ZPAROI, YPERIO, ZPERIO, YSYM, ZSYM
      DOUBLE PRECISION PHIDAT
      DOUBLE PRECISION ALPHA, EK_VOR, EE_VOR, V_VOR, W_VOR
      DOUBLE PRECISION LT, LK, DEUXPI
C
C=======================================================================
C 1. CALCUL DE GAMMA
C=======================================================================
C
      ALPHA = 4.D0 * SQRT(PI * XSURFV(IENT)/
     &       (3.D0 * NVOR*(2.D0*LOG(3.D0)-3.D0*LOG(2.D0))))
C
      DO II = 1, NVOR
        YY = YZVOR(II,1)
        ZZ = YZVOR(II,2)
        III = 0
        EK_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &       YDAT(1,IENT),ZDAT(1,IENT),KDAT(1,IENT),III)
        GAMMA(II,1) = ALPHA*SQRT(EK_VOR)*SIGNV(II)
        GAMMA(II,2) = GAMMA(II,1)
      ENDDO
C
C=======================================================================
C 2. CALCUL DE SIGMA
C=======================================================================
C
      IF(ISGMVO(IENT).EQ.1) THEN
        DO II = 1, NVOR
          SIGMA(II) = XSGMVO(IENT)
        ENDDO
      ELSEIF (ISGMVO(IENT).EQ.2) THEN
        DO II = 1, NVOR
          YY = YZVOR(II,1)
          ZZ = YZVOR(II,2)
          III = 0
          EK_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &         YDAT(1,IENT),ZDAT(1,IENT),KDAT(1,IENT),III)
          EE_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &         YDAT(1,IENT),ZDAT(1,IENT),EPSDAT(1,IENT),III)
          SIGMA(II) = (CMU**(0.75D0))*(EK_VOR**(1.5D0))/EE_VOR
        ENDDO
      ELSEIF (ISGMVO(IENT).EQ.3) THEN
        DO II = 1, NVOR
          YY = YZVOR(II,1)
          ZZ = YZVOR(II,2)
          III = 0
          EK_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &         YDAT(1,IENT),ZDAT(1,IENT),KDAT(1,IENT),III)
          EE_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &         YDAT(1,IENT),ZDAT(1,IENT),EPSDAT(1,IENT),III)
          XVISC  =  VISCO(IVORCE(II))
          LT = SQRT(5.D0*XVISC*EK_VOR/EE_VOR)
          LK = 200.D0*(XVISC**3/EE_VOR)**(0.25D0)
          SIGMA(II) = MAX(LT,LK)
        ENDDO
      ENDIF
C
C=======================================================================
C 3. CALCUL DU CHAMP DE VITESSE INDUIT (AU CENTRE DES CELLULES)
C=======================================================================
C
      DEUXPI = 2.D0*PI
C
      DO II = 1, NCEVOR
        VV = 0.D0
        WW = 0.D0
        DO JJ = 1, NVOR
          YY = YZVOR(JJ,1)-YZCEL(II,1)
          ZZ = YZVOR(JJ,2)-YZCEL(II,2)
          NORME = YY**2+ZZ**2
          ALPHA = SIGMA(JJ)
          IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
            THETA = -NORME/(2.D0*ALPHA**2)
            THETA = EXP(THETA)
            DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
            DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
            VV = VV - DV
            WW = WW + DW
          ENDIF
C
C=======================================================================
C 4. TRAITEMENT DES CONDITIONS AUX LIMITES
C=======================================================================
C suite des boucles en DO
C -----------------------
C Conduite rectangulaire
C -----------------------
          IF(ICAS(IENT).EQ.1) THEN
C
C Periodicites
C
            IF(ICLVOR(1,IENT).EQ.3) THEN
              IF(YZVOR(JJ,1).GT.0.D0) THEN
                YPERIO = YZVOR(JJ,1) - LLY(IENT)
                ZPERIO = YZVOR(JJ,2)
              ELSE
                YPERIO = YZVOR(JJ,1) + LLY(IENT)
                ZPERIO = YZVOR(JJ,2)
              ENDIF
              YY = YPERIO - YZCEL(II,1)
              ZZ = ZPERIO - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
                DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
                VV = VV - DV
                WW = WW + DW
              ENDIF
            ENDIF
C
            IF(ICLVOR(2,IENT).EQ.3) THEN
              IF(YZVOR(JJ,2).GT.0.D0) THEN
                YPERIO = YZVOR(JJ,1)
                ZPERIO = YZVOR(JJ,2) - LLZ(IENT)
              ELSE
                YPERIO = YZVOR(JJ,1)
                ZPERIO = YZVOR(JJ,2) + LLZ(IENT)
              ENDIF
              YY = YPERIO - YZCEL(II,1)
              ZZ = ZPERIO - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
                DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
                VV = VV - DV
                WW = WW + DW
              ENDIF
            ENDIF
C
C Parois
C
            IF(ICLVOR(1,IENT).EQ.1) THEN
              YPAROI = LLY(IENT)/2.D0
              ZPAROI = YZCEL(II,2)
              YPAROI = 2.D0*YPAROI - YZVOR(JJ,1)
              ZPAROI = 2.D0*ZPAROI - YZVOR(JJ,2)
              YY = YPAROI - YZCEL(II,1)
              ZZ = ZPAROI - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
                DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
                VV = VV - DV
                WW = WW + DW
              ENDIF
            ENDIF
C
            IF(ICLVOR(3,IENT).EQ.1) THEN
              YPAROI = - LLY(IENT)/2.D0
              ZPAROI = YZCEL(II,2)
              YPAROI = 2.D0*YPAROI - YZVOR(JJ,1)
              ZPAROI = 2.D0*ZPAROI - YZVOR(JJ,2)
              YY = YPAROI - YZCEL(II,1)
              ZZ = ZPAROI - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
                DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
                VV = VV - DV
                WW = WW + DW
              ENDIF
            ENDIF
C
            IF(ICLVOR(2,IENT).EQ.1) THEN
              YPAROI = YZCEL(II,1)
              ZPAROI = LLZ(IENT)/2.D0
              YPAROI = 2.D0*YPAROI - YZVOR(JJ,1)
              ZPAROI = 2.D0*ZPAROI - YZVOR(JJ,2)
              YY = YPAROI - YZCEL(II,1)
              ZZ = ZPAROI - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
                DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
                VV = VV - DV
                WW = WW + DW
              ENDIF
            ENDIF
C
            IF(ICLVOR(4,IENT).EQ.1) THEN
              YPAROI = YZCEL(II,1)
              ZPAROI = -LLZ(IENT)/2.D0
              YPAROI = 2.D0*YPAROI - YZVOR(JJ,1)
              ZPAROI = 2.D0*ZPAROI - YZVOR(JJ,2)
              YY = YPAROI - YZCEL(II,1)
              ZZ = ZPAROI - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
                DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
                VV = VV - DV
                WW = WW + DW
              ENDIF
            ENDIF
C
C Symetries
C
            IF(ICLVOR(1,IENT).EQ.2) THEN
              YSYM = LLY(IENT)/2.D0
              ZSYM = YZCEL(II,2)
              YSYM = 2.D0*YSYM - YZVOR(JJ,1)
              ZSYM = 2.D0*ZSYM - YZVOR(JJ,2)
              YY = YSYM - YZCEL(II,1)
              ZZ = ZSYM - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
                VV = VV - DV
              ENDIF
            ENDIF
C
            IF(ICLVOR(3,IENT).EQ.2) THEN
              YSYM = - LLY(IENT)/2.D0
              ZSYM = YZCEL(II,2)
              YSYM = 2.D0*YSYM - YZVOR(JJ,1)
              ZSYM = 2.D0*ZSYM - YZVOR(JJ,2)
              YY = YSYM - YZCEL(II,1)
              ZZ = ZSYM - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
                VV = VV - DV
              ENDIF
            ENDIF
C
            IF(ICLVOR(2,IENT).EQ.2) THEN
              YSYM = YZCEL(II,1)
              ZSYM = LLZ(IENT)/2.D0
              YSYM = 2.D0*YSYM - YZVOR(JJ,1)
              ZSYM = 2.D0*ZSYM - YZVOR(JJ,2)
              YY = YSYM - YZCEL(II,1)
              ZZ = ZSYM - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
                WW = WW + DW
              ENDIF
            ENDIF
C
            IF(ICLVOR(4,IENT).EQ.2) THEN
              YSYM = YZCEL(II,1)
              ZSYM = -LLZ(IENT)/2.D0
              YSYM = 2.D0*YSYM - YZVOR(JJ,1)
              ZSYM = 2.D0*ZSYM - YZVOR(JJ,2)
              YY = YSYM - YZCEL(II,1)
              ZZ = ZSYM - YZCEL(II,2)
              NORME = YY**2+ZZ**2
              ALPHA = SIGMA(JJ)
              IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
                THETA = -NORME/(2.D0*ALPHA**2)
                THETA = EXP(THETA)
                DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
                WW = WW + DW
              ENDIF
            ENDIF
C
C --------------------
C Conduite circulaire
C --------------------
          ELSEIF(ICAS(IENT).EQ.2) THEN
C
            YPAROI = YZCEL(II,1)*((LLD(IENT)/2.0D0)
     &           /SQRT(YZCEL(II,1)**2 + YZCEL(II,2)**2))
            ZPAROI = YZCEL(II,2)*((LLD(IENT)/2.0D0)
     &           /SQRT(YZCEL(II,1)**2 + YZCEL(II,2)**2))
            YPAROI = 2.D0*YPAROI - YZVOR(JJ,1)
            ZPAROI = 2.D0*ZPAROI - YZVOR(JJ,2)
C
            YY = YPAROI - YZCEL(II,1)
            ZZ = ZPAROI - YZCEL(II,2)
            NORME = YY**2+ZZ**2
            ALPHA = SIGMA(JJ)
C
            IF(NORME.GT.EPZERO.AND.NORME.LT.4.D0*ALPHA) THEN
              THETA = -NORME/(2.D0*ALPHA**2)
              THETA = EXP(THETA)
              DV = ZZ/NORME*(1.D0-THETA)*GAMMA(JJ,1)*THETA/DEUXPI
              DW = YY/NORME*(1.D0-THETA)*GAMMA(JJ,2)*THETA/DEUXPI
              VV = VV - DV
              WW = WW + DW
            ENDIF
C
          ENDIF
C
        ENDDO
        XV(II) = VV
        XW(II) = WW
      ENDDO
C
C=======================================================================
C 5. AJOUT DE LA VITESSE MOYENNE POUR V ET W
C=======================================================================
      IF(ICAS(IENT).EQ.1.OR.ICAS(IENT).EQ.2.OR.ICAS(IENT).EQ.3) THEN
        DO II = 1, NCEVOR
          YY = YZCEL(II,1)
          ZZ = YZCEL(II,2)
          III = 0
          V_VOR = PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &         YDAT(1,IENT),ZDAT(1,IENT),VDAT(1,IENT),III)
          W_VOR = PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &         YDAT(1,IENT),ZDAT(1,IENT),WDAT(1,IENT),III)
          XV(II) = XV(II) + V_VOR
          XW(II) = XW(II) + W_VOR
        ENDDO
      ENDIF
C
C----
C FIN
C----
C
      RETURN
      END
c@z
