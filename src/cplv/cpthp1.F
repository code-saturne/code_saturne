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
                        SUBROUTINE CPTHP1
C                       *****************
C     -------------------------------------------------------------
     & ( MODE  , EH     , XESP   , F1MC   , F2MC   ,
     &   TP     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC CALCUL DE LA TEMPERATURE DU GAZ
CFONC   EN FONCTION DE L'ENTHALPIE DU GAZ ET DES CONCENTRATIONS
CFONC   SI MODE = 1
CFONC CALCUL DE L'ENTHALPIE DU GAZ
CFONC   EN FONCTION DE LA TEMPERATURE DU GAZ ET DES CONCENTRATIONS
CFONC   SI MODE = -1
CFONC
c@fonce
C                             ARGUMENTS
c@argub
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! EH           ! TR !  -> ! ENTHALPIE DU GAZ                     !
CARGU !              !    !     ! (J/kg de melange gazeux)             !
CARGU ! XESP         ! TR !  -> ! FRACTION MASSIQUE DES ESPECES        !
CARGU ! F1MC         ! TR !  -> ! f1 moyen                             !
CARGU ! F2MC         ! TR !  -> ! f2 moyen                             !
CARGU ! TP           ! TR ! <-  ! TEMPERATURE DU GAZ (Kelvin)          !
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
C     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C**********************************************************************
C     DONNEES EN COMMON
C**********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          MODE
      DOUBLE PRECISION EH,TP
      DOUBLE PRECISION XESP(NGAZEM)
      DOUBLE PRECISION F1MC(NCHARM),F2MC(NCHARM)
C
C VARIABLES LOCALES
C
      INTEGER          I , ICHA
C
      DOUBLE PRECISION YCHX10 , YCHX20 , EHCHX1 , EHCHX2
      DOUBLE PRECISION DEN1   , DEN2
      DOUBLE PRECISION EH0 , EH1
C
C***********************************************************************
C=======================================================================
C 1. CALCUL DE LA TEMPERATURE A PARTIR DE l'ENTHALPIE
C=======================================================================
C
      IF ( MODE .EQ. 1 ) THEN
C
        I = NPO
C
C --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
C                                            et CHx2m a TH(NPO)
        EHCHX1 = ZERO
        EHCHX2 = ZERO
        YCHX10 = ZERO
        YCHX20 = ZERO
        DO ICHA = 1, NCHARB
          DEN1   = 1.D0
     &         / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
          YCHX10 = YCHX10 + DEN1 *
     &         ( F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          EHCHX1 = EHCHX1 + DEN1 *
     &         ( EHGAZE(ICHX1C(ICHA),I)*
     &         F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          DEN2   = 1.D0
     &         / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
          YCHX20 = YCHX20 + DEN2 *
     &         ( F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          EHCHX2 = EHCHX2 + DEN2 *
     &         ( EHGAZE(ICHX2C(ICHA),I)*
     &         F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
        ENDDO
        IF ( YCHX10.GT.EPZERO ) THEN
          EHCHX1 = EHCHX1 / YCHX10
        ELSE
          EHCHX1 = EHGAZE(ICHX1,I)
        ENDIF
        IF ( YCHX20.GT.EPZERO ) THEN
          EHCHX2 = EHCHX2 / YCHX20
        ELSE
          EHCHX2 = EHGAZE(ICHX2,I)
        ENDIF
C
C --- Clipping eventuel de TP a TH(NPO) si EH > EH1
C
        EH1 = XESP(ICHX1)*EHCHX1
     &      + XESP(ICHX2)*EHCHX2
     &      + XESP(ICO  )*EHGAZE(ICO ,I)
     &      + XESP(IO2  )*EHGAZE(IO2 ,I)
     &      + XESP(ICO2 )*EHGAZE(ICO2,I)
     &      + XESP(IH2O )*EHGAZE(IH2O,I)
     &      + XESP(IN2  )*EHGAZE(IN2 ,I)
C
        IF ( EH .GE. EH1 ) THEN
          TP = TH(I)
          GOTO 501
        ENDIF
C
        I = 1
C
C --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
C                                            et CHx2m a TH(1)
        EHCHX1 = ZERO
        EHCHX2 = ZERO
        YCHX10 = ZERO
        YCHX20 = ZERO
        DO ICHA = 1, NCHARB
          DEN1   = 1.D0
     &         / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
          YCHX10 = YCHX10 + DEN1 *
     &     ( F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          EHCHX1 = EHCHX1 + DEN1 *
     &     ( EHGAZE(ICHX1C(ICHA),I)*
     &       F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          DEN2   = 1.D0
     &         / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
          YCHX20 = YCHX20 + DEN2 *
     &     ( F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          EHCHX2 = EHCHX2 + DEN2 *
     &     ( EHGAZE(ICHX2C(ICHA),I)*
     &       F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
        ENDDO
        IF ( YCHX10.GT.EPZERO ) THEN
          EHCHX1 = EHCHX1 / YCHX10
        ELSE
          EHCHX1 = EHGAZE(ICHX1,I)
        ENDIF
        IF ( YCHX20.GT.EPZERO ) THEN
          EHCHX2 = EHCHX2 / YCHX20
        ELSE
          EHCHX2 = EHGAZE(ICHX2,I)
        ENDIF
C
C --- Clipping eventuel de TP a TH(1) si EH < EH0
C
        EH0 = XESP(ICHX1)*EHCHX1
     &      + XESP(ICHX2)*EHCHX2
     &      + XESP(ICO  )*EHGAZE(ICO ,I)
     &      + XESP(IO2  )*EHGAZE(IO2 ,I)
     &      + XESP(ICO2 )*EHGAZE(ICO2,I)
     &      + XESP(IH2O )*EHGAZE(IH2O,I)
     &      + XESP(IN2  )*EHGAZE(IN2 ,I)
C
        IF ( EH .LE. EH0 ) THEN
          TP= TH(I)
          GOTO 501
        ENDIF
C
C
 500    CONTINUE
        I = I + 1
C
C --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
C                                            et CHx2m pour TH(I-1)
        EHCHX1 = ZERO
        EHCHX2 = ZERO
        YCHX10 = ZERO
        YCHX20 = ZERO
        DO ICHA = 1, NCHARB
          DEN1   = 1.D0
     &           / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
          YCHX10 = YCHX10 + DEN1 *
     &      ( F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          EHCHX1 = EHCHX1 + DEN1 *
     &      ( EHGAZE(ICHX1C(ICHA),I-1)*
     &        F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          DEN2   = 1.D0
     &           / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
          YCHX20 = YCHX20 + DEN2 *
     &       ( F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          EHCHX2 = EHCHX2 + DEN2 *
     &       ( EHGAZE(ICHX2C(ICHA),I-1)*
     &         F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
        ENDDO
        IF ( YCHX10.GT.EPZERO ) THEN
          EHCHX1 = EHCHX1 / YCHX10
        ELSE
          EHCHX1 = EHGAZE(ICHX1,I-1)
        ENDIF
        IF ( YCHX20.GT.EPZERO ) THEN
          EHCHX2 = EHCHX2 / YCHX20
        ELSE
          EHCHX2 = EHGAZE(ICHX2,I-1)
        ENDIF
        EH0 = XESP(ICHX1)*EHCHX1
     &      + XESP(ICHX2)*EHCHX2
     &      + XESP(ICO  )*EHGAZE(ICO ,I-1)
     &      + XESP(IO2  )*EHGAZE(IO2 ,I-1)
     &      + XESP(ICO2 )*EHGAZE(ICO2,I-1)
     &      + XESP(IH2O )*EHGAZE(IH2O,I-1)
     &      + XESP(IN2  )*EHGAZE(IN2 ,I-1)
C
C --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
C                                            et CHx2m pour TH(I)
        EHCHX1 = ZERO
        EHCHX2 = ZERO
        YCHX10 = ZERO
        YCHX20 = ZERO
        DO ICHA = 1, NCHARB
          DEN1   = 1.D0
     &           / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
          YCHX10 = YCHX10 + DEN1 *
     &      ( F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          EHCHX1 = EHCHX1 + DEN1 *
     &      ( EHGAZE(ICHX1C(ICHA),I)*
     &        F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          DEN2   = 1.D0
     &           / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
          YCHX20 = YCHX20 + DEN2 *
     &       ( F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          EHCHX2 = EHCHX2 + DEN2 *
     &       ( EHGAZE(ICHX2C(ICHA),I)*
     &         F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
        ENDDO
        IF ( YCHX10.GT.EPZERO ) THEN
          EHCHX1 = EHCHX1 / YCHX10
        ELSE
          EHCHX1 = EHGAZE(ICHX1,I)
        ENDIF
        IF ( YCHX20.GT.EPZERO ) THEN
          EHCHX2 = EHCHX2 / YCHX20
        ELSE
          EHCHX2 = EHGAZE(ICHX2,I)
        ENDIF
C
        EH1 = XESP(ICHX1)*EHCHX1
     &      + XESP(ICHX2)*EHCHX2
     &      + XESP(ICO  )*EHGAZE(ICO ,I)
     &      + XESP(IO2  )*EHGAZE(IO2 ,I)
     &      + XESP(ICO2 )*EHGAZE(ICO2,I)
     &      + XESP(IH2O )*EHGAZE(IH2O,I)
     &      + XESP(IN2  )*EHGAZE(IN2 ,I)
C
        IF ( EH .GE. EH0  .AND. EH .LE. EH1  ) THEN
          TP = TH(I-1) + (EH-EH0) *
     &                 (TH(I)-TH(I-1))/(EH1-EH0)
          GOTO 501
        ENDIF
        GOTO 500
 501    CONTINUE
C
C=======================================================================
C 1. CALCUL DE L'ENTHALPIE A PARTIR DE LA TEMPERATURE
C=======================================================================
C
      ELSE IF ( MODE .EQ. -1 ) THEN
C
        I = NPO
C
C --- Calcul en Max
C
        IF ( TP .GE. TH(I) ) THEN
          EHCHX1 = ZERO
          EHCHX2 = ZERO
          YCHX10 = ZERO
          YCHX20 = ZERO
          DO ICHA = 1, NCHARB
            DEN1   = 1.D0
     &           / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
            YCHX10 = YCHX10 + DEN1 *
     &           ( F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            EHCHX1 = EHCHX1 + DEN1 *
     &           ( EHGAZE(ICHX1C(ICHA),I)*
     &           F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            DEN2   = 1.D0
     &           / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
            YCHX20 = YCHX20 + DEN2 *
     &           ( F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
            EHCHX2 = EHCHX2 + DEN2 *
     &           ( EHGAZE(ICHX2C(ICHA),I)*
     &           F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          ENDDO
          IF ( YCHX10.GT.EPZERO ) THEN
            EHCHX1 = EHCHX1 / YCHX10
          ELSE
            EHCHX1 = EHGAZE(ICHX1,I)
          ENDIF
          IF ( YCHX20.GT.EPZERO ) THEN
            EHCHX2 = EHCHX2 / YCHX20
          ELSE
            EHCHX2 = EHGAZE(ICHX2,I)
          ENDIF
          EH = XESP(ICHX1)*EHCHX1
     &       + XESP(ICHX2)*EHCHX2
     &       + XESP(ICO  )*EHGAZE(ICO ,I)
     &       + XESP(IO2  )*EHGAZE(IO2 ,I)
     &       + XESP(ICO2 )*EHGAZE(ICO2,I)
     &       + XESP(IH2O )*EHGAZE(IH2O,I)
     &       + XESP(IN2  )*EHGAZE(IN2 ,I)
          GOTO 601
        ENDIF
C
C Clipping en Min
C
        I = 1
C
        IF ( TP .LE. TH(I) ) THEN
          EHCHX1 = ZERO
          EHCHX2 = ZERO
          YCHX10 = ZERO
          YCHX20 = ZERO
          DO ICHA = 1, NCHARB
            DEN1   = 1.D0
     &           / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
            YCHX10 = YCHX10 + DEN1 *
     &           ( F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            EHCHX1 = EHCHX1 + DEN1 *
     &           ( EHGAZE(ICHX1C(ICHA),I)*
     &           F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            DEN2   = 1.D0
     &           / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
            YCHX20 = YCHX20 + DEN2 *
     &           ( F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
            EHCHX2 = EHCHX2 + DEN2 *
     &           ( EHGAZE(ICHX2C(ICHA),I)*
     &           F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          ENDDO
          IF ( YCHX10.GT.EPZERO ) THEN
            EHCHX1 = EHCHX1 / YCHX10
          ELSE
            EHCHX1 = EHGAZE(ICHX1,I)
          ENDIF
          IF ( YCHX20.GT.EPZERO ) THEN
            EHCHX2 = EHCHX2 / YCHX20
          ELSE
            EHCHX2 = EHGAZE(ICHX2,I)
          ENDIF
C
C --- Clipping eventuel de TP a TH(1) si EH < EH0
C
          EH = XESP(ICHX1)*EHCHX1
     &       + XESP(ICHX2)*EHCHX2
     &       + XESP(ICO  )*EHGAZE(ICO ,I)
     &       + XESP(IO2  )*EHGAZE(IO2 ,I)
     &       + XESP(ICO2 )*EHGAZE(ICO2,I)
     &       + XESP(IH2O )*EHGAZE(IH2O,I)
     &       + XESP(IN2  )*EHGAZE(IN2 ,I)
          GOTO 601
        ENDIF
C
C Interpolation dans la table
C
        I = 1
 600    CONTINUE
C
        I = I + 1
        IF ( TP .LE. TH(I) ) THEN
C
C --- Calcul de l'enthalpie de l'espece gazeuse TH(I-1)
C
          EHCHX1 = ZERO
          EHCHX2 = ZERO
          YCHX10 = ZERO
          YCHX20 = ZERO
          DO ICHA = 1, NCHARB
            DEN1   = 1.D0
     &             / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
            YCHX10 = YCHX10 + DEN1 *
     &        ( F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            EHCHX1 = EHCHX1 + DEN1 *
     &        ( EHGAZE(ICHX1C(ICHA),I-1)*
     &          F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            DEN2   = 1.D0
     &             / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
            YCHX20 = YCHX20 + DEN2 *
     &         ( F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
            EHCHX2 = EHCHX2 + DEN2 *
     &         ( EHGAZE(ICHX2C(ICHA),I-1)*
     &           F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          ENDDO
          IF ( YCHX10.GT.EPZERO ) THEN
            EHCHX1 = EHCHX1 / YCHX10
          ELSE
            EHCHX1 = EHGAZE(ICHX1,I-1)
          ENDIF
          IF ( YCHX20.GT.EPZERO ) THEN
            EHCHX2 = EHCHX2 / YCHX20
          ELSE
            EHCHX2 = EHGAZE(ICHX2,I-1)
          ENDIF
C
          EH0 = XESP(ICHX1)*EHCHX1
     &        + XESP(ICHX2)*EHCHX2
     &        + XESP(ICO  )*EHGAZE(ICO ,I-1)
     &        + XESP(IO2  )*EHGAZE(IO2 ,I-1)
     &        + XESP(ICO2 )*EHGAZE(ICO2,I-1)
     &        + XESP(IH2O )*EHGAZE(IH2O,I-1)
     &        + XESP(IN2  )*EHGAZE(IN2 ,I-1)
C
C --- Calcul de l'enthalpie de l'espece gazeuse TH(I)
C
          EHCHX1 = ZERO
          EHCHX2 = ZERO
          YCHX10 = ZERO
          YCHX20 = ZERO
          DO ICHA = 1, NCHARB
            DEN1   = 1.D0
     &             / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
            YCHX10 = YCHX10 + DEN1 *
     &        ( F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            EHCHX1 = EHCHX1 + DEN1 *
     &        ( EHGAZE(ICHX1C(ICHA),I)*
     &          F1MC(ICHA)*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            DEN2   = 1.D0
     &             / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
            YCHX20 = YCHX20 + DEN2 *
     &         ( F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
            EHCHX2 = EHCHX2 + DEN2 *
     &         ( EHGAZE(ICHX2C(ICHA),I)*
     &           F2MC(ICHA)*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          ENDDO
          IF ( YCHX10.GT.EPZERO ) THEN
            EHCHX1 = EHCHX1 / YCHX10
          ELSE
            EHCHX1 = EHGAZE(ICHX1,I)
          ENDIF
          IF ( YCHX20.GT.EPZERO ) THEN
            EHCHX2 = EHCHX2 / YCHX20
          ELSE
            EHCHX2 = EHGAZE(ICHX2,I)
          ENDIF
C
          EH1 = XESP(ICHX1)*EHCHX1
     &        + XESP(ICHX2)*EHCHX2
     &        + XESP(ICO  )*EHGAZE(ICO ,I)
     &        + XESP(IO2  )*EHGAZE(IO2 ,I)
     &        + XESP(ICO2 )*EHGAZE(ICO2,I)
     &        + XESP(IH2O )*EHGAZE(IH2O,I)
     &        + XESP(IN2  )*EHGAZE(IN2 ,I)
C
          EH =  EH0
     &           + (EH1-EH0)*(TP-TH(I-1))/(TH(I)-TH(I-1))
          GOTO 601
        ENDIF
        GOTO 600
C
 601    CONTINUE
C
      ELSE
        WRITE(NFECRA,1000) MODE
        CALL CSEXIT(1)
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR DANS CPTHP1                          ',/,
     &'@    *********                                               ',/,
     &'@    VALEUR INCORRECTE DE L''ARGUMENT MODE                   ',/,
     &'@    CE DOIT ETRE UN ENTIER EGAL A 1 OU -1                   ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
