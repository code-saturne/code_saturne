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
                        SUBROUTINE CPTEHM
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NTBMCI , NTBMCR ,
     &   EH     , RTP    ,
     &   FUEL1  , FUEL2  , FUEL3  , OXYD   , PROD1  , PROD2  ,
     &   XINER   ,
     &   TP     ,
     &   TBMCI  , TBMCR  ,
     &   EH0    , EH1    , X1     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC CALCUL DE LA TEMPERATURE DU MELANGE
CFONC   EN FONCTION DE L'ENTHALPIE DU MELANGE ET DES CONCENTRATIONS
CFONC
c@fonce
C                             ARGUMENTS
c@argub
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NTBMCI       ! E  !  -> ! TAILLE DU MACRO TABLEAU MC ENTIERS   !
CARGU ! NTBMCR       ! E  !  -> ! TAILLE DU MACRO TABLEAU MC REELS     !
CARGU ! EH           ! TR !  -> ! ENTHALPIE DU GAZ                     !
CARGU !              !    !     ! (J/kg de melange gazeux)             !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! FUEL1        ! TR !  -> ! FRACTION MASSIQUE CHx1               !
CARGU ! FUEL2        ! TR !  -> ! FRACTION MASSIQUE CHx2               !
CARGU ! FUEL3        ! TR !  -> ! FRACTION MASSIQUE CO                 !
CARGU ! OXYD         ! TR !  -> ! FRACTION MASSIQUE O2                 !
CARGU ! PROD1        ! TR !  -> ! FRACTION MASSIQUE CO2                !
CARGU ! PROD2        ! TR !  -> ! FRACTION MASSIQUE H2O                !
CARGU ! XINER        ! TR !  -> ! FRACTION MASSIQUE N2                 !
CARGU ! TP           ! TR ! <-  ! TEMPERATURE DU MELANGE (Kelvin)      !
CARGU ! TBMCI        ! TR ! <-> ! MACRO TABLEAU ENTIER MC TRAVAIL      !
CARGU ! TBMCR        ! TR ! <-> ! MACRO TABLEAU REEL   MC TRAVAIL      !
CARGU ! EH0          ! TR ! <-> ! TABLEAU REEL DE TRAVAIL              !
CARGU ! EH1          ! TR ! <-> ! TABLEAU REEL DE TRAVAIL              !
CARGU ! X1           ! TR ! <-> ! TABLEAU REEL DE TRAVAIL              !
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
      INTEGER          NCELET , NCEL
      INTEGER          NTBMCI , NTBMCR
      INTEGER          TBMCI(NCELET,NTBMCI)
C
      DOUBLE PRECISION EH(NCELET) , RTP(NCELET,*)
      DOUBLE PRECISION FUEL1(NCELET), FUEL2(NCELET), FUEL3(NCELET)
      DOUBLE PRECISION OXYD(NCELET), XINER(NCELET)
      DOUBLE PRECISION PROD1(NCELET),PROD2(NCELET)
      DOUBLE PRECISION TP(NCELET)
      DOUBLE PRECISION TBMCR(NCELET,NTBMCR)
      DOUBLE PRECISION EH0(NCELET), EH1(NCELET), X1(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          I, ICEL, ICHA , ICLA
C
      DOUBLE PRECISION YCHX10 , YCHX20 , EHCHX1 , EHCHX2
      DOUBLE PRECISION DEN1   , DEN2
      DOUBLE PRECISION XCH , XCK , XASH , X2 , X2H2
C
C***********************************************************************
C
C=======================================================================
C 1. CALCULS PRELIMINAIRES
C=======================================================================
C
C
C RQ IMPORTANTE : Pour la definition des pointeurs pour TBMC
C                 VOIR PPINI1.F
C
C --- Calcul de X1
C
      DO ICEL = 1, NCEL
        X1(ICEL) = 1.D0
      ENDDO
C
      DO ICLA = 1, NCLACP
        DO ICEL = 1, NCEL
          XCK  = RTP(ICEL,ISCA(IXCK(ICLA)))
          XCH  = RTP(ICEL,ISCA(IXCH(ICLA)))
          XASH = RTP(ICEL,ISCA(INP(ICLA)))*XMASH(ICLA)
          X2  = XCH + XCK + XASH
          X1(ICEL) = X1(ICEL) - X2
        ENDDO
      ENDDO
C
C
C=======================================================================
C 2. CALCUL DE LA TEMPERATURE DU MELANGE GAZEUX
C=======================================================================
C
      I = NPO-1
C
      DO 100 ICEL = 1, NCEL
C
C --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
C                                            et CHx2m a TH(NPO)
        EHCHX1 = ZERO
        EHCHX2 = ZERO
        YCHX10 = ZERO
        YCHX20 = ZERO
        DO ICHA = 1, NCHARB
          DEN1   = 1.D0
     &           / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
          YCHX10 = YCHX10 + DEN1 *
     &       ( TBMCR(ICEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          EHCHX1 = EHCHX1 + DEN1 *
     &       ( EHGAZE(ICHX1C(ICHA),I+1)*
     &         TBMCR(ICEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          DEN2   = 1.D0
     &           / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
          YCHX20 = YCHX20 + DEN2 *
     &       ( TBMCR(ICEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          EHCHX2 = EHCHX2 + DEN2 *
     &       ( EHGAZE(ICHX2C(ICHA),I+1)*
     &         TBMCR(ICEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
        ENDDO
        IF ( YCHX10.GT.EPZERO ) THEN
          EHCHX1 = EHCHX1 / YCHX10
        ELSE
          EHCHX1 = EHGAZE(ICHX1,I+1)
        ENDIF
        IF ( YCHX20.GT.EPZERO ) THEN
          EHCHX2 = EHCHX2 / YCHX20
        ELSE
          EHCHX2 = EHGAZE(ICHX2,I+1)
        ENDIF
C
C --- Calcul de Somme (X2.H2) a TH(NPO)
C
        X2H2 = ZERO
        DO ICLA = 1, NCLACP
          XCK  = RTP(ICEL,ISCA(IXCK(ICLA)))
          XCH  = RTP(ICEL,ISCA(IXCH(ICLA)))
          XASH = RTP(ICEL,ISCA(INP(ICLA)))*XMASH(ICLA)
          X2H2 = X2H2
     &         + XCH  * EHSOLI(ICH(ICHCOR(ICLA) ),I+1)
     &         + XCK  * EHSOLI(ICK(ICHCOR(ICLA) ),I+1)
     &         + XASH * EHSOLI(IASH(ICHCOR(ICLA)),I+1)
        ENDDO
C
C --- Clipping eventuel de TP a TH(NPO) si EH > EH1
C
        EH1(ICEL) =  (  FUEL1(ICEL)*EHCHX1
     &                + FUEL2(ICEL)*EHCHX2
     &                + FUEL3(ICEL)*EHGAZE(ICO ,I+1)
     &                + OXYD(ICEL) *EHGAZE(IO2 ,I+1)
     &                + PROD1(ICEL)*EHGAZE(ICO2,I+1)
     &                + PROD2(ICEL)*EHGAZE(IH2O,I+1)
     &                + XINER(ICEL)*EHGAZE(IN2 ,I+1) ) * X1(ICEL)
     &             + X2H2
C
        IF ( EH(ICEL).GE.EH1(ICEL) ) TP(ICEL)= TH(I+1)
 100  CONTINUE
C
C
C
      I = 1
      DO 200 ICEL = 1, NCEL
C
C --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
C                                            et CHx2m a TH(1)
        EHCHX1 = ZERO
        EHCHX2 = ZERO
        YCHX10 = ZERO
        YCHX20 = ZERO
        DO ICHA = 1, NCHARB
          DEN1   = 1.D0
     &           / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
          YCHX10 = YCHX10 + DEN1 *
     &       ( TBMCR(ICEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          EHCHX1 = EHCHX1 + DEN1 *
     &       ( EHGAZE(ICHX1C(ICHA),I)*
     &         TBMCR(ICEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          DEN2   = 1.D0
     &           / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
          YCHX20 = YCHX20 + DEN2 *
     &       ( TBMCR(ICEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          EHCHX2 = EHCHX2 + DEN2 *
     &       ( EHGAZE(ICHX2C(ICHA),I)*
     &         TBMCR(ICEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
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
C --- Calcul de Somme (X2.H2) a TH(1)
C
        X2H2 = ZERO
        DO ICLA = 1, NCLACP
          XCK  = RTP(ICEL,ISCA(IXCK(ICLA)))
          XCH  = RTP(ICEL,ISCA(IXCH(ICLA)))
          XASH = RTP(ICEL,ISCA(INP(ICLA)))*XMASH(ICLA)
          X2H2 = X2H2
     &         + XCH  * EHSOLI(ICH(ICHCOR(ICLA) ),I)
     &         + XCK  * EHSOLI(ICK(ICHCOR(ICLA) ),I)
     &         + XASH * EHSOLI(IASH(ICHCOR(ICLA)),I)
        ENDDO
C
C --- Clipping eventuel de TP a TH(1) si EH < EH0
C
        EH0(ICEL) = (  FUEL1(ICEL)*EHCHX1
     &               + FUEL2(ICEL)*EHCHX2
     &               + FUEL3(ICEL)*EHGAZE(ICO ,I)
     &               + OXYD(ICEL) *EHGAZE(IO2 ,I)
     &               + PROD1(ICEL)*EHGAZE(ICO2,I)
     &               + PROD2(ICEL)*EHGAZE(IH2O,I)
     &               + XINER(ICEL)*EHGAZE(IN2 ,I) ) * X1(ICEL)
     &             + X2H2
C
        IF ( EH(ICEL).LE.EH0(ICEL) ) TP(ICEL)= TH(I)
 200  CONTINUE
C
C
      DO I = 1, NPO-1
        DO ICEL = 1, NCEL
C
C --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
C                                            et CHx2m pour TH(I)
          EHCHX1 = ZERO
          EHCHX2 = ZERO
          YCHX10 = ZERO
          YCHX20 = ZERO
          DO ICHA = 1, NCHARB
            DEN1   = 1.D0
     &             / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
            YCHX10 = YCHX10 + DEN1 *
     &        ( TBMCR(ICEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            EHCHX1 = EHCHX1 + DEN1 *
     &        ( EHGAZE(ICHX1C(ICHA),I)*
     &          TBMCR(ICEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            DEN2   = 1.D0
     &             / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
            YCHX20 = YCHX20 + DEN2 *
     &         ( TBMCR(ICEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
            EHCHX2 = EHCHX2 + DEN2 *
     &         ( EHGAZE(ICHX2C(ICHA),I)*
     &           TBMCR(ICEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
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
C --- Calcul de Somme (X2.H2) pour TH(I)
C
          X2H2 = ZERO
          DO ICLA = 1, NCLACP
            XCK  = RTP(ICEL,ISCA(IXCK(ICLA)))
            XCH  = RTP(ICEL,ISCA(IXCH(ICLA)))
            XASH = RTP(ICEL,ISCA(INP(ICLA)))*XMASH(ICLA)
            X2H2 = X2H2
     &           + XCH  * EHSOLI(ICH(ICHCOR(ICLA) ),I  )
     &           + XCK  * EHSOLI(ICK(ICHCOR(ICLA) ),I  )
     &           + XASH * EHSOLI(IASH(ICHCOR(ICLA)),I  )
          ENDDO
C
          EH0(ICEL) = (  FUEL1(ICEL)*EHCHX1
     &                 + FUEL2(ICEL)*EHCHX2
     &                 + FUEL3(ICEL)*EHGAZE(ICO ,I  )
     &                 + OXYD(ICEL) *EHGAZE(IO2 ,I  )
     &                 + PROD1(ICEL)*EHGAZE(ICO2,I  )
     &                 + PROD2(ICEL)*EHGAZE(IH2O,I  )
     &                 + XINER(ICEL)*EHGAZE(IN2 ,I  ) ) * X1(ICEL)
     &               + X2H2
C
C --- Calcul de l'enthalpie de l'espece gazeuse CHx1m
C                                            et CHx2m pour TH(I+1)
          EHCHX1 = ZERO
          EHCHX2 = ZERO
          YCHX10 = ZERO
          YCHX20 = ZERO
          DO ICHA = 1, NCHARB
            DEN1   = 1.D0
     &             / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
            YCHX10 = YCHX10 + DEN1 *
     &        ( TBMCR(ICEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            EHCHX1 = EHCHX1 + DEN1 *
     &        ( EHGAZE(ICHX1C(ICHA),I+1)*
     &          TBMCR(ICEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
            DEN2   = 1.D0
     &             / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
            YCHX20 = YCHX20 + DEN2 *
     &         ( TBMCR(ICEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
            EHCHX2 = EHCHX2 + DEN2 *
     &         ( EHGAZE(ICHX2C(ICHA),I+1)*
     &           TBMCR(ICEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          ENDDO
          IF ( YCHX10.GT.EPZERO ) THEN
            EHCHX1 = EHCHX1 / YCHX10
          ELSE
            EHCHX1 = EHGAZE(ICHX1,I+1)
          ENDIF
          IF ( YCHX20.GT.EPZERO ) THEN
            EHCHX2 = EHCHX2 / YCHX20
          ELSE
            EHCHX2 = EHGAZE(ICHX2,I+1)
          ENDIF
C
C --- Calcul de Somme (X2.H2) pour TH(I+1)
C
          X2H2 = ZERO
          DO ICLA = 1, NCLACP
            XCK  = RTP(ICEL,ISCA(IXCK(ICLA)))
            XCH  = RTP(ICEL,ISCA(IXCH(ICLA)))
            XASH = RTP(ICEL,ISCA(INP(ICLA)))*XMASH(ICLA)
            X2H2 = X2H2
     &           + XCH  * EHSOLI(ICH(ICHCOR(ICLA) ),I+1)
     &           + XCK  * EHSOLI(ICK(ICHCOR(ICLA) ),I+1)
     &           + XASH * EHSOLI(IASH(ICHCOR(ICLA)),I+1)
          ENDDO
C
          EH1(ICEL) = (  FUEL1(ICEL)*EHCHX1
     &                 + FUEL2(ICEL)*EHCHX2
     &                 + FUEL3(ICEL)*EHGAZE(ICO ,I+1)
     &                 + OXYD(ICEL) *EHGAZE(IO2 ,I+1)
     &                 + PROD1(ICEL)*EHGAZE(ICO2,I+1)
     &                 + PROD2(ICEL)*EHGAZE(IH2O,I+1)
     &                 + XINER(ICEL)*EHGAZE(IN2 ,I+1) ) * X1(ICEL)
     &               + X2H2
C
          IF ( EH(ICEL).GE.EH0(ICEL) .AND. EH(ICEL).LE.EH1(ICEL) ) THEN
            TP(ICEL)= TH(I) + (EH(ICEL)-EH0(ICEL)) *
     &                        (TH(I+1)-TH(I))/(EH1(ICEL)-EH0(ICEL))
          ENDIF
C
        ENDDO
      ENDDO
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
