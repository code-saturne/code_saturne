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
                        SUBROUTINE CPCYM2
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NRTBMC ,
     &   INDPDF ,
     &   RTP    ,
     &   F1M    , F2M    , F3M   , F4M   , F5M  , F6M , F7M ,
     &   F3MAX  ,
     &   F1CL   , F2CL   , F3CL  , F4CL  , F5CL , F6CL , F7CL ,
     &   F4M1   , F4M2   , D4CL  , D4F4  , HREC ,
     &   RTBMC  , W1     ,
     &   FUEL1  , FUEL2  , FUEL3 , OXYD  , PROD1 , PROD2  ,
     &   XINER  )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CALCUL DES CONCENTRATIONS GAZEUSES MOYENNES
CFONC   VALEURS CELLULES
CFONC   ----------------
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! FVAP         ! TR !  -> ! MOYENNE DU TRACEUR 1 MVl [FOVm+CO]   !
CARGU ! FHET         ! TR !  -> ! MOYENNE DU TRACEUR 3 C héterogène    !
CARGU ! F4P2M        ! TR !  -> ! VARIANCE DU TRACEUR 4 (air)          !
CARGU ! INDPDF       ! TE !  -> ! PASSAGE PAR LES PDF                  !
CARGU ! F1M          ! TR !  -> ! MOYENNE DU TRACEUR 1 MVl [CHx1m+CO]  !
CARGU ! F2M          ! TR !  -> ! MOYENNE DU TRACEUR 2 MVL [CHx2m+CO]  !
CARGU ! F3M          ! TR !  -> ! MOYENNE DU TRACEUR 3 (CO C.HET)      !
CARGU ! F4M          ! TR !  -> ! MOYENNE DU TRACEUR 4 (AIR)           !
CARGU ! F5M          ! TR !  -> ! MOYENNE DU TRACEUR 5 (H2O)           !
CARGU ! FUEL1        ! TR !  <- ! FRACTION MASSIQUE CHx1m              !
CARGU ! FUEL2        ! TR !  <- ! FRACTION MASSIQUE CHx2m              !
CARGU ! FUEL3        ! TR !  <- ! FRACTION MASSIQUE CO                 !
CARGU ! OXYD         ! TR !  <- ! FRACTION MASSIQUE O2                 !
CARGU ! PROD1        ! TR !  <- ! FRACTION MASSIQUE CO2                !
CARGU ! PROD2        ! TR !  <- ! FRACTION MASSIQUE H2O                !
CARGU ! XINER        ! TR !  <- ! FRACTION MASSIQUE N2                 !
CARGU ! F4M1         ! TR ! <-> ! BORNE MINIMUM                        !
CARGU ! F4M2         ! TR ! <-> ! BORNE MAX                            !
CARGU ! D4CL         ! TR ! <-> ! AMPLITUDE DU PIC DE DIRAC EN F4CL    !
CARGU ! D4F4         ! TR ! <-> ! AMPLITUDE DU PIC DE DIRAC EN 1       !
CARGU !              !    !     ! (AIR PUR)                            !
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
      INCLUDE "parall.h"
      INCLUDE "pointe.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
      INCLUDE "ppcpfu.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL , NRTBMC
      INTEGER          INDPDF(NCELET)
C
      DOUBLE PRECISION RTP(NCELET,*)
      DOUBLE PRECISION F1M(NCELET)  , F2M(NCELET)
      DOUBLE PRECISION F3M(NCELET)  , F4M(NCELET) , F5M(NCELET)
      DOUBLE PRECISION F6M(NCELET)  , F7M(NCELET) , F3MAX(NCELET)
C
      DOUBLE PRECISION F1CL(NCELET) , F2CL(NCELET)
      DOUBLE PRECISION F3CL(NCELET) , F4CL(NCELET), F5CL(NCELET)
      DOUBLE PRECISION F6CL(NCELET) , F7CL(NCELET)
C
      DOUBLE PRECISION F4M1(NCELET) , F4M2(NCELET) , D4CL(NCELET)
      DOUBLE PRECISION D4F4(NCELET) , HREC(NCELET)
C
      DOUBLE PRECISION FUEL1(NCELET), FUEL2(NCELET) , FUEL3(NCELET)
      DOUBLE PRECISION OXYD(NCELET) , PROD1(NCELET), PROD2(NCELET)
      DOUBLE PRECISION XINER(NCELET)
C
      DOUBLE PRECISION RTBMC(NCELET,NRTBMC)
      DOUBLE PRECISION W1(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL , ICHA , II
C
      INTEGER          N1 , N2 , N3 , N4 , N5 , N6 , N7 , N8 , N9
      INTEGER          N10
      INTEGER          ITBR , ICLA
      INTEGER          NCLO2,NCLCO,NCLO21,NCLCO1
      INTEGER          NBRINT
      PARAMETER       (NBRINT = 10)
      INTEGER          INTTMP(NBRINT)
C
      DOUBLE PRECISION ANU1 , ANU2 , ANU3
      DOUBLE PRECISION SOMM
C
      DOUBLE PRECISION ZCHX10 , ZCHX20 , ZCO0   , ZCO20
      DOUBLE PRECISION ZN20   , ZH2O0  , ZO20
      DOUBLE PRECISION YCHX10 , YCHX20
      DOUBLE PRECISION ZCHX11 , ZCHX21 , ZCO1   , ZCO21
      DOUBLE PRECISION ZN21   , ZH2O1  , ZO21
      DOUBLE PRECISION ZCHX12 , ZCHX22 , ZCO2   , ZCO22
      DOUBLE PRECISION ZN22   , ZH2O2  , ZO22
      DOUBLE PRECISION ZCHX13 , ZCHX23 , ZCO3   , ZCO23
      DOUBLE PRECISION ZN23   , ZH2O3  , ZO23
      DOUBLE PRECISION REAC1  , REAC2  , REAC3
      DOUBLE PRECISION X2     , XCH    , XCK    , XASH , XWAT
C
      DOUBLE PRECISION ZCOF10 , ZCOF20 , ZH2O10 , ZH2O20
      DOUBLE PRECISION CX1M   , CX2M
      DOUBLE PRECISION WMCHX1 , WMCHX2 , WMCO   , WMCO2 , WMH2O
      DOUBLE PRECISION WMO2   , WMN2   , WMC
      DOUBLE PRECISION ACHX1F1, ACOF1  , ACHX2F2, ACOF2
      DOUBLE PRECISION YO2OX
      DOUBLE PRECISION AH2OF1 , AH2OF2
C
C     Variables pour le support de la Pdf et sa description
C
      DOUBLE PRECISION F42M
      DOUBLE PRECISION ZZCL(7),ZZS1(7),ZZS2(7),ZZS3(7)
      DOUBLE PRECISION ZZ(7)  ,ZZOXY(7)
C
      DOUBLE PRECISION F4S1 , F4S2 , F4S3
      DOUBLE PRECISION BB1  , BB2  , DEN1 , DEN2
      DOUBLE PRECISION FP1 , FP2 , FP3 , FP4 , FP5 , ZCO2T , ZCOT
C
      DOUBLE PRECISION F4LOC,F6LOC,F7LOC,S4LOC,S6LOC,S7LOC
      DOUBLE PRECISION FOXY ,FOXYCL
C
      DOUBLE PRECISION SOMMIN,SOMMAX,AO2MIN,AO2MAX
C
      DOUBLE PRECISION CMAXS1,CMINS1,CMAXA1,CMINA1
      DOUBLE PRECISION CMAXS2,CMINS2,CMAXA2,CMINA2
C
C***********************************************************************
C
C=======================================================================
C 0. INITIALISATIONS
C=======================================================================
C
C --- Initialisation tableau de travail
C
      DO IEL = 1, NCEL
        W1(IEL) = ZERO
        DO ITBR = 1, NRTBMC
          RTBMC(IEL,ITBR) = ZERO
        ENDDO
      ENDDO
C
C=======================================================================
C 1. CALCULS PRELIMINAIRES
C=======================================================================
C
C --- Masses molaires
C
      WMCO   = WMOLE(ICO  )
      WMO2   = WMOLE(IO2  )
      WMCO2  = WMOLE(ICO2 )
      WMH2O  = WMOLE(IH2O )
      WMN2   = WMOLE(IN2)
      WMC    = WMOLAT(IATC)
C
C --- Calcul de F1M de chaque charbon dans RTBMC
C            de F2M de chaque charbon dans RTBMC
C
C ------ W1 = - Somme des X2(icla)
C
      DO ICLA = 1, NCLACP
        DO IEL = 1, NCEL
          XCK  = RTP(IEL,ISCA(IXCK(ICLA)))
          XCH  = RTP(IEL,ISCA(IXCH(ICLA)))
          XASH = RTP(IEL,ISCA(INP (ICLA)))*XMASH(ICLA)
          IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
            XWAT = RTP(IEL,ISCA(IXWT(ICLA)))
          ELSE
            XWAT = 0.D0
          ENDIF
          X2   = XCH + XCK + XASH + XWAT
          W1(IEL) =  W1(IEL) - X2
        ENDDO
      ENDDO
C
C ------ RTBMC(IF1MC(ICHA)) = F1M(ICHA)
C        RTBMC(IF2MC(ICHA)) = F2M(ICHA)
C
      DO ICHA = 1, NCHARB
        DO IEL = 1, NCEL
          RTBMC(IEL,IF1MC(ICHA)) = RTP(IEL,ISCA(IF1M(ICHA)))
     &                           / (1.D0+W1(IEL))
          RTBMC(IEL,IF2MC(ICHA)) = RTP(IEL,ISCA(IF2M(ICHA)))
     &                           / (1.D0+W1(IEL))
        ENDDO
      ENDDO
C
C
C --- Calcul de CX1M, CX2M
C     et des coefficient relatifs a l'expression de
C        ZCHx1m0, ZCHx2m0, ZCO0, ZO20 et ZN20
C        introduction de ZH2O0
C        fonctions de F1, F2, F3 et F4
C
      DO IEL = 1, NCEL
C
C ---- YCHX1,20 en kg/kg de melange et
C      ZCHX1,20 en nb de moles/ kg de melange
C
        YCHX10 = ZERO
        YCHX20 = ZERO
        ZCHX10 = ZERO
        ZCHX20 = ZERO
        ZCOF10 = ZERO
        ZCOF20 = ZERO
        ZH2O10 = ZERO
        ZH2O20 = ZERO
        DO ICHA = 1, NCHARB
          DEN1   = 1.D0
     &           / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))
     &              +B1(ICHA)*WMOLE(ICO)
     &              +C1(ICHA)*WMOLE(IH2O) )
          YCHX10 = YCHX10 + DEN1 *
     &           ( RTBMC(IEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          ZCHX10 = ZCHX10 + DEN1 *
     &           ( RTBMC(IEL,IF1MC(ICHA))*A1(ICHA) )
          ZCOF10 = ZCOF10 + DEN1 *
     &           ( RTBMC(IEL,IF1MC(ICHA))*B1(ICHA) )
          ZH2O10 = ZH2O10 + DEN1 *
     &           ( RTBMC(IEL,IF1MC(ICHA))*C1(ICHA) )
          DEN2   = 1.D0
     &           / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))
     &              +B2(ICHA)*WMOLE(ICO)
     &              +C2(ICHA)*WMOLE(IH2O) )
          YCHX20 = YCHX20 + DEN2 *
     &           ( RTBMC(IEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          ZCHX20 = ZCHX20 + DEN2 *
     &           ( RTBMC(IEL,IF2MC(ICHA))*A2(ICHA) )
          ZCOF20 = ZCOF20 + DEN2 *
     &           ( RTBMC(IEL,IF2MC(ICHA))*B2(ICHA) )
          ZH2O20 = ZH2O20 + DEN2 *
     &           ( RTBMC(IEL,IF2MC(ICHA))*C2(ICHA) )
        ENDDO
        RTBMC(IEL,IX1MC) = 4.D0
        IF ( ZCHX10.GT.EPZERO )
     &    RTBMC(IEL,IX1MC) = ( YCHX10/ZCHX10-WMOLAT(IATC) )
     &                     / WMOLAT(IATH)
        RTBMC(IEL,IX2MC) = 2.D0
        IF ( ZCHX20.GT.EPZERO )
     &    RTBMC(IEL,IX2MC) = ( YCHX20/ZCHX20-WMOLAT(IATC) )
     &                     / WMOLAT(IATH)
C
        RTBMC(IEL,ICHX1F1) = 0.D0
        RTBMC(IEL,ICOF1  ) = 0.D0
        RTBMC(IEL,IH2OF1 ) = 0.D0
        IF ( F1M(IEL).GT.EPZERO ) THEN
          RTBMC(IEL,ICHX1F1) = ZCHX10 / F1M(IEL)
          RTBMC(IEL,ICOF1  ) = ZCOF10 / F1M(IEL)
          RTBMC(IEL,IH2OF1 ) = ZH2O10 / F1M(IEL)
        ENDIF
        RTBMC(IEL,ICHX2F2) = 0.D0
        RTBMC(IEL,ICOF2  ) = 0.D0
        RTBMC(IEL,IH2OF2 ) = 0.D0
        IF ( F2M(IEL).GT.EPZERO ) THEN
          RTBMC(IEL,ICHX2F2) = ZCHX20 / F2M(IEL)
          RTBMC(IEL,ICOF2  ) = ZCOF20 / F2M(IEL)
          RTBMC(IEL,IH2OF2 ) = ZH2O20 / F2M(IEL)
        ENDIF
C
      ENDDO
C
C ---- Initialisation des fractions massiques avec de l'air
C
      DO IEL = 1, NCEL
        FUEL1(IEL)  = ZERO
        FUEL2(IEL)  = ZERO
        FUEL3(IEL)  = ZERO
        OXYD(IEL)   = WMO2*AO2F4
        PROD1(IEL)  = ZERO
        PROD2(IEL)  = ZERO
        XINER(IEL)  = WMN2*AN2F4
      ENDDO
C
      NCLO2  = 0
      NCLCO  = 0
      NCLO21 = 0
      NCLCO1 = 0
C
      CMINA1 =  1.D+20
      CMAXA1 = -1.D+20
      CMINA2 =  1.D+20
      CMAXA2 = -1.D+20
      CMINS1 =  1.D+20
      CMAXS1 = -1.D+20
      CMINS2 =  1.D+20
      CMAXS2 = -1.D+20
C
C=======================================================================
C 3. CALCUL DE LA COMPOSITION DU MELANGE SANS LES PDF
C    SI LES FLUCTUATIONS DE S SONT TROP FAIBLES
C=======================================================================
C
      DO IEL = 1, NCEL
C
        IF ( INDPDF(IEL).EQ.0 ) THEN
C
C --> Calculs preliminaires
C
          IF ( IEQCO2 .EQ. 1 ) THEN
            ZCO2T = (RTP(IEL,ISCA(IYCO2))/ (1.D0+W1(IEL)))/WMCO2
          ELSE IF ( IEQCO2 .EQ. 2 ) THEN
            ZCOT  = (RTP(IEL,ISCA(IYCO2))/ (1.D0+W1(IEL)))/WMCO
          ENDIF
C
          CX1M    = RTBMC(IEL,IX1MC)
          CX2M    = RTBMC(IEL,IX2MC)
          WMCHX1  = WMOLAT(IATC) + CX1M*WMOLAT(IATH)
          WMCHX2  = WMOLAT(IATC) + CX2M*WMOLAT(IATH)
          ACHX1F1 = RTBMC(IEL,ICHX1F1)
          ACHX2F2 = RTBMC(IEL,ICHX2F2)
          ACOF1   = RTBMC(IEL,ICOF1)
          ACOF2   = RTBMC(IEL,ICOF2)
          AH2OF1  = RTBMC(IEL,IH2OF1)
          AH2OF2  = RTBMC(IEL,IH2OF2)
C
C --> Composition de la phase gazeuse avant combustion (indice 0)
C
C          ZCHX10 = ACHX1F1*F1
C          ZCHX20 = ACHX2F2*F2
C          ZCO0   = ACOF1*F1 + ACOF2*F2 + ACOF3*F3
C          ZCO20  = ZERO
C          ZH2O0  = ZERO
C          ZO20   = AO2F4*F4 + AO2F3*F3
C          ZN20   = AN2F4*F4
C
          ZCHX10  = ACHX1F1*F1M(IEL)
          ZCHX20  = ACHX2F2*F2M(IEL)
          ZCO0    = (ACOF1*F1M(IEL)+ACOF2*F2M(IEL) +
     &               ACOF3*F3M(IEL))
          ZO20    = AO2F4*F4M(IEL)+AO2F6*F6M(IEL)+AO2F7*F7M(IEL)
     &             +AO2F3*F3M(IEL)
          ZN20    = AN2F4*F4M(IEL)+AN2F6*F6M(IEL)+AN2F7*F7M(IEL)
          ZCO20   = ACO2F4*F4M(IEL)+ACO2F6*F6M(IEL)+ACO2F7*F7M(IEL)
          ZH2O0   = AH2OF4*F4M(IEL)+AH2OF6*F6M(IEL)+AH2OF7*F7M(IEL)
     &             +AH2OF5*F5M(IEL)+AH2OF1*F1M(IEL)+AH2OF2*F2M(IEL)
C
C --> Calcul de la composition du melange
C
          ANU1   = 0.25D0*(CX1M-CX2M)
          REAC1  = MIN(ZCHX10, (ZO20/ANU1))
          ZCHX11 = ZCHX10 -           REAC1
          ZO21   = ZO20   -      ANU1*REAC1
          ZCHX21 = ZCHX20 +           REAC1
          ZCO1   = ZCO0
          ZCO21  = ZCO20
          ZH2O1  = ZH2O0  + 2.D0*ANU1*REAC1
          ZN21   = ZN20
C
          ANU2   = 0.25D0*(2.D0 + CX2M)
          REAC2  = MIN(ZCHX21, (ZO21/ANU2))
          ZCHX12 = ZCHX11
          ZCHX22 = ZCHX21 -           REAC2
          ZO22   = ZO21   -      ANU2*REAC2
          ZCO2   = ZCO1   +           REAC2
          ZCO22  = ZCO21
          ZH2O2  = ZH2O1  + .5D0*CX2M*REAC2
          ZN22   = ZN21
C
          ANU3   = 0.5D0
          IF ( IEQCO2 .EQ. 0 ) THEN
            REAC3  = MIN(ZCO2, (ZO22/ANU3))
            ZCHX13 = ZCHX12
            ZCHX23 = ZCHX22
            ZO23   = ZO22  - ANU3*REAC3
            ZCO3   = ZCO2  -      REAC3
            ZCO23  = ZCO22 +      REAC3
            ZH2O3  = ZH2O2
            ZN23   = ZN22
C
          ELSE IF ( IEQCO2 .EQ. 1 ) THEN
C
C      Transport de CO2
C
            ZCHX13 = ZCHX12
            ZCHX23 = ZCHX22
C
            IF ( (ZCO2-(ZCO2T-ZCO22)) .LT. -EPZERO ) THEN
              NCLCO = NCLCO +1
              CMINS1 = MIN((ZCO2-ZCO2T),CMINS1)
              CMAXS1 = MAX((ZCO2-ZCO2T),CMINS1)
            ENDIF
            IF ( (ZO22-ANU3*(ZCO2T-ZCO22)) .LT. -EPZERO ) THEN
              NCLO2 = NCLO2 +1
              CMINS2 = MIN((ZO22-ANU3*ZCO2T),CMINS2)
              CMAXS2 = MAX((ZO22-ANU3*ZCO2T),CMINS2)
            ENDIF
            REAC3 = MIN(MAX(ZCO2T-ZCO22,0.D0),ZCO2,ZO22/ANU3)
C
C            RTP(IEL,ISCA(IYCO2))= ZCO2T*(1.D0+W1(IEL))*WMCO2
C
            ZO23   = ZO22  - ANU3*REAC3
            ZCO3   = ZCO2  -      REAC3
            ZCO23  = ZCO22 +      REAC3
C
            ZH2O3  = ZH2O2
            ZN23   = ZN22
C
          ELSE
C
            WRITE(NFECRA,3000) IEQCO2
            CALL CSEXIT(1)
C
          ENDIF
C
          FUEL1(IEL) = ZCHX13 * WMCHX1
          FUEL2(IEL) = ZCHX23 * WMCHX2
          FUEL3(IEL) = ZCO3   * WMCO
          PROD1(IEL) = ZCO23  * WMCO2
          PROD2(IEL) = ZH2O3  * WMH2O
          OXYD(IEL)  = ZO23   * WMO2
          XINER(IEL) = ZN23   * WMN2
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 4. CALCUL DE LA COMPOSITION DU MELANGE AVEC LES PDF
C=======================================================================
C
      DO IEL = 1, NCEL
C
        IF ( INDPDF(IEL).EQ.1 ) THEN
C
C --> Calculs preliminaires
C
          IF ( IEQCO2 .EQ. 1 ) THEN
            ZCO2T = (RTP(IEL,ISCA(IYCO2))/(1.D0+W1(IEL)))/WMCO2
          ELSE IF ( IEQCO2 .EQ. 2 ) THEN
            ZCOT  = (RTP(IEL,ISCA(IYCO2))/(1.D0+W1(IEL)))/WMCO
          ENDIF
C
          CX1M    = RTBMC(IEL,IX1MC)
          CX2M    = RTBMC(IEL,IX2MC)
          WMCHX1  = WMOLAT(IATC) + CX1M*WMOLAT(IATH)
          WMCHX2  = WMOLAT(IATC) + CX2M*WMOLAT(IATH)
          ACHX1F1 = RTBMC(IEL,ICHX1F1)
          ACHX2F2 = RTBMC(IEL,ICHX2F2)
          ACOF1   = RTBMC(IEL,ICOF1)
          ACOF2   = RTBMC(IEL,ICOF2)
          AH2OF1  = RTBMC(IEL,IH2OF1)
          AH2OF2  = RTBMC(IEL,IH2OF2)
C
          F4LOC = F4M(IEL)/(F4M(IEL)+F6M(IEL)+F7M(IEL))
          F6LOC = F6M(IEL)/(F4M(IEL)+F6M(IEL)+F7M(IEL))
          F7LOC = F7M(IEL)/(F4M(IEL)+F6M(IEL)+F7M(IEL))
          S4LOC = OXYO2 (1)*WMO2 +OXYN2 (1)*WMN2
     &           +OXYH2O(1)*WMH2O+OXYCO2(1)*WMCO2
          S6LOC = OXYO2 (2)*WMO2 +OXYN2 (2)*WMN2
     &           +OXYH2O(2)*WMH2O+OXYCO2(2)*WMCO2
          S7LOC = OXYO2 (3)*WMO2 +OXYN2 (3)*WMN2
     &           +OXYH2O(3)*WMH2O+OXYCO2(3)*WMCO2
C
          YO2OX  =  WMO2
     &            /(F4LOC*S4LOC+F6LOC*S6LOC+F7LOC*S7LOC)
          AN2F3  = (1.D0-YO2OX)/WMOLE(IN2) * (1.D0-F3MAX(IEL))
C
          ANU1   = 0.25D0*(CX1M-CX2M)
          ANU2   = 0.25D0*(2.D0 + CX2M)
          ANU3   = 0.5D0
C
          FOXY   = F4M(IEL)+F5M(IEL)+F6M(IEL)+F7M(IEL)
          FOXYCL = F4CL(IEL)+F5CL(IEL)+F6CL(IEL)+F7CL(IEL)
C
C         Concentrations aux extrémités
C
          DO II = 1,7
            ZZCL(II)  = ZERO
            ZZOXY(II) = ZERO
          ENDDO
C
          ZZCL(ICHX1)= ACHX1F1*F1CL(IEL)
          ZZCL(ICHX2)= ACHX2F2*F2CL(IEL)
          ZZCL(ICO)  = ACOF1 *F1CL(IEL)+ACOF2*F2CL(IEL)
     &                                 +ACOF3*F3CL(IEL)
          ZZCL(IO2)  = AO2F4*F4CL(IEL)+AO2F6*F6CL(IEL)+AO2F7*F7CL(IEL)
     &                +AO2F3*F3CL(IEL)
          ZZCL(IN2)  = AN2F4*F4CL(IEL)+AN2F6*F6CL(IEL)+AN2F7*F7CL(IEL)
          ZZCL(IH2O) = AH2OF4*F4CL(IEL)+AH2OF6*F6CL(IEL)
     &                +AH2OF7*F7CL(IEL)
     &                +AH2OF1*F1CL(IEL)+AH2OF2*F2CL(IEL)
     &                +AH2OF5*F5CL(IEL)
          ZZCL(ICO2) = ACO2F4*F4CL(IEL)+ACO2F6*F6CL(IEL)
     &                +ACO2F7*F7CL(IEL)
C
          ZZOXY(IO2)  = (AO2F4*F4M(IEL)+AO2F6*F6M(IEL)+AO2F7*F7M(IEL))
     &                 /FOXY
          ZZOXY(IN2)  = (AN2F4*F4M(IEL)+AN2F6*F6M(IEL)+AN2F7*F7M(IEL))
     &                 /FOXY
          ZZOXY(IH2O) = ( AH2OF4*F4M(IEL)+AH2OF6*F6M(IEL)
     &                   +AH2OF7*F7M(IEL)
     &                   +AH2OF5*F5M(IEL) )
     &                 /FOXY
          ZZOXY(ICO2) = (ACO2F4*F4M(IEL)+ACO2F6*F6M(IEL)
     &                                  +ACO2F7*F7M(IEL))
     &                 /FOXY
C
C         Calcul de la stoechiométrie de la première réaction
C
          IF ( ZZCL(ICHX1)  .GT. 0.D0 .OR.
     &         ZZOXY(IO2)   .GT. 0.D0         ) THEN
            F4S1 = ( ANU1 * ZZCL(ICHX1) + FOXYCL * ZZOXY(IO2) )
     &           / ( ANU1 * ZZCL(ICHX1) +          ZZOXY(IO2) )
          ELSE
            F4S1 = 1.D0
          ENDIF
C
C         Calcul des concentrations au premier point stoechio
C         avant réaction
C
          DO II = 1,7
            ZZS1(II) = ZZCL(II)
     &                + (F4S1-FOXYCL)/(1.D0-FOXYCL)
     &           *(ZZOXY(II)-ZZCL(II))
          ENDDO
C
          ZZS1(ICHX2) = ZZS1(ICHX2) + ZZS1(ICHX1)
          ZZS1(IO2)   = ZERO
          ZZS1(IH2O)  = ZZS1(IH2O) + 2.D0*ANU1*ZZS1(ICHX1)
          ZZS1(ICHX1) = ZERO
C
C         Calcul de la stoechiometrie de la deuxième réaction
C
          IF ( ZZS1(ICHX2) .GT. 0.D0 .OR.
     &         ZZOXY(IO2)  .GT. 0.D0         ) THEN
            F4S2 = ( ANU2 * ZZS1(ICHX2) + F4S1 * ZZOXY(IO2) )
     &            /( ANU2 * ZZS1(ICHX2) +        ZZOXY(IO2) )
          ELSE
            F4S2 = 1.D0
          ENDIF
C
C         Calcul des concentrations au deuxième point stoechio
C         avant réaction
C
          IF ( ABS(1.D0-F4S1).GT. 0.D0 ) THEN
            DO II = 1,7
               ZZS2(II) = ZZS1(II)
     &           + (F4S2-F4S1) / (1.D0-F4S1) * (ZZOXY(II) - ZZS1(II))
            ENDDO
          ELSE
            DO II = 1,7
               ZZS2(II) = ZZS1(II)
            ENDDO
          ENDIF
C
          ZZS2(ICO)   = ZZS2(ICO)  + ZZS2(ICHX2)
          ZZS2(IH2O)  = ZZS2(IH2O) + (2.D0*ANU2-1.D0)*ZZS2(ICHX2)
          ZZS2(IO2)   = ZERO
          ZZS2(ICHX2) = ZERO
C
C         Calcul de la stoechiometrie de la troisième réaction
C         CO + O2 => CO2
C         Les concentrations sont désormais linéaires entre F4S2 et 1
C         ZZ(F4,II)  = ZZS2(II) + (F4-F4S2)/(1-F4S2)*(ZZF4(II)-ZZS2(II))
C         On cherche le point tel que
C         ZZ(f4s3,ICO) = ZZ(F4S3,IO2)/ANU3
C         La formule est semblable à la précedente en substituant
C         ANU2 par ANU3
C
          IF ( ZZS2(ICO) .GT. 0.D0 .OR.
     &         ZZS2(IO2) .GT. 0.D0         ) THEN
            F4S3 = ( ANU3 * ZZS2(ICO ) + F4S2 * ZZOXY(IO2) )
     &           / ( ANU3 * ZZS2(ICO ) +        ZZOXY(IO2) )
          ELSE
            F4S3 = 1.D0
          ENDIF
C
C         Calcul des concentrations au troisième point stoechio
C         avant réaction
C
          IF ( ABS(1.D0-F4S2).GT. 0.D0 ) THEN
            DO II = 1,7
               ZZS3(II) = ZZS2(II)
     &           + (F4S3-F4S2) / (1.D0-F4S2) * (ZZOXY(II) - ZZS2(II))
            ENDDO
          ELSE
            DO II = 1,7
               ZZS3(II) = ZZS2(II)
            ENDDO
          ENDIF
C
C         Le nombre de moles de réaction est le nombre de moles de CO
C
          IF ( IEQCO2 .EQ. 0 ) THEN
            ZZS3(ICO2) = ZZS3(ICO2) + ZZS3(ICO)
            ZZS3(IO2)  = ZERO
            ZZS3(ICO)  = ZERO
          ENDIF
C
C        Désormais on connait les concentrations en 5 points
C        cl,s1,s2,s3,f4
C        les concentrations intermédiaires sont linéaires par morceaux
C        et les paramêtres de la pdf D4cl,D4f4, F4m1,F4m2,HREC
C
C        On initialise par les concentrations aux extrémités
         DO II = 1,7
           ZZ(II) = D4CL(IEL)*ZZCL(II)+D4F4(IEL)*ZZOXY(II)
         ENDDO
C
C        Intégration sur le premier intervalle de richesse (entre cl et s1)
C
         BB1 = MAX(F4M1(IEL),FOXYCL)
         BB2 = MIN(F4M2(IEL),F4S1)
         IF( BB2.GT.BB1 ) THEN
           DO II = 1,7
             ZZ(II) =  ZZ(II) + HREC(IEL)*(BB2-BB1)/(F4S1-FOXYCL)
     &               *( (ZZCL(II)*F4S1-ZZS1(II)*FOXYCL)
     &              +   (ZZS1(II)-ZZCL(II))*(BB1+BB2)*0.5D0)
           ENDDO
         ENDIF
C
C        Intégration sur le deuxième intervalle de (entre s1 et s2)
C
         BB1 = MAX(F4M1(IEL),F4S1)
         BB2 = MIN(F4M2(IEL),F4S2)
         IF( BB2.GT.BB1 ) THEN
           DO II = 1,7
             ZZ(II) = ZZ(II) + HREC(IEL)*(BB2-BB1)/(F4S2-F4S1)
     &             * ( (ZZS1(II)*F4S2-ZZS2(II)*F4S1)
     &           +   (ZZS2(II)-ZZS1(II))*(BB1+BB2)*0.5D0)
           ENDDO
         ENDIF
C
C        Intégration de s2 à s3
C
         BB1 = Max(F4M1(IEL),F4S2)
         BB2 = Min(F4M2(IEL),F4S3)
         IF( BB2.GT.BB1 ) THEN
           DO II = 1,7
             ZZ(II) = ZZ(II) + HREC(IEL)*(BB2-BB1)/(F4S3-F4S2)
     &           * ( (ZZS2(II)*F4S3-ZZS3(II)*F4S2)
     &           +   (ZZS3(II)-ZZS2(II))*(BB1+BB2)*0.5D0)
           ENDDO
         ENDIF
C
C       Intégration de s3 à f4
C
         BB1 = Max(F4M1(IEL),F4S3)
         BB2 = Min(F4M2(IEL),1.D0)
         IF ( BB2.GT.BB1 ) THEN
           DO II = 1,7
            ZZ(II) = ZZ(II) + HREC(IEL)*(BB2-BB1)/(1.D0-F4S3)
     &             * ( (ZZS3(II)*1.D0-ZZOXY(II)*F4S3)
     &                +(ZZOXY(II)-ZZS3(II))*(BB1+BB2)*0.5D0)
           ENDDO
         ENDIF
C
         IF ( IEQCO2 .EQ. 1 ) THEN
C
C  Transport de CO2
C
           IF ( (ZZ(ICO) - (ZCO2T-ZZ(ICO2))) .LT. -EPZERO ) THEN
             NCLCO1 = NCLCO1 +1
             CMINA1 = MIN((ZZ(ICO)- (ZCO2T-ZZ(ICO2))),CMINA1)
             CMAXA1 = MAX((ZZ(ICO)- (ZCO2T-ZZ(ICO2))),CMAXA1)
           ENDIF
           IF ( (ZZ(IO2)-ANU3*(ZCO2T-ZZ(ICO2))) .LT. -EPZERO ) THEN
             NCLO21 = NCLO21 +1
             CMINA2 = MIN((ZZ(IO2) -ANU3*(ZCO2T-ZZ(ICO2))),CMINA2)
             CMAXA2 = MAX((ZZ(IO2) -ANU3*(ZCO2T-ZZ(ICO2))),CMAXA2)
           ENDIF

           REAC3 = MIN(MAX(ZCO2T-ZZ(ICO2),0.D0),ZZ(ICO),
     &                                          ZZ(IO2)/ANU3)
C
           ZZ(ICO ) = ZZ(ICO) -      REAC3
           ZZ(IO2 ) = ZZ(IO2) - ANU3*REAC3
           ZZ(ICO2) = ZZ(ICO2)+      REAC3
C
C           RTP(IEL,ISCA(IYCO2))= ZCO2T*(1.D0+W1(IEL))*WMCO2
C
         ELSE IF ( IEQCO2 .NE. 0 ) THEN
C
           WRITE(NFECRA,3000) IEQCO2
           CALL CSEXIT(1)
C
         ENDIF
C
         FUEL1(IEL) = ZZ(ICHX1) * WMCHX1
         FUEL2(IEL) = ZZ(ICHX2) * WMCHX2
         FUEL3(IEL) = ZZ(ICO  ) * WMCO
         PROD1(IEL) = ZZ(ICO2 ) * WMCO2
         PROD2(IEL) = ZZ(IH2O ) * WMH2O
         OXYD(IEL)  = ZZ(IO2  ) * WMO2
         XINER(IEL) = ZZ(IN2  ) * WMN2
C
        ENDIF
C
      ENDDO
C
      N2 = 0
      DO IEL = 1, NCEL
        IF ( FUEL3(IEL) .LT. (-EPZERO) ) THEN
          N2= N2 + 1
        ENDIF
      ENDDO
C
      IF ( IEQCO2 .EQ. 1 .OR. IEQCO2 .EQ.2 ) THEN
C

        IF (IRANGP.GE.0) THEN
C
          CALL PARCPT(N2)
          CALL PARCPT(NCLCO)
          CALL PARCPT(NCLCO1)
          CALL PARCPT(NCLO2)
          CALL PARCPT(NCLO21)
C
          CALL PARMIN(CMINS1)
          CALL PARMIN(CMINA1)
          CALL PARMIN(CMINS2)
          CALL PARMIN(CMINA2)
C
          CALL PARMAX(CMAXS1)
          CALL PARMAX(CMAXA1)
          CALL PARMAX(CMAXS2)
          CALL PARMAX(CMAXA2)
C
        ENDIF
C
        WRITE(NFECRA,*) ' Point a CO < 0 ',N2
        WRITE(NFECRA,*) ' Clipping CO2 - CO sans avec : ',NCLCO,NCLCO1
        WRITE(NFECRA,*) ' Min Max sans                : ',CMINS1,CMAXS1
        WRITE(NFECRA,*) ' Min Max avec                : ',CMINA1,CMAXA1
        WRITE(NFECRA,*) ' Clipping CO2 - O2 sans avec : ',NCLO2,NCLO21
        WRITE(NFECRA,*) ' Min Max sans                : ',CMINS2,CMAXS2
        WRITE(NFECRA,*) ' Min Max avec                : ',CMINA2,CMAXA2
C
      ENDIF
C
C=======================================================================
C 5.  IMPRESSION
C=======================================================================
C
      N2  = 0
      N3  = 0
      N4  = 0
      N5  = 0
      N6  = 0
      N7  = 0
      N8  = 0
      N9  = 0
      N10 = 0
C
C --> Controle des parametres de la pdf
C
      N1 = NCEL
C
      DO IEL = 1, NCEL
        IF ( INDPDF(IEL).NE.0 ) THEN
          N2 = N2 +1
        ENDIF
      ENDDO
C
C --> Controle des differentes valeurs des fractions massiques
C
      SOMMIN = 1.D+20
      SOMMAX =-1.D+20
      DO IEL = 1, NCEL
C
        SOMM = FUEL1(IEL) + FUEL2(IEL) + FUEL3(IEL)
     &       + OXYD(IEL)
     &       + PROD1(IEL) + PROD2(IEL) + XINER(IEL)
C
        SOMMIN = MIN(SOMMIN,SOMM)
        SOMMAX = MAX(SOMMAX,SOMM)
C
        IF ( ABS(SOMM-1.D0).LT.EPSICP )   THEN
          N3  = N3 +1
        ENDIF
C
        IF ( FUEL1(IEL).LT.(-EPZERO) .OR.
     &       FUEL1(IEL).GT.(1.D0+EPZERO) ) N4 = N4+1
        IF ( FUEL2(IEL).LT.(-EPZERO) .OR.
     &       FUEL2(IEL).GT.(1.D0+EPZERO) ) N5 = N5+1
        IF ( FUEL3(IEL).LT.(-EPZERO) .OR.
     &       FUEL3(IEL).GT.(1.D0+EPZERO) ) N6 = N6+1
        IF ( OXYD(IEL).LT.(-EPZERO) .OR.
     &       OXYD(IEL).GT.(1.D0+EPZERO)  ) N7 = N7+1
        IF ( XINER(IEL).LT.(-EPZERO) .OR.
     &       XINER(IEL).GT.(1.D0+EPZERO) ) N8 = N8+1
        IF ( PROD1(IEL).LT.(-EPZERO) .OR.
     &       PROD1(IEL).GT.(1.D0+EPZERO) ) N9 = N9+1
        IF ( PROD2(IEL).LT.(-EPZERO) .OR.
     &       PROD2(IEL).GT.(1.D0+EPZERO) ) N10 = N10+1
C
      ENDDO
C
      N1 = NCEL
C
      IF (IRANGP.GE.0) THEN
        INTTMP( 1) = N1
        INTTMP( 2) = N2
        INTTMP( 3) = N3
        INTTMP( 4) = N4
        INTTMP( 5) = N5
        INTTMP( 6) = N6
        INTTMP( 7) = N7
        INTTMP( 8) = N8
        INTTMP( 9) = N9
        INTTMP(10) = N10
        CALL PARISM(NBRINT,INTTMP)
C       ===========
      ENDIF
C
      WRITE(NFECRA,1000) N1 , N2
C
      WRITE(NFECRA,2200) N3 , N4, N5, N6, N7, N8,
     &                        N9, N10
C
      IF ( IRANGP .GE. 0 ) THEN
        CALL PARMIN(SOMMIN)
        CALL PARMAX(SOMMAX)
      ENDIF
C
      write(NFECRA,*) ' Somme Min MAX = ',SOMMIN,SOMMAX
C
C-------
C FORMAT
C-------
C
 1000 FORMAT (/,
     &'MODELISATION DE LA COMBUSTION AVEC LE MODELE DE DIFFUSION ',
     &'TURBULENTE (CPCYM2)',/,
     &'CHIMIE RAPIDE A 3 CORPS - EXTENSION A 3 COMBUSTIBLES ',
     &'(Application au FUEL)',/,
     &'==========================================================',
     &'==================',/,
     &' Nb de points de calculs                                     = ',
     &   I9,/,
     &' Nb de points turbulents (passage par les PDF)               = ',
     &   I9)
C
 2200 FORMAT(/,
     &'CONTROLE DES VALEURS DES FRACTIONS MASSIQUES',/,
     &' Nb de points de calculs qui respectent somme des Yi = 1    = ',
     &   I9,/,
     &' Nb de points YCHX1 ,YCHX2    < 0 ou > 1 = ',I9,I9,/,
     &' Nb de points YC0             < 0 ou > 1 = ',I9,/,
     &' Nb de points YO2  ,YN2       < 0 ou > 1 = ',I9,I9,/,
     &' Nb de points YCO2 ,YH2O      < 0 ou > 1 = ',I9,I9,/)
C
 3000 FORMAT(/,
     &'MODELE DE CO CHARBON ACTIVE : ',/,
     &'       AVEC  IEQCO2 = ',I2,/,
     &'HORS SEUL LES OPTIONS 0 , 1 et 2 SONT DISPONIBLES',/,
     &'    ARRET DU CALCUL : VERIFIER USPPMO')
C
C
      RETURN
      END
c@z
