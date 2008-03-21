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
     &   F1M    , F2M    , F3M   , F4M   , F5M  ,
     &   F1CL   , F2CL   , F3CL  , F4CL  , F5CL ,
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
      DOUBLE PRECISION F1CL(NCELET) , F2CL(NCELET)
      DOUBLE PRECISION F3CL(NCELET) , F4CL(NCELET), F5CL(NCELET)
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
      INTEGER          NBRINT
      PARAMETER        (NBRINT = 11)
      INTEGER          INTTMP(NBRINT)
C
      DOUBLE PRECISION ANU1 , ANU2 , ANU3
      DOUBLE PRECISION SOMM
C
      DOUBLE PRECISION ZCHX10 , ZCHX20 , ZCO0   , ZCO20
      DOUBLE PRECISION ZN20   , ZH2O0  , ZO20
      DOUBLE PRECISION YCHX10 , YCHX20 , YCO0   , YCO20
      DOUBLE PRECISION YN20   , YH2O0  , YO20
      DOUBLE PRECISION ZCHX11 , ZCHX21 , ZCO1   , ZCO21
      DOUBLE PRECISION ZN21   , ZH2O1  , ZO21
      DOUBLE PRECISION ZCHX12 , ZCHX22 , ZCO2   , ZCO22
      DOUBLE PRECISION ZN22   , ZH2O2  , ZO22
      DOUBLE PRECISION ZCHX13 , ZCHX23 , ZCO3   , ZCO23
      DOUBLE PRECISION ZN23   , ZH2O3  , ZO23
      DOUBLE PRECISION REAC1  , REAC2  , REAC3
      DOUBLE PRECISION X2     , XCH    , XCK    , XASH , XWAT
C
      DOUBLE PRECISION ZCOF10 , ZCOF20
      DOUBLE PRECISION CX1M   , CX2M
      DOUBLE PRECISION WMCHX1 , WMCHX2 , WMCO   , WMCO2 , WMH2O
      DOUBLE PRECISION WMO2   , WMN2   , WMC
      DOUBLE PRECISION ACHX1F1, ACOF1  , ACHX2F2, ACOF2
      DOUBLE PRECISION ACOF3  , AO2F3  , AO2F4  , AN2F4
      DOUBLE PRECISION YO2OX  , AN2F3  , AH2OF5
C
C     Variables pour le support de la Pdf et sa description
C
      DOUBLE PRECISION ZZCL(7),ZZS1(7),ZZS2(7),ZZS3(7),ZZF4(7)
      DOUBLE PRECISION ZZ(7)
C
      DOUBLE PRECISION F4S1 , F4S2 , F4S3
      DOUBLE PRECISION BB1  , BB2  , DEN1 , DEN2
      DOUBLE PRECISION ZCO2T
C
      DOUBLE PRECISION SOMMIN,SOMMAX,AO2MIN,AO2MAX
      INTEGER          NCLO2,NCLCO
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
C        fonctions de F1, F2, F3 et F4
C
      ACOF3   =  1.D0/WMC
      AO2F3   = -.5D0/WMC
      AO2F4   = 1.D0 / (WMO2+XSI*WMN2)
      AN2F4   =  XSI / (WMO2+XSI*WMN2)
      AH2OF5  = 1.D0/WMH2O
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
        DO ICHA = 1, NCHARB
          DEN1   = 1.D0
     &           / ( A1(ICHA)*WMOLE(ICHX1C(ICHA))+B1(ICHA)*WMOLE(ICO))
          YCHX10 = YCHX10 + DEN1 *
     &           ( RTBMC(IEL,IF1MC(ICHA))*A1(ICHA)*WMOLE(ICHX1C(ICHA)) )
          ZCHX10 = ZCHX10 + DEN1 *
     &           ( RTBMC(IEL,IF1MC(ICHA))*A1(ICHA) )
          ZCOF10 = ZCOF10 + DEN1 *
     &           ( RTBMC(IEL,IF1MC(ICHA))*B1(ICHA) )
          DEN2   = 1.D0
     &           / ( A2(ICHA)*WMOLE(ICHX2C(ICHA))+B2(ICHA)*WMOLE(ICO))
          YCHX20 = YCHX20 + DEN2 *
     &           ( RTBMC(IEL,IF2MC(ICHA))*A2(ICHA)*WMOLE(ICHX2C(ICHA)) )
          ZCHX20 = ZCHX20 + DEN2 *
     &           ( RTBMC(IEL,IF2MC(ICHA))*A2(ICHA) )
          ZCOF20 = ZCOF20 + DEN2 *
     &           ( RTBMC(IEL,IF2MC(ICHA))*B2(ICHA) )
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
        IF ( F1M(IEL).GT.EPZERO ) THEN
          RTBMC(IEL,ICHX1F1) = ZCHX10 / F1M(IEL)
          RTBMC(IEL,ICOF1  ) = ZCOF10 / F1M(IEL)
        ENDIF
        RTBMC(IEL,ICHX2F2) = 0.D0
        RTBMC(IEL,ICOF2  ) = 0.D0
        IF ( F2M(IEL).GT.EPZERO ) THEN
          RTBMC(IEL,ICHX2F2) = ZCHX20 / F2M(IEL)
          RTBMC(IEL,ICOF2  ) = ZCOF20 / F2M(IEL)
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
      NCLO2 = 0
      NCLCO = 0
C
C=======================================================================
C 3. CALCUL DE LA COMPOSITION DU MELANGE SANS LES PDF
C    SI LES FLUCTUATIONS DE S SONT TROP FAIBLES
C=======================================================================
C
C
      DO IEL = 1, NCEL
C
        ZCO2T = (RTP(IEL,ISCA(IYCO2))/ (1.D0+W1(IEL)))/WMCO2
C
        IF ( INDPDF(IEL).EQ.0 ) THEN
C
C --> Calculs preliminaires
C
          CX1M    = RTBMC(IEL,IX1MC)
          CX2M    = RTBMC(IEL,IX2MC)
          WMCHX1  = WMOLAT(IATC) + CX1M*WMOLAT(IATH)
          WMCHX2  = WMOLAT(IATC) + CX2M*WMOLAT(IATH)
          ACHX1F1 = RTBMC(IEL,ICHX1F1)
          ACHX2F2 = RTBMC(IEL,ICHX2F2)
          ACOF1   = RTBMC(IEL,ICOF1)
          ACOF2   = RTBMC(IEL,ICOF2)
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
          YCHX10  = ZCHX10*WMCHX1
          ZCHX20  = ACHX2F2*F2M(IEL)
          YCHX20  = ZCHX20*WMCHX2
          ZCO0    = (ACOF1*F1M(IEL)+ACOF2*F2M(IEL) +
     &               ACOF3*F3M(IEL))
          YCO0    = ZCO0 * WMCO
          ZO20    = (AO2F4*F4M(IEL)+AO2F3*F3M(IEL))
          YO20    = ZO20*WMO2
          ZN20    = (AN2F4*F4M(IEL))
          YN20    = ZN20*WMN2
          ZCO20   = ZERO
          YCO20   = ZCO20*WMCO2
          ZH2O0   = ZERO + AH2OF5*F5M(IEL)
          YH2O0   = ZH2O0*WMH2O
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
          ELSE
            ZCHX13 = ZCHX12
            ZCHX23 = ZCHX22
C
            IF ( (ZCO2-ZCO2T) .LT. -EPZERO ) THEN
              NCLCO = NCLCO +1
            ENDIF
            IF ( (ZO22-ANU3*ZCO2T) .LT. -EPZERO ) THEN
              NCLO2 = NCLO2 +1
            ENDIF
            ZCO2T = MIN(ZCO2T,ZCO2,ZO22/ANU3)
            RTP(IEL,ISCA(IYCO2))= ZCO2T*(1.D0+W1(IEL))*WMCO2
C
            ZO23   = ZO22  - ANU3*ZCO2T
            ZCO3   = ZCO2  -      ZCO2T
            ZCO23  = ZCO2T
C
            ZH2O3  = ZH2O2
            ZN23   = ZN22
          ENDIF
C
          FUEL1(IEL) = ZCHX13 * WMCHX1
          FUEL2(IEL) = ZCHX23 * WMCHX2
          FUEL3(IEL) = ZCO3  * WMCO
          PROD1(IEL) = ZCO23 * WMCO2
          PROD2(IEL) = ZH2O3 * WMH2O
          OXYD(IEL)  = ZO23  * WMO2
          XINER(IEL) = ZN23  * WMN2
C
        ENDIF
C
      ENDDO
C
C=======================================================================
C 4. CALCUL DE LA COMPOSITION DU MELANGE AVEC LES PDF
C=======================================================================
C
      AO2MIN = 1.D+20
      AO2MAX = -1.D+20
      DO IEL = 1, NCEL
C
        ZCO2T = (RTP(IEL,ISCA(IYCO2))/ (1.D0+W1(IEL)))/WMCO2
C
        IF ( INDPDF(IEL).EQ.1 ) THEN
C
C --> Calculs preliminaires
C
          CX1M    = RTBMC(IEL,IX1MC)
          CX2M    = RTBMC(IEL,IX2MC)
          WMCHX1  = WMOLAT(IATC) + CX1M*WMOLAT(IATH)
          WMCHX2  = WMOLAT(IATC) + CX2M*WMOLAT(IATH)
          ACHX1F1 = RTBMC(IEL,ICHX1F1)
          ACHX2F2 = RTBMC(IEL,ICHX2F2)
          ACOF1   = RTBMC(IEL,ICOF1)
          ACOF2   = RTBMC(IEL,ICOF2)
C
          YO2OX  = WMO2 / (WMO2+XSI*WMN2)
          AN2F3  = (1.D0-YO2OX)/WMOLE(IN2) * (1.D0-F3MAX)
C
          AH2OF5 = 1.D0/WMOLE(IH2O)
C
          ANU1   = 0.25D0*(CX1M-CX2M)
          ANU2   = 0.25D0*(2.D0 + CX2M)
          ANU3   = 0.5D0
C
C         Concentrations aux extrémités
C
          DO II = 1,7
            ZZCL(II) = ZERO
            ZZF4(II) = ZERO
          ENDDO
C
          ZZCL(ICHX1) = ACHX1F1*F1CL(IEL)
          ZZCL(ICHX2) = ACHX2F2*F2CL(IEL)
          ZZCL(ICO)   = ACOF1 *F1CL(IEL)+ACOF2 *F2CL(IEL)
     &                                  +ACOF3 *F3CL(IEL)
          ZZCL(IH2O)  =                  AH2OF5*F5CL(IEL)
          ZZCL(IN2)   =                  AN2F4 *F4CL(IEL)
          ZZCL(IO2)   = AO2F4*F4CL(IEL)+AO2F3*F3CL(IEL)
C
          AO2MIN = mIN(AO2MIN, ZZCL(IO2))
          AO2MAX = MAX(AO2MAX, ZZCL(IO2))
C
          ZZF4(IO2)   = AO2f4*F4M(IEL)/(F4M(IEL)+F5M(IEL))
          ZZF4(IN2)   = AN2f4*F4M(IEL)/(F4M(IEL)+F5M(IEL))
          ZZF4(IH2O)  = AH2OF5*F5M(IEL)/(F4M(IEL)+F5M(IEL))
C
C         Calcul de la stoechiométrie de la première réaction
C
          F4S1 = ( ANU1 * ZZCL(ICHX1) + F4CL(IEL) * ZZF4(IO2) )
     &         / ( ANU1 * ZZCL(ICHX1) +             ZZF4(IO2) )
C
C         Calcul des concentrations au premier point stoechio
C         avant réaction
C
          DO II = 1,7
            ZZS1(II) = ZZCL(II)
     &          + (F4S1-F4CL(IEL))/(1.D0-F4CL(IEL))
     &           *(ZZF4(II)-ZZCL(II))
          ENDDO
C
          ZZS1(ICHX2) = ZZS1(ICHX2) + ZZS1(ICHX1)
          ZZS1(IO2)   = ZERO
          ZZS1(IH2O)  = ZZS1(IH2O) + 2.D0*ANU1*ZZS1(ICHX1)
          ZZS1(ICHX1) = ZERO
C
C         Calcul de la stoechiometrie de la deuxième réaction
C
          F4S2 = ( ANU2 * ZZS1(ICHX2) + F4S1 * ZZF4(IO2) )
     &          /( ANU2 * ZZS1(ICHX2) +        ZZF4(IO2) )
C
C         Calcul des concentrations au deuxième point stoechio
C         avant réaction
C
          DO II = 1,7
             ZZS2(II) = ZZS1(II)
     &          + (F4S2-F4S1) / (1.D0-F4S1) * (ZZF4(II) - ZZS1(II))
          ENDDO
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
          F4S3 = ( ANU3 * ZZS2(ICO ) + F4S2 * ZZF4(IO2) )
     &         / ( ANU3 * ZZS2(ICO ) +        ZZF4(IO2) )
C
C         Calcul des concentrations au troisième point stoechio
C         avant réaction
C
          DO II = 1,7
             ZZS3(II) = ZZS2(II)
     &          + (F4S3-F4S2) / (1.D0-F4S2) * (ZZF4(II) - ZZS2(II))
          ENDDO
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
           ZZ(II) = D4CL(IEL)*ZZcl(II)+D4F4(IEL)*ZZF4(II)
         ENDDO
C
C        Intégration sur le premier intervalle de richesse (entre cl et s1)
C
         BB1 = MAX(F4M1(IEL),F4CL(IEL))
         BB2 = MIN(F4M2(IEL),F4S1)
         IF( BB2.GT.BB1 ) THEN
           DO II = 1,7
             ZZ(II) =  ZZ(II) + HREC(IEL)*(BB2-BB1)/(F4S1-F4CL(IEL))
     &               *( (ZZCL(II)*F4S1-ZZS1(II)*F4CL(IEL))
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
     &             * ( (ZZS3(II)*1.D0-ZZF4(II)*F4S3)
     &                +(ZZF4(II)-ZZS3(II))*(BB1+BB2)*0.5D0)
           ENDDO
         ENDIF
C
         IF ( IEQCO2 .EQ. 1 ) THEN
C
           IF ( (ZZ(ICO) - ZCO2T) .LT. -EPZERO ) THEN
             NCLCO = NCLCO +1
           ENDIF
           IF ( (ZZ(IO2) -0.5D0*ZCO2T) .LT. -EPZERO ) THEN
             NCLO2 = NCLO2 +1
           ENDIF
           ZCO2T = MIN(ZCO2T,ZZ(ICO),2.D0*ZZ(IO2))
C
           ZZ(ICO ) = ZZ(ICO) - ZCO2T
           ZZ(IO2 ) = ZZ(IO2) -0.5D0*ZCO2T
           ZZ(ICO2) = ZCO2T
C
           RTP(IEL,ISCA(IYCO2))= ZCO2T*(1.D0+W1(IEL))*WMCO2
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
     &'TURBULENTE (CPCYM1)',/,
     &'CHIMIE RAPIDE A 3 CORPS - EXTENSION A 3 COMBUSTIBLES ',
     &'(Application au FUEL)',/,
     &'==========================================================',
     &'==================',/,
     &' Nb de points de calculs                                     = ',
     &   I6,/,
     &' Nb de points turbulents (passage par les PDF)               = ',
     &   I6)
C
 2200 FORMAT(/,
     &'CONTROLE DES VALEURS DES FRACTIONS MASSIQUES',/,
     &' Nb de points de calculs qui respectent somme des Yi = 1    = ',
     &   I6,/,
     &' Nb de points YCHX1 ,YCHX2    < 0 ou > 1 = ',I6,I6,/,
     &' Nb de points YC0             < 0 ou > 1 = ',I6,/,
     &' Nb de points YO2  ,YN2       < 0 ou > 1 = ',I6,I6,/,
     &' Nb de points YCO2 ,YH2O      < 0 ou > 1 = ',I6,I6,/)
C
C
      RETURN
      END
c@z
