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
                        SUBROUTINE CPLYM1
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NITBMC , NRTBMC ,
     &   RTP    ,
     &   F1M    , F2M    , F3M    , F4M    ,
     &   INDPDF , SI7    , SI8    , SP2M   , F4I7   ,
     &   DSI7   , DSI8   , SDEB   , SFIN   , HAUT   ,
     &   FUEL1  , FUEL2  , FUEL3  , OXYD   , PROD1  , PROD2  ,
     &   XINER  ,
     &   ITBMC  , RTBMC  ,
     &   IPI7I8 ,
     &   SP     , SI     , SR     , W1     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
CFONC    --------------------------------------------------------------
CFONC
CFONC     ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE
CFONC
CFONC       COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
CFONC       TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON
CFONC
CFONC       CALCUL DES CONCENTRATIONS GAZEUSES MOYENNES
CFONC        VALEURS CELLULES
CFONC        ----------------
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NITBMC       ! E  !  -> ! TAILLE DU MACRO TABLEAU MC ENTIERS   !
CARGU ! NRTBMC       ! E  !  -> ! TAILLE DU MACRO TABLEAU MC REELS     !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! F1M          ! TR !  -> ! MOYENNE DU TRACEUR 1 MVl [CHx1m+CO]  !
CARGU ! F2M          ! TR !  -> ! MOYENNE DU TRACEUR 2 MVL [CHx2m+CO]  !
CARGU ! F3M          ! TR !  -> ! MOYENNE DU TRACEUR 3 (CO C.HET)      !
CARGU ! F4M          ! TR !  -> ! MOYENNE DU TRACEUR 4 (AIR)           !
CARGU ! INDPDF       ! TE !  -> ! PASSAGE PAR LES PDF                  !
CARGU ! SI7          ! TR !  -> ! ABSCISSE CURVILIGNE AU POINT I7      !
CARGU ! SI8          ! TR !  -> ! ABSCISSE CURVILIGNE AU POINT I8      !
CARGU ! SP2M         ! TR !  -> ! VARIANCE DROITE SUPPORT              !
CARGU ! F4I7         ! TR !  -> ! F4 AU POINT I7                       !
CARGU ! DSI7         ! TR !  -> ! DIRAC AU POINT I7                    !
CARGU ! DSI8         ! TR !  -> ! DIRAC AU POINT I8                    !
CARGU ! SDEB         ! TR !  -> ! ABSCISSE DEBUT RECTANGLE             !
CARGU ! SFIN         ! TR !  -> ! ABSCISSE FIN RECTANGLE               !
CARGU ! HREC         ! TR !  -> ! HAUTEUR RECTANGLE                    !
CARGU ! FUEL1        ! TR !  <- ! FRACTION MASSIQUE CHx1m              !
CARGU ! FUEL2        ! TR !  <- ! FRACTION MASSIQUE CHx2m              !
CARGU ! FUEL3        ! TR !  <- ! FRACTION MASSIQUE CO                 !
CARGU ! OXYD         ! TR !  <- ! FRACTION MASSIQUE O2                 !
CARGU ! PROD1        ! TR !  <- ! FRACTION MASSIQUE CO2                !
CARGU ! PROD2        ! TR !  <- ! FRACTION MASSIQUE H2O                !
CARGU ! XINER        ! TR !  <- ! FRACTION MASSIQUE N2                 !
CARGU ! ITBMC        ! TR ! <-> ! MACRO TABLEAU ENTIER MC TRAVAIL      !
CARGU ! RTBMC        ! TR ! <-> ! MACRO TABLEAU REEL   MC TRAVAIL      !
CARGU ! IPI7I8       ! TE ! <-> ! POSITION DES POINTS I7 ET I8 DANS    !
CARGU !              !    !     ! LES DIFFERENTS DOMAINES DE RICHESSE  !
CARGU ! SP           ! TR ! <-> ! DELIMITATION SUR LE SUPPORT DES      !
CARGU !              !    !     ! DOMAINES PAUVRE ET INTERMEDIAIRE     !
CARGU ! SI           ! TR ! <-> ! DELIMITATION SUR LE SUPPORT DES      !
CARGU !              !    !     ! DOMAINES INTERMEDIAIRE ET RICHE      !
CARGU ! SR           ! TR ! <-> ! DELIMITATION SUR LE SUPPORT DES      !
CARGU !              !    !     ! DOMAINES RICHE ET HYPER-RICHE        !
CARGU ! W1           ! TR ! <-> ! TABLEAU DE TRAVAIL                   !
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
      INTEGER          NCELET , NCEL
      INTEGER          NITBMC , NRTBMC
      INTEGER          ITBMC(NCELET,NITBMC)
      INTEGER          INDPDF(NCELET) , IPI7I8(NCELET)
C
      DOUBLE PRECISION RTP(NCELET,*)
      DOUBLE PRECISION F1M(NCELET)  , F2M(NCELET)
      DOUBLE PRECISION F3M(NCELET)  , F4M(NCELET)
      DOUBLE PRECISION DSI7(NCELET) , DSI8(NCELET) , SDEB(NCELET)
      DOUBLE PRECISION SFIN(NCELET) , HAUT(NCELET)
      DOUBLE PRECISION SI7(NCELET)  , SI8(NCELET)  , SP2M(NCELET)
      DOUBLE PRECISION F4I7(NCELET)
      DOUBLE PRECISION FUEL1(NCELET), FUEL2(NCELET) , FUEL3(NCELET)
      DOUBLE PRECISION OXYD(NCELET) , PROD1(NCELET) , PROD2(NCELET)
      DOUBLE PRECISION XINER(NCELET)
      DOUBLE PRECISION RTBMC(NCELET,NRTBMC)
      DOUBLE PRECISION SP(NCELET) , SI(NCELET) , SR(NCELET)
      DOUBLE PRECISION W1(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL, ICHA, ITBR, ITBI, IGE
C
      INTEGER               N1 , N2 , N3 , N4 , N5 , N6 , N7 , N8 , N9
      INTEGER          N10, N11, N12, N13, N14, N15, N16, N17, N18, N19
      INTEGER          N20, N21, N22, N23, N24, N25, N26, N27, N28, N29
      INTEGER          N30, N31, N32
      INTEGER          NBRINT
      PARAMETER        (NBRINT = 32)
      INTEGER          INTTMP(NBRINT)
C
      DOUBLE PRECISION ZCOF10 , ZCOF20
      DOUBLE PRECISION CX1M   , CX2M
      DOUBLE PRECISION WMCHX1 , WMCHX2 , WMCO   , WMCO2 , WMH2O
      DOUBLE PRECISION WMO2   , WMN2   , WMC
      DOUBLE PRECISION ACHX1F1, ACOF1  , ACHX2F2, ACOF2
      DOUBLE PRECISION CC1    , CC2    , CC3    , CC4
      DOUBLE PRECISION F1I7   , F2I7   , F3I7
      DOUBLE PRECISION VARS
      DOUBLE PRECISION ZC(NGAZEM)
      DOUBLE PRECISION AC0(NGAZEM) , BC0(NGAZEM)
      DOUBLE PRECISION AC1(NGAZEM) , BC1(NGAZEM)
      DOUBLE PRECISION AC2(NGAZEM) , BC2(NGAZEM)
      DOUBLE PRECISION AC3(NGAZEM) , BC3(NGAZEM)
      DOUBLE PRECISION AI(NGAZEM)  , BI(NGAZEM)
      DOUBLE PRECISION ANU1 , ANU2 , ANU3
      DOUBLE PRECISION S1, S2, S3, SM1, SM2
      DOUBLE PRECISION DEN1   , DEN2   , SOMM
#ifdef DEBUG
      DOUBLE PRECISION SOMMF  , MOM0   , MOM1   , MOM2
#endif
      DOUBLE PRECISION ZCHX10 , ZCHX20 , ZCO0   , ZCO20
      DOUBLE PRECISION ZN20   , ZH2O0  , ZO20
      DOUBLE PRECISION YCHX10 , YCHX20
#ifdef DEBUG
      DOUBLE PRECISION YCO0   , YCO20
      DOUBLE PRECISION YN20   , YH2O0  , YO20
#endif
      DOUBLE PRECISION ZCHX11 , ZCHX21 , ZCO1   , ZCO21
      DOUBLE PRECISION ZN21   , ZH2O1  , ZO21
      DOUBLE PRECISION ZCHX12 , ZCHX22 , ZCO2   , ZCO22
      DOUBLE PRECISION ZN22   , ZH2O2  , ZO22
      DOUBLE PRECISION ZCHX13 , ZCHX23 , ZCO3   , ZCO23
      DOUBLE PRECISION ZN23   , ZH2O3  , ZO23
      DOUBLE PRECISION REAC1  , REAC2  , REAC3
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
C --- Initialisation tableau de travail
C
      DO IEL = 1, NCEL
        W1(IEL) = ZERO
        DO ITBR = 1, NRTBMC
          RTBMC(IEL,ITBR) = ZERO
        ENDDO
        DO ITBI = 1, NITBMC
          ITBMC(IEL,ITBI) = 0
        ENDDO
      ENDDO
C
C=======================================================================
C 2. CALCULS PRELIMINAIRES
C=======================================================================
C
C --- Calcul de F1M de chaque charbon dans RTBMC
C            de F2M de chaque charbon dans RTBMC
C
C ------ RTBMC(IF1MC(ICHA)) = F1M(ICHA)
C        RTBMC(IF2MC(ICHA)) = F2M(ICHA)
C
      DO ICHA = 1, NCHARB
        DO IEL = 1, NCEL
          RTBMC(IEL,IF1MC(ICHA)) = RTP(IEL,ISCA(IF1M(ICHA)))
          RTBMC(IEL,IF2MC(ICHA)) = RTP(IEL,ISCA(IF2M(ICHA)))
        ENDDO
      ENDDO
C
C
C --- Masses molaires
C
      WMCO   = WMOLE(ICO  )
      WMO2   = WMOLE(IO2  )
      WMCO2  = WMOLE(ICO2 )
      WMH2O  = WMOLE(IH2O )
      WMN2   = WMOLE(IN2  )
      WMC    = WMOLAT(IATC)
C
C
C --- Calcul de CX1M, CX2M
C     et des coefficient relatifs a l'expression de
C        ZCHx1m0, ZCHx2m0, ZCO0, ZO20 et ZN20
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
        IPI7I8(IEL) = 0
        SP(IEL)     = ZERO
        SI(IEL)     = ZERO
        SR(IEL)     = ZERO
      ENDDO
C
C
C
C=======================================================================
C 3. CALCUL DE LA COMPOSITION DU MELANGE SANS LES PDF
C    SI LES FLUCTUATIONS DE S SONT TROP FAIBLES
C    SI SI7 ET SI8 SONT TROP PROCHE DE 0
C    ALORS F1 = F1M, F2 = F2M , F3 = F3M ET F4 = F4M
C=======================================================================
C
      DO IEL = 1, NCEL
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
          ZCHX20  = ACHX2F2*F2M(IEL)
          ZCO0    = (ACOF1*F1M(IEL)+ACOF2*F2M(IEL) +
     &               ACOF3*F3M(IEL))
          ZO20    = (AO2F4*F4M(IEL)+AO2F3*F3M(IEL))
          ZN20    = (AN2F4*F4M(IEL))
          ZCO20   = ZERO
          ZH2O0   = ZERO
C ---- Test
#ifdef DEBUG
          YCHX10  = ZCHX10*WMCHX1
          YCHX20  = ZCHX20*WMCHX2
          YCO0    = ZCO0 * WMCO
          YO20    = ZO20*WMO2
          YN20    = ZN20*WMN2
          YCO20   = ZCO20*WMCO2
          YH2O0   = ZH2O0*WMH2O
          SOMM    = YCHX10 + YCHX20 + YCO0 + YO20 + YN20 +
     &              YCO20  + YH2O0
          IF (ABS(SOMM-1.D0).GT.EPSICP) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'SOMME0', SOMM
            write(NFECRA,*) 'F1M',F1M(IEL)
            write(NFECRA,*) 'F2M',F2M(IEL)
            write(NFECRA,*) 'F3M',F3M(IEL)
            write(NFECRA,*) 'F4M',F4M(IEL)
          ENDIF
#endif
C ---- Fin test
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
          REAC3  = MIN(ZCO2, (ZO22/ANU3))
          ZCHX13 = ZCHX12
          ZCHX23 = ZCHX22
          ZO23   = ZO22  - ANU3*REAC3
          ZCO3   = ZCO2  -      REAC3
          ZCO23  = ZCO22 +      REAC3
          ZH2O3  = ZH2O2
          ZN23   = ZN22
C
          FUEL1(IEL) = ZCHX13 * WMCHX1
          FUEL2(IEL) = ZCHX23 * WMCHX2
          FUEL3(IEL) = ZCO3   * WMCO
          PROD1(IEL) = ZCO23  * WMCO2
          PROD2(IEL) = ZH2O3  * WMH2O
          OXYD(IEL)  = ZO23   * WMO2
          XINER(IEL) = ZN23   * WMN2
C
C ---- Test
#ifdef DEBUG
          SOMM = FUEL1(IEL) + FUEL2(IEL) + FUEL3(IEL)
     &         + PROD1(IEL) + PROD2(IEL)
     &         + OXYD(IEL)  + XINER(IEL)
c          IF (ABS(SOMM-1.D0).GT.EPSICP) THEN
          IF ( (ABS(SOMM-1.D0).GT.EPSICP) .OR.
     &         (PROD2(IEL).LT.(-EPZERO) .OR.
     &          PROD2(IEL).GT.(1.D0+EPZERO)) ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'INDPDF0',INDPDF(IEL)
            write(NFECRA,*) 'F1M',F1M(IEL)
            write(NFECRA,*) 'F2M',F2M(IEL)
            write(NFECRA,*) 'F3M',F3M(IEL)
            write(NFECRA,*) 'F4M',F4M(IEL)
            SOMMF = F1M(IEL) + F2M(IEL) + F3M(IEL) + F4M(IEL)
            write(NFECRA,*) 'SOMME F', SOMMF
            write(NFECRA,*) 'AO2F4, AO2F3 ', AO2F4, AO2F3
            write(NFECRA,*) 'ZCHX1, 0,1,2,3',
     &                       ZCHX10,ZCHX11,ZCHX12,ZCHX13
            write(NFECRA,*) 'ZCHX2, 0,1,2,3',
     &                       ZCHX20,ZCHX21,ZCHX22,ZCHX23
            write(NFECRA,*) 'ZCO, 0,1,2,3',
     &                       ZCO0,ZCO1,ZCO2,ZCO3
            write(NFECRA,*) 'ZCO2, 0,1,2,3',
     &                       ZCO20,ZCO21,ZCO22,ZCO23
            write(NFECRA,*) 'ZH2O, 0,1,2,3',
     &                       ZH2O0,ZH2O1,ZH2O2,ZH2O3
            write(NFECRA,*) 'ZO2, 0,1,2,3',
     &                       ZO20,ZO21,ZO22,ZO23
            write(NFECRA,*) 'ZN2, 0,1,2,3',
     &                       ZN20,ZN21,ZN22,ZN23
            write(NFECRA,*) 'Fuel1',FUEL1(IEL)
            write(NFECRA,*) 'Fuel2',FUEL2(IEL)
            write(NFECRA,*) 'Fuel3',FUEL3(IEL)
            write(NFECRA,*) 'Oxyd',OXYD(IEL)
            write(NFECRA,*) 'Iner',XINER(IEL)
            write(NFECRA,*) 'Prod1',PROD1(IEL)
            write(NFECRA,*) 'Prod2',PROD2(IEL)
            write(NFECRA,*) 'SOMME Y', SOMM
          ENDIF
#endif
C ---- Fin Test
C
        ENDIF
C
      ENDDO
C
C
C
C=======================================================================
C 4. CALCUL DE LA COMPOSITION DU MELANGE AVEC LES PDF
C    ATTENTION : ON SE LIMITE POUR L'INSTANT A INDPDF = 3 SUPPORT (I4,M)
C=======================================================================
C
C
C -->  I7 = I4 : F1I7 = F2I7 = F3I7 = 0
C                F4I7 = 1
C
C -->  RELATIONS SUR LE SUPPORT DE LA PDF
C        FJ = FJM + CJ*S avec CJ = (FJM-FJI7)/SI7 (J = 1 a 4)
C
C Devectorisation de la boucle pour VPP 5000
C
!OCL SCALAR
C
      DO IEL = 1, NCEL
C
        IF ( INDPDF(IEL).EQ.3 ) THEN
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
C ---- Calcul de F1I7, F2I7, F3I7
C
          F1I7   = ZERO
          F2I7   = ZERO
          F3I7   = 1.D0 - F4I7(IEL)
C
C ----- Calcul de C1, C2, C3 et C4
C
          CC1 = (F1I7     -F1M(IEL)) / SI7(IEL)
          CC2 = (F2I7     -F2M(IEL)) / SI7(IEL)
          CC3 = (F3I7     -F3M(IEL)) / SI7(IEL)
          CC4 = (F4I7(IEL)-F4M(IEL)) / SI7(IEL)
C
C ---- Rq : les autres pointeurs sont deja definis dans cplecd.F
C
C --> Initialisation
C
          DO IGE = 1, (NGAZE-2*NCHARB)
            ZC(IGE) = ZERO
          ENDDO
C
C --> Calculs des coefficients avant combustion phase gaz
C       ZC(IGE) = AC0(IGE) + BC0(IGE) * S
C
          AC0(ICHX1) = ACHX1F1*F1M(IEL)
          BC0(ICHX1) = ACHX1F1*CC1
          AC0(ICHX2) = ACHX2F2*F2M(IEL)
          BC0(ICHX2) = ACHX2F2*CC2
          AC0(ICO  ) = ACOF1*F1M(IEL)+ACOF2*F2M(IEL) +
     &                 ACOF3*F3M(IEL)
          BC0(ICO  ) = ACOF1*CC1     +ACOF2*CC2      +
     &                 ACOF3*CC3
          AC0(ICO2 ) = ZERO
          BC0(ICO2 ) = ZERO
          AC0(IH2O ) = ZERO
          BC0(IH2O ) = ZERO
          AC0(IO2  ) = AO2F3*F3M(IEL)+AO2F4*F4M(IEL)
          BC0(IO2  ) = AO2F3*CC3     +AO2F4*CC4
          AC0(IN2  ) = AN2F4*F4M(IEL)
          BC0(IN2  ) = AN2F4*CC4
C
C ---- Test
#ifdef DEBUG
          WMOLE(ICHX1) = WMCHX1
          WMOLE(ICHX2) = WMCHX2
          SOMM  = 0.D0
          SOMMF = 0.D0
          DO IGE = 1, (NGAZE-2*NCHARB)
            SOMM  = SOMM  + AC0(IGE)*WMOLE(IGE)
            SOMMF = SOMMF + BC0(IGE)*WMOLE(IGE)
          ENDDO
          IF (ABS(SOMM-1.D0).GT.EPSICP) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'F4I7 ', F4I7(IEL)
            write(NFECRA,*) 'SOMM AC0', SOMM
          ENDIF
          IF (ABS(SOMMF).GT.EPSICP) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'SOMM BC0', SOMMF
          ENDIF
#endif
C ---- Fin Test
C
C --> Les pics de Dirac sont situes dans des domaines
C     ou la reaction n'est pas possible (points froids)
C
          DO IGE = 1, (NGAZE-2*NCHARB)
            ZC(IGE) = DSI7(IEL)*(AC0(IGE)+BC0(IGE)*SI7(IEL))
     &              + DSI8(IEL)*(AC0(IGE)+BC0(IGE)*SI8(IEL))
          ENDDO
C
C --> Initialisation des coefficients avant integration
C
          DO IGE = 1, (NGAZE-2*NCHARB)
            AC1(IGE) = AC0(IGE)
            BC1(IGE) = BC0(IGE)
            AC2(IGE) = AC0(IGE)
            BC2(IGE) = BC0(IGE)
            AC3(IGE) = AC0(IGE)
            BC3(IGE) = BC0(IGE)
          ENDDO
C
C --> Integration dans le domaine hyper riche
C
C Rq : On rajoute le test suivant CAR
C      si F1M = 0 alors le domaine hyper riche n'existe pas
C
          IF ( F1M(IEL).GT.EPSICP ) THEN
C
          ANU1 = 0.25D0*(CX1M-CX2M)
          S1   = - (AC0(IO2)-AC0(ICHX1)*ANU1)
     &           / (BC0(IO2)-BC0(ICHX1)*ANU1)
          SM1  = MAX(SDEB(IEL),S1)
          SM2  = SFIN(IEL)
C
          DO IGE = 1, (NGAZE-2*NCHARB)
            AI(IGE) = AC0(IGE)
            BI(IGE) = BC0(IGE)
          ENDDO
          BI(ICHX1) = (AC0(ICHX1)+BC0(ICHX1)*SI8(IEL))/(SI8(IEL)-S1)
          AI(ICHX1) = -S1*BI(ICHX1)
          AI(ICHX2) = AC0(ICHX2) +           (AC0(ICHX1)-AI(ICHX1))
          BI(ICHX2) = BC0(ICHX2) +           (BC0(ICHX1)-BI(ICHX1))
          AI(IO2  ) = AC0(IO2  ) -      ANU1*(AC0(ICHX1)-AI(ICHX1))
          BI(IO2  ) = BC0(IO2  ) -      ANU1*(BC0(ICHX1)-BI(ICHX1))
          AI(IH2O ) = AC0(IH2O ) + 2.D0*ANU1*(AC0(ICHX1)-AI(ICHX1))
          BI(IH2O ) = BC0(IH2O ) + 2.D0*ANU1*(BC0(ICHX1)-BI(ICHX1))
C
          IF ( SM1.LT.SM2 ) THEN
C
            DO IGE = 1, (NGAZE-2*NCHARB)
              ZC(IGE) = ZC(IGE)
     &          + HAUT(IEL)*(SM2-SM1)*(AI(IGE)+.5D0*BI(IGE)*(SM1+SM2))
            ENDDO
C
          ENDIF
C
          DO IGE = 1, (NGAZE-2*NCHARB)
            AC1(IGE) = ( AC0(IGE)*S1-AI(IGE)*SI7(IEL) +
     &                     S1*SI7(IEL)*(BC0(IGE)-BI(IGE)) )
     &                 / ( S1-SI7(IEL) )
            BC1(IGE) = (AI(IGE)+BI(IGE)*S1-AC1(IGE)) / S1
          ENDDO
C
C ---- Apres la premiere reaction, on se trouve hors du domaine hyper riche
C      CHX1 ne doit plus exister AC1(ICHX1) = BC1(ICHX1) = 0
C
C ---- Test
#ifdef DEBUG
          SOMM  = 0.D0
          SOMMF = 0.D0
          DO IGE = 1, (NGAZE-2*NCHARB)
            SOMM  = SOMM  + AC1(IGE)*WMOLE(IGE)
            SOMMF = SOMMF + BC1(IGE)*WMOLE(IGE)
          ENDDO
            IF ( ABS(AC1(ICHX1)).GT.EPSICP .OR.
     &         ABS(BC1(ICHX1)).GT.EPSICP ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'AC1(ICHX1)', AC1(ICHX1)
            write(NFECRA,*) 'BC1(ICHX1)', BC1(ICHX1)
          ENDIF
            IF ( ABS(SOMM-1.D0).GT.EPSICP ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'SOMM AC1', SOMM
          ENDIF
            IF ( ABS(SOMMF).GT.EPSICP ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'SOMM BC1', SOMMF
          ENDIF
            IF ( ABS(AC1(IO2)+BC1(IO2)*S1).GT.EPSICP ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'AC01(IO2)', AC0(IO2),AC1(IO2)
            write(NFECRA,*) 'BC01(IO2)', BC0(IO2),BC1(IO2)
          ENDIF
#endif
C ---- Fin Test
C
          ENDIF
C
C --> Integration dans le domaine riche
C
          ANU2 = 0.25D0*(2.D0 + CX2M)
          S2   = - (AC1(IO2)-AC1(ICHX2)*ANU2)
     &           / (BC1(IO2)-BC1(ICHX2)*ANU2)
C
          SM1  = MAX(SDEB(IEL),S2)
          SM2  = MIN(SFIN(IEL),S1)
C
          DO IGE = 1, (NGAZE-2*NCHARB)
            AI(IGE) = AC1(IGE)
            BI(IGE) = BC1(IGE)
          ENDDO
C
          BI(ICHX2) = (AC1(ICHX2)+BC1(ICHX2)*S1)/(S1-S2)
          AI(ICHX2) = -S2*BI(ICHX2)
          AI(IO2  ) = AC1(IO2 ) -      ANU2*(AC1(ICHX2)-AI(ICHX2))
          BI(IO2  ) = BC1(IO2 ) -      ANU2*(BC1(ICHX2)-BI(ICHX2))
          AI(ICO  ) = AC1(ICO ) +           (AC1(ICHX2)-AI(ICHX2))
          BI(ICO  ) = BC1(ICO ) +           (BC1(ICHX2)-BI(ICHX2))
          AI(IH2O ) = AC1(IH2O) + .5D0*CX2M*(AC1(ICHX2)-AI(ICHX2))
          BI(IH2O ) = BC1(IH2O) + .5D0*CX2M*(BC1(ICHX2)-BI(ICHX2))
C
          IF ( SM1.LT.SM2 ) THEN
C
            DO IGE = 1, (NGAZE-2*NCHARB)
              ZC(IGE) = ZC(IGE)
     &          + HAUT(IEL)*(SM2-SM1)*(AI(IGE)+.5D0*BI(IGE)*(SM1+SM2))
            ENDDO
C
          ENDIF
C
          DO IGE = 1, (NGAZE-2*NCHARB)
c              AC2(IGE) = (AC1(IGE)*S2-AI(IGE)*SI7(IEL)) +
c     &                    S2*SI7(IEL)*(BC1(IGE)-BI(IGE))/(S2-SI7(IEL))
c              BC2(IGE) = (AI(IGE)+BI(IGE)*S2-AC2(IGE))/(S2-SI7(IEL))
            AC2(IGE) = ( AC1(IGE)*S2-AI(IGE)*SI7(IEL) +
     &                     S2*SI7(IEL)*(BC1(IGE)-BI(IGE)) )
     &                 / ( S2-SI7(IEL) )
            BC2(IGE) = (AI(IGE)+BI(IGE)*S2-AC2(IGE)) / S2
          ENDDO
C
C ---- Apres la seconde reaction, on se trouve hors du domaine riche
C      CHX1 et CHX2 ne doivent plus exister
C      AC2(ICHX1) = BC2(ICHX1) = AC2(ICHX2) = BC2(ICHX2) = 0
C
C ---- Test
#ifdef DEBUG
          SOMM  = 0.D0
          SOMMF = 0.D0
          DO IGE = 1, (NGAZE-2*NCHARB)
            SOMM  = SOMM  + AC2(IGE)*WMOLE(IGE)
            SOMMF = SOMMF + BC2(IGE)*WMOLE(IGE)
          ENDDO
          IF (ABS(AC2(ICHX1)).GT.EPSICP .OR.
     &        ABS(BC2(ICHX1)).GT.EPSICP .OR.
     &        ABS(AC2(ICHX2)).GT.EPSICP .OR.
     &        ABS(BC2(ICHX2)).GT.EPSICP      ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'AC2(ICHX1)', AC2(ICHX1)
            write(NFECRA,*) 'BC2(ICHX1)', BC2(ICHX1)
            write(NFECRA,*) 'AC2(ICHX2)', AC2(ICHX2)
            write(NFECRA,*) 'BC2(ICHX2)', BC2(ICHX2)
          ENDIF
          IF (ABS(SOMM-1.D0).GT.EPSICP) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'SOMM AC2', SOMM
          ENDIF
          IF (ABS(SOMMF).GT.EPSICP) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'SOMM BC2', SOMMF
          ENDIF
          IF (ABS(AC2(IO2)+BC2(IO2)*S2).GT.EPSICP) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'AC12(IO2)', AC1(IO2),AC2(IO2)
            write(NFECRA,*) 'BC12(IO2)', BC1(IO2),BC2(IO2)
          ENDIF
#endif
C ---- Fin Test
C
C --> Integration dans le domaine intermediaire
C
          ANU3 = 0.5D0
          S3   = - (AC2(IO2)-AC2(ICO)*ANU3)
     &         / (BC2(IO2)-BC2(ICO)*ANU3)
C
          SM1  = MAX(SDEB(IEL),S3)
          SM2  = MIN(SFIN(IEL),S2)
C
          DO IGE = 1, (NGAZE-2*NCHARB)
            AI(IGE) = AC2(IGE)
            BI(IGE) = BC2(IGE)
          ENDDO
C
          BI(ICO)   = (AC2(ICO)+BC2(ICO)*S2)/(S2-S3)
          AI(ICO)   = -S3*BI(ICO)
          AI(IO2)   = AC2(IO2 ) - ANU3*(AC2(ICO)-AI(ICO))
          BI(IO2)   = BC2(IO2 ) - ANU3*(BC2(ICO)-BI(ICO))
          AI(ICO2)  = AC2(ICO2) +      (AC2(ICO)-AI(ICO))
          BI(ICO2)  = BC2(ICO2) +      (BC2(ICO)-BI(ICO))
C
          IF ( SM1.LT.SM2 ) THEN
C
            DO IGE = 1, (NGAZE-2*NCHARB)
              ZC(IGE) = ZC(IGE)
     &          + HAUT(IEL)*(SM2-SM1)*(AI(IGE)+.5D0*BI(IGE)*(SM1+SM2))
            ENDDO
C
          ENDIF
C
          DO IGE = 1, (NGAZE-2*NCHARB)
            AC3(IGE) = ( AC2(IGE)*S3-AI(IGE)*SI7(IEL) +
     &                     S3*SI7(IEL)*(BC2(IGE)-BI(IGE)) )
     &                 / ( S3-SI7(IEL) )
            BC3(IGE) = (AI(IGE)+BI(IGE)*S3-AC3(IGE)) / S3
          ENDDO
C
C ---- Apres la troisieme reaction, om se trouve dans le domaine pauvre
C      CHX1,CHX2 et CO doivent etre nuls
C      AC3(ICHX1) = BC3(ICHX1) = 0
C      AC3(ICHX2) = BC3(ICHX2) = 0
C      AC3(ICO  ) = BC3(ICO  ) = 0
C
C ---- Test
#ifdef DEBUG
          SOMM  = 0.D0
          SOMMF = 0.D0
          DO IGE = 1, (NGAZE-2*NCHARB)
            SOMM  = SOMM  + AC3(IGE)*WMOLE(IGE)
            SOMMF = SOMMF + BC3(IGE)*WMOLE(IGE)
          ENDDO
          IF (ABS(AC3(ICHX1)).GT.EPSICP .OR.
     &        ABS(BC3(ICHX1)).GT.EPSICP .OR.
     &        ABS(AC3(ICHX2)).GT.EPSICP .OR.
     &        ABS(BC3(ICHX2)).GT.EPSICP .OR.
     &        ABS(AC3(ICO)).GT.EPSICP .OR.
     &        ABS(BC3(ICO)).GT.EPSICP     ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'AC3(ICHX1)', AC3(ICHX1)
            write(NFECRA,*) 'BC3(ICHX1)', BC3(ICHX1)
            write(NFECRA,*) 'AC3(ICHX2)', AC3(ICHX2)
            write(NFECRA,*) 'BC3(ICHX2)', BC3(ICHX2)
            write(NFECRA,*) 'AC3(ICO)', AC3(ICO)
            write(NFECRA,*) 'BC3(ICO)', BC3(ICO)
          ENDIF
          IF (ABS(SOMM-1.D0).GT.EPSICP) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'SOMM AC3', SOMM
          ENDIF
          IF (ABS(SOMMF).GT.EPSICP) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'SOMM BC3', SOMMF
          ENDIF
          IF (ABS(AC3(IO2)+BC3(IO2)*S3).GT.EPSICP) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'AC23(IO2)', AC2(IO2),AC3(IO2)
            write(NFECRA,*) 'BC23(IO2)', BC2(IO2),BC3(IO2)
          ENDIF
#endif
C ---- Fin Test
C
C --> Integration dans le domaine pauvre
C
          SM1 = SDEB(IEL)
          SM2  = MIN(SFIN(IEL),S3)
C
          DO IGE = 1, (NGAZE-2*NCHARB)
            AI(IGE) = AC3(IGE)
            BI(IGE) = BC3(IGE)
          ENDDO
C
          IF ( SM1.LT.SM2 ) THEN
C
            DO IGE = 1, (NGAZE-2*NCHARB)
              ZC(IGE) = ZC(IGE)
     &          + HAUT(IEL)*(SM2-SM1)*(AI(IGE)+.5D0*BI(IGE)*(SM1+SM2))
            ENDDO
C
          ENDIF
C
          FUEL1(IEL) = ZC(ICHX1) * WMCHX1
          FUEL2(IEL) = ZC(ICHX2) * WMCHX2
          FUEL3(IEL) = ZC(ICO  ) * WMCO
          PROD1(IEL) = ZC(ICO2 ) * WMCO2
          PROD2(IEL) = ZC(IH2O ) * WMH2O
          OXYD(IEL)  = ZC(IO2  ) * WMO2
          XINER(IEL) = ZC(IN2  ) * WMN2
C
C ---- Test
#ifdef DEBUG
          SOMM = FUEL1(IEL) + FUEL2(IEL) + FUEL3(IEL)
     &         + PROD1(IEL) + PROD2(IEL)
     &         + OXYD(IEL)  + XINER(IEL)
          IF ( (ABS(SOMM-1.D0).GT.EPSICP) .OR.
     &         (PROD2(IEL).LT.(-EPZERO) .OR.
     &          PROD2(IEL).GT.(1.D0+EPZERO)) ) THEN
            write(NFECRA,*) 'PB CELLULE ',IEL
            write(NFECRA,*) 'INDPDF3',INDPDF(IEL)
            write(NFECRA,*) 'ACHX1F1',ACHX1F1
            write(NFECRA,*) 'ACHX2F2',ACHX2F2
            write(NFECRA,*) 'ACOF1',ACOF1
            write(NFECRA,*) 'ACOF2',ACOF2
            write(NFECRA,*) 'AO2F4',AO2F4
            write(NFECRA,*) 'AO2F3',AO2F3
            write(NFECRA,*) 'AN2F4',AN2F4
            write(NFECRA,*) 'F1M',F1M(IEL)
            write(NFECRA,*) 'F2M',F2M(IEL)
            write(NFECRA,*) 'F3M',F3M(IEL)
            write(NFECRA,*) 'F4M',F4M(IEL)
            SOMMF = F1M(IEL) + F2M(IEL) + F3M(IEL) + F4M(IEL)
            write(NFECRA,*) 'SOMME F', SOMMF
            write(NFECRA,*) 'SI7  ',SI7(IEL)
            write(NFECRA,*) 'SI8  ',SI8(IEL)
            write(NFECRA,*) 'SDEB  ',SDEB(IEL)
            write(NFECRA,*) 'SFIN  ',SFIN(IEL)
            write(NFECRA,*) 'HAUT  ',HAUT(IEL)
            write(NFECRA,*) 'DSI7  ',DSI7(IEL)
            write(NFECRA,*) 'DSI8  ',DSI8(IEL)
            ANU1 = 0.25D0*(CX1M-CX2M)
            S1   = - (ANU1*AC0(IO2)-AC0(ICHX1))
     &           / (ANU1*BC0(IO2)-BC0(ICHX1))
            write(NFECRA,*) 'S1 ', S1
            ANU2 = 0.25D0*(2.D0 + CX2M)
            S2   = - (ANU2*AC1(IO2)-AC1(ICHX2))
     &           / (ANU2*BC1(IO2)-BC1(ICHX2))
            write(NFECRA,*) 'S2 ', S2
            ANU3 = 0.5D0
            S3   = - (ANU3*AC2(IO2)-AC2(ICO))
     &           / (ANU3*BC2(IO2)-BC2(ICO))
            write(NFECRA,*) 'S3 ', S3
            MOM0 = HAUT(IEL)*(SFIN(IEL)-SDEB(IEL))
     &           + DSI8(IEL)+DSI7(IEL)
            write(NFECRA,*) 'MOM0 ', MOM0, 1.D0
            MOM1 = HAUT(IEL)*(SFIN(IEL)**2-SDEB(IEL)**2)/2.D0
     &           + DSI7(IEL)*SI7(IEL)+DSI8(IEL)*SI8(IEL)
            write(NFECRA,*) 'MOM1 ', MOM1, 0.D0
            MOM2 = HAUT(IEL)*(SFIN(IEL)**3-SDEB(IEL)**3)/3.D0
     &           + DSI7(IEL)*SI7(IEL)**2+DSI8(IEL)*SI8(IEL)**2
            write(NFECRA,*) 'MOM2 ', MOM2, SP2M(IEL)
            write(NFECRA,*) 'Fuel1',FUEL1(IEL)
            write(NFECRA,*) 'Fuel2',FUEL2(IEL)
            write(NFECRA,*) 'Fuel3',FUEL3(IEL)
            write(NFECRA,*) 'Oxyd',OXYD(IEL)
            write(NFECRA,*) 'Iner',XINER(IEL)
            write(NFECRA,*) 'Prod1',PROD1(IEL)
            write(NFECRA,*) 'Prod2',PROD2(IEL)
            write(NFECRA,*) 'SOMME', SOMM
          ENDIF
#endif
C ---- Fin Test
C
        ENDIF
C
      ENDDO
C
C
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
      N11 = 0
      N12 = 0
      N13 = 0
      N14 = 0
      N15 = 0
      N16 = 0
      N17 = 0
      N18 = 0
      N19 = 0
      N20 = 0
      N21 = 0
      N22 = 0
      N23 = 0
      N24 = 0
      N25 = 0
      N26 = 0
      N27 = 0
      N28 = 0
      N29 = 0
      N30 = 0
      N31 = 0
      N32 = 0
C
C --> Controle des differentes valeurs des fractions massiques
C
      DO IEL = 1, NCEL
        SOMM = FUEL1(IEL) + FUEL2(IEL) + FUEL3(IEL) + OXYD(IEL)
     &       + PROD1(IEL) + PROD2(IEL) + XINER(IEL)
        IF ( FUEL1(IEL).LT.(-EPZERO) .OR.
     &       FUEL1(IEL).GT.(1.D0+EPZERO) ) N14 = N14+1
        IF ( FUEL2(IEL).LT.(-EPZERO) .OR.
     &       FUEL2(IEL).GT.(1.D0+EPZERO) ) N15 = N15+1
        IF ( FUEL3(IEL).LT.(-EPZERO) .OR.
     &       FUEL3(IEL).GT.(1.D0+EPZERO) ) N16 = N16+1
        IF ( OXYD(IEL).LT.(-EPZERO) .OR.
     &       OXYD(IEL).GT.(1.D0+EPZERO)  ) N17 = N17+1
        IF ( XINER(IEL).LT.(-EPZERO) .OR.
     &       XINER(IEL).GT.(1.D0+EPZERO) ) N18 = N18+1
        IF ( PROD1(IEL).LT.(-EPZERO) .OR.
     &       PROD1(IEL).GT.(1.D0+EPZERO) ) N19 = N19+1
        IF ( PROD2(IEL).LT.(-EPZERO) .OR.
     &       PROD2(IEL).GT.(1.D0+EPZERO) ) N20 = N20+1
C
        IF ( ABS(SOMM-1.D0).LT.EPSICP )    N3  = N3 +1
      ENDDO
C
C --> Controle des parametres de la pdf
C
      DO IEL = 1, NCEL
        IF ( INDPDF(IEL).NE.0 ) THEN
          VARS = SP2M(IEL) / (-SI7(IEL)*SI8(IEL))
          IF ( VARS.GT.EPZERO .AND.
     &      VARS.LE.(1.-EPZERO) ) N2 = N2+1
          IF ( DSI7(IEL).LT.EPZERO .AND. DSI8(IEL).LT.EPZERO )
     &      N4 = N4+1
          IF ( DSI7(IEL).GT.EPZERO .AND. DSI8(IEL).LT.EPZERO )
     &      N5 = N5+1
          IF ( DSI7(IEL).LT.EPZERO .AND. DSI8(IEL).GT.EPZERO )
     &      N6 = N6+1
          IF ( DSI7(IEL).GT.EPZERO .AND. DSI8(IEL).GT.EPZERO )
     &      N7 = N7+1
          IF ( DSI7(IEL).LT.(-EPZERO) )   N8 = N8+1
          IF ( DSI8(IEL).LT.(-EPZERO) )   N9 = N9+1
          IF ( HAUT(IEL).LT.(-EPZERO) )  N10 = N10+1
          N27 = N27+1
          IF ( SDEB(IEL).LT.SI7(IEL) .OR. SDEB(IEL).GT.SI8(IEL) )
     &      N11 = N11+1
          IF ( SFIN(IEL).LT.SI7(IEL) .OR. SFIN(IEL).GT.SI8(IEL) )
     &      N12 = N12+1
          IF ( (SFIN(IEL)-SDEB(IEL)).LT.ZERO ) N13 = N13+1
          IF ( INDPDF(IEL).EQ.12 ) N28 = N28+1
          IF ( INDPDF(IEL).EQ.2  ) N29 = N29+1
          IF ( INDPDF(IEL).EQ.13 ) N31 = N31+1
          IF ( INDPDF(IEL).EQ.11 ) N30 = N30+1
          IF ( INDPDF(IEL).EQ.3  ) N32 = N32+1
          IF ( IPI7I8(IEL).GT.10 .AND. IPI7I8(IEL).LT.15 ) N21 = N21+1
          IF ( IPI7I8(IEL).GT.20 .AND. IPI7I8(IEL).LT.25 ) N22 = N22+1
          IF ( IPI7I8(IEL).EQ.11 .OR.  IPI7I8(IEL).EQ.21 ) N23 = N23+1
          IF ( IPI7I8(IEL).EQ.12 .OR.  IPI7I8(IEL).EQ.22 ) N24 = N24+1
          IF ( IPI7I8(IEL).EQ.13 .OR.  IPI7I8(IEL).EQ.23 ) N25 = N25+1
          IF ( IPI7I8(IEL).EQ.14 .OR.  IPI7I8(IEL).EQ.24 ) N26 = N26+1
        ENDIF
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
        INTTMP(11) = N11
        INTTMP(12) = N12
        INTTMP(13) = N13
        INTTMP(14) = N14
        INTTMP(15) = N15
        INTTMP(16) = N16
        INTTMP(17) = N17
        INTTMP(18) = N18
        INTTMP(19) = N19
        INTTMP(20) = N20
        INTTMP(21) = N21
        INTTMP(22) = N22
        INTTMP(23) = N23
        INTTMP(24) = N24
        INTTMP(25) = N25
        INTTMP(26) = N26
        INTTMP(27) = N27
        INTTMP(28) = N28
        INTTMP(29) = N29
        INTTMP(30) = N30
        INTTMP(31) = N31
        INTTMP(32) = N32
        CALL PARISM(NBRINT,INTTMP)
C       ===========
      ENDIF
C
      WRITE(NFECRA,1000) N1 , N27
      WRITE(NFECRA,1100) N28, N29, N31, N30, N32, N21, N22, N23,
     &                   N24, N25, N26
      WRITE(NFECRA,2000) N4 , N5 , N6 , N7
      WRITE(NFECRA,2100) N27, N2 , N8 , N9 , N10, N11, N12, N13
      WRITE(NFECRA,2200) N1 , N3 , N14, N15, N16, N17, N18, N19, N20
C
C-------
C FORMAT
C-------
C
 1000 FORMAT (/,
     &'MODELISATION DE LA COMBUSTION AVEC LE MODELE DE DIFFUSION ',
     &'TURBULENTE (CPLYM1)',/,
     &'CHIMIE RAPIDE A 3 CORPS - EXTENSION A 3 COMBUSTIBLES ',
     &'(Application au CP)',/,
     &'==========================================================',
     &'==================',/,
     &' Nb de points de calculs                                     = ',
     &   I9,/,
     &' Nb de points turbulents (passage par les PDF)               = ',
     &   I9)
C
 1100 FORMAT(
C ..v.7..1....v....2....v....3....v....4....v....5....v....6....v....7.I
     &' Nb de points turbulents pour lesquels I7 app. [I4,L3] T12   = ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I7 app. [I4,L5] T2    = ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I7 app. [I4,L5] T13   = ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I7 app. [L5,I3max] T11= ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I7 = I4  T3           = ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I7 domaine P     = ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I7 domaine I     = ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I8 domaine P     = ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I8 domaine I     = ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I8 domaine R     = ',
     &   I9,/,
     &' - - - - - - - - - - - - pour lesquels I8 domaine HR    = ',
     &   I9)
C
 2000 FORMAT(
     &'PDF CONJOINTE DEGENEREE EN PDF MONODIMENSIONNELLE',/,
     &' Nb de points PDF rectangle sans Dirac            = ',I9,/,
     &' - - - - - - - - - - - - -  et Dirac en I7        = ',I9,/,
     &' - - - - - - - - - - - - -  et Dirac en I8        = ',I9,/,
     &' - - - - - - - - - - - - -  et Diracs en I8 et I7 = ',I9)
C
 2100 FORMAT(
     &'CONTROLE DES PARAMETRES DE LA PDF',/,
     &' Nb de points turbulents (passage par les PDF)              = ',
     &   I9,/,
     &' Nb de points turbulents pour lesquels VARS est physique    = ',
     &   I9,/
     &' Nb de points Dirac_I7, Dirac_I8, HAUT < 0  = ',I9,I9,I9,/,
     &' Nb de points Abs.Deb < SI7 ou > SI8        = ', I9,/,
     &' Nb de points Abs.Fin < SI7 ou > SI8        = ', I9,/,
     &' Nb de points (Abs.Fin-Abs.Deb) < 0         = ', I9)
C
 2200 FORMAT(
     &'CONTROLE DES VALEURS DES FRACTIONS MASSIQUES',/,
     &' Nb de points de calculs                                    = ',
     &   I9,/,
     &' Nb de points de calculs qui respectent somme des Yi = 1    = ',
     &   I9,/,
     &' Nb de points YCHx1m,YCHx2m,YCO < 0 ou > 1 = ',I9,I9,I9,/,
     &' Nb de points YO2,YN2        < 0 ou > 1    = ',I9,I9,/,
     &' Nb de points YCO2,YH2O      < 0 ou > 1    = ',I9,I9)
C
C
      RETURN
      END
c@z
