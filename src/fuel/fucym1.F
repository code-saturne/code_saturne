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
                        SUBROUTINE FUCYM1
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   INDPDF ,
     &   F1M    , F3M    , F4M   , F1CL  , F3CL , F4CL ,
     &   F4M1   , F4M2   , D4CL  , D4F4  , HREC ,
     &   FUEL1  , FUEL3  , OXYD  , PROD1 , PROD2  ,
     &   XINER  , XH2S   , XSO2   )
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
CARGU ! FUEL1        ! TR !  <- ! FRACTION MASSIQUE FOVm               !
CARGU ! FUEL3        ! TR !  <- ! FRACTION MASSIQUE CO                 !
CARGU ! OXYD         ! TR !  <- ! FRACTION MASSIQUE O2                 !
CARGU ! PROD1        ! TR !  <- ! FRACTION MASSIQUE CO2                !
CARGU ! PROD2        ! TR !  <- ! FRACTION MASSIQUE H2O                !
CARGU ! XINER        ! TR !  <- ! FRACTION MASSIQUE N2                 !
CARGU ! XH2S         ! TR !  <- ! FRACTION MASSIQUE H2S                !
CARGU ! XSO2         ! TR !  <- ! FRACTION MASSIQUE SO2                !
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
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL
      INTEGER          INDPDF(NCELET)
C
      DOUBLE PRECISION F1M(NCELET) , F3M(NCELET)  , F4M(NCELET)
      DOUBLE PRECISION F1CL(NCELET), F3CL(NCELET) , F4CL(NCELET)
C
      DOUBLE PRECISION F4M1(NCELET) , F4M2(NCELET) , D4CL(NCELET)
      DOUBLE PRECISION D4F4(NCELET) , HREC(NCELET)
C
      DOUBLE PRECISION FUEL1(NCELET), FUEL3(NCELET)
      DOUBLE PRECISION OXYD(NCELET) , PROD1(NCELET), PROD2(NCELET)
      DOUBLE PRECISION XINER(NCELET)
      DOUBLE PRECISION XH2S(NCELET) , XSO2(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL
C
      INTEGER          N1 , N2 , N3 , N4 , N5 , N6 , N7 , N8 , N9
      INTEGER          N10, N11
      INTEGER          NBRINT
      PARAMETER        (NBRINT = 11)
      INTEGER          INTTMP(NBRINT)
C
      DOUBLE PRECISION WMFOV , WMCO   , WMCO2 , WMH2O
      DOUBLE PRECISION WMO2   , WMN2   , WMC
      DOUBLE PRECISION WMH2S , WMSO2
C
      DOUBLE PRECISION ANU1 , ANU2 , ANU3
      DOUBLE PRECISION SOMM
C
      DOUBLE PRECISION ZFOV0 , ZCO0   , ZCO20
      DOUBLE PRECISION ZN20  , ZH2O0  , ZO20
      DOUBLE PRECISION YFOV0 , YCO0   , YCO20
      DOUBLE PRECISION YN20  , YH2O0  , YO20
      DOUBLE PRECISION ZFOV1 , ZCO1   , ZCO21
      DOUBLE PRECISION ZN21  , ZH2O1  , ZO21
      DOUBLE PRECISION ZFOV2 , ZCO2   , ZCO22
      DOUBLE PRECISION ZN22  , ZH2O2  , ZO22
      DOUBLE PRECISION ZFOV3 , ZCO3   , ZCO23
      DOUBLE PRECISION ZN23  , ZH2O3  , ZO23
      DOUBLE PRECISION REAC1 , REAC2  , REAC3
C
C     Variables pour le support de la Pdf et sa description
C
      DOUBLE PRECISION F42M
      DOUBLE PRECISION ZZCL(8),ZZS1(8),ZZS2(8),ZZS3(8),ZZF4(8)
      DOUBLE PRECISION ZZ(8)
C
      INTEGER II
      DOUBLE PRECISION BB1,BB2
      DOUBLE PRECISION ZH2S0,ZH2S1,ZH2S2,ZH2S3
      DOUBLE PRECISION ZSO20,ZSO21,ZSO22,ZSO23
      DOUBLE PRECISION F4S1,F4S2,F4S3
C
C***********************************************************************
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
      WMFOV  = WMOLE(IFOV)
      WMH2S  = WMOLE(IH2S)
      WMSO2  = WMOLE(ISO2)
C
      WMC    = WMOLAT(IATC)
C
C
C --- Calcul des coefficient relatifs a l'expression de
C        ZCHx1m0,  ZCO0, ZO20 et ZN20
C        fonctions de F1,  F3 et F4
C
C     Rappel : 1 = vapeur de fuel, 2 = non utilisé
C              3 = C sous forme CO provenant de l'oxydation du coke
C              4 = air
C
C ---- Initialisation des fractions massiques avec de l'air
C
      DO IEL = 1, NCEL
        FUEL1(IEL)  = ZERO
        FUEL3(IEL)  = ZERO
        OXYD(IEL)   = WMO2*AO2F4
        PROD1(IEL)  = ZERO
        PROD2(IEL)  = ZERO
        XINER(IEL)  = WMN2*AN2F4
        XH2S (IEL)  = ZERO
        XSO2 (IEL)  = ZERO
      ENDDO
C
C
C  ANU1 : nb de moles d'oxygène pour convertir
C         une mole de FOV en CO et H2O
C
      ANU1   = 0.5D0+0.25D0*nHCFOV
C
C  ANU2 nb de moles d'oxygènes pour la réaction H2S -> SO2 + H2O
C
      ANU2 = 1.5D0
C
C  ANU3 nb de mole d'oxygene pour la reaction CO -> CO2
C
      ANU3   = 0.5D0
C
C
C=======================================================================
C 3. CALCUL DE LA COMPOSITION DU MELANGE SANS LES PDF
C    SI LES FLUCTUATIONS DE S SONT TROP FAIBLES
C=======================================================================
C
C
      DO IEL = 1, NCEL
C
        IF ( INDPDF(IEL).EQ.0 ) THEN
C
C --> Composition de la phase gazeuse avant combustion (indice 0)
C
C          ZFOV0  = AFOVF1*F1
C          ZCO0   = ACOF1*F1 +  ACOF3*F3
C          ZCO20  = ZERO
C          ZH2S0  = AH2SF1*F1 + AH2SF3*F3
C          ZSO20  = ZERO
C          ZH2O0  = AH2OF3*F3
C          ZO20   = AO2F4*F4 + AO2F3*F3
C          ZN20   = AN2F4*F4
C
C ..v.7..x....v....x....v....x....v....x....v....L....v....x....v....x.I
          ZFOV0  = AFOVF1* F1M(IEL)
          ZCO0   = ACOF1 * F1M(IEL) + ACOF3 * F3M(IEL)
          ZH2S0  = AH2SF1* F1M(IEL) + AH2SF3* F3M(IEL)
          ZO20   = AO2F4 * F4M(IEL) + AO2F3 * F3M(IEL)
          ZN20   = AN2F4 * F4M(IEL) + AN2F3 * F3M(IEL)
          ZCO20  = ZERO
          ZSO20  = ZERO
          ZH2O0  = AH2OF3 * F3M(IEL)
C
C --> Calcul de la composition du melange
C
          REAC1  = MIN(ZFOV0, (ZO20/ANU1))
          ZFOV1 = ZFOV0  -           REAC1
          ZO21   = ZO20  -      ANU1*REAC1
          ZCO1   = ZCO0  +           REAC1
          ZCO21  = ZCO20
          ZH2O1  = ZH2O0 + 0.5D0*nHCFOV*REAC1
          ZN21   = ZN20
          ZH2S1  = ZH2S0
          ZSO21  = ZSO20
C
          REAC2 = MIN(ZH2S1,(ZO21/ANU2))
          ZFOV2 = ZFOV1
          ZCO2 = ZCO1
          ZO22  = ZO21  - ANU2 * REAC2
          ZH2S2 = ZH2S1 -        REAC2
          ZH2O2 = ZH2O1 +        REAC2
          ZSO22 = ZSO21 +        REAC2
          ZCO22 = ZCO21
          ZN22  = ZN21
C
          REAC3 = MIN(ZCO2,(ZO22/ANU3))
          ZFOV3 = ZFOV2
          ZO23  = ZO22  - ANU3*REAC3
          ZCO3  = ZCO2  -      REAC3
          ZCO23 = ZCO22 +      REAC3
          ZH2O3 = ZH2O2
          ZN23  = ZN22
          ZH2S3 = ZH2S2
          ZSO23 = ZSO22
C
          FUEL1(IEL) = ZFOV3 * WMFOV
          FUEL3(IEL) = ZCO3  * WMCO
          PROD1(IEL) = ZCO23 * WMCO2
          PROD2(IEL) = ZH2O3 * WMH2O
          OXYD(IEL)  = ZO23  * WMO2
          XINER(IEL) = ZN23  * WMN2
          XH2S (IEL) = ZH2S3 * WMH2S
          XSO2 (IEL) = ZSO23 * WMSO2
C
        ENDIF
C
      ENDDO
C
C
C
C=======================================================================
C 4. CALCUL DE LA COMPOSITION DU MELANGE AVEC LES PDF
C=======================================================================
C
C  14/08/2006
C  En toute rigueur, on pourrait calculer l'enthalpie de mélange aux 5 pts
C  en déduire la température et intégrer l'inverse de la masse volumique
C  mais ... les flammes fuligineuses sont rarement adiabatiques.
C  Il faudrait donc rechercher l'enthalpie aux points stoechio comme variable
C  auxilliaire pour retrouver l'enthalpie moyenne (qui tient compte des pertes)
C  En mettant les choses au mieux (i.e. on saurait formuler une hypothèse de
C  répartition des pertes entre les 3 points) on en déduirait une température
C  linéaire par morceaux, les fractions mmassiques étant linéaires,
C  1/Rho serait quadratique.
C  Déjà fait une fois (flamme de diff) mais pénible.
C
C  Pour l'instant on conserve le calcul de la température et de la masse
C  volumique d'aprés les concentrations moyennes dans FUPHY1
C  comme pour le cp
C
C
      DO IEL = 1, NCEL
C
        IF ( INDPDF(IEL).EQ.1 ) THEN
C
C         Concentrations aux extrémités
C
          DO II = 1,8
            ZZCL(II) = ZERO
            ZZF4(II) = ZERO
          ENDDO
C
          ZZCL(IFOV) = AFOVF1*F1CL(IEL)
          ZZCL(ICO)  = ACOF1 *F1CL(IEL)+ACOF3 *F3CL(IEL)
          ZZCL(IH2S) = AH2SF1*F1CL(IEL)+AH2SF3*F3CL(IEL)
          ZZCL(IH2O) =                  AH2OF3*F3CL(IEL)
          ZZCL(IN2)  =                  AN2F3 *F3CL(IEL)
          ZZF4(IO2)  = AO2F4
          ZZF4(IN2)  = AN2F4
C
C         Calcul de la stoechiométrie de la première réaction
C         FOV + O2 => CO + H2O
C         FOV et O2 sont, avant la première réaction linéaire en F4
C         ZZ(F4,II) = ZZcl(II) + (F4-F4cl)/(1-F4cl) * (ZZF4(II)-ZZcl(II))
C         On cherche le point tel que
C         ZZ(F4s1,IFOV) = ZZ(F4s1,IO2) / (0.5+0.25*nHCFOV)
C         ANU1*ZZ(F4s1,IFOV) = ZZ(F4s1,IO2)
C         ANU1*(1-F4cl)*ZZcl(IFOV) + ANU1*(F4s1-F4cl)*(ZZF4(IFOV)-ZZcl(IFOV)) =
C              (1-F4cl)*ZZcl(IO2)  +      (F4s1-F4cl)*(ZZF4(IO2) -ZZcl(IO2))
C         on remarque que ZZcl(IO2) = ZZF4(IFOV) = 0
C
          F4S1 = ( ANU1 * ZZCL(IFOV) + F4CL(IEL) * ZZF4(IO2) )
     &         / ( ANU1 * ZZCL(IFOV) +             ZZF4(IO2) )
C
C         Calcul des concentrations au premier point stoechio
C         avant réaction
C
          DO II = 1,8
            ZZS1(II) = ZZCL(II)
     &           + (F4S1-F4CL(IEL))/(1.D0-F4CL(IEL))
     &           *(ZZF4(II)-ZZCL(II))
          ENDDO
C
C         Le nombre de moles de réaction est le nombre de moles de fov
C
          ZZS1(ICO) = ZZS1(ICO) + ZZS1(IFOV)
          ZZS1(IO2) = ZERO
          ZZS1(IH2O)= ZZS1(IH2O) + 0.5D0*NHCFOV*ZZS1(IFOV)
          ZZS1(IFOV)= ZERO
C
C         Calcul de la stoechiometrie de la deuxième réaction
C         H2S + O2 => SO2 + H2O
C         Les concentrations sont désormais linéaires entre F4s1 et 1
C         ZZ(F4,II)  = ZZs1(II) + (F4-F4s1)/(1-F4s1)*(ZZF4(II)-ZZs1(II))
C         On cherche le point tel que
C         ZZ(f4s2,IH2S) = ZZ(F4S2,IO2)/ANU2
C         a nouveau ZZs1(IO2) = ZZF4(IH2S) = 0
C         La formule est semblable à la précedente en substituant
C         ANU1 par ANU2
C         IFOV par IH2S
C         cl par s1
C         s1 par s2
C
          F4S2 = ( ANU2 * ZZS1(IH2S) + F4S1 * ZZF4(IO2) )
     &         /( ANU2 * ZZS1(IH2S) +        ZZF4(IO2) )
C
C     Calcul des concentrations au deuxième point stoechio
C     avant réaction
C
          DO II = 1,8
            ZZS2(II) = ZZS1(II)
     &           + (F4S2-F4S1) / (1.D0-F4S1) * (ZZF4(II) - ZZS1(II))
          ENDDO
C
C         Le nombre de moles de réaction est le nombre de moles de H2S
C
          ZZS2(ISO2) = ZZS2(ISO2) + ZZS2(IH2S)
          ZZS2(IH2O) = ZZS2(IH2O) + ZZS2(IH2S)
          ZZS2(IO2)  = ZERO
          ZZS2(IH2S) = ZERO
C
C         Calcul de la stoechiometrie de la troisième réaction
C         CO + O2 => CO2
C         Les concentrations sont désormais linéaires entre F4S2 et 1
C         ZZ(F4,II)  = ZZS2(II) + (F4-F4S2)/(1-F4S2)*(ZZF4(II)-ZZS2(II))
C         On cherche le point tel que
C         ZZ(f4s3,ICO) = ZZ(F4S3,IO2)/ANU3
C         La formule est semblable à la précedente en substituant
C         ANU2 par ANU3
C         IH2S par ICO
C         s1 par s2
C         s2 par s3
C
          F4S3 = ( ANU3 * ZZS2(ICO ) + F4S2 * ZZF4(IO2) )
     &         / ( ANU3 * ZZS2(ICO ) +        ZZF4(IO2) )
C
C         Calcul des concentrations au troisième point stoechio
C         avant réaction
C
          DO II = 1,8
            ZZS3(II) = ZZS2(II)
     &           + (F4S3-F4S2) / (1.D0-F4S2) * (ZZF4(II) - ZZS2(II))
          ENDDO
C
C         Le nombre de moles de réaction est le nombre de moles de CO
C
          ZZS3(ICO2) = ZZS3(ICO2) + ZZS3(ICO)
          ZZS3(IO2)  = ZERO
          ZZS3(ICO)  = ZERO
C
C        Désormais on connait les concentrations en 5 points
C        cl,s1,s2,s3,f4
C        les concentrations intermédiaires sont linéaires par morceaux
C        et les paramêtres de la pdf D4cl,D4f4, F4m1,F4m2,HREC
C
C        On initialise par les concentrations aux extrémités
          DO II = 1,8
            ZZ(II) = D4CL(IEL)*ZZcl(II)+D4F4(IEL)*ZZF4(II)
          ENDDO
C
C        Intégration sur le premier intervalle de richesse (entre cl et s1)
C
          BB1 = MAX(F4M1(IEL),F4CL(IEL))
          BB2 = MIN(F4M2(IEL),F4S1)
          IF( BB2.GT.BB1 ) THEN
            DO II = 1,8
              ZZ(II) =  ZZ(II) + HREC(IEL)*(BB2-BB1)/(F4S1-F4CL(IEL))
     &             *( (ZZCL(II)*F4S1-ZZS1(II)*F4CL(IEL))
     &             +   (ZZS1(II)-ZZCL(II))*(BB1+BB2)*0.5D0)
            ENDDO
          ENDIF
C
C        Intégration sur le deuxième intervalle de (entre s1 et s2)
C
          BB1 = MAX(F4M1(IEL),F4S1)
          BB2 = MIN(F4M2(IEL),F4S2)
          IF( BB2.GT.BB1 ) THEN
            DO II = 1,8
              ZZ(II) = ZZ(II) + HREC(IEL)*(BB2-BB1)/(F4S2-F4S1)
     &             * ( (ZZS1(II)*F4S2-ZZS2(II)*F4S1)
     &             +   (ZZS2(II)-ZZS1(II))*(BB1+BB2)*0.5D0)
            ENDDO
          ENDIF
C
C        Intégration de s2 à s3
C
          BB1 = MAX(F4M1(IEL),F4S2)
          BB2 = MIN(F4M2(IEL),F4S3)
          IF( BB2.GT.BB1 ) THEN
            DO II = 1,8
              ZZ(II) = ZZ(II) + HREC(IEL)*(BB2-BB1)/(F4S3-F4S2)
     &             * ( (ZZS2(II)*F4S3-ZZS3(II)*F4S2)
     &             +   (ZZS3(II)-ZZS2(II))*(BB1+BB2)*0.5D0)
            ENDDO
          ENDIF
C
C       Intégration de s3 à f4
C
          BB1 = MAX(F4M1(IEL),F4S3)
          BB2 = MIN(F4M2(IEL),1.D0)
          IF ( BB2.GT.BB1 ) THEN
            DO II = 1,8
              ZZ(II) = ZZ(II) + HREC(IEL)*(BB2-BB1)/(1.D0-F4S3)
     &             * ( (ZZS3(II)*1.D0-ZZF4(II)*F4S3)
     &             +(ZZF4(II)-ZZS3(II))*(BB1+BB2)*0.5D0)
            ENDDO
          ENDIF
C
          FUEL1(IEL) = ZZ(IFOV ) * WMFOV
          FUEL3(IEL) = ZZ(ICO  ) * WMCO
          PROD1(IEL) = ZZ(ICO2 ) * WMCO2
          PROD2(IEL) = ZZ(IH2O ) * WMH2O
          OXYD(IEL)  = ZZ(IO2  ) * WMO2
          XINER(IEL) = ZZ(IN2  ) * WMN2
          XH2S(IEL)  = ZZ(IH2S ) * WMH2S
          XSO2(IEL)  = ZZ(ISO2 ) * WMSO2
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
      N11 = 0
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
      DO IEL = 1, NCEL
C
        SOMM = FUEL1(IEL) + FUEL3(IEL) + OXYD(IEL)
     &       + PROD1(IEL) + PROD2(IEL) + XINER(IEL)
     &       + XH2S(IEL) + XSO2(IEL)
C
        IF ( ABS(SOMM-1.D0).LT.EPSICP )    N3  = N3 +1
        IF ( FUEL1(IEL).LT.(-EPZERO) .OR.
     &       FUEL1(IEL).GT.(1.D0+EPZERO) ) N4 = N4+1
        IF ( FUEL3(IEL).LT.(-EPZERO) .OR.
     &       FUEL3(IEL).GT.(1.D0+EPZERO) ) N5 = N5+1
        IF ( OXYD(IEL).LT.(-EPZERO) .OR.
     &       OXYD(IEL).GT.(1.D0+EPZERO)  ) N6 = N6+1
        IF ( XINER(IEL).LT.(-EPZERO) .OR.
     &       XINER(IEL).GT.(1.D0+EPZERO) ) N7 = N7+1
        IF ( PROD1(IEL).LT.(-EPZERO) .OR.
     &       PROD1(IEL).GT.(1.D0+EPZERO) ) N8 = N8+1
        IF ( PROD2(IEL).LT.(-EPZERO) .OR.
     &       PROD2(IEL).GT.(1.D0+EPZERO) ) N9 = N9+1
        IF ( XH2S(IEL).LT.(-EPZERO) .OR.
     &       XH2S(IEL).GT.(1.D0+EPZERO) ) N10 = N10+1
        IF ( XSO2(IEL).LT.(-EPZERO) .OR.
     &       XSO2(IEL).GT.(1.D0+EPZERO) ) N11 = N11+1
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
        INTTMP(11) = N11
        CALL PARISM(NBRINT,INTTMP)
C       ===========
      ENDIF
C
      WRITE(NFECRA,1000) N1 , N2
C
      WRITE(NFECRA,2200) N3 , N4, N5, N6, N7, N8,
     &                        N9, N10, N11
C
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
     &' Nb de points YFOV ,YCO    < 0 ou > 1 = ',I6,I6,/,
     &' Nb de points YO2  ,YN2    < 0 ou > 1 = ',I6,I6,/,
     &' Nb de points YCO2 ,YH2O   < 0 ou > 1 = ',I6,I6,/,
     &' Nb de points YH2S ,YSO2   < 0 ou > 1 = ',I6,I6,/)
C
C
      RETURN
      END
c@z
