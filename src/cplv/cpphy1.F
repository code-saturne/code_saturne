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
                        SUBROUTINE CPPHY1
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCELET , NCEL   ,
     &   NITBCP , NRTBCP , NITBMC , NRTBMC , NITBWO , NRTBWO ,
     &   F1M    , F2M    , F3M    , F4M    , F5M    ,
     &   F6M    , F7M    , F4P2M  , F3MAX  ,
     &   ENTH   ,
     &   RTP    , PROPCE , ROM1   ,
     &   ITBCP  , RTBCP  ,
     &   ITBMC  , RTBMC  ,
     &   ITBWO  , RTBWO   )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE GAZEUSE
CFONC   VALEURS CELLULES
CFONC   ----------------
CFONC   TEMPERATURE, MASSE VOLUMIQUE ET CONCENTRATIONS MOYENNES
CFONC   (UTILISATION D'UNE PDF RECTANGLE-DIRAC)
CFONC
CFONC  ==> CHIMIE RAPIDE MODELE EN 3 POINTS
CFONC      EXTENSION A TROIS COMBUSTIBLES POUR LE CHARBON PULVERISE
CFONC                                          --------------------
CFONC
CFONC  REACTIONS HETEROGENES
CFONC    - Pyrolyse
CFONC      Composition elementaire de la mole de matieres volatiles
CFONC      Le charbon reactif s'ecrit C(1)H(ALPHA)O(BETA)
CFONC
CFONC        -(k1)-> ALPHA/4 CH4  + BETA CO + (1-ALPHA/4-BETA)    Coke
CFONC      Charbon reactif
CFONC        -(k2)-> ALPHA/Y CXHY + BETA CO + (1-ALPHA/RYSX-BETA) Coke
CFONC
CFONC        Avec RYSX = Y/X
CFONC
CFONC    - Combustion heterogene
CFONC
CFONC      Coke + 1/2 (O2 + XSI N2) -> CO + XSI/2 N2
CFONC
CFONC    - Reactions en phase gaz
CFONC
CFONC (4/(4-RYSX)) CH4 + (O2 + XSI N2)   -(1)->  4/X/(4-RYSX)*CXHY + 2 H2O
CFONC                                            + XSI N2
CFONC CXHY + X/4*(2+RYSX) (O2 + XSI N2)  -(2)->  X CO + Y/2 H2O
CFONC                                            + X/4*(2+RYSX)*XSI N2
CFONC            CO + 1/2 (O2 + XSI N2)  -(3)->  CO2 + XSI/2 N2
CFONC
CFONC CHOIX DES VARIABLES
CFONC
CFONC   F1 est la fractions massique des matieres volatiles : CH4  + CO
CFONC   F2 est la fractions massique des matieres volatiles : CXHY + CO
CFONC   F3 est la fraction massique de carbone venant de la combustion
CFONC     heterogene
CFONC
CFONC   Soit Y les fractions massiques et Z les concentrations (moles/kg)
CFONC     indice f avant reaction, b final
CFONC
CFONC PDF CONJOINTE DEGENERE EN UNE PDF 1D DE TYPE RECTANGLE - DIRAC
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NITBCP       ! E  !  -> ! TAILLE DU MACRO TABLEAU CP ENTIERS   !
CARGU ! NRTBCP       ! E  !  -> ! TAILLE DU MACRO TABLEAU CP REELS     !
CARGU ! NITBMC       ! E  !  -> ! TAILLE DU MACRO TABLEAU MC ENTIERS   !
CARGU ! NRTBMC       ! E  !  -> ! TAILLE DU MACRO TABLEAU MC REELS     !
CARGU ! NITBWO       ! E  !  -> ! TAILLE DU MACRO TABLEAU WORK ENTIERS !
CARGU ! NRTBWO       ! E  !  -> ! TAILLE DU MACRO TABLEAU WORK REELS   !
CARGU ! PA           ! TR !  -> ! PRESSION ABSOLUE EN PASCALS          !
CARGU ! F1M          ! TR !  -> ! MOYENNE DU TRACEUR 1 MVl [CHx1+CO]   !
CARGU ! F2M          ! TR !  -> ! MOYENNE DU TRACEUR 2 MVL [CHx2+CO]   !
CARGU ! F3M          ! TR !  -> ! MOYENNE DU TRACEUR 3 (CO C.HET)      !
CARGU ! F4M          ! TR !  -> ! MOYENNE DU TRACEUR 4 (AIR)           !
CARGU ! F4M          ! TR !  -> ! MOYENNE DU TRACEUR 5 (H2O)           !
CARGU ! F4P2M        ! TR !  -> ! VARIANCE DU TRACEUR 4 (AIR)          !
CARGU ! ENTH         ! TR !  -> ! ENTHALPIE EN J/KG  SOIT DU GAZ       !
CARGU !              !    !     !                    SOIT DU MELANGE   !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! ITBCP        ! TR ! <-> ! MACRO TABLEAU ENTIER CP TRAVAIL      !
CARGU ! RTBCP        ! TR ! <-> ! MACRO TABLEAU REEL   CP TRAVAIL      !
CARGU ! ITBMC        ! TR ! <-> ! MACRO TABLEAU ENTIER MC TRAVAIL      !
CARGU ! RTBMC        ! TR ! <-> ! MACRO TABLEAU REEL   MC TRAVAIL      !
CARGU ! ITBWO        ! TR ! <-> ! MACRO TABLEAU ENTIER TRAVAIL         !
CARGU ! RTBWO        ! TR ! <-> ! MACRO TABLEAU REEL   TRAVAIL         !
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
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
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
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NCELET , NCEL
      INTEGER          NITBCP , NRTBCP
      INTEGER          NITBMC , NRTBMC
      INTEGER          NITBWO , NRTBWO
      INTEGER          ITBCP(NCELET,NITBCP)
      INTEGER          ITBMC(NCELET,NITBMC)
      INTEGER          ITBWO(NCELET,NITBWO)
C
      DOUBLE PRECISION F1M(NCELET), F2M(NCELET)
      DOUBLE PRECISION F3M(NCELET), F4M(NCELET) , F5M(NCELET)
      DOUBLE PRECISION F6M(NCELET), F7M(NCELET)
      DOUBLE PRECISION F4P2M(NCELET)
      DOUBLE PRECISION F3MAX(NCELET)
      DOUBLE PRECISION ENTH(NCELET)
      DOUBLE PRECISION RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION ROM1(NCELET)
      DOUBLE PRECISION RTBCP(NCELET,NRTBCP)
      DOUBLE PRECISION RTBMC(NCELET,NRTBMC)
      DOUBLE PRECISION RTBWO(NCELET,NRTBWO)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL    , IPHAS  , ICE
      INTEGER          IITBCP , IITBMC , IITBWO
      INTEGER          IPCTE1
      INTEGER          IPCYF1 , IPCYF2 , IPCYF3 , IPCYOX
      INTEGER          IPCYP1 , IPCYP2 , IPCYIN , IPCYCE
      INTEGER          IPDF   , NBF3
C
      DOUBLE PRECISION WMOLME , WMCHX1 , WMCHX2
      DOUBLE PRECISION F1CL , F2CL , F3CL , F4CL , F5CL
      DOUBLE PRECISION F6CL , F7CL
      DOUBLE PRECISION FP1  , FP2  , FP4  , FP5, FP6, FP7
      DOUBLE PRECISION XALPHA , XBETA , XGAM
      DOUBLE PRECISION F4CMIN,F4CMAX,XDEN
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
C --- Initialisation memoire
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C --- Initialisation des tableaux d'entiers de travail
C
      DO IEL = 1, NCEL
        DO IITBCP = 1, NITBCP
          ITBCP(IEL,IITBCP) = 0
        ENDDO
        DO IITBMC = 1, NITBMC
          ITBMC(IEL,IITBMC) = 0
        ENDDO
        DO IITBWO = 1, NITBWO
          ITBWO(IEL,IITBWO) = 0
        ENDDO
      ENDDO
C
C --- Initialisation des tableaux de reels de travail
C
      DO IEL = 1, NCEL
        DO IITBCP = 1, NRTBCP
          RTBCP(IEL,IITBCP) = ZERO
        ENDDO
        DO IITBMC = 1, NRTBMC
          RTBMC(IEL,IITBMC) = ZERO
        ENDDO
        DO IITBWO = 1, NRTBWO
          RTBWO(IEL,IITBWO) = ZERO
        ENDDO
      ENDDO
C
C=======================================================================
C 2. DETERMINATION DU TYPE DE PDF
C=======================================================================
C
      NBF3 = 0
      DO IEL = 1, NCEL
C
C       Calcul de Fcl
C
        IF ( F1M(IEL)+F2M(IEL) .GT. 0.D0 ) THEN
          FP1 = F1M(IEL)/(F1M(IEL)+F2M(IEL))
          FP2 = F2M(IEL)/(F1M(IEL)+F2M(IEL))
        ELSE
          FP1 = 0.D0
          FP2 = 0.D0
        ENDIF
        IF (F4M(IEL)+F5M(IEL)+F6M(IEL)+F7M(IEL) .GT. 0.D0 ) THEN
          FP4 = F4M(IEL)/(F4M(IEL)+F5M(IEL)+F6M(IEL)+F7M(IEL))
          FP5 = F5M(IEL)/(F4M(IEL)+F5M(IEL)+F6M(IEL)+F7M(IEL))
          FP6 = F6M(IEL)/(F4M(IEL)+F5M(IEL)+F6M(IEL)+F7M(IEL))
          FP7 = F7M(IEL)/(F4M(IEL)+F5M(IEL)+F6M(IEL)+F7M(IEL))
        ELSE
          FP4 = 0.D0
          FP5 = 0.D0
          FP6 = 0.D0
          FP7 = 0.D0
        ENDIF
C
C      Calcul de F3MAX
C
        IF ( ( AO2F4*F4M(IEL)+AO2F6*F6M(IEL)
     &        +AO2F7*F7M(IEL) ) .GT. 0.D0     ) THEN
          F3MAX(IEL) = 1.D0
     &     /(1.D0+0.5D0*(1.D0/WMOLAT(IATC))
     &                 *(F4M(IEL)+F5M(IEL)+F6M(IEL)+F7M(IEL))
     &                 /(AO2F4*F4M(IEL)+AO2F6*F6M(IEL)+AO2F7*F7M(IEL)) )
        ELSE
          F3MAX(IEL) = 1.D0
        ENDIF
C
        IF ( F3M(IEL) .GT. F3MAX(IEL) ) THEN
          F3M(IEL) = F3MAX(IEL)
          NBF3 = NBF3 + 1
        ENDIF
C
        XDEN = F3M(IEL)/F3MAX(IEL)+F1M(IEL)+F2M(IEL)
C
        IF ( XDEN .GT. 1.D-15 ) THEN
C
          XALPHA = 1.D0-1.D0/XDEN
C
          F1CL = (1.D0-XALPHA)*F1M(IEL)
          F2CL = (1.D0-XALPHA)*F2M(IEL)
          F3CL = (1.D0-XALPHA)*F3M(IEL)
          F4CL = XALPHA*FP4 + (1.D0-XALPHA)*F4M(IEL)
          F5CL = XALPHA*FP5 + (1.D0-XALPHA)*F5M(IEL)
          F6CL = XALPHA*FP6 + (1.D0-XALPHA)*F6M(IEL)
          F7CL = XALPHA*FP7 + (1.D0-XALPHA)*F7M(IEL)
C
        ELSE
          F1CL = FP1
          F2CL = FP2
          F3CL = 0.D0
          F4CL = 0.D0
          F5CL = 0.D0
          F6CL = 0.D0
          F7CL = 0.D0
        ENDIF
C
        RTBCP(IEL,1) = F1CL
        RTBCP(IEL,2) = F2CL
        RTBCP(IEL,3) = F3CL
        RTBCP(IEL,4) = F4CL
        RTBCP(IEL,5) = F5CL
        RTBCP(IEL,6) = F6CL
        RTBCP(IEL,7) = F7CL
C
C       bornes min et max de la pdf : F4CL a 1
C
        RTBCP(IEL,8) = 1.D0
C
C Somme de F4+F5
C
        RTBCP(IEL,14) = F4M(IEL)+F5M(IEL)+F6M(IEL)+F7M(IEL)
C
      ENDDO
C
      IF ( IRANGP .GE. 0 ) THEN
        CALL PARCPT(NBF3)
      ENDIF
      WRITE(NFECRA,*) ' Nombre de clipping sur F3 : ',NBF3
C
      CALL PPPDFR
C     ===========
     & ( NCELET,NCEL,
     &   ITBCP(1,1) ,
     &   RTBCP(1,14), F4P2M ,
C          F4+F5
     &   RTBCP(1,4), RTBCP(1,8),
C           FMINI        FMAXI
     &   RTBCP(1,9) , RTBCP(1,10) , RTBCP(1,11) , RTBCP(1,12),
C           D4CL         D4F4          F4M1        F4M2
     &    RTBCP(1,13) )
C           HREC
C
C=======================================================================
C 2.CALCUL DES CONCENTRATIONS MOYENNES
C=======================================================================
C
C
      IPCYF1 = IPPROC(IYM1(ICHX1))
      IPCYF2 = IPPROC(IYM1(ICHX2))
      IPCYF3 = IPPROC(IYM1(ICO  ))
      IPCYOX = IPPROC(IYM1(IO2  ))
      IPCYP1 = IPPROC(IYM1(ICO2 ))
      IPCYP2 = IPPROC(IYM1(IH2O ))
      IPCYIN = IPPROC(IYM1(IN2  ))
C
      CALL CPCYM2
C     ===========
     & ( NCELET , NCEL   , NRTBMC ,
     &   ITBCP(1,1) ,
C         INTPDF
     &   RTP    ,
     &   F1M    , F2M , F3M , F4M , F5M , F6M , F7M ,
     &   F3MAX  ,
     &   RTBCP(1,1),RTBCP(1,2),RTBCP(1,3),RTBCP(1,4),RTBCP(1,5),
C         F1CL        F2CL       F3CL      F4CL        F5CL
     &  RTBCP(1,6) ,  RTBCP(1,7) ,
C         F6CL         F7CL
     &   RTBCP(1,11) , RTBCP(1,12) , RTBCP(1,9) ,
C           F4M1         F4M2        D4CL
     &   RTBCP(1,10) ,RTBCP(1,13) ,
C           D4F4         HREC
     &   RTBMC      , RTBWO(1,1)  ,
     &   PROPCE(1,IPCYF1) , PROPCE(1,IPCYF2) , PROPCE(1,IPCYF3) ,
     &   PROPCE(1,IPCYOX) , PROPCE(1,IPCYP1) , PROPCE(1,IPCYP2) ,
     &   PROPCE(1,IPCYIN)  )
C
C --> Clipping eventuel des fractions massiques
C
      DO IEL = 1, NCEL
        DO ICE = 1, NGAZE
          IPCYCE = IPPROC(IYM1(ICE))
          IF ( ABS(PROPCE(IEL,IPCYCE)) .LT. EPSICP )
     &         PROPCE(IEL,IPCYCE) = ZERO
        ENDDO
      ENDDO
C
C=======================================================================
C 4. CALCUL DE LA TEMPERATURE ET DE LA MASSE VOLUMIQUE
C=======================================================================
C
      IPCTE1 = IPPROC(ITEMP1)
C
C --- Transport d'H2
C
      CALL CPTEH1
C     ===========
     & ( NCELET , NCEL   , NITBMC , NRTBMC ,
     &   ENTH,
     &   PROPCE(1,IPCYF1), PROPCE(1,IPCYF2), PROPCE(1,IPCYF3),
     &   PROPCE(1,IPCYOX), PROPCE(1,IPCYP1), PROPCE(1,IPCYP2),
     &   PROPCE(1,IPCYIN),
     &   PROPCE(1,IPCTE1),
     &   ITBMC      , RTBMC      ,
C          MACRO TABLEAU MULTI CHARBONS ENTIERS REELS
     &   RTBWO(1,1) , RTBWO(1,2) )
C          TABLEAUX DE TRAVAIL
C
      IPHAS  = 1
      IPCTE1 = IPPROC(ITEMP1)
      DO IEL = 1, NCEL
        WMCHX1 = WMOLAT(IATC)+RTBMC(IEL,IX1MC)*WMOLAT(IATH)
        WMCHX2 = WMOLAT(IATC)+RTBMC(IEL,IX2MC)*WMOLAT(IATH)
        WMOLME = PROPCE(IEL,IPCYF1)/WMCHX1
     &         + PROPCE(IEL,IPCYF2)/WMCHX2
     &         + PROPCE(IEL,IPCYF3)/WMOLE(ICO )
     &         + PROPCE(IEL,IPCYOX)/WMOLE(IO2 )
     &         + PROPCE(IEL,IPCYP1)/WMOLE(ICO2)
     &         + PROPCE(IEL,IPCYP2)/WMOLE(IH2O)
     &         + PROPCE(IEL,IPCYIN)/WMOLE(IN2 )
C
C stockage de la masse molaire du melange
C
        PROPCE(IEL,IPPROC(IMMEL)) = 1.D0 / WMOLME
C
C ---- On ne met pas la pression mecanique RTP(IEL,IPR(IPHAS))
C      mais P0(IPHAS)
C
        ROM1(IEL) = P0(IPHAS) / (WMOLME*RR*PROPCE(IEL,IPCTE1))
      ENDDO
C
C=======================================================================
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
