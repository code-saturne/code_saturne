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
                        SUBROUTINE FUPHY1
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCELET , NCEL   ,
     &   NITBFU , NRTBFU , NITBWO , NRTBWO ,
     &   FVAP   , FHTF   , F4P2M  ,
     &   ENTH   ,
     &   RTP    , PROPCE , ROM1   ,
     &   ITBFU  , RTBFU  ,
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
CFONC    - Evaporation
CFONC      Composition de la vapeur (FOV pour Fuel Oil Vapor)
CFONC
CFONC         Le FOV est supposé être un mélange de H2S, CO, CHn
C          Les fractions massiques sont HSFOV pour H2S
C                                       COFOV      CO
C                                       CHFOV      CHn
C          l'hydrocarbure moyen est détermine par nHCFOV
CFONC
CFONC    - Combustion heterogene
CFONC      La composition massque élémentairee du coke est donnée par
C          CKF, HKF, OKF, SKF
C           (et InKF inertes qui resteront dans l'inclusion)
CFONC       lors de la réaction héterognène, on dégaze H2S, H2O, CO
CFONC
C           Attention, ceci signifie qu'en presence de FHET il y a
C           eut prélevement d'O2 dans l'air environnant
C           (avant les réaction homogènes).
C
CFONC    - Reactions en phase gaz
CFONC
CFONC      Avec l'O2 restant dans l'air (après dilution et oxydation hétérogène)
CFONC      on considère des réaction s successives dans leur ordre
C          de priorité pour l'accès à l'O2
CFONC            CHn + (1/2+n/4)O2 -(1)-> CO  + n/2 H2O
C                H2S + 3/2 O2      -(2)-> SO2 + H2O
CFONC            CO + 1/2 O2       -(3)-> CO2
CFONC
CFONC CHOIX DES VARIABLES
CFONC
CFONC
CFONC   Soit Y les fractions massiques et Z les concentrations (moles/kg)
CFONC     indice f avant reaction, b final
CFONC
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
CARGU ! NITBFU       ! E  !  -> ! TAILLE DU MACRO TABLEAU FUEL ENTIERS !
CARGU ! NRTBFU       ! E  !  -> ! TAILLE DU MACRO TABLEAU FUEL REELS   !
CARGU ! NITBWO       ! E  !  -> ! TAILLE DU MACRO TABLEAU WORK ENTIERS !
CARGU ! NRTBWO       ! E  !  -> ! TAILLE DU MACRO TABLEAU WORK REELS   !
CARGU ! PA           ! TR !  -> ! PRESSION ABSOLUE EN PASCALS          !
CARGU ! FVAP         ! TR !  -> ! MOYENNE DU TRACEUR 1 FOV [CHn+CO]    !
CARGU ! FHTF         ! TR !  -> ! MOYENNE DU TRACEUR 3 (CO C.HET)      !
CARGU ! F4P2M        ! TR !  -> ! VARIANCE DU TRACEUR 4 (AIR)          !
CARGU ! ENTH         ! TR !  -> ! ENTHALPIE EN J/KG  SOIT DU GAZ       !
CARGU !              !    !     !                    SOIT DU MELANGE   !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! ITBFU        ! TR ! <-> ! MACRO TABLEAU ENTIER FUEL TRAVAIL    !
CARGU ! RTBFU        ! TR ! <-> ! MACRO TABLEAU REEL   FUEL TRAVAIL    !
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
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NCELET , NCEL
      INTEGER          NITBFU , NRTBFU
      INTEGER          NITBWO , NRTBWO
      INTEGER          ITBFU(NCELET,NITBFU)
      INTEGER          ITBWO(NCELET,NITBWO)
C
      DOUBLE PRECISION FVAP(NCELET), FHTF(NCELET)
      DOUBLE PRECISION F4P2M(NCELET), ENTH(NCELET)
      DOUBLE PRECISION RTP(NCELET,*), PROPCE(NCELET,*)
      DOUBLE PRECISION ROM1(NCELET)
      DOUBLE PRECISION RTBFU(NCELET,NRTBFU)
      DOUBLE PRECISION RTBWO(NCELET,NRTBWO)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
      INTEGER          IEL    , IPHAS  , ICE
      INTEGER          IPCTE1
      INTEGER          IPCYF1 , IPCYF2 , IPCYF3 , IPCYOX
      INTEGER          IPCYP1 , IPCYP2 , IPCYIN , IPCYCE
      INTEGER          IPCY2S , IPCYSO
      DOUBLE PRECISION WMOLME
      DOUBLE PRECISION F1M,F3M,F4M,F1CL,F3CL,F4CL
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
C=======================================================================
C 2. DETERMINATION DU TYPE DE PDF
C=======================================================================
C
      DO IEL = 1, NCEL
C
C       Traceur virtuel au point moyen
C
        F1m =  FVAP(IEL)
        F3m =  FHTF(IEL)/FF3MAX
        F4m = 1.D0 - F1M - F3M
C
C       Calcul des caractéristiques du point correspondant au
C       combustible moyen
C       F3cl : fraction de masse provenant de F3max et non de F3
C
        F1CL = F1M*FF3MAX/(F3M+F1M*FF3MAX)
        F3CL = 1.D0-F1CL
        F4CL = (1.D0-FF3MAX)*F3CL
C
        RTBFU(IEL,1) = F3M
        RTBFU(IEL,2) = F4M
        RTBFU(IEL,3) = F1CL
        RTBFU(IEL,4) = F3CL
        RTBFU(IEL,5) = F4CL
C
C       bornes min et max de la pdf : F4CL a 1
C
        RTBFU(IEL,6) = 1.D0
C
      ENDDO
C
      CALL PPPDFR
C     ===========
     & ( NCELET,NCEL,
     &   ITBFU(1,1) ,
     &   RTBFU(1,2), RTP(1,ISCA(IF4P2M)),
C           F4M
     &   RTBFU(1,5), RTBFU(1,6),
C           FMINI        FMAXI
     &   RTBFU(1,7) , RTBFU(1,8) , RTBFU(1,9) , RTBFU(1,10),
C           D4CL         D4F4          F4M1        F4M2
     &    RTBFU(1,11) )
C           HREC
C
C=======================================================================
C 2.CALCUL DES CONCENTRATIONS MOYENNES
C=======================================================================
C
C
      IPCYF1 = IPPROC(IYM1(IFOV))
      IPCYF3 = IPPROC(IYM1(ICO  ))
      IPCYOX = IPPROC(IYM1(IO2  ))
      IPCYP1 = IPPROC(IYM1(ICO2 ))
      IPCYP2 = IPPROC(IYM1(IH2O ))
      IPCYIN = IPPROC(IYM1(IN2  ))
      IPCY2S = IPPROC(IYM1(IH2S ))
      IPCYSO = IPPROC(IYM1(ISO2 ))
C
       CALL FUCYM1
C      ===========
     & ( NCELET , NCEL   ,
     &   ITBFU(1,1) ,
C         INTPDF
     &   FVAP   ,   RTBFU(1,1) , RTBFU(1,2) ,
C         F1M           F3M         F4M
     &  RTBFU(1,3) , RTBFU(1,4) ,RTBFU(1,5) ,
C         F1CL         F3CL         F4CL

     &   RTBFU(1,9) , RTBFU(1,10) , RTBFU(1,7) ,
C           F4M1         F4M2        D4CL
     &   RTBFU(1,8) ,RTBFU(1,11) ,
C           D4F4         HREC
     &   PROPCE(1,IPCYF1) , PROPCE(1,IPCYF3) ,
     &   PROPCE(1,IPCYOX) , PROPCE(1,IPCYP1) , PROPCE(1,IPCYP2) ,
     &   PROPCE(1,IPCYIN) ,
     &   PROPCE(1,IPCY2S) , PROPCE(1,IPCYSO) )
C
C --> Clipping eventuel des fractions massiques
C
      DO IEL = 1, NCEL
        DO ICE = 1, NGAZE
          IPCYCE = IPPROC(IYM1(ICE))
          IF ( ABS(PROPCE(IEL,IPCYCE)) .LT. EPSIFL )
     &         PROPCE(IEL,IPCYCE) = ZERO
        ENDDO
      ENDDO
C
C
C=======================================================================
C 3. CALCUL DE LA TEMPERATURE ET DE LA MASSE VOLUMIQUE
C=======================================================================
C
      IPCTE1 = IPPROC(ITEMP1)
C
      CALL FUTEH1
C     ===========
     & ( NCELET , NCEL   ,
     &   ENTH,
     &   PROPCE(1,IPCYF1), PROPCE(1,IPCYF3),
     &   PROPCE(1,IPCYOX), PROPCE(1,IPCYP1), PROPCE(1,IPCYP2),
     &   PROPCE(1,IPCYIN), PROPCE(1,IPCY2S), PROPCE(1,IPCYSO),
     &   PROPCE(1,IPCTE1),
     &   RTBWO(1,1) , RTBWO(1,2) )
C
C          TABLEAUX DE TRAVAIL
C
      IPHAS  = 1
      IPCTE1 = IPPROC(ITEMP1)
      DO IEL = 1, NCEL
        WMOLME = PROPCE(IEL,IPCYF1) / WMOLE(IFOV)
     &         + PROPCE(IEL,IPCYF3) / WMOLE(ICO )
     &         + PROPCE(IEL,IPCYOX) / WMOLE(IO2 )
     &         + PROPCE(IEL,IPCYP1) / WMOLE(ICO2)
     &         + PROPCE(IEL,IPCYP2) / WMOLE(IH2O)
     &         + PROPCE(IEL,IPCYIN) / WMOLE(IN2 )
     &         + PROPCE(IEL,IPCY2S) / WMOLE(IH2S)
     &         + PROPCE(IEL,IPCYSO) / WMOLE(ISO2)
C
C stockage de la masse molaire du melange
C
        PROPCE(IEL,IPPROC(IMMEL)) = 1.D0 / WMOLME
C
C ---- On ne met pas la pression mecanique RTP(IEL,IPR(IPHAS))
C      mais P0(IPHAS)
C
        ROM1(IEL) = P0(IPHAS) / (WMOLME * RR * PROPCE(IEL,IPCTE1) )
      ENDDO
C
C=======================================================================
C FORMATS
C----
C
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
