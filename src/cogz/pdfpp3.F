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
                         SUBROUTINE PDFPP3
C                        *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL  ,
     &   FM     , FP2M  , YFM    , YFP2M , COYFP  ,
     &   PROPCE )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL DES PARAMETRES DE LA PDF
CFONC  PDF LIBBY - WILLIAMS 3 POINTS AVEC HYPOTHESE DE CURL
CFONC        PDFPP3 RIVISE
CFONC
CFONC  COMMENTAIRES :
CFONC  ------------
CFONC
CFONC     Dans un diagramme (F, Yf), on construit deux droites:
CFONC          - La droite de combustion complete
CFONC          - La droite de melange
CFONC
CFONC     Dans ce domaine, nous allons trouver deux pics qui
CFONC     definiront une troisieme droite sur laquelle on definit
CFONC     une abscisse curviligne G.
CFONC
CFONC
CFONC  LE RESULTAT EST :
CFONC  ---------------
CFONC
CFONC     CALCUL DES PARAMETRES ASSOCIES AUX FONCTIONS DIRAC
CFONC
CFONC       Les Diracs sont en position [F(.,1),Y(.,1)] et [F(.,2),Y(.,2)]
CFONC       Leurs amplitudes respectives sont D(.,1) et D(.,2)
CFONC       Pour chaque dirac,
CFONC           on calcule la temperature [T(.,1), T(.,2)]
CFONC                   la masse volumique [RHO(.,1), RHO(.,2)]
CFONC                   le terme source chimique [W(.,1),W(.,2)]
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! FM           ! TR !  -> ! MOYENNE DE LA FRACTION DE MELANGE    !
CARGU ! FP2M         ! TR !  -> ! VARIANCE DE LA FRACTION DE MELANGE   !
CARGU ! YFM          ! TR !  -> ! MOYENNE DE LA FRACTION MASSIQUE      !
CARGU ! YFP2M        ! TR !  -> ! VARIANCE DE LA FRACTION MASSIQUE     !
CARGU ! COYFP        ! TR !  -> ! COVARIANCE
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
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
      INCLUDE "pointe.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "coincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET, NCEL
C
      DOUBLE PRECISION FM(NCELET)   , FP2M(NCELET)
      DOUBLE PRECISION YFM(NCELET)  , YFP2M(NCELET)
      DOUBLE PRECISION COYFP(NCELET)
      DOUBLE PRECISION PROPCE(NCELET,*)
C
C VARIABLES LOCALES
C
      INTEGER          IEL, IGG, IDIRAC, IPHAS
      INTEGER          MODE
      INTEGER          CLICOY, CLIY1, CLIY2 , CLIY2P
      INTEGER          CLIF  , CLIY , CLIFP2, CLYFP2
C
      DOUBLE PRECISION COEFG(NGAZGM), EPSI
      DOUBLE PRECISION YFUEL
      DOUBLE PRECISION YOXYD, YO2
      DOUBLE PRECISION YPROD,FMP,FP2MP,YFMP,YFP2MP,COYFPP
      DOUBLE PRECISION Y1, Y2, F1, F2
      DOUBLE PRECISION CSTFA1, CSTFA2, YFMPMX, CST
      DOUBLE PRECISION DYFMP
      DOUBLE PRECISION CLIMAX, CLIMIN
C
C --- Position des Diracs
C
      DOUBLE PRECISION F(NDRACM)    , Y(NDRACM)     , D(NDRACM)
      DOUBLE PRECISION H(NDRACM)    , TEML(NDRACM)
      DOUBLE PRECISION MAML(NDRACM) , W(NDRACM)
      DOUBLE PRECISION RHOL(NDRACM) , THETA(NDRACM)
      DOUBLE PRECISION YMIN(NDRACM)  , YMAX(NDRACM)
      DOUBLE PRECISION Y2P(NDRACM)
C
C --- Pointeurs & autres
C
      INTEGER          IPCTEM, IPCMAM,  IPCROM
      INTEGER          IPAMPL(NDRACM), IPFMEL(NDRACM)
      INTEGER          IPFMAL(NDRACM), IPTEML(NDRACM)
      INTEGER          IPMAML(NDRACM)
      INTEGER          IPRHOL(NDRACM)
      INTEGER          IPTSCL(NDRACM)
      INTEGER          IPCFUE, IPCOXY, IPCPRO, IPCTSC
      DOUBLE PRECISION NBMOL,  TEMSMM
      DOUBLE PRECISION SUM7, SUM8, SUM9, SUM10, SUM11, SUM12, SUM17
      DOUBLE PRECISION SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, SUM16, SUM15
C
      INTEGER IPASS
      DATA    IPASS /0/
      SAVE    IPASS
C
C***********************************************************************
C
C ---> Position des variables
C
      DO IDIRAC = 1, NDIRAC
        IPAMPL(IDIRAC) = IPPROC(IAMPL(IDIRAC))
        IPFMEL(IDIRAC) = IPPROC(IFMEL(IDIRAC))
        IPFMAL(IDIRAC) = IPPROC(IFMAL(IDIRAC))
        IPTEML(IDIRAC) = IPPROC(ITEML(IDIRAC))
        IPMAML(IDIRAC) = IPPROC(IMAML(IDIRAC))
        IPRHOL(IDIRAC) = IPPROC(IRHOL(IDIRAC))
        IPTSCL(IDIRAC) = IPPROC(ITSCL(IDIRAC))
      ENDDO
      IPCFUE = IPPROC(IYM(1))
      IPCOXY = IPPROC(IYM(2))
      IPCPRO = IPPROC(IYM(3))
      IPCTSC = IPPROC(ITSC)
      IPCTEM = IPPROC(ITEMP)
      IPHAS = 1
      IPCROM = IPPROC(IROM(IPHAS))
      IPCMAM = IPPROC(IMAM)
C
C ---> Initialisation
C
      DO IGG = 1, NGAZGM
        COEFG(IGG) = ZERO
      ENDDO
      EPSI = 1.D-6
C
C     Compteur pour clipping
C
      CLICOY  = 0
      CLIY1   = 0
      CLIY2   = 0
      CLIY2P  = 0
      CLIF    = 0
      CLIY    = 0
      CLIFP2  = 0
      CLYFP2 = 0
C
C=======================================================================
C 0.  CALCULS PRELIMINAIRES
C=======================================================================
C
      IPASS = IPASS + 1
C
C Controle de l intialisation de FMIN et FMAX
C
      IF ( (IPASS.LE.1 .AND. ISUITE.EQ.0 ) .OR.
     &     (IPASS.LE.1 .AND. ISUITE.EQ.1
     &                 .AND. INITRO(IPHAS).NE.1) ) THEN
        FMIN = 4.405286343612334E-02
        FMAX = 5.506607929515418E-02
      ENDIF
C
      DO IEL =1, NCEL
C
C-- F
        IF ((FM(IEL).LE.FMIN).OR.(FM(IEL).GE.FMAX)) THEN
          FMP    =  MAX (MIN(FMAX,FM(IEL)),FMIN)
          CLIF = CLIF +1
        ELSE
          FMP = FM(IEL)
        ENDIF
C---Y
        CLIMAX =  (FMAX-FMP)*FMIN/(FMAX-FMIN)
        CLIMIN = MAX(ZERO,(FMP-FS(1))/(1.D0-FS(1)))
        IF (( YFM(IEL).GE.CLIMAX).OR.(YFM(IEL).LT.CLIMIN)) THEN
          YFMP = MAX(CLIMIN,MIN(YFM(IEL),CLIMAX))
          CLIY = CLIY + 1
        ELSE
          YFMP = YFM(IEL)
        ENDIF
C-- FP2M
        CLIMAX = (FMAX -FMP)*(FMP-FMIN)
        CLIMIN = ZERO
C
        IF ((FP2M(IEL).GE.CLIMAX).OR.(FP2M(IEL).LT.CLIMIN)) THEN
          FP2MP = MAX(CLIMIN,MIN(FP2M(IEL),CLIMAX))
          CLIFP2 = CLIFP2 + 1
        ELSE
          FP2MP = FP2M(IEL)
        ENDIF
C -- YFP2M
C    YFMAX = FMIN dans le cas Moreau
        CLIMAX = (FMIN-YFMP)*YFMP
        CLIMIN = ZERO
        IF ((YFP2M(IEL).GE.CLIMAX).OR.(YFP2M(IEL).LT.CLIMIN)) THEN
          YFP2MP = MAX(CLIMIN,MIN(YFP2M(IEL),CLIMAX))
          CLYFP2 = CLYFP2 + 1
        ELSE
          YFP2MP = YFP2M(IEL)
        ENDIF

C --> Clip pour la covariance
C
        CLIMAX = SQRT(FP2MP*YFP2MP)
        CLIMIN = -SQRT(FP2MP*YFP2MP)
        IF (COYFP(IEL).GE.CLIMAX) THEN
          COYFPP = CLIMAX
          CLICOY = CLICOY + 1
        ELSEIF (COYFP(IEL).LE.CLIMIN) THEN
          COYFPP = CLIMIN
          CLICOY = CLICOY + 1
        ELSE
          COYFPP = COYFP(IEL)
        ENDIF
c
        YFMPMX = (FMAX - FMP)*FMIN/(FMAX - FMIN)
        DYFMP   = (YFMPMX - YFMP)
C
        IF (((FP2MP.LT.EPSI).AND.(YFP2MP.LT.EPSI))
     &       .OR.
     &       ((YFMP.LT.EPSI ).OR.(DYFMP.LT.EPSI))
     &       .OR.
     &       ((FMP -FMIN).LT.EPSI)
     &       .OR.
     &       ( FP2MP.LT.EPSI     )
     &       .OR.
     &       (((FMAX -FMP).LT.EPSI)             )) THEN
C
c======================================================================
C   1.    NON PASSAGE PAR PDF
C======================================================================
C
          SUM1 = ZERO
          SUM2 = ZERO
          SUM3 = ZERO
          SUM4 = ZERO
          SUM5 = ZERO
          SUM6 = ZERO
          SUM15 = ZERO
          SUM16 = ZERO
C
          DO IDIRAC =1, NDIRAC
            D(IDIRAC) = 1.D0 / NDIRAC
            F(IDIRAC) = FMP
            Y(IDIRAC) = YFMP
C
C---> Calcul de l'enthalpie
C
            H(IDIRAC) = ((HMAX-HMIN)*F(IDIRAC) + HMIN*FMIN - HMAX*FMAX)
     &          / (FMIN-FMAX)
C
C ---> Calcul de la fraction massique des gaz (F, O et P)
C
            YFUEL = Y(IDIRAC)
            YOXYD = 1.D0 - (COEFF3+1.0D0)*F(IDIRAC) + COEFF3*Y(IDIRAC)
            YPROD = 1.D0 - YFUEL - YOXYD
            YO2   = COEFF1 - (COEFF1 + COEFF2) * F(IDIRAC)
     &           + COEFF2 * Y(IDIRAC)
C
C ---> Calcul de la masse molaire et de la temperature
C
            COEFG(1) = YFUEL
            COEFG(2) = YOXYD
            COEFG(3) = YPROD
C
C ------ Masse molaire
C
            NBMOL = 0.D0
            DO IGG = 1, NGAZG
              NBMOL = NBMOL + COEFG(IGG)/WMOLG(IGG)
            ENDDO
            MAML(IDIRAC) = 1.D0/NBMOL
C
C ------ Calcul de la temperature pour le pic 1 et 2
C
            MODE    = 1
            CALL COTHHT
C           ===========
     &           ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &           NPO    , NPOT  , TH      , EHGAZG ,
     &           H(IDIRAC)     , TEML(IDIRAC)    )
C
C     ---> Calcul de la masse volumique en 1 et 2
C
            IF ( IPASS.GT.1 .OR.
     &          (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
              RHOL(IDIRAC) = P0(IPHAS) * MAML(IDIRAC)
     &             / (RR*TEML(IDIRAC))
            ELSE
              RHOL(IDIRAC) = RO0(IPHAS)
            ENDIF
C
C ---> Calcul du terme source en 1 et 2 du scalaire YFM
C
            THETA(IDIRAC) = TA / TEML(IDIRAC)
     &           * (1.D0 - TEML(IDIRAC) / TSTAR)
C
            W(IDIRAC) = VREF / LREF * (- D(IDIRAC)*RHOL(IDIRAC)
     &           * YFUEL*YO2
     &           * EXP( -THETA(IDIRAC) ))
C     -> BO 27/06 Controle du signe de W
C
            W(IDIRAC) = MIN( W(IDIRAC), ZERO)
C
C ---> Masse molaire du melange
C
            SUM1 = SUM1 + D(IDIRAC)*MAML(IDIRAC)
C
C ---> Temperature du melange
C
            SUM2 = SUM2 + D(IDIRAC)*TEML(IDIRAC)
C
C ---> Temperature / Masse molaire
C
            SUM3 = SUM3 + D(IDIRAC)*TEML(IDIRAC)/MAML(IDIRAC)
C
C ---> Fractions massiques des especes globales
C
            SUM4 = SUM4 + YFUEL*D(IDIRAC)
C
            SUM5 = SUM5 + YOXYD*D(IDIRAC)
C
            SUM6 = SUM6 + YPROD*D(IDIRAC)
C
            SUM15 = SUM15 +RHOL(IDIRAC)*D(IDIRAC)
C
            SUM16 = SUM16 +W(IDIRAC)
C
C ---> Stockage des proprietes via PROPCE
C
            PROPCE(IEL,IPAMPL(IDIRAC)) = D(IDIRAC)
            PROPCE(IEL,IPFMEL(IDIRAC)) = F(IDIRAC)
            PROPCE(IEL,IPFMAL(IDIRAC)) = Y(IDIRAC)
            PROPCE(IEL,IPTEML(IDIRAC)) = TEML(IDIRAC)
            PROPCE(IEL,IPMAML(IDIRAC)) = MAML(IDIRAC)
            PROPCE(IEL,IPRHOL(IDIRAC)) = RHOL(IDIRAC)
            PROPCE(IEL,IPTSCL(IDIRAC)) = W(IDIRAC)
C
          ENDDO
C
          PROPCE(IEL,IPCMAM) = SUM1
          PROPCE(IEL,IPCTEM) = SUM2
          TEMSMM             = SUM3
          PROPCE(IEL,IPCFUE) = SUM4
          PROPCE(IEL,IPCOXY) = SUM5
          PROPCE(IEL,IPCPRO) = SUM6
          PROPCE(IEL,IPCTSC) = SUM16
C
C---> Masse volumique du melange
C
          IF ( IPASS.GT.1 .OR.
     &        (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
            PROPCE(IEL,IPCROM) = SRROM*PROPCE(IEL,IPCROM)
     &                          +(1.D0-SRROM)*(P0(IPHAS)/(RR*TEMSMM))
          ENDIF
C
        ELSE
C
C==========================================================================
C    2.    PASSAGE PAR PDF
C==========================================================================
C
C -------> Constantes
C
          CST = 1.D0
C
C-------->Calcul de F1 et F2 avec Curl em F
C
          CALL LWCURL
C         =========
     &         ( CST  , FMP , FP2MP  ,
     &         FMIN , FMAX ,
     &         F1 , F2 , CSTFA1 , CSTFA2 )
C
C ------>  On calcul les Moyennes conditionnelles Y1, Y2
C
          Y2 = ((FMP*YFMP + COYFPP) - F1*YFMP)
     &         /(CSTFA1*(F2 - F1))
          Y1 = (YFMP - CSTFA2*Y2)/CSTFA1
C
          YMIN(1) = MAX(ZERO , ((F1- FS(1))/(1D0-FS(1))))
          YMAX(1) = (FMAX - F1)*FMIN/(FMAX - FMIN)
          YMIN(2) = MAX(ZERO , ((F2- FS(1))/(1D0-FS(1))))
          YMAX(2) = (FMAX - F2)*FMIN/(FMAX - FMIN)
C
C clipping pour les moyennes conditionnelles
C
c Y1 = MAX(YMIN(1),MIN(Y1,YMAX(1)))
C  ===> compteur
C
          IF (Y1.GE.YMAX(1)) THEN
            Y1 = YMAX(1)
            CLIY1 = CLIY1 + 1
          ELSEIF (Y1.LE.YMIN(1)) THEN
            Y1 = YMIN(1)
            CLIY1 = CLIY1 + 1
          ENDIF
C == fin
C Y2 = MAX(YMIN(2),MIN(Y2,YMAX(2)))
C
          IF (Y2.GE.YMAX(2)) THEN
            Y2 = YMAX(2)
            CLIY2 = CLIY2 + 1
          ELSEIF (Y2.LE.YMIN(2)) THEN
            Y2 = YMIN(2)
            CLIY2 = CLIY2 + 1
          ENDIF
C
          Y2P(1)  = ((YFMP**2 + YFP2MP) - CSTFA2*(Y2**2))
     &         /CSTFA1 - Y1**2
c          WRITE(NFECRA,*) ' Y2P(1)=',Y2P(1)
C
C clipping pour variance conditionnelles
C
c          Y2P(1) = MAX(ZERO,
c     &           MIN(Y2P(1),((Y1-YMIN(1))*(YMAX(1) - Y1))))
C  ====> Compteur
C
          CLIMAX =( (Y1-YMIN(1))*(YMAX(1) - Y1))
          CLIMIN = ZERO
          IF (Y2P(1).GE.CLIMAX) THEN
            Y2P(1) = CLIMAX
            CLIY2P = CLIY2P + 1
          ELSEIF (Y2P(1).LE.CLIMIN) THEN
            Y2P(1) = CLIMIN
            CLIY2P = CLIY2P + 1
          ENDIF
C==== fin
          CALL LWCURL
C         =========
     &         ( CSTFA1  , Y1   , Y2P(1)  ,
     &         YMIN(1) , YMAX(1) ,
     &         Y(1) , Y(2) , D(1) , D(2) )
C
C ---------> Parametres des dirac en F1
C
          F(1) = F1
          F(2) = F1
C
C ---------> Parametres du dirac en F2
C
          F(3) = F2
          Y(3) = Y2
          D(3) = CSTFA2
C
C=======================================================================
C 3.  DETERMINATION DES GRANDEURS THERMOCHIMIQUES DES DEUX PICS
C=======================================================================
C
C ---> Calcul de l'enthalpies en 1 et 2
C
          SUM7  = ZERO
          SUM8  = ZERO
          SUM9  = ZERO
          SUM10 = ZERO
          SUM11 = ZERO
          SUM12 = ZERO
          SUM17 = ZERO
C
          DO IDIRAC = 1, NDIRAC
            H(IDIRAC) = ( (HMAX-HMIN)*F(IDIRAC)
     &                   + HMIN*FMIN - HMAX*FMAX) / (FMIN-FMAX)
C
C ---> Calcul de la fraction massique des gaz (F, O et P) en 1 et 2
C
            YFUEL = Y(IDIRAC)
            YOXYD = 1.D0 - (COEFF3+1.0D0)*F(IDIRAC)
     &           + COEFF3*Y(IDIRAC)
            YPROD = 1.D0 - YFUEL - YOXYD
            YO2   = COEFF1 - (COEFF1 + COEFF2) * F(IDIRAC)
     &                     + COEFF2 * Y(IDIRAC)
C
C ---> Calcul de la masse molaire et de la temperature en 1 et 2
C
            COEFG(1) = YFUEL
            COEFG(2) = YOXYD
            COEFG(3) = YPROD
C
C ------ Masse molaire pour le pic 1 et 2
C
            NBMOL = 0.D0
            DO IGG = 1, NGAZG
              NBMOL = NBMOL + COEFG(IGG)/WMOLG(IGG)
            ENDDO
            MAML(IDIRAC) = 1.D0/NBMOL
C
C ------ Calcul de la temperature pour le pic 1 et 2
C
            MODE    = 1
            CALL COTHHT
C           ===========
     &           ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &           NPO    , NPOT  , TH      , EHGAZG ,
     &           H(IDIRAC)      , TEML(IDIRAC)    )
C
C ---> Calcul de la masse volumique en 1 et 2
C
            IF ( IPASS.GT.1 .OR.
     &          (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
              RHOL(IDIRAC) = P0(IPHAS) * MAML(IDIRAC)
     &             /(RR*TEML(IDIRAC))
            ELSE
              RHOL(IDIRAC) = RO0(IPHAS)
            ENDIF
C
C ---> Calcul du terme source en 1 et 2 du scalaire YFM
C
            THETA(IDIRAC) = TA / TEML(IDIRAC)
     &           *(1.D0 - TEML(IDIRAC)/TSTAR)
C
            W(IDIRAC) = VREF / LREF
     &           *(- D(IDIRAC)*RHOL(IDIRAC)
     &           * YFUEL*YO2
     &           * EXP( -THETA(IDIRAC) ))
C ---> Controle du signe de W
C
            W(IDIRAC) = MIN( W(IDIRAC), ZERO)
C
C
C ---> Masse molaire du melange
C
            SUM7 = SUM7 + D(IDIRAC)*MAML(IDIRAC)
C
C ---> Temperature du melange
C
            SUM8 = SUM8 + D(IDIRAC)*TEML(IDIRAC)
C
C ---> Temperature / Masse molaire
C
            SUM9 = SUM9 + D(IDIRAC)*TEML(IDIRAC)
     &           /MAML(IDIRAC)
C
C ---> Fractions massiques des especes globales
C
            SUM10 = SUM10 + YFUEL*D(IDIRAC)
C
            SUM11 = SUM11 + YOXYD*D(IDIRAC)
C
            SUM12 = SUM12 + YPROD*D(IDIRAC)
C
            SUM17 = SUM17 + W(IDIRAC)
C
C ---> Stockage des proprietes via PROPCE
C
            PROPCE(IEL,IPAMPL(IDIRAC)) = D(IDIRAC)
            PROPCE(IEL,IPFMEL(IDIRAC)) = F(IDIRAC)
            PROPCE(IEL,IPFMAL(IDIRAC)) = Y(IDIRAC)
            PROPCE(IEL,IPMAML(IDIRAC)) = MAML(IDIRAC)
            PROPCE(IEL,IPTEML(IDIRAC)) = TEML(IDIRAC)
            PROPCE(IEL,IPRHOL(IDIRAC)) = RHOL(IDIRAC)
            PROPCE(IEL,IPTSCL(IDIRAC)) = W(IDIRAC)
C
C
            IF ((F(IDIRAC).NE.ZERO).AND.(Y(IDIRAC).NE.ZERO)) THEN
              IF ((F(IDIRAC).GT.FS(1)).OR.
     &             (F(IDIRAC).LT.0.8*FS(1))) THEN
                WRITE(NFECRA,*)'**************F OUT**************',
     &               IDIRAC
                WRITE(NFECRA,*)'F',F(IDIRAC)
                WRITE(NFECRA,*)' IEL =', IEL
              ENDIF
C
              IF ((Y(IDIRAC).GT.F(IDIRAC))
     &             .OR.(Y(IDIRAC).LT.-EPSI)) THEN
                WRITE(NFECRA,*)'*************Y OUT*****************',
     &               IDIRAC
                WRITE(NFECRA,*)'Y',Y(IDIRAC)
                WRITE(NFECRA,*)' IEL =', IEL
              ENDIF
            ENDIF
          ENDDO
C
          PROPCE(IEL,IPCMAM) = SUM7
          PROPCE(IEL,IPCTEM) = SUM8
          TEMSMM             = SUM9
          PROPCE(IEL,IPCFUE) = SUM10
          PROPCE(IEL,IPCOXY) = SUM11
          PROPCE(IEL,IPCPRO) = SUM12
          PROPCE(IEL,IPCTSC) = SUM17
C
C ---> Masse volumique du melange
C
          IF ( IPASS.GT.1 .OR.
     &        (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
            PROPCE(IEL,IPCROM) = SRROM * PROPCE(IEL,IPCROM)
     &           + (1.D0-SRROM) * (P0(IPHAS)/(RR*TEMSMM))
          ENDIF
C
        ENDIF
      ENDDO
C
      END
