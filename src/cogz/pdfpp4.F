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
                         SUBROUTINE PDFPP4
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
CFONC  PDF LIBBY - WILLIAMS 4  POINTS AVEC HYPOTHESE DE CURL
CFONC              OU CURL MODIFIE
CFONC
CFONC
CFONC  COMMENTAIRES :
CFONC  ------------
CFONC
CFONC     Dans un diagramme (F, Yf), on construit deux droites:
CFONC          - La droite de combustion complete
CFONC          - La droite de melange
CFONC
CFONC     Dans ce domaine, nous allons construire  deux pics sur F qui
CFONC     seront dedoubler chaque un en deux avec un Curl sur Yf
CFONC
CFONC
CFONC  LE RESULTAT EST :
CFONC  ---------------
CFONC
CFONC     CALCUL DES PARAMETRES ASSOCIES AUX FONCTIONS DIRAC
CFONC
CFONC       Les positions des Pics sont :
CFONC          [F(1),Y1(1)] et [F(1),Y1(2)]
CFONC          [F(2),Y2(1)] et [F(2),Y2(2)]
CFONC       Leurs amplitudes respectives sont :
CFONC                D1(1) et D1(2)
CFONC                D2(1) et D2(2)
CFONC       Pour chaque dirac, on calcule :
CFONC       la temperature Ti(j),
CFONC       la masse volumique RHOi(j),
CFONC       le terme source chimique Wi(j),
CFONC            i etant la positiion sur F du pic de Dirac
CFONC            j etant la positiion sur Yf du pic de Dirac
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
CARGU ! COYFP        ! TR !  -> ! COVARIANCE                           !
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
      INCLUDE "parall.h"
C
C***********************************************************************
C ARGUMENTS
C***********************************************************************
C
      INTEGER          NCELET, NCEL
C
      DOUBLE PRECISION FM(NCELET)   , FP2M(NCELET)
      DOUBLE PRECISION YFM(NCELET)  , YFP2M(NCELET)
      DOUBLE PRECISION COYFP(NCELET)
      DOUBLE PRECISION PROPCE(NCELET,*)
C
C***********************************************************************
C VARIABLES LOCALES
C***********************************************************************
C
      INTEGER          IEL, IGG, IDIRAC
      INTEGER          MODE, IPHAS
C
      INTEGER          IPCTEM, IPCMAM,  IPCROM
      INTEGER          IPAMPL(NDRACM), IPFMEL(NDRACM)
      INTEGER          IPFMAL(NDRACM), IPTEML(NDRACM)
      INTEGER          IPMAML(NDRACM)
      INTEGER          IPRHOL(NDRACM)
      INTEGER          IPTSCL(NDRACM)
      INTEGER          IPCFUE, IPCOXY, IPCPRO, IPCTSC
      INTEGER          IPCKAB, IPT4 , IPT3
C
      DOUBLE PRECISION F(NDRACM), Y(NDRACM), D(NDRACM)
      DOUBLE PRECISION H(NDRACM), TEML(NDRACM)
      DOUBLE PRECISION MAML(NDRACM), W(NDRACM)
      DOUBLE PRECISION RHOL(NDRACM), THETA(NDRACM)
C
      DOUBLE PRECISION COEFG(NGAZGM), EPSI
      DOUBLE PRECISION YFUEL
      DOUBLE PRECISION YOXYD, YO2
      DOUBLE PRECISION YPROD,FMP,FP2MP,YFMP,YFP2MP,COYFPP
      DOUBLE PRECISION YFP2MAX
C
      DOUBLE PRECISION NBMOL,  TEMSMM
      DOUBLE PRECISION SUM7, SUM8, SUM9, SUM10, SUM11, SUM12, SUM17
      DOUBLE PRECISION SUM1, SUM2, SUM3, SUM4, SUM5, SUM6, SUM16, SUM15
      DOUBLE PRECISION WMIN, WMAX, TETMIN, TETMAX, O2MIN, O2MAX
C
      DOUBLE PRECISION Y1, Y2, F1, F2
      DOUBLE PRECISION CSTFA1, CSTFA2, CST
      DOUBLE PRECISION CSTVAR, CSTVA2
C
      DOUBLE PRECISION YMIN(2), YMAX(2)
      DOUBLE PRECISION Y2P(2)
      DOUBLE PRECISION Y2PMIN(2), Y2PMAX(2)
C
C ---> variables pour clipping
C
      INTEGER          ICPT1,ICPT2
      INTEGER          CLY2P1, CLY2P2, CLY2P3,CLY2P4
      INTEGER          CLFP2M, CLYF21, CLYF22
      INTEGER          CLCYF1, CLCYF2
C
      DOUBLE PRECISION VYMX, VFMX, VYMN
      DOUBLE PRECISION VCYFMX, VCYFMN
      DOUBLE PRECISION MXCYFP,MXCFP, MNCYFP
      DOUBLE PRECISION MXCCYF, MNCCYF
      DOUBLE PRECISION MXYFP2,MAXFP2, MNYFP2
      DOUBLE PRECISION MXCOYF, MNCOYF
      DOUBLE PRECISION MCY2P1, MCY2P3
      DOUBLE PRECISION MCY2P2, MCY2P4
      DOUBLE PRECISION MY2P1 , MY2P3, MY2P2, MY2P4
C
      INTEGER IPASS
      DATA    IPASS /0/
      SAVE    IPASS
C
C***********************************************************************
C
C=======================================================================
C     0. POSITION ET INITIALISATION DES VARIABLES
C=======================================================================
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
C
      IPCFUE = IPPROC(IYM(1))
      IPCOXY = IPPROC(IYM(2))
      IPCPRO = IPPROC(IYM(3))
      IPCTSC = IPPROC(ITSC)
      IPCTEM = IPPROC(ITEMP)
      IPHAS = 1
      IPCROM = IPPROC(IROM(IPHAS))
      IPCMAM = IPPROC(IMAM)
C
C      IF ( IRAYPP.GT.0 ) THEN
C        IPCKAB = IPPROC(ICKABS)
C        IPT4   = IPPROC(IT4M)
C        IPT3   = IPPROC(IT3M)
C      ENDIF
C
C ---> Initialisation des variables et compteurs
C
      DO IGG = 1, NGAZGM
        COEFG(IGG) = ZERO
      ENDDO
C
      EPSI = 1.D-10
C
      CLFP2M  = 0
      CLYF21 = 0
      CLYF22 = 0
      CLCYF1 = 0
      CLCYF2 = 0
C
      CLY2P1 = 0
      CLY2P2 = 0
      CLY2P3 = 0
      CLY2P4 = 0
C
      WMAX = -1.D+10
      WMIN =  1.D+10
      TETMIN = 1.D+10
      TETMAX = -1.D+10
      O2MAX = -1.D+10
      O2MIN =  1.D+10
C
      MXCFP = -1.D+10
      MXCYFP = -1.D+10
      MNCYFP = -1.D+10
      MXCCYF = -1.D+10
      MNCCYF = -1.D+10
C
      MAXFP2= -1.D+10
      MXYFP2= -1.D+10
      MNYFP2= -1.D+10
      MXCOYF= -1.D+10
      MNCOYF= -1.D+10
C
      MCY2P1=-1.D+10
      MCY2P3=-1.D+10
      MCY2P2=-1.D+10
      MCY2P4=-1.D+10
C
      MY2P1=-1.D+10
      MY2P3=-1.D+10
      MY2P2=-1.D+10
      MY2P4=-1.D+10
C
      ICPT1 = 0
      ICPT2 = 0
C
C=======================================================================
C     1. BOUCLE SUR LES CELLULES
C=======================================================================
C
      IPASS = IPASS + 1
C
      DO IEL = 1, NCEL
C
C ---> position des variables
C
      FMP=FM(IEL)
      YFMP=MAX(YFM(IEL), ZERO)
      FP2MP=FP2M(IEL)
      YFP2MP=YFP2M(IEL)
      COYFPP=COYFP(IEL)
C
C=======================================================================
C    TEST DE PASSAGE OU NON PAR LA PDF
C    (non passage par pdf si aucune variance ni en f ni en y)
C    (on peut si on veut toujours passer par la pdf,
C     le cas des variances nulles y est traite)
C=======================================================================
C
       IF (((FP2MP.LT.EPSI).AND.(YFP2MP.LT.EPSI))) THEN
C
c======================================================================
C    2.    NON PASSAGE PAR PDF
C======================================================================
C
          ICPT1 = ICPT1 +1
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
C ---> boucle sur chaque DIRAC
C
          DO IDIRAC =1, NDIRAC
C
C ---> calcul de f, Y, Amplitude de chaque DIRAC == Val moy
C
            D(IDIRAC) = 1.D0 / NDIRAC
            F(IDIRAC) = FMP
            Y(IDIRAC) = YFMP
C
C ---> Calcul de l'enthalpie
C
            H(IDIRAC) = ((HMAX-HMIN)*F(IDIRAC) + HMIN*FMAX - HMAX*FMIN)
     &                / (FMAX-FMIN)
C
C ---> Calcul de la fraction massique des gaz (F, O et P)
C
            YFUEL = Y(IDIRAC)
            YOXYD = 1.D0 - (COEFF3+1.0D0)*F(IDIRAC) + COEFF3*Y(IDIRAC)
            YPROD = 1.D0 - YFUEL - YOXYD
            YO2   = COEFF1 - (COEFF1 + COEFF2) * F(IDIRAC)
     &                      + COEFF2 * Y(IDIRAC)
C
C ---> Calcul de la masse molaire et de la temperature
C
            COEFG(1) = YFUEL
            COEFG(2) = YOXYD
            COEFG(3) = YPROD
C
C ---> Masse molaire
C
            NBMOL = ZERO
            DO IGG = 1, NGAZG
              NBMOL = NBMOL + COEFG(IGG)/WMOLG(IGG)
            ENDDO
            MAML(IDIRAC) = 1.D0/NBMOL
C
C ---> Calcul de la temperature pour chaque pic
C
            MODE    = 1
            CALL COTHHT
C           ===========
     &      ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &        NPO    , NPOT  , TH      , EHGAZG ,
     &        H(IDIRAC)     , TEML(IDIRAC)    )
C
C ---> Calcul de la masse volumique pour chaque pic
C
            IF ( IPASS.GT.1 .OR.
     &          (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
              RHOL(IDIRAC) = P0(IPHAS) * MAML(IDIRAC)
     &             / (RR*TEML(IDIRAC))
            ELSE
              RHOL(IDIRAC) = RO0(IPHAS)
            ENDIF
C
C ---> Calcul du terme source du scalaire YFM pour chaque pic
C
            THETA(IDIRAC) = TA / TEML(IDIRAC)
     &           * (1.D0 - TEML(IDIRAC) / TSTAR)
C
            W(IDIRAC) = VREF / LREF * (- D(IDIRAC)
     &           * YFUEL*YO2
     &           * EXP( -THETA(IDIRAC) ))
C
C ---> Controle du signe de W
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
            SUM5 = SUM5 + YOXYD*D(IDIRAC)
            SUM6 = SUM6 + YPROD*D(IDIRAC)
            SUM15 = SUM15 +RHOL(IDIRAC)*D(IDIRAC)
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
Cfin de boucle sur chaque DIRAC
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
C ---> Masse volumique du melange
C
          IF ( IPASS.GT.1.OR.
     &        (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
            PROPCE(IEL,IPCROM) = SRROM*PROPCE(IEL,IPCROM)
     &           +(1.D0-SRROM)*(P0(IPHAS)/(RR*TEMSMM))
          ENDIF
C
        ELSE
C
C==========================================================================
C    3.    PASSAGE PAR PDF
C==========================================================================
C
          ICPT2 = ICPT2 + 1
C
C ---> Clipping sur la variance en f
C
          VFMX = (FMAX-FMP)*(FMP-FMIN)
C
          IF (FP2MP.GT.VFMX) THEN
            IF ((FP2MP-VFMX).GT.MXCFP) THEN
              MXCFP=(FP2MP-VFMX)
              MAXFP2=VFMX
            ENDIF
            FP2MP=VFMX
            CLFP2M=CLFP2M+1
          ENDIF
C
C ---> Calcul des positions et amplitudes en F des Pics avec lWCURL en F
C
          CST = 1.D0
C
          CALL LWCURL
C         =========
     &   ( CST  , FMP , FP2MP  ,
     &       FMIN , FMAX ,
     &       F1 , F2 , CSTFA1 , CSTFA2 )
C
          F(1) = F1
          F(2) = F1
          F(3) = F2
          F(4) = F2
C
C ---> Determination des valeur max et min de Yf en F1 et F2
C
          YMIN(1) = MAX(ZERO , ((F1- FS(1))/(1D0-FS(1))))
          YMAX(1) = F1
          YMIN(2) = MAX(ZERO , ((F2- FS(1))/(1D0-FS(1))))
          YMAX(2) = F2
C
C ---> clipping COVARIANCE en fonction des bornes des moyennes conditionnelles
C
          VCYFMX=MIN(YMAX(2)*CSTFA2*(F2-F1)-YFMP*(FMP-F1)
     &         ,YFMP*(F2-FMP)-YMIN(1)*CSTFA1*(F2-F1))
C
          VCYFMN=MAX(YMIN(2)*CSTFA2*(F2-F1)-YFMP*(FMP-F1)
     &         ,YFMP*(F2-FMP)-YMAX(1)*CSTFA1*(F2-F1))
C
          IF (COYFPP.GT.VCYFMX) THEN
            IF ((COYFPP-VCYFMX).GT.MXCCYF) THEN
              MXCCYF=(COYFPP-VCYFMX)
              MXCOYF=VCYFMX
            ENDIF
            COYFPP=VCYFMX
            CLCYF1=CLCYF1+1
          ELSEIF (COYFPP.LT.VCYFMN) THEN
            IF ((VCYFMN-COYFPP).GT.MNCCYF) THEN
              MNCCYF=(VCYFMN-COYFPP)
              MNCOYF=VCYFMN
            ENDIF
            COYFPP=VCYFMN
            CLCYF2=CLCYF2+1
          ENDIF
C
C --->  On calcul les Moyennes conditionnelles Y1, Y2
C
          IF ((F2-F1).GT.EPSI) THEN
            Y2 = (YFMP*(FMP- F1) + COYFPP)
     &           /(CSTFA2*(F2 - F1))
            Y1 = (YFMP*(F2-FMP) - COYFPP)
     &           /(CSTFA1*(F2 - F1))
          ELSE
            Y2 = YFMP
            Y1 = YFMP
          ENDIF
          IF ((FMP-YFMP).LT.EPSI) THEN
            Y2=F2
            Y1=F1
          ENDIF
          IF ((YFMP-MAX(ZERO,(FMP -FS(1))
     &         /(1.D0 - FS(1)))).LT.EPSI) THEN
            Y2=MAX(ZERO,(F2 -FS(1))/(1.D0 - FS(1)))
            Y1=MAX(ZERO,(F1 -FS(1))/(1.D0 - FS(1)))
          ENDIF
C
C ---> Determination des valeurs max et min des variances de y en F1 et F2
C
          Y2PMAX(1) =(Y1-YMIN(1))*(YMAX(1) - Y1)
          Y2PMIN(1) = ZERO
          Y2PMAX(2) = (Y2-YMIN(2))*(YMAX(2) - Y2)
          Y2PMIN(2) = ZERO
C
C ---> calcul de la variance en Y max avec Y2PMAX 1 et 2
C
          YFP2MAX=CSTFA1*(Y1**2+Y2PMAX(1))
     &         +CSTFA2*(Y2**2+Y2PMAX(2))-YFMP**2
          VYMX=YFP2MAX
C
C ---> rapport des variances conditionnelles
C
          IF (((YMAX(2)-Y2).GT.EPSI).AND.((Y2-YMIN(2)).GT.EPSI)
     &         .AND.((YMAX(1)-Y1).GT.EPSI).AND.((Y1-YMIN(1)).GT.EPSI))
     &         THEN
            CSTVAR = ((YMAX(2)-Y2)*(Y2-YMIN(2)))
     &           /((YMAX(1)-Y1)*(Y1 -YMIN(1)))
          ENDIF
C
C
C ---> clipping VARIANCE Y
C      (on peut soit clipper la variance en y en fonction
C       des valeurs extremes des variances conditionnelles
C       ou clipper directement les variances conditionnelles)
C
C         CSTVA2 =  (CSTFA1*Y1**2+CSTFA2*Y2**2)-YFMP**2
C
C         IF (((YMAX(2)-Y2).GT.EPSI).AND.((Y2-YMIN(2)).GT.EPSI)
C     &  .AND.((YMAX(1)-Y1).GT.EPSI).AND.((Y1-YMIN(1)).GT.EPSI))
C     &       THEN
C
C      VYMX=MIN(Y2PMAX(1)*(CSTFA1+CSTFA2*CSTVAR)+CSTVA2
C     &         ,Y2PMAX(2)*(CSTFA1/CSTVAR+CSTFA2)+CSTVA2)
C       VYMN=MAX(Y2PMIN(1)*(CSTFA1+CSTFA2*CSTVAR)+CSTVA2
C     &           ,Y2PMIN(2)*(CSTFA1/CSTVAR+CSTFA2)+CSTVA2)
C
C         ELSEIF (((YMAX(2)-Y2).GT.EPSI)
C     &         .AND.((Y2-YMIN(2)).GT.EPSI)) THEN
C  VYMX=Y2PMAX(2)*CSTFA2+CSTVA2
C  VYMN=Y2PMIN(2)*CSTFA2+CSTVA2
C
C         ELSEIF (((YMAX(1)-Y1).GT.EPSI)
C     &         .AND.((Y1-YMIN(1)).GT.EPSI)) THEN
C  VYMX=Y2PMAX(1)*CSTFA1+CSTVA2
C  VYMN=Y2PMIN(1)*CSTFA1+CSTVA2
C  ENDIF
C
          IF (YFP2MP.GT.VYMX) THEN
            IF ((YFP2MP-VYMX).GT.MXCYFP) THEN
              MXCYFP=(YFP2MP-VYMX)
              MXYFP2=VYMX
            ENDIF
            YFP2MP=VYMX
            CLYF21=CLYF21+1
C   ELSEIF (YFP2M(IEL).LT.VYMN) THEN
C   IF ((VYMN-YFP2MP).GT.MNCYFP) THEN
C      MNCYFP=(VYMN-YFP2MP)
C      MNYFP2=VYMN
C   ENDIF
C   YFP2MP=VYMN
C   CLYF22=CLYF22+1
          ENDIF
C
C ---> calcul des variances conditionnelles
C
          IF (((YMAX(2)-Y2).GT.EPSI).AND.((Y2-YMIN(2)).GT.EPSI)
     &         .AND.((YMAX(1)-Y1).GT.EPSI).AND.((Y1-YMIN(1)).GT.EPSI))
     &         THEN
C
            Y2P(1) = ((YFMP**2)+YFP2MP
     &           -CSTFA2*(Y2**2)-CSTFA1*(Y1**2))
     &           /(CSTFA1 + CSTFA2*CSTVAR)
C
            Y2P(2) = ((YFMP**2)+YFP2MP
     &           -CSTFA2*(Y2**2)-CSTFA1*(Y1**2))
     &           /(CSTFA1/CSTVAR + CSTFA2)
C
          ELSEIF (((YMAX(2)-Y2).GT.EPSI)
     &           .AND.((Y2-YMIN(2)).GT.EPSI)) THEN
C
            Y2P(1) = ZERO
C
            Y2P(2) = ((YFMP**2)+YFP2MP
     &           -CSTFA2*(Y2**2)-CSTFA1*(Y1**2))
     &           /CSTFA2
C
          ELSEIF (((YMAX(1)-Y1).GT.EPSI)
     &           .AND.((Y1-YMIN(1)).GT.EPSI)) THEN
C
            Y2P(2) = ZERO
C
            Y2P(1) = ((YFMP**2)+YFP2MP
     &           -CSTFA2*(Y2**2)-CSTFA1*(Y1**2))
     &           /CSTFA1
C
          ELSE
            Y2P(1) = ZERO
            Y2P(2) = ZERO
          ENDIF
C
C ---> clipping pour variances conditionnelles
C
          IF (Y2P(1).GT.Y2PMAX(1)) THEN
            IF ((Y2P(1)-Y2PMAX(1)).GT.MCY2P1) THEN
              MCY2P1=(Y2P(1)-Y2PMAX(1))
              MY2P1=Y2PMAX(1)
            ENDIF
            Y2P(1) = Y2PMAX(1)
            Y2P(2) = (((YFMP**2)+YFP2MP-CSTFA1*
     &           ((Y1**2)+Y2P(1)))/CSTFA2)-(Y2**2)
            CLY2P1 =  CLY2P1 + 1
          ELSEIF (Y2P(1).LT.Y2PMIN(1)) THEN
            IF ((Y2PMIN(1)-Y2P(1)).GT.MCY2P3) THEN
              MCY2P3=(Y2PMIN(1)-Y2P(1))
              MY2P3=Y2PMIN(1)
            ENDIF
            Y2P(1) = Y2PMIN(1)
            Y2P(2) = (((YFMP**2)+YFP2MP-CSTFA1*
     &           ((Y1**2)+Y2P(1)))/CSTFA2)-(Y2**2)
            CLY2P3 = CLY2P3 + 1
          ENDIF
          IF (Y2P(2).GT.Y2PMAX(2)) THEN
            IF ((Y2P(2)-Y2PMAX(2)).GT.MCY2P2) THEN
              MCY2P2=(Y2P(2)-Y2PMAX(2))
              MY2P2=Y2PMAX(2)
            ENDIF
            Y2P(2) = Y2PMAX(2)
            Y2P(1) = (((YFMP**2)+YFP2MP-CSTFA2*
     &           ((Y2**2)+Y2P(2)))/CSTFA1)-(Y1**2)
            CLY2P2 =  CLY2P2 + 1
          ELSEIF (Y2P(2).LT.Y2PMIN(2)) THEN
            IF ((Y2PMIN(2)-Y2P(2)).GT.MCY2P4) THEN
              MCY2P4=(Y2PMIN(2)-Y2P(2))
              MY2P4=Y2PMIN(2)
            ENDIF
            Y2P(2) = Y2PMIN(2)
            Y2P(1) = (((YFMP**2)+YFP2MP-CSTFA2*
     &           ((Y2**2)+Y2P(2)))/CSTFA1)-(Y1**2)
            CLY2P4 = CLY2P4 + 1
          ENDIF
C
C ---> calcul des positions et amplitudes des pics en Y sur F1
C
          CALL LWCURL
C         =========
     &  ( CSTFA1  , Y1   , Y2P(1)  ,
     &    YMIN(1) , YMAX(1) ,
     &    Y(1) , Y(2) , D(1) , D(2) )
C
C ---> calcul des positions et amplitudes des pics en Y sur F2
C
          CALL LWCURL
C         =========
     &  ( CSTFA2  , Y2   , Y2P(2)  ,
     &    YMIN(2) , YMAX(2) ,
     &    Y(3) , Y(4) , D(3) , D(4) )
C
C
C=======================================================================
C 2.  DETERMINATION DES GRANDEURS THERMOCHIMIQUES DES DEUX PICS
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
C ---> boucle sur chaque DIRAC
C
        DO IDIRAC = 1, NDIRAC
C
C ---> Calcul de l'enthalpie
C
          H(IDIRAC) = ( (HMAX-HMIN)*F(IDIRAC)
     &         + HMIN*FMAX - HMAX*FMIN) / (FMAX-FMIN)
C
C ---> Calcul de la fraction massique des gaz (F, O et P) pour chaque pic
C
          YFUEL = Y(IDIRAC)
          YOXYD = 1.D0 - (COEFF3+1.0D0)*F(IDIRAC)
     &         + COEFF3*Y(IDIRAC)
          YPROD = 1.D0 - YFUEL - YOXYD
          YO2   = COEFF1 - (COEFF1 + COEFF2) * F(IDIRAC)
     &         + COEFF2 * Y(IDIRAC)
C
C ---> Coefficients d'absorption pour chaque pic
C
c            IF ( IRAYPP .GT. 0  ) THEN
c              KABSGF = YFUEGF(IEL)*KABSG(1) + YOXYGF(IEL)*KABSG(2)
c     &                +YPROGF(IEL)*KABSG(3)
c              KABSGB = YFUEGB(IEL)*KABSG(1) + YOXYGB(IEL)*KABSG(2)
c     &                +YPROGB(IEL)*KABSG(3)
c            ENDIF
C
C
C ---> Calcul de la masse molaire et de la temperature pour chaque pic
C
          COEFG(1) = YFUEL
          COEFG(2) = YOXYD
          COEFG(3) = YPROD
C
C ---> Masse molaire pour chaque pic
C
          NBMOL = ZERO
          DO IGG = 1, NGAZG
            NBMOL = NBMOL + COEFG(IGG)/WMOLG(IGG)
          ENDDO
          MAML(IDIRAC) = 1.D0/NBMOL
C
C --->Calcul de la temperature pour chaque pic
C
          MODE    = 1
          CALL COTHHT
C         ===========
     & ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &   NPO    , NPOT  , TH      , EHGAZG ,
     &   H(IDIRAC)      , TEML(IDIRAC)    )
C
C ---> Calcul de la masse volumique  pour chaque pic
C
          IF ( IPASS.GT.1.OR.
     &        (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1)) THEN
            RHOL(IDIRAC) = P0(IPHAS) * MAML(IDIRAC)
     &           /(RR*TEML(IDIRAC))
          ELSE
            RHOL(IDIRAC) = RO0(IPHAS)
          ENDIF
C
C ---> Calcul du terme source du scalaire YFM  pour chaque pic
C
          THETA(IDIRAC) = TA / TEML(IDIRAC)
     &         *(1.D0 - TEML(IDIRAC)/TSTAR)
C
          TETMAX = MAX(THETA(IDIRAC),TETMAX)
          TETMIN = MIN(THETA(IDIRAC),TETMIN)
C
          W(IDIRAC) = VREF / LREF
     &         *(- D(IDIRAC)
     &         * YFUEL*YO2
     &         * EXP( -THETA(IDIRAC) ))
C
          WMAX = MAX(W(IDIRAC),WMAX)
          WMIN = MIN(W(IDIRAC),WMIN)
          O2MAX = MAX(YO2,O2MAX)
          O2MIN = MIN(YO2,O2MIN)
C
C ---> Controle du signe de W
C
          W(IDIRAC) = MIN( W(IDIRAC), ZERO)
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
     &         /MAML(IDIRAC)
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
C ---> Grandeurs relatives au rayonnement
C
c           IF ( IRAYPP.GT.0 ) THEN
c             PROPCE(IEL,IPCKAB) = YGFM*KABSGF + YGBM*KABSGB
c             PROPCE(IEL,IPT4)   = YGFM*TGF**3+YGBM*TGB**3
c             PROPCE(IEL,IPT3)   = YGFM*TGF**3+YGBM*TGB**3
c           ENDIF
C
C fin de boucle sur chaque DIRAC
        ENDDO
C
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
     &      (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
          PROPCE(IEL,IPCROM) = SRROM * PROPCE(IEL,IPCROM)
     &         + (1.D0-SRROM) * (P0(IPHAS)/(RR*TEMSMM))
        ENDIF
C
C de passage ou non par la PDF
      ENDIF
C
C fin de boucle sur les cellules
      ENDDO
C
C ---> impression clipping
C
      IF(IRANGP.GE.0) THEN
        CALL PARCPT(CLFP2M)
        CALL PARMAX(MXCFP)
        CALL PARMAX(MAXFP2)
      ENDIF
C
      WRITE(NFECRA,*)' nombre de clip haut sur
     & la variance en f =', CLFP2M
      WRITE(NFECRA,*)' ecart maximum
     & (valeur atteinte-valeur max) =', MXCFP
      WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MAXFP2
      WRITE(NFECRA,*)'     '
C
        IF(IRANGP.GE.0) THEN
        CALL PARCPT(CLYF21)
        CALL PARMAX(MXCYFP)
        CALL PARMAX(MXYFP2)
        ENDIF
C
      WRITE(NFECRA,*)' nombre de clip haut sur
     & la variance en y =', CLYF21
      WRITE(NFECRA,*)' ecart maximum
     & (valeur atteinte-valeur max)  =', MXCYFP
      WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MXYFP2
      WRITE(NFECRA,*)'     '
C
C     WRITE(NFECRA,*)' nombre de clip bas sur
C    & la variance en y =', CLYF22
C     WRITE(NFECRA,*)' ecart maximum
C    & (valeur min-valeur atteinte) =', MNCYFP
C     WRITE(NFECRA,*)' valeur min (pour l ecart max) =', MNYFP2
C     WRITE(NFECRA,*)'     '
C
      IF(IRANGP.GE.0) THEN
        CALL PARCPT(CLCYF1)
        CALL PARMAX(MXCCYF)
        CALL PARMAX(MXCOYF)
      ENDIF
C
      WRITE(NFECRA,*)' nombre de clip haut sur
     & la covariance =', CLCYF1
      WRITE(NFECRA,*)' ecart maximum
     & (valeur atteinte-valeur max)F =', MXCCYF
      WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MXCOYF
      WRITE(NFECRA,*)'     '
C
      IF(IRANGP.GE.0) THEN
        CALL PARCPT(CLCYF2)
        CALL PARMAX(MNCCYF)
        CALL PARMIN(MNCOYF)
      ENDIF
C
      WRITE(NFECRA,*)' nombre de clip bas sur
     & la covariance =', CLCYF2
      WRITE(NFECRA,*)' ecart maximum
     & (valeur min-valeur atteinte) =', MNCCYF
      WRITE(NFECRA,*)' valeur min (pour l ecart max) =', MNCOYF
      WRITE(NFECRA,*)'     '
C
      IF(IRANGP.GE.0) THEN
        CALL PARCPT(CLY2P1)
        CALL PARMAX(MCY2P1)
        CALL PARMAX(MY2P1)
      ENDIF
C
      WRITE(NFECRA,*)' nombre de clip haut sur
     & la variance conditionnelle 1 =', CLY2P1
      WRITE(NFECRA,*)' ecart maximum
     & (valeur atteinte-valeur max) =', MCY2P1
      WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MY2P1
      WRITE(NFECRA,*)'     '
C
      IF(IRANGP.GE.0) THEN
        CALL PARCPT(CLY2P3)
        CALL PARMAX(MCY2P3)
        CALL PARMIN(MY2P3)
      ENDIF
C
      WRITE(NFECRA,*)' nombre de clip bas sur
     & la variance conditionnelle 1 =', CLY2P3
      WRITE(NFECRA,*)' ecart maximum
     & (valeur min-valeur atteinte) =', MCY2P3
      WRITE(NFECRA,*)' valeur min (pour l ecart max) =', MY2P3
      WRITE(NFECRA,*)'     '
C
      IF(IRANGP.GE.0) THEN
        CALL PARCPT(CLY2P2)
        CALL PARMAX(MCY2P2)
        CALL PARMAX(MY2P2)
      ENDIF
C
      WRITE(NFECRA,*)' nombre de clip haut sur
     & la variance conditionnelle 2 =', CLY2P2
      WRITE(NFECRA,*)' ecart maximum
     & (valeur atteinte-valeur max) =', MCY2P2
      WRITE(NFECRA,*)' valeur max (pour l ecart max) =', MY2P2
      WRITE(NFECRA,*)'     '
C
      IF(IRANGP.GE.0) THEN
        CALL PARCPT(CLY2P4)
        CALL PARMAX(MCY2P4)
        CALL PARMIN(MY2P4)
      ENDIF
C
      WRITE(NFECRA,*)' nombre de clip bas sur
     & la variance conditionnelle 2=', CLY2P4
      WRITE(NFECRA,*)' ecart maximum
     & (valeur min-valeur atteinte) =', MCY2P4
      WRITE(NFECRA,*)' valeur min (pour l ecart max) =', MY2P4
      WRITE(NFECRA,*)'     '
C
      IF(IRANGP.GE.0) THEN
        CALL PARCPT(ICPT1)
        CALL PARCPT(ICPT2)
        CALL PARMAX(O2MAX)
        CALL PARMIN(O2MIN)
        CALL PARMAX(TETMAX)
        CALL PARMIN(TETMIN)
        CALL PARMAX(WMAX)
        CALL PARMIN(WMIN)
      ENDIF
C
      WRITE(NFECRA,*) ' POINT NON  PDF = ',ICPT1
      WRITE(NFECRA,*) ' POINT AVEC PDF = ',ICPT2
      WRITE(NFECRA,*) ' MIN MAX O2     = ',O2MIN,O2MAX
      WRITE(NFECRA,*) ' MIN MAX THETA  = ',TETMIN,TETMAX
      WRITE(NFECRA,*) ' MIN MAX W      = ',WMIN,WMAX
C
       END
