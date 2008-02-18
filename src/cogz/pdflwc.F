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
                         SUBROUTINE PDFLWC
C                        *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL  ,
     &   FM     , FP2M  , YFM    , YFP2M ,
     &   PROPCE )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL DES PARAMETRES DE LA PDF
CFONC  PDF LIBBY - WILLIAMS 2 POINTS AVEC HYPOTHESE DE CURL
CFONC
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
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET, NCEL
C
      DOUBLE PRECISION FM(NCELET)   , FP2M(NCELET)
      DOUBLE PRECISION YFM(NCELET)  , YFP2M(NCELET)
      DOUBLE PRECISION PROPCE(NCELET,*)
C
C VARIABLES LOCALES
C
      INTEGER          IEL, IGG, IDIRAC, IPHAS
      INTEGER          MODE
C
      DOUBLE PRECISION COEFG(NGAZGM), EPSI
      DOUBLE PRECISION YFUEL
      DOUBLE PRECISION YOXYD, YO2
      DOUBLE PRECISION YPROD,FMP,FP2MP,YFMP,YFP2MP
C
C --- Expression des droites
C
      DOUBLE PRECISION AA1, AA2
      DOUBLE PRECISION AA3, AA6
      DOUBLE PRECISION BB1, BB2
      DOUBLE PRECISION BB3, BB6
      DOUBLE PRECISION F12, F13, F14
      DOUBLE PRECISION F15, F16
C
C --- Extrema
C
      DOUBLE PRECISION GMIN, GMAX, GP2M
      DOUBLE PRECISION G1MAX, G2MAX, G3MAX
      DOUBLE PRECISION G1MIN, G2MIN
      DOUBLE PRECISION G3MIN, G4MIN
C
C --- Position des Diracs
C
      DOUBLE PRECISION F(NDRACM)    , Y(NDRACM) , D(NDRACM)
      DOUBLE PRECISION G(NDRACM)    , H(NDRACM)
      DOUBLE PRECISION TEML(NDRACM) , MAML(NDRACM)
      DOUBLE PRECISION W(NDRACM)    , RHOL(NDRACM)
      DOUBLE PRECISION THETA(NDRACM)
C
C --- Pointeurs & autres
C
      INTEGER IPCTEM, IPCMAM,  IPCROM
      INTEGER IPAMPL(NDRACM), IPFMEL(NDRACM)
      INTEGER IPFMAL(NDRACM), IPTEML(NDRACM)
      INTEGER IPMAML(NDRACM)
      INTEGER IPRHOL(NDRACM), IPTSCL(NDRACM)
      INTEGER IPCFUE, IPCOXY, IPCPRO, IPCTSC
C      INTEGER IPCKAB, IPT4 , IPT3
      DOUBLE PRECISION NBMOL,  TEMSMM
      DOUBLE PRECISION SUM1, SUM2, SUM3, SUM4,  SUM5,  SUM6 , SUM16
      DOUBLE PRECISION SUM7, SUM8, SUM9, SUM10, SUM11, SUM12, SUM15
      DOUBLE PRECISION SUM17
C
      INTEGER IPASS
      DATA    IPASS /0/
      SAVE    IPASS
C
C***********************************************************************
C
C=======================================================================
C 0.  POSITION DES VARIABLES
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
C        IPT4  = IPPROC(IT4M)
C        IPT3  = IPPROC(IT3M)
C      ENDIF
C
C ---> Initialisation
C
      DO IGG = 1, NGAZGM
        COEFG(IGG) = ZERO
      ENDDO
      EPSI = 1.D-09
C
C=======================================================================
C 1.  CALCULS PRELIMINAIRES
C=======================================================================
C
      IPASS = IPASS + 1
C
      DO IEL = 1, NCEL
C

        FMP    =  MAX(MIN(FMAX,FM(IEL)),FMIN)
        YFMP = FMP-((FMP-MAX(ZERO,(FMP-FS(1))/(1.D0-FS(1))))
     &                  /(FM(IEL)-MAX(ZERO,
     &       (FM(IEL)-FS(1))/(1.D0-FS(1))))*(FM(IEL)-YFM(IEL)))
        YFMP   =  MAX(MIN(YFM(IEL),FMP),
     &                    ZERO,(FMP-FS(1))/(1.D0-FS(1)))
        FP2MP  =  MAX(MIN(FP2M(IEL),
     &               (MIN(FMAX,(1-FS(1))*YFMP+FS(1))-FMP)
     &             *      (FMP-MAX(FMIN,YFMP))),ZERO)
        YFP2MP =  MAX(MIN(YFP2M(IEL),(FMP-YFMP)
     &             *      (YFMP-MAX(ZERO,
     &                    (FMP-FS(1))/(1.D0-FS(1))))),ZERO)
C
C
C======================================================================
C 2. NON PASSAGE PAR LA PDF
C======================================================================
C
        IF ( (FP2MP.LE.EPSI).AND.(YFP2MP.LE.EPSI )) THEN
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
C ---> Calcul de l'enthalpie
C
            H(IDIRAC) = ((HMAX-HMIN)*F(IDIRAC) + HMIN*FMIN - HMAX*FMAX)
     &                / (FMIN-FMAX)
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
     &      ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &        NPO    , NPOT  , TH      , EHGAZG ,
     &        H(IDIRAC)     , TEML(IDIRAC)    )
C
C ---> Calcul de la masse volumique en 1 et 2
C
            IF ( IPASS.GT.1.OR.
     &          (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
              RHOL(IDIRAC) = P0(IPHAS) * MAML(IDIRAC)
     &                     / (RR*TEML(IDIRAC))
            ELSE
              RHOL(IDIRAC) = RO0(IPHAS)
            ENDIF
C
C ---> Calcul du terme source en 1 et 2 du scalaire YFM
C
            THETA(IDIRAC) = TA / TEML(IDIRAC)
     &                    * (1.D0 - TEML(IDIRAC) / TSTAR)
C
            W(IDIRAC) = VREF / LREF * (- D(IDIRAC)*RHOL(IDIRAC)
     &                  * YFUEL*YO2
     &                  * DEXP( -THETA(IDIRAC) ))
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
C ---> Masse volumique du melange
C
          IF ( IPASS.GT.1 .OR.
     &        (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
            PROPCE(IEL,IPCROM) = SRROM * PROPCE(IEL,IPCROM)
     &                         + (1.D0-SRROM) * (P0(IPHAS)/(RR*TEMSMM))
          ENDIF
C
        ELSE
C
C======================================================================
C 3.  PASSAGE PAR LA PDF
C======================================================================
C
          IF  ( ( (FP2MP.GT.EPSI) .AND. (YFP2MP.LT.EPSI) )
     &      .OR.( (FP2MP.LT.EPSI) .AND. (YFP2MP.GT.EPSI) ) )THEN
C
C -- Choix suivant les differents cas
C
C --- CAS 1 : Droite Verticale
C     =====
            IF (FP2MP .LT. EPSI) THEN
C
C ---> Choix des extrema
C
              GMIN =   MAX ( -  YFMP,
     &                    -  (YFMP-(FMP-FS(1))/(1.D0-FS(1))))
              GMAX = FMP - YFMP
C
C ---> Calcul des amplitudes des Diracs
C
              D(1) = GMAX / ( GMAX-GMIN )
              D(2) = (1.D0 - D(1))
C
C ---> Calcul de la variance de l'abscisse curviligne (GP2M)
C
              GP2M = FP2MP + YFP2MP
C
C ---> Test sur GP2M
C
              G(1) = -DSQRT( -GMIN/GMAX*GP2M )
              G(2) =  -D(1) * G(1) / D(2)
C
              DO IDIRAC = 1, NDIRAC
                F(IDIRAC) =
     &                    FMP + G(IDIRAC) * DSQRT( FP2MP/GP2M)
                Y(IDIRAC) =
     &                    YFMP + G(IDIRAC) * DSQRT( YFP2MP/GP2M)
              ENDDO
C
C --- CAS 2 : Droite Horizontale
C     =====
C
            ELSEIF    (YFP2MP .LT. EPSI) THEN
C
C ---> Calcul des differentes intersections des droites
C
              F12 = YFMP
              F13 = FS(1) + YFMP * (1-FS(1))
              F14 = FMIN
              F15 = FMAX
C
C ---> Calcul des extrema de l'abscisse curviligne
C      Appel de la function G
C
              CALL LWCGFU(G2MAX, F15, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G3MAX, F13, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G1MIN, F12, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G2MIN, F14, FMP, YFP2MP, FP2MP)
C
C ---> Choix des extrema
C
              GMIN =   MAX(  G1MIN,  G2MIN)
              GMAX =   MIN(  G2MAX,  G3MAX)
C ---> Calcul des amplitudes des Diracs
C
              D(1) = GMAX / ( GMAX-GMIN )
              D(2) = (1.D0 - D(1))
C
C ---> Calcul de la variance de l'abscisse curviligne (GP2M)
C
              GP2M = FP2MP + YFP2MP
C
              G(1) = -DSQRT( -GMIN/GMAX*GP2M )
              G(2) =  -D(1) * G(1) / D(2)
C
              DO IDIRAC = 1, NDIRAC
                F(IDIRAC) =
     &              FMP + G(IDIRAC) * DSQRT( FP2MP/GP2M)
                Y(IDIRAC) =
     &              YFMP + G(IDIRAC) * DSQRT( YFP2MP/GP2M)
              ENDDO
            ENDIF
C
          ELSEIF ( (FP2MP.GT.EPSI) .AND. (YFP2MP.GT.EPSI) ) THEN
C
C --- CAS 3 : Parallelisme avec la Droite de Melange
C     =====
C
            IF ( ((YFP2MP / FP2MP) .LT. 1.D0 + EPSI)
     &        .AND. ((YFP2MP / FP2MP) .GT. 1.D0 - EPSI) ) THEN
C
              AA1 = 1.D0
              BB1 = YFMP - FMP
C
              AA3 = 1.D0 / (1.D0-FS(1))
              BB3 = -FS(1) / (1.D0-FS(1))
C
              AA6 = ZERO
              BB6 = ZERO
C
C ---> Calcul des differentes intersections de ces droites
C
              F13 = (BB3-BB1) / (AA1-AA3)
              F14 = FMIN
              F15 = FMAX
              F16 = (BB6-BB1) / (AA1-AA6)
C
C ---> Calcul des extrema de l'abscisse curviligne
C      Appel de la function G
C
              CALL LWCGFU(G2MAX,F15, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G3MAX,F13, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G2MIN,F14, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G3MIN,F16, FMP, YFP2MP, FP2MP)
C
              GMIN =   MAX(  G2MIN,  G3MIN)
              GMAX =   MIN(  G2MAX,  G3MAX)
C
C ---> Calcul des amplitudes des Diracs
C
              D(1) = GMAX / ( GMAX-GMIN )
              D(2) = (1.D0 - D(1))
C
C ---> Calcul de la variance de l'abscisse curviligne (GP2M)
C
              GP2M = FP2MP + YFP2MP
C
C ---> Test sur GP2M
C
              G(1) = -DSQRT( -GMIN/GMAX*GP2M )
              G(2) =  -D(1) * G(1) / D(2)
C
              DO IDIRAC = 1, NDIRAC
                F(IDIRAC) =
     &              FMP + G(IDIRAC) * DSQRT( FP2MP/GP2M)
                Y(IDIRAC) =
     &              YFMP + G(IDIRAC) * DSQRT( YFP2MP/GP2M)
              ENDDO
C
C --- CAS 4 : Parallelisme avec la Droite de Combustion Complete
C     =====
            ELSEIF ( ((SQRT(YFP2MP/FP2MP) * (1.D0-FS(1)))
     &                                        .LT. 1.D0 + EPSI)
     &        .AND.  ((SQRT(YFP2MP/FP2MP) * (1.D0-FS(1)))
     &                                        .GT. 1.D0 - EPSI) ) THEN
C
              AA1 = DSQRT( YFP2MP/FP2MP )
              BB1 = YFMP - DSQRT( YFP2MP/FP2MP ) * FMP
C
              AA2 = 1.D0
              BB2 = ZERO
C
              AA6 = ZERO
              BB6 = ZERO
C
C ---> Calcul des differentes intersections de ces droites
C
              F12 = (BB2-BB1) / (AA1-AA2)
              F14 = FMIN
              F15 = FMAX
              F16 = (BB6-BB1) / (AA1-AA6)
C
C ---> Calcul des extrema de l'abscisse curviligne
C      Appel de la function G
C
              CALL LWCGFU(G2MAX,F15, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G1MAX,F12, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G2MIN,F14, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G3MIN,F16, FMP, YFP2MP, FP2MP)
C
              GMIN =   MAX(  G2MIN,  G3MIN)
              GMAX =   MIN(  G1MAX,  G2MAX)
C
C ---> Calcul des amplitudes des Diracs
C
              D(1) = GMAX / ( GMAX-GMIN )
              D(2) = (1.D0 - D(1))
C
C ---> Calcul de la variance de l'abscisse curviligne (GP2M)
C
              GP2M = FP2MP + YFP2MP
C
C ---> Test sur GP2M
C
              G(1) = -DSQRT( -GMIN/GMAX*GP2M )
              G(2) =  -D(1) * G(1) / D(2)
C
              DO IDIRAC = 1, NDIRAC
                F(IDIRAC) =
     &              FMP + G(IDIRAC) * DSQRT( FP2MP/GP2M)
                Y(IDIRAC) =
     &             YFMP + G(IDIRAC) * DSQRT( YFP2MP/GP2M)
              ENDDO
C
C --- CAS GENERAL
C     ===========
            ELSE
C
C ---> Expression des differentes droites: Y=AX+B
C
C
              AA1 = -DSQRT( YFP2MP/FP2MP )
              BB1 = YFMP - AA1 * FMP
C
              AA2 = 1.D0
              BB2 = ZERO
C
              AA3 = 1.D0 / (1.D0-FS(1))
              BB3 = -FS(1) / (1.D0-FS(1))
C
              AA6 = ZERO
              BB6 = ZERO
C
C ---> Calcul des differentes intersections de ces droites
C
              F12 = (BB2-BB1) / (AA1-AA2)
              F13 = (BB3-BB1) / (AA1-AA3)
              F14 = FMIN
              F15 = FMAX
              F16 = (BB6-BB1) / (AA1-AA6)
C
C ---> Calcul des extrema de l'abscisse curviligne
C      Appel de la function G
C
              CALL LWCGFU(G1MAX, F12, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G2MAX, F15, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G3MAX, F13, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G1MIN, F12, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G2MIN, F14, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G3MIN, F16, FMP, YFP2MP, FP2MP)
              CALL LWCGFU(G4MIN, F13, FMP, YFP2MP, FP2MP)
C
C=======================================================================
C 4.  CALCUL DES PARAMETRES DES DEUX PICS DE DIRAC
C=======================================================================
C
C ---> Choix des extrema suivant les pentes des droites
C
              IF ( AA1 .GT. AA3 ) THEN
                GMIN =   MAX (  G2MIN,  G3MIN,  G4MIN )
                GMAX =   MIN ( G1MAX, G2MAX )
C
              ELSEIF (( AA1 .GT. AA2 ).AND. ( AA1 .LE. AA3 )) THEN
                GMIN =   MAX (  G2MIN,  G3MIN )
                GMAX =   MIN ( G1MAX, G2MAX, G3MAX )
C
              ELSEIF ( AA1 .LE. AA2 ) THEN
                GMIN =   MAX (  G1MIN,  G2MIN,  G3MIN )
                GMAX =   MIN ( G2MAX, G3MAX )
              ENDIF
C
C ---> Calcul des amplitudes des Diracs
C
              D(1) = GMAX / ( GMAX-GMIN )
              D(2) = (1.D0 - D(1))
C
C ---> Calcul de la variance de l'abscisse curviligne (GP2M)
C
              GP2M = FP2MP + YFP2MP
C
              G(1) = -DSQRT( GMIN/GMAX*GP2M )
              G(2) =  -D(1) * G(1) / D(2)
C
              DO IDIRAC = 1, NDIRAC
                F(IDIRAC) =
     &              FMP + G(IDIRAC) * DSQRT( FP2MP/GP2M)
                Y(IDIRAC) =
     &              YFMP + G(IDIRAC) * DSQRT( YFP2MP/GP2M)
              ENDDO

            ENDIF
          ENDIF
C
C=======================================================================
C 5.  DETERMINATION DES GRANDEURS THERMOCHIMIQUES DES DEUX PICS
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
            H(IDIRAC) = ((HMAX-HMIN)*F(IDIRAC) + HMIN*FMIN - HMAX*FMAX)
     &                / (FMIN-FMAX)
C
C ---> Calcul de la fraction massique des gaz (F, O et P) en 1 et 2
C
            YFUEL = Y(IDIRAC)
            YOXYD = 1.D0 - (COEFF3+1.0D0)*F(IDIRAC) + COEFF3*Y(IDIRAC)
            YPROD = 1.D0 - YFUEL - YOXYD
            YO2   = COEFF1 - (COEFF1 + COEFF2) * F(IDIRAC)
     &                      + COEFF2 * Y(IDIRAC)
C
C ---> Coefficients d'absorption pour les pics 1 et 2
C
c            IF ( IRAYPP .GT. 0  ) THEN
c              KABSGF = YFUEGF(IEL)*KABSG(1) + YOXYGF(IEL)*KABSG(2)
c     &                                      + YPROGF(IEL)*KABSG(3)
c              KABSGB = YFUEGB(IEL)*KABSG(1) + YOXYGB(IEL)*KABSG(2)
c     &                                      + YPROGB(IEL)*KABSG(3)
c            ENDIF
C
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
     &      ( MODE   , NGAZG , NGAZGM  , COEFG  ,
     &        NPO    , NPOT  , TH      , EHGAZG ,
     &        H(IDIRAC)     , TEML(IDIRAC)    )
C
C ---> Calcul de la masse volumique en 1 et 2
C
            IF ( IPASS.GT.1 .OR.
     &          (ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1) ) THEN
              RHOL(IDIRAC) = P0(IPHAS) * MAML(IDIRAC)
     &                     / (RR*TEML(IDIRAC))
            ELSE
              IPHAS = 1
              RHOL(IDIRAC) = RO0(IPHAS)
            ENDIF
C
C ---> Calcul du terme source en 1 et 2 du scalaire YFM
C
            THETA(IDIRAC) = TA / TEML(IDIRAC)
     &                    * (1.D0 - TEML(IDIRAC) / TSTAR)
C
            W(IDIRAC) = VREF / LREF * (- D(IDIRAC)*RHOL(IDIRAC)
     &                  * YFUEL*YO2
     &                  * DEXP( -THETA(IDIRAC) ))
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
            SUM9 = SUM9 + D(IDIRAC)*TEML(IDIRAC)/MAML(IDIRAC)
C
C ---> Fractions massiques des especes globales
C
            SUM10 = SUM10 + YFUEL*D(IDIRAC)
C
            SUM11 = SUM11 + YOXYD*D(IDIRAC)
C
            SUM12 = SUM12 + YPROD*D(IDIRAC)
C
            SUM17 = SUM17 +W(IDIRAC)
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
c            IF ( IRAYPP .GT. 0  ) THEN
c              PROPCE(IEL,IPCKAB) = YGFM*KABSGF + YGBM*KABSGB
c              PROPCE(IEL,IPT4)    = YGFM*TGF**4+YGBM*TGB**4
c              PROPCE(IEL,IPT3)    = YGFM*TGF**3+YGBM*TGB**3
c            ENDIF
C
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
        IF (IPASS.GT.1.OR.(ISUITE.EQ.1.AND.INITRO(IPHAS).EQ.1)) THEN
            PROPCE(IEL,IPCROM) = SRROM * PROPCE(IEL,IPCROM)
     &                       + (1.D0-SRROM) * (P0(IPHAS)/(RR*TEMSMM))
          ENDIF
C
        ENDIF
C
      ENDDO
C
C----
C FIN
C----
C
      END
c@z
