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
                        SUBROUTINE FULECD
C                       *****************
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  LECTURE DU FICHIER DE DONNEES PHYSIQUE PARTICULIERE
CFONC       RELATIF A LA COMBUSTION FUEL
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
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
      INCLUDE "pointe.h"
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
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
C VARIABLES LOCALES
C
      CHARACTER *150 CHAIN1,CHAIN2
      CHARACTER *12   NOMCOE(NGAZEM)
C
      INTEGER          IT     , ICE    , IAT    , IOS
      INTEGER          NCOEL  , INICOE
      INTEGER          IDEBCH , IFINCH , LONCH  , ICHAI  , ICHCOE
      INTEGER          ATCOEL(NGAZEM,NATOM), INICHA
C
      DOUBLE PRECISION TMIN   , TMAX
      DOUBLE PRECISION WMOLCE(NGAZEM), EHCOEL(NGAZEM,NPOT)
      DOUBLE PRECISION CPCOEL(NGAZEM,NPOT)
      DOUBLE PRECISION nCFOV,nHFOV,nOFOV,nSFOV
      DOUBLE PRECISION mHSFOV,mCOFOV,mCHFOV,mtoFOV
      DOUBLE PRECISION nHSFOV,nCOFOV,nCMV,nHMV
      DOUBLE PRECISION CH2FV,CH4FV,H02FOV,YO2Ox
C
C***********************************************************************
C==========================================
C 1. LECTURE DU FICHIER DONNEES SPECIFIQUES
C==========================================
C
C --> Ouverture du fichier
C
      OPEN ( UNIT=IMPFPP, FILE=FICFPP,
     &        STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL',
     &                                        IOSTAT=IOS, ERR=99 )
      REWIND (UNIT=IMPFPP,ERR=99 )
C
C --> Lecture thermochimie
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
C ---- Nb de constituants elementaires (gazeux,liquide et solide)
C
      READ ( IMPFPP,*,ERR=999,END=999 ) NCOEL
      IF ( NCOEL.GT.NGAZGM ) THEN
        WRITE(NFECRA,9991) NGAZGM,NCOEL
        CALL CSEXIT (1)
      ENDIF
C
C ---- Nb de points de tabulation ENTH-TEMP
C
      READ ( IMPFPP,*,ERR=999,END=999 ) NPO
      IF ( NPO.GT.NPOT ) THEN
        WRITE(NFECRA,9992) NPOT,NPO
        CALL CSEXIT (1)
      ENDIF
C
C --- Lecture des noms des constituants elementaires
C
      DO ICE=1,NCOEL
        DO INICOE=1,LEN(NOMCOE(ICE))
          NOMCOE(ICE)(INICOE:INICOE)=' '
        ENDDO
      ENDDO
C
      DO INICHA=1,LEN(CHAIN1)
        CHAIN1(INICHA:INICHA)=' '
      ENDDO
C
      DO INICHA=1,LEN(CHAIN2)
        CHAIN2(INICHA:INICHA)=' '
      ENDDO
C
      READ (IMPFPP,*,ERR=999,END=999)
c     READ (IMPFPP,*,ERR=999,END=999) CHAIN2
      READ (IMPFPP,1010,ERR=999,END=999 ) CHAIN1
      CALL VERLON (CHAIN1, IDEBCH, IFINCH, LONCH)
      CHAIN2(1:LONCH)=CHAIN1(IDEBCH:IFINCH)
C
      ICE=1
      ICHCOE=0
      DO ICHAI = 1, LONCH
        IF (CHAIN2(ICHAI:ICHAI).NE.' ') THEN
          ICHCOE=ICHCOE+1
          NOMCOE(ICE)(ICHCOE:ICHCOE) =CHAIN2(ICHAI:ICHAI)
        ELSE
          IF (ICHCOE.NE.0) THEN
            ICE=ICE+1
            ICHCOE=0
          ENDIF
        ENDIF
      ENDDO
C
 1010 FORMAT(A150)
C
C --- Temperature Min et Max
C
      READ (IMPFPP,*,ERR=999,END=999) TMIN
      READ (IMPFPP,*,ERR=999,END=999) TMAX
C
C
C ---- Nb especes atomiques (C, H, O, N, ...)
C
      READ (IMPFPP,*,ERR=999,END=999 ) NATO
      IF ( NATO.GT.NATOM ) THEN
        WRITE(NFECRA,9993) NATOM,NATO
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C ---- Masse molaire especes atomiques
C      Composition des constituants elementaires en fonction
C        des especes elementaires
C
      DO IAT = 1, NATO
        READ (IMPFPP,*,ERR=999,END=999 ) WMOLAT(IAT),
     &                      ( ATCOEL(ICE,IAT),ICE=1,NCOEL )
      ENDDO
C
C ---- Calcul des masses molaires des constituants elementaires
C
      DO ICE = 1, NCOEL
        WMOLCE(ICE) = 0.D0
        DO IAT = 1, NATO
          WMOLCE(ICE)= WMOLCE(ICE) + ATCOEL(ICE,IAT)*WMOLAT(IAT)
        ENDDO
      ENDDO
C
C
C --> Lecture rayonnement
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
C      IRAYPP =  0 : pas de rayonnement,
C             =  1 : DOM + calcul a partir des CABS des especes
C             =  2 : DOM + CABS par Modak
C             =  3 : P-1 + calcul a partir des CABS des especes
C             =  4 : P-1 + CABS par Modak
C
      READ ( IMPFPP,*,ERR=999,END=999 ) IRAYPP
      IF ( IRAYPP.GT.4 .OR. IRAYPP.LT.0 ) THEN
        WRITE(NFECRA,9994) IRAYPP
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C ---- Coefficient d'absorption du melange gazeux
C
      READ (IMPFPP,*,ERR=999,END=999 ) CKABS1
C
C
C --> Lecture caracteristiques fuel
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
C --> Diametre initial  (mm)
C
      READ (IMPFPP,*,ERR=999,END=999 ) DINIFL
C
C --> Composition elementaire en C, H, O, S, In (% en masse)
C     In désigne les inertes (métaux, etc.) qui resteront
C        dans le résidu solide
C
      READ (IMPFPP,*,ERR=999,END=999 ) CFOL
      READ (IMPFPP,*,ERR=999,END=999 ) HFOL
      READ (IMPFPP,*,ERR=999,END=999 ) OFOL
      READ (IMPFPP,*,ERR=999,END=999 ) SFOL
C
      CFOL = 1.D-2 * CFOL
      HFOL = 1.D-2 * HFOL
      OFOL = 1.D-2 * OFOL
      SFOL = 1.D-2 * SFOL
      XINFOL = 1.D0-CFOL-HFOL-OFOL-SFOL
      if (XINFOL .LT. ZERO) then
         WRITE(NFECRA,*)'Erreur dans les fractions massiques du FOL'
C         STOP
      endif
      WRITE (NFECRA,*) 'Fractions massiques elementaires / FOL  '
      WRITE (NFECRA,*) ' C = ',CFOL
      WRITE (NFECRA,*) ' H = ',HFOL
      WRITE (NFECRA,*) ' O = ',OFOL
      WRITE (NFECRA,*) ' S = ',SFOL
      WRITE (NFECRA,*) ' In= ',XINFOL
C
C
C --> PCI
C
      READ (IMPFPP,*,ERR=999,END=999 ) PCIFOL
C
C --> CP moyen du fuel sec (J/kg/K)
C
      READ (IMPFPP,*,ERR=999,END=999 ) CP2FOL
C
C --> Masse volumique initiale (kg/m3)
C
      READ (IMPFPP,*,ERR=999,END=999 ) RHO0FL
C
C --> Caracteristiques du coke
C
      READ (IMPFPP,*,ERR=999,END=999)
C
C ------- Composition elementaire en C, H, O, S (% / pur)
C
      READ (IMPFPP,*,ERR=999,END=999 ) CKF
      READ (IMPFPP,*,ERR=999,END=999 ) HKF
      READ (IMPFPP,*,ERR=999,END=999 ) OKF
      READ (IMPFPP,*,ERR=999,END=999 ) SKF
C
      CKF = 1.D-2 * CKF
      HKF = 1.D-2 * HKF
      OKF = 1.D-2 * OKF
      SKF = 1.D-2 * SKF
C
      IF ( ABS(CKF+HKF+OKF+SKF-1.D0) .GT. 1.D-15 ) THEN
        WRITE(NFECRA,9990) CKF+HKF+OKF+SKF
        CALL CSEXIT(1)
      ENDIF
C
C ------ PCI
C
      READ (IMPFPP,*,ERR=999,END=999 ) PCIKF
C
C ---- Fraction de coke dans le fuel
C
      READ (IMPFPP,*,ERR=999,END=999) FKC
      WRITE (NFECRA,*)' Fraction massique de coke / FOL'
      WRITE (NFECRA,*) FKC
C
C     Les inertes restent dans le coke
      XINKF = ZERO
      IF ( FKC .GT. ZERO) XINKF = XINFOL/FKC
      IF ( (CKF+HKF+OKF+SKF) .GT. 1.D0) THEN
         WRITE(NFECRA,*)'Erreur dans les fractions massiques du KF'
C         STOP
      ENDIF
C
      WRITE (NFECRA,*) 'Fractions massiques elementaires / coke '
      WRITE (NFECRA,*) ' C = ',CKF*(1.D0-XINKF)
      WRITE (NFECRA,*) ' H = ',HKF*(1.D0-XINKF)
      WRITE (NFECRA,*) ' O = ',OKF*(1.D0-XINKF)
      WRITE (NFECRA,*) ' S = ',SKF*(1.D0-XINKF)
      WRITE (NFECRA,*) ' In= ',XInKF
C
C     Compatibilite des fractions massiques et des formules moléculaires
C     masses élémentaires dans le fuel, le coke, les vapeurs
C        F      K        MV
C   C    CFOL   CKF*FKC  CFOL-CKF*FKC
C   H    HFOL   HKF*FKC  HFOL-HKF*FKC
C   O    OFOL   OKF*FKC  OFOL-OKF*FKC
C   S    SFOL   SKF*FKC  SFOL-SKF*FKC
C   In   XInFOL  XInFOL    0
C      elements dans les vapeurs
      NCFOV  = (CFOL-CKF*FKC*(1.D0-XINKF))/WMOLAT(IATC)/(1.D0-FKC)
      NHFOV  = (HFOL-HKF*FKC*(1.D0-XINKF))/WMOLAT(IATH)/(1.D0-FKC)
      NOFOV  = (OFOL-OKF*FKC*(1.D0-XINKF))/WMOLAT(IATO)/(1.D0-FKC)
      NSFOV  = (SFOL-SKF*FKC*(1.D0-XINKF))/WMOLAT(IATS)/(1.D0-FKC)
C       on considére que S se dégage sous forme H2S
C                    que O                      CO
      NHSFOV = NSFOV
      NCOFOV = NOFOV
      NCMV   = NCFOV - NCOFOV
      NHMV   = NHFOV - 2.D0*NHSFOV
C
C   Les vapeurs sont alors constituées de nHSFOV moles de H2S
C                                         nCOFOV          CO
C                                         nCMV            CHn
C   où CHn est un hydrocarbure modèle de formule moyenne avec
      NHCFOV  = NHMV/NCMV
      WRITE(NFECRA,*) ' nHCFOV = ',NHCFOV ,NHMV,NCMV
C
C   Les masses dans les vapeurs sont
      MHSFOV = (WMOLAT(IATS)+2.D0*WMOLAT(IATH))*NHSFOV
      MCOFOV = (WMOLAT(IATC)+WMOLAT(IATO))*NCOFOV
      MCHFOV = WMOLAT(IATC)*NCMV+WMOLAT(IATH)*NHMV
      MTOFOV = MHSFOV+MCOFOV+MCHFOV
C
      WRITE(NFECRA,*) ' mtoFOV = ',MTOFOV
C
C   Les fractions massiques dans les vapeurs sont
      HSFOV = MHSFOV / MTOFOV
      COFOV = MCOFOV / MTOFOV
      CHFOV = MCHFOV / MTOFOV
      WRITE (NFECRA,*) 'Fractions massiques spécifiques / FOV '
      WRITE (NFECRA,*) ' H2S = ',HSFOV
      WRITE (NFECRA,*) ' CO  = ',COFOV
      WRITE (NFECRA,*) ' CHn = ',CHFOV
      WRITE (NFECRA,*) ' ..n = ',nHCFOV
      CH4FV = ZERO
      CH2FV = CHFOV
      IF ( NHCFOV.GE.2.D0 .AND. NHCFOV.LE.4.D0 ) THEN
        WRITE(NFECRA,*) 'Le FOV est equivalent a un melange '
        CH2FV = 2.D0-0.5D0*NHCFOV
         CH4FV = (1-CH2FV)*16.D0/(12.D0+NHCFOV)
         CH2FV = CH2FV*14.D0/(12.D0+NHCFOV)
         CH4FV = CH4FV * CHFOV
         CH2FV = CH2FV * CHFOV
         WRITE (NFECRA,*) ' H2S = ',HSFOV
         WRITE (NFECRA,*) ' CO  = ',COFOV
         WRITE (NFECRA,*) ' CH4 = ',CH4FV
         WRITE (NFECRA,*) 'C2H4 = ',CH2FV
      ENDIF
      WRITE(NFECRA,*) ' nHCFOV 2 = ',NHCFOV
C
C ---- Parametre d'evaporation
C
       READ (IMPFPP,*,ERR=999,END=999) TEVAP1
       READ (IMPFPP,*,ERR=999,END=999) TEVAP2
C
C ---- Parametres combustion heterogene (modele a sphere retrecissante)
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
      READ (IMPFPP,*,ERR=999,END=999 ) AHETFL
      READ (IMPFPP,*,ERR=999,END=999 ) EHETFL
      READ (IMPFPP,*,ERR=999,END=999 ) IOFHET
C
C --> Fermeture du fichier (ne pas oublier, car l'unite sert pour janaf)
C
      CLOSE(IMPFPP)
C
C======================================
C 2.
C======================================
C
C
C --> Discretisation de la temperature
C
      DO IT = 1, NPO
        TH(IT) = DBLE(IT-1)*(TMAX-TMIN)/DBLE(NPO-1) + TMIN
      ENDDO
C
C --> Calcul des enthalpies pour les differentes especes courantes
C
      CALL PPTBHT
C     ===========
     & ( NCOEL  ,
     &   NOMCOE , EHCOEL , CPCOEL , WMOLCE )
C
C --> Calcul tabulation enthalpie - temperature pour le melange gazeux
C
C ---- Nb de constituants gazeux
C     ATTENTION ON COMPTE EGALEMENT H2S et le monomere SO2
C
      NGAZE = 8
C
C ---- Definition des pointeurs pour les tableaux WMOLE et EHGAZE
C      REMARQUE : Cette position de pointeurs va egalement servir
C                 pour le tableau de pointeurs IYM1 relatif aux
C                 tableaux PROPCE et PROPFB
C                 ON BALAYE JUSTE DE 1 A NGAZE
C
C     ATTENTION : ordre des espèces dans EHCOEL, WMOLCE
C                 vient du fichier data_FUE
C     Actuellement 1   CH4
C                  2   C2H4
C                  3   CO
C                  4   O2
C                  5   CO2
C                  6   H2O
C                  7   N2
C                  8   C (solide)
C                  9   H2S
C                 10   SO2
      IFOV = 1
      ICO  = 2
      IH2S = 3
      IH2O = 4
      ICO2 = 5
      ISO2 = 6
      IO2  = 7
      IN2  = 8
C
C ---- Remplissage de EHGAZE et WMOLE
C         a partir de EHCOEL et WMOLCE
C
      DO IT = 1, NPO
        EHGAZE(IFOV ,IT) = ( CH4FV*EHCOEL(1,IT) + CH2FV*EHCOEL(2,IT) )
        EHGAZE(ICO  ,IT) = EHCOEL( 3,IT)
        EHGAZE(IO2  ,IT) = EHCOEL( 4,IT)
        EHGAZE(ICO2 ,IT) = EHCOEL( 5,IT)
        EHGAZE(IH2O ,IT) = EHCOEL( 6,IT)
        EHGAZE(IN2  ,IT) = EHCOEL( 7,IT)
        EHGAZE(IH2S ,IT) = EHCOEL( 9,IT)
        EHGAZE(ISO2 ,IT) = EHCOEL(10,IT)
      ENDDO
      WMOLE(IFOV ) = (CH4FV+CH2FV)/(CH4FV/WMOLCE(1)+CH2FV/WMOLCE(2))
C
      WMOLE(IFOV ) = (1.D0*0.012D0 + nHCFOV *0.001D0 )
      WRITE(NFECRA,*) ' Wmole IFOV = ',WMOLE(IFOV ),CH4FV,CH2FV
      WMOLE(ICO  ) = WMOLCE( 3)
      WMOLE(IO2  ) = WMOLCE( 4)
      WMOLE(ICO2 ) = WMOLCE( 5)
      WMOLE(IH2O ) = WMOLCE( 6)
      WMOLE(IN2  ) = WMOLCE( 7)
      WMOLE(IH2S ) = WMOLCE( 9)
      WMOLE(ISO2 ) = WMOLCE(10)
C
C     Concentrations dans les vapeurs
      AFOVF1 = CHFOV / WMOLE(IFOV)
      ACOF1  = COFOV / WMOLE(ICO)
      AH2SF1  = HSFOV / WMOLE(IH2S)
C      Caractérisation de l'oxydant
      YO2Ox = WMOLE(IO2) / (WMOLE(IO2)+XSI*WMOLE(IN2))
C     Caractéristion des gaz issus de la combustion hétérogène
C      Max est le point où F3 est maximal ; il correspond à un
C      mélange stoechiometrique de coke et d'oxydant
C      On suppose, pour l'instant, que l'oxydant est un melange O2, N2
C      à modifier si recyclage de fumées
C      FF3MAX(CKF*C+OKF*O+SKF*S+HKF*H) + (1-F3max)*(YO2Ox*O2+(1-YO2Ox)*N2) =>
C            CO + H2O + H2S + N2
C      On suppose que S est prioritaire pour H
C      Masse de C dans les réactants FF3MAX*CKF
C      Masse de O                    FF3MAX*OKF + (1-F3max)*YO2Ox
C      Masse de S                    FF3MAX*SKF
C      Masse de H                    FF3MAX*HKF
C      Masse de N                    (1-FF3MAX)*(1-YO2Ox)
C      Nombre de C                   FF3MAX*CKF / WMOLAT(IATC)
C      Nombre de O                   (FF3MAX*OKF+(1-F3max*YO2Ox)/WMOLAT(IATO)
C      Nombre de S                   FF3MAX*SKF / WMOLAT(IATS)
C      Nombre de H                   FF3MAX*HKF / WMOLAT(IATH)
C      Nombre de moles de H2S dans les produits = nombre de moles de S
C      Nombre de moles de H2O dans les produits = 1/2 nombre H-nb moles H2S
C      Nombre de moles de CO = nombre de moles de C
C      Nombre de moles de O dans les produits = nb moles CO + H2O
C      Et il vient :
C      FF3MAX * ( CKF/WMOLAT(IATC) + 0.5*HKF/WOLAT(IATH) - SKF/WMOLAT(IATS)
C                -OKF/WMOALT(IATO) - YO2Ox/WOLAT(IATO) )
C      =     YO2Ox/WMOLAT(IATO)
C
C      FF3MAX = YO2Ox /(0.016*(CKF/.012 + 0.5*HKF/.001 -SKF/.032)-OKF+YO2Ox)
C
       FF3MAX = YO2OX/(WMOLAT(IATO)*(CKF/WMOLAT(IATC)
     &                       +0.5D0*HKF/WMOLAT(IATH)
     &                             -SKF/WMOLAT(IATS))
     &                -OKF + YO2OX)
C
C      AXXF3 nb de moles de l'espèce XX en kilog de F3 calcules en FF3MAX
       AH2SF3 = FF3MAX*SKF/WMOLAT(IATS)
       AH2OF3 = 0.5d0*FF3MAX*HKF/WMOLAT(IATH)-AH2SF3
       ACOF3  = FF3MAX*CKF/WMOLAT(IATC)
C
C       AO2F3  = 0.5d0*FF3MAX*(     OKF/WMOLAT(IATO)-CKF/WMOLAT(IATC)
C     &                       -0.5*HKF/WMOLAT(IATH)+SKF/WMOLAT(IATS))
       AO2F3  = 0.D0
C
       AN2F3  = (1.D0-YO2OX)/WMOLE(IN2) * (1.D0-FF3MAX)
C
       AO2F4 = YO2OX/WMOLE(IO2)
       AN2F4 = (1.D0-YO2OX)/WMOLE(IN2)
C
C      Avec cette convention la concentration en O2 est négative en FF3MAX
C      ceci correspond à l'oxygène absorbé par la réaction hétérogène.
C
C --> Calcul tabulation enthalpie - temperature pour la phase dispersee
C     Fuel Oil Liquid et  Coke
C
C ---- Nb de constituants solide
C
      NSOLID = 2
C
C ---- Definition des pointeurs IFOL et IKF
C
      IFOL = 1
      IKF = 2
C
C ------ Calcul de H02FOL
C
C       H0, EH & PCI en J/kg
C       CFOL, HFOL sont des fractions massiques élémentaires
C       rapports des masses molaires des produits aux éléments du
C                combustible (le comburant est dans l'état de réf.)
C
C ------ Calcul de HRFVAP
C
C  L'enthalpie de formation du fuel gazeux est connue (mélange CH4, C2H4),
C  Le PCI du fuel liquide est connu , on peut donc reconstituer son
C   enthalpie de formation (on néglige l'effet de H2S => SO2)
C   on introduit les enthalpies de formation massique du CO2 et de H2O
C
        H02FOL = PCIFOL
     &         + CFOL * 44.D0/12.D0 * EHCOEL(5,1)
     &         + HFOL * 18.D0/2.D0  * EHCOEL(6,1)
C       H02FOL en J/kg (de fol)
C       L'enthalpie de formation de la vapeur de fuel
C       est supposée etre des celle des seuls hydrocarbures
C       (i.e. on néglige, pour l'instant, CO et H2S)
        H02FOV = CH4FV * EHCOEL(1,1) + CH2FV * EHCOEL(2,1)
C  L'enthalpie de formation du coke peut-être considérée nulle
C  (pas loin du graphite)
C
C  L'enthalpie de changement de phase est donc celle de la réaction
C  Fuel_Liquide => FKC*Coke + (1-FKC)*Fuel_Vapeur
        HRFVAP =  (1.D0-FKC)*H02FOV-H02FOL
        WRITE(NFECRA,*) 'Donnees thermo pour le fuel'
        WRITE(NFECRA,*) 'PCIFOL ',PCIFOL
        WRITE(NFECRA,*) 'H02FOL ',H02FOL
        WRITE(NFECRA,*) 'CP2FOL ',CP2FOL
        WRITE(NFECRA,*) 'HRFVAP ',HRFVAP
        WRITE(NFECRA,*) 'H02FOV ',H02FOV
C
C  L'enthalpie de la réaction hétérogène est directement celle de la
C  formation d'une  mole de CO à partir de carbone à l'état de référence
C  il est d'usage d'ajouter cette enthalpie à celle de la phase
C  dispersée
C
C ------ Calcul de EHSOLI pour le fuel
C        Si CP2FOL > 0 : HFOL = H02FOL + CP2FOL(T2-TREFTH)
C
          DO IT = 1, NPO
            EHSOLI(IFOL,IT) = H02FOL
     &                            + CP2FOL * ( TH(IT) - TREFTH )
          ENDDO
C
C ---- Calcul relatif au coke
C
C ------ Coke = CH(GAMMA)O(DELTA)
C
C        On considere le PCI constant qqs T
C
C          Soit le PCI est connu et fourni dans le fichier
C          soit on considère qu'il est entierement fourni
C          par la combustion de la fraction carbone
C          supposée à l'état de référence
       DO IT = 1, NPO
          EHSOLI(IKF,IT) = CP2FOL * ( TH(IT) - TREFTH )
       ENDDO
C
      WRITE(NFECRA,*) ' Verification des enthalpies de formation'
      WRITE(NFECRA,*) ' CH4  ',EHCOEL(1,1)
      WRITE(NFECRA,*) ' C2H4 ',EHCOEL(2,1)
      WRITE(NFECRA,*) ' FOV  ',EHGAZE(IFOV,1)
      WRITE(NFECRA,*) ' FOL  ',EHSOLI(IFOL,1)
      WRITE(NFECRA,*) ' KF   ',EHSOLI(IKF,1)
C
C     Masse Vol + Diametre (en milimetres)
C        on suppose que les masse vol sont les memes
C        pour le fuel, coke et residu
C
      RHOKF  = RHO0FL
      DINIKF = DINIFL*(FKC*RHO0FL/RHOKF)**(1.D0/3.D0)
      DINIIN = DINIFL*(XINFOL*RHO0FL/RHO0FL)**(1.D0/3.D0)
      WRITE(NFECRA,*) ' D = ',DINIFL,DINIKF,DINIIN
C
      RETURN
C
C
C====================
C 3. SORTIE EN ERREUR
C====================
C
  99  CONTINUE
      WRITE ( NFECRA,9998 )
      CALL CSEXIT (1)
C     ===========
C
  999 CONTINUE
      WRITE ( NFECRA,9999 )
      CALL CSEXIT (1)
C     ===========
C--------
C FORMATS
C--------
C
C
 9990 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (FULECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (FUEL)                          ',/,
     &'@                                                            ',/,
     &'@  Erreur sur la composition du Coke :                       ',/,
     &'@   la somme des compositions elementaires doit etre egal    ',/, 
     &'@   a 1, elle vaut ici : ',G15.7,'                           ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9991 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (FULECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (FUEL)                          ',/,
     &'@                                                            ',/,
     &'@  Le nombre d''especes courantes doit etre inferieur        ',/,
     &'@                                  ou egal a',I10             ,/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9992 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (FULECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (FUEL)                          ',/,
     &'@                                                            ',/,
     &'@  Le nombre de points de tabulation est limite a ',I10       ,/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9993 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (FULECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (FUEL)                          ',/,
     &'@                                                            ',/,
     &'@  Le nombre d''especes elementaires est limite a ',I10       ,/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9994 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (FULECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (FUEL)                          ',/,
     &'@                                                            ',/,
     &'@  IRAYPP doit etre un entier egal a 0, 1 2 3 ou 4.          ',/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9998 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (FULECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (FUEL)                          ',/,
     &'@                                                            ',/,
     &'@  Erreur a l''ouverture du fichier parametrique.            ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (FULECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (FUEL)                          ',/,
     &'@                                                            ',/,
     &'@  Erreur a la lecture du fichier parametrique.              ',/,
     &'@    Le fichier a ete ouvert mais est peut etre incomplet    ',/,
     &'@    ou son format inadapte.                                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      END
C
C
c@z
