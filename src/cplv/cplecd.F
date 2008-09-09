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
                        SUBROUTINE CPLECD
C                       *****************
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  LECTURE DU FICHIER DE DONNEES PHYSIQUE PARTICULIERE
CFONC       RELATIF A LA COMBUSTION CHARBON PULVERISE
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
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
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
C VARIABLES LOCALES
C
      CHARACTER *150 CHAIN1,CHAIN2
      CHARACTER *12   NOMCOE(NGAZEM)
C
      INTEGER          IT     , ICE    , IAT    , IOS , II , IOXY
      INTEGER          NCOEL  , INICOE , INICHA , IERROR
      INTEGER          IDECAL , ICLA   , ICLAPC , ICHA   , IS
      INTEGER          IDEBCH , IFINCH , LONCH  , ICHAI  , ICHCOE
      INTEGER          ATCOEL(NGAZEM,NATOM)
C
      DOUBLE PRECISION FCOR   , PCISEC , PCIPUR , XASHPC
      DOUBLE PRECISION HCO20  , HO20   , HH2O0
      DOUBLE PRECISION DEN1   , DEN2
      DOUBLE PRECISION TMIN   , TMAX
      DOUBLE PRECISION IPCI(NCHARM)
      DOUBLE PRECISION WMOLCE(NGAZEM), EHCOEL(NGAZEM,NPOT)
      DOUBLE PRECISION CPCOEL(NGAZEM,NPOT),DET,MATDET
      DOUBLE PRECISION WMV1,WMVCH1,WMV2,WMVCH2
      DOUBLE PRECISION A11,A12,A13,A21,A22,A23,A31,A32,A33
      DOUBLE PRECISION DHVOL1,DHVOL2,HCOKE,HCHAR,EHVOL1,EHVOL2
      DOUBLE PRECISION WMCO,WMO2,WMCO2,WMH2O,WMN2,WMC
      DOUBLE PRECISION DMF4,DMF6,DMF7
C
C
C***********************************************************************
C
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
C ---- Nb de constituants elementaires (gazeux et solide)
C
      READ ( IMPFPP,*,ERR=999,END=999 ) NCOEL
      IF ( NCOEL.GT.NGAZEM ) THEN
        WRITE(NFECRA,9991) NGAZEM,NCOEL
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
      READ (IMPFPP,*,ERR=999,END=999 )
      READ (IMPFPP,1010,ERR=999,END=999 ) CHAIN1
      CALL VERLON (CHAIN1, IDEBCH, IFINCH, LONCH)
      CHAIN2(1:LONCH)=CHAIN1(IDEBCH:IFINCH)
C
      ICE=1
      ICHCOE=0
      DO ICHAI=1,LONCH
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
C --> Lecture caracteristiques charbon
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
C ---- Nb de charbons
C
      READ ( IMPFPP,*,ERR=999,END=999 ) NCHARB
      IF ( NCHARB.GT.NCHARM ) THEN
        WRITE(NFECRA,9995) NCHARM,NCHARB
        CALL CSEXIT (1)
      ENDIF
C
C ---- Nb de classes par charbon
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &                 ( NCLPCH(ICHA),ICHA=1,NCHARB )
      DO ICHA = 1, NCHARB
        IF ( NCLPCH(ICHA).GT.NCLCPM ) THEN
          WRITE(NFECRA,9996) NCLCPM,NCLPCH(ICHA),ICHA
          CALL CSEXIT (1)
        ENDIF
      ENDDO
C
C ---- Calcul du nb de classes et remplissage de ICHCOR
C
      NCLACP = 0
      DO ICHA = 1, NCHARB
        NCLACP = NCLACP + NCLPCH(ICHA)
      ENDDO
      IDECAL = 0
      DO ICHA = 1, NCHARB
        DO ICLAPC = 1, NCLPCH(ICHA)
          ICLA = ICLAPC+IDECAL
          ICHCOR(ICLA) = ICHA
        ENDDO
        IDECAL = NCLPCH(ICHA)
      ENDDO
C
C ---- Diametre initial par classe (m)
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( DIAM20(ICLA),ICLA=1,NCLACP )
C
C ---- Composition elementaire en C, H, O sur sec (% en masse)
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( CCH(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( HCH(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( OCH(ICHA),ICHA=1,NCHARB )
C
C ---- PCI sur charbon sec ou pur suivant la valeur de IPCI
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( IPCI(ICHA),PCICH(ICHA),ICHA=1,NCHARB )
C
C ---- CP moyen du charbon sec (J/kg/K)
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( CP2CH(ICHA),ICHA=1,NCHARB )
C
C ---- Masse volumique initiale (kg/m3)
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( RHO0CH(ICHA),ICHA=1,NCHARB )
C
C ---- Caracteristiques coke
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
C ------- Composition elementaire en C, H, O sur sec (%)
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( CCK(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( HCK(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( OCK(ICHA),ICHA=1,NCHARB )
C
C ------ PCI sur charbon sec
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( PCICK(ICHA),ICHA=1,NCHARB )
C
C ---- Caracteristiques cendres
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
C ------ Taux de cendre (kg/kg) en %
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( XASHCH(ICHA),ICHA=1,NCHARB )
C        Transformation en kg/kg
      DO ICHA = 1, NCHARB
        XASHCH(ICHA) = XASHCH(ICHA)/100.D0
      ENDDO
C
C ------ Enthalpie de formation des cendres (J/kg)
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( H0ASHC(ICHA),ICHA=1,NCHARB )
C
C ------ CP des cendres (J/kg/K)
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( CPASHC(ICHA),ICHA=1,NCHARB )
C
C ------ Taux d'humidite (kg/kg) en %
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( XWATCH(ICHA),ICHA=1,NCHARB )
C
C        Transformation en kg/kg
      DO ICHA = 1, NCHARB
        XWATCH(ICHA) = XWATCH(ICHA)/100.D0
      ENDDO
C
C      Transformation du taux de cendre de sec
C      sur humide en kg/kg
C
      DO ICHA = 1, NCHARB
        XASHCH(ICHA) = XASHCH(ICHA)*(1.D0-XWATCH(ICHA))
      ENDDO
C
C ---- Parametres de devolatilisation (modele de Kobayashi)
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( IY1CH(ICHA),Y1CH(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( IY2CH(ICHA),Y2CH(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( A1CH(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( A2CH(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( E1CH(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( E2CH(ICHA),ICHA=1,NCHARB )
C
C ---- Parametres combustion heterogene pour O2(modele a sphere retrecissante)
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( AHETCH(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( EHETCH(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( IOCHET(ICHA),ICHA=1,NCHARB)
C
C ---- Parametres combustion heterogene pour CO2(modele a sphere retrecissante)
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( AHETC2(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( EHETC2(ICHA),ICHA=1,NCHARB )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( IOETC2(ICHA),ICHA=1,NCHARB)
C
C --> Lecture caracteristiques Oxydants
C
      READ (IMPFPP,*,ERR=999,END=999 )
C
C ---- Nb d'oxydants
C
      READ ( IMPFPP,*,ERR=999,END=999 ) NOXYD
      IF ( NOXYD.LT.1 .OR. NOXYD .GT. 3 ) THEN
        WRITE(NFECRA,9895) NOXYD
        CALL CSEXIT (1)
      ENDIF
C
C ---- Composition en O2,N2,H2O,N2
C
      DO IOXY=1,3
        OXYO2 (IOXY) = 0.D0
        OXYN2 (IOXY) = 0.D0
        OXYH2O(IOXY) = 0.D0
        OXYCO2(IOXY) = 0.D0
      ENDDO
C
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( OXYO2(IOXY),IOXY=1,NOXYD )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( OXYN2(IOXY),IOXY=1,NOXYD )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( OXYH2O(IOXY),IOXY=1,NOXYD )
      READ (IMPFPP,*,ERR=999,END=999 )
     &     ( OXYCO2(IOXY),IOXY=1,NOXYD )
C
C --> Fermeture du fichier (ne pas oublier, car l'unite sert pour janaf)
C
      CLOSE(IMPFPP)
C
C
C======================================
C 2. CALCULS DE DONNEES COMPLEMENTAIRES
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
C     ATTENTION ON COMPTE EGALEMENT CH4 et le monomere CH2
C
      NGAZE = 2 + 5 + 2*NCHARB
C
C ---- Definition des pointeurs pour les tableaux WMOLE et EHGAZE
C      REMARQUE : Cette position de pointeurs va egalement servir
C                 pour le tableau de pointeurs IYM1 relatif aux
C                 tableaux PROPCE et PROPFB
C                 ON BALAYE JUSTE DE 1 A (NGAZE-2*NCHARB)
C
      IS    = 0
      IS    = IS + 1
      ICHX1 = IS
      IS    = IS + 1
      ICHX2 = IS
      IS   = IS + 1
      ICO  = IS
      IS   = IS + 1
      IO2  = IS
      IS   = IS + 1
      ICO2 = IS
      IS   = IS + 1
      IH2O = IS
      IS   = IS + 1
      IN2  = IS
      DO ICHA = 1, NCHARB
        ICHX1C(ICHA) = IS + ICHA
        ICHX2C(ICHA) = IS + ICHA + NCHARB
      ENDDO
C
C ---- Remplissage de EHGAZE et WMOLE
C         a partir de EHCOEL et WMOLCE
C
C ------ ATTENTION :
C        On prend par defaut
C          du CH4 pour CHX1m
C          et le monomere CH2 (jouant le role du C2H4) pour CHX2m
C
      DO IT = 1, NPO
        EHGAZE(ICHX1,IT) = EHCOEL(1,IT)
        EHGAZE(ICHX2,IT) = EHCOEL(2,IT)
        EHGAZE(ICO  ,IT) = EHCOEL(3,IT)
        EHGAZE(IO2  ,IT) = EHCOEL(4,IT)
        EHGAZE(ICO2 ,IT) = EHCOEL(5,IT)
        EHGAZE(IH2O ,IT) = EHCOEL(6,IT)
        EHGAZE(IN2  ,IT) = EHCOEL(7,IT)
      ENDDO
      WMOLE(ICHX1) = WMOLCE(1)
C ------ ATTENTION : Il faut prendre la masse molaire du monomere CH2
C                    et non celle du C2H4
      WMOLE(ICHX2) = WMOLAT(IATC)+WMOLAT(IATH)*2.D0
      WMOLE(ICO  ) = WMOLCE(3)
      WMOLE(IO2  ) = WMOLCE(4)
      WMOLE(ICO2 ) = WMOLCE(5)
      WMOLE(IH2O ) = WMOLCE(6)
      WMOLE(IN2  ) = WMOLCE(7)
C
C
C --> Calcul tabulation enthalpie - temperature pour le solide
C     Charbon reactif, Coke et Cendres
C
C ---- Nb de constituants solide
C
      NSOLID = 4*NCHARB
C
C ---- Definition des pointeurs ICH, ICK et IASH
C
      IS = 0
      DO ICHA = 1, NCHARB
        IS         = IS + 1
        ICH(ICHA)  = IS
        IS         = IS + 1
        ICK(ICHA)  = IS
        IS         = IS + 1
        IASH(ICHA) = IS
C
        IS         = IS + 1
        IWAT(ICHA) = IS
C
      ENDDO
C
C ---- Calcul relatif au charbon reactif
C
      ALPHAM = ZERO
      BETAM  = ZERO
C
      DO ICHA = 1, NCHARB
C
C ------ Composition du charbon reactif
C        sous la forme CH(ALPHA)O(BETA)
C        et calcul de la masse molaire
C
        ALPHA(ICHA) = (HCH(ICHA)/WMOLAT(IATH))
     &              / (CCH(ICHA)/WMOLAT(IATC))
        BETA(ICHA)  = (OCH(ICHA)/WMOLAT(IATO))
     &              / (CCH(ICHA)/WMOLAT(IATC))
        ALPHAM      = ALPHAM + ALPHA(ICHA)
        BETAM       = BETAM  + BETA(ICHA)
        WMOLS(ICH(ICHA)) = WMOLAT(IATC) + ALPHA(ICHA)*WMOLAT(IATH)
     &                   + BETA(ICHA)*WMOLAT(IATO)
      ENDDO
      ALPHAM = ALPHAM / DBLE(NCHARB)
      BETAM  = BETAM  / DBLE(NCHARB)
C
C ------ Transformation par la formule de Schaff du
C        PCI sur sec en PCI sur pur LORSQUE IPCI = 1
C
      FCOR = 1.2
      DO ICHA = 1, NCHARB
C
        IF (IPCI(ICHA).EQ.1) THEN
C    On donne directement le PCI sur sec J/kg ---> kcal/kg
          PCISEC = PCICH(ICHA)/1000.D0/XCAL2J
          XASHPC = XASHCH(ICHA)*100.D0
          PCIPUR = ( PCISEC+5.95D0*(ZERO+(1.D0-FCOR)*XASHPC) )*100.D0
     &           / ( 100.D0-ZERO-FCOR*XASHPC )
C    Conversion de PCI sur pur de kcal/kg ----> J/kg
          PCIPUR = PCIPUR*1000.D0*XCAL2J
C    On ecrase PCICH sur sec par PCIPUR
          PCICH(ICHA) = PCIPUR
        ENDIF
C
      ENDDO
C
C ------ Calcul de : H02 pour CH , CK , ASH et WT
C
      DO ICHA = 1, NCHARB
        HCO20 = EHGAZE(ICO2,1)*WMOLE(ICO2)
        HH2O0 = EHGAZE(IH2O,1)*WMOLE(IH2O)
        HO20  = EHGAZE(IO2 ,1)*WMOLE(IO2 )
C
        H02CH(ICHA) = PCICH(ICHA) +
     &     ( HCO20 + ALPHA(ICHA)/2.D0*HH2O0 -
     &       (1.D0+ALPHA(ICHA)/4.D0-BETA(ICHA)/2.D0)*HO20 )
     &     / WMOLS(ICH(ICHA))
C
        EH0SOL(ICH(ICHA)) = H02CH(ICHA)
C
      ENDDO
C
C   Pour l'instant H02CK et H02ASH sont nulles
C
      DO ICHA = 1, NCHARB
C
        EH0SOL(ICK (ICHA)) = 0.D0
        EH0SOL(IASH(ICHA)) = 0.D0
C
      ENDDO
C
C   Calcul de H02 pour l'eau liquide : d'apres JANAF
C                   H0 = -15871802.2
C                   CP = 4181.35
C
      DO ICHA = 1, NCHARB
        EH0SOL(IWAT(ICHA)) = -15871802.2D0
      ENDDO
C
C Construction table enthalpie/temperature Charbon
C   ---> Pour l'instant meme table que le Gaz
C
C  Nombre de tabulation : NPOC
C
      NPOC = NPO
C
      DO II = 1, NPOC
        THC(II) = TH(II)
      ENDDO
C
C
C ------ Calcul de EHSOLI pour le charbon reactif
C        Si CP2CH > 0 : HCH = H02CH + CP2CH(T2-TREFTH)
C        Sinon       : On considere le PCI constant qqs T2
C
      DO ICHA = 1, NCHARB
        IF ( CP2CH(ICHA).GT.EPSICP ) THEN
          DO IT = 1, NPOC
            EHSOLI(ICH(ICHA),IT) = EH0SOL(ICH(ICHA))
     &                            + CP2CH(ICHA)*(THC(IT)-TREFTH)
          ENDDO
        ELSE
          DO IT = 1, NPOC
            HCO20 = EHGAZE(ICO2,IT)*WMOLE(ICO2)
            HH2O0 = EHGAZE(IH2O,IT)*WMOLE(IH2O)
            HO20  = EHGAZE(IO2,IT) *WMOLE(IO2)
            EHSOLI(ICH(ICHA),IT) = PCICH(ICHA) +
     &       ( HCO20 + ALPHA(ICHA)/2.D0*HH2O0 -
     &         (1.D0+ALPHA(ICHA)/4.D0-BETA(ICHA)/2.D0)*HO20 )
     &       / WMOLS(ICH(ICHA))
          ENDDO
        ENDIF
      ENDDO
C
C ---- Calcul relatif au coke
C
C ------ Par defaut Coke = Carbone solide
C                   GAMMA = 0 et DELTA = 0
C            Si CP2CH > 0 : HCK = H02CH + CP2CH(T2-TREFTH)
C            Sinon       : HCK = Enthalpie du carbone pur
C
      GAMMA(ICHA) = ZERO
      DELTA(ICHA) = ZERO
      DO ICHA = 1, NCHARB
        WMOLS(ICK(ICHA)) = WMOLCE(8)
        IF (CP2CH(ICHA).GT.EPSICP) THEN
          DO IT = 1, NPOC
            EHSOLI(ICK(ICHA),IT) = EH0SOL(ICK(ICHA))
     &                           + CP2CH(ICHA)*(THC(IT)-TREFTH)
          ENDDO
        ELSE
          DO IT = 1, NPOC
            EHSOLI(ICK(ICHA),IT) = EHCOEL(8,IT)
          ENDDO
        ENDIF
      ENDDO
C
C ------ Coke = CH(GAMMA)O(DELTA)
C
      DO ICHA = 1, NCHARB
        IF ( PCICK(ICHA).GT.EPSICP ) THEN
          GAMMA(ICHA) = HCK(ICHA)/CCK(ICHA)
          DELTA(ICHA) = OCK(ICHA)/CCK(ICHA)
          WMOLS(ICK(ICHA)) = WMOLAT(IATC)+GAMMA(ICHA)*WMOLAT(IATH)
     &                     + DELTA(ICHA)*WMOLAT(IATO)
C        On considere le PCI constant qqs T
          DO IT = 1, NPOC
            HCO20 = EHGAZE(ICO2,IT)*WMOLE(ICO2)
            HH2O0 = EHGAZE(IH2O,IT)*WMOLE(IH2O)
            HO20  = EHGAZE(IO2,IT) *WMOLE(IO2)
            EHSOLI(ICK(ICHA),IT) = PCICK(ICHA) +
     &       ( HCO20 + GAMMA(ICHA)/2.D0*HH2O0 -
     &        (1.D0+GAMMA(ICHA)/4.D0-DELTA(ICHA)/2.D0)*HO20 )
     &       / WMOLS(ICH(ICHA))
          ENDDO
        ENDIF
      ENDDO
C
C ---- Calcul relatif aux cendres
C
      DO ICHA = 1, NCHARB
        IF (CP2CH(ICHA).GT.EPSICP) THEN
          DO IT = 1, NPOC
            EHSOLI(IASH(ICHA),IT) = EH0SOL(IASH(ICHA))
     &                            + CP2CH(ICHA)*(THC(IT)-TREFTH)
          ENDDO
        ELSE
          DO IT = 1, NPOC
            EHSOLI(IASH(ICHA),IT) = H0ASHC(ICHA)
     &                            + CPASHC(ICHA)*(THC(IT)-TREFTH)
          ENDDO
        ENDIF
        WMOLS(IASH(ICHA)) = ZERO
      ENDDO
C
C ---- Calcul relatif a l'eau
C
      DO ICHA = 1, NCHARB
        CP2WAT(ICHA) = 4181.35D0
C
        DO IT = 1, NPOC
          EHSOLI(IWAT(ICHA),IT) = EH0SOL(IWAT(ICHA))
     &                           +CP2WAT(ICHA)*(THC(IT)-TREFTH)
        ENDDO
      ENDDO
C
C --> Calcul relatifs aux matieres volatiles
C
C  On test que l'on utilise la meme option pour Y1 et Y2
C  pour chaque charbon
C
      DO ICHA = 1, NCHARB
        IF ( IY1CH(ICHA).NE. IY2CH(ICHA) ) THEN
C
          WRITE(NFECRA,9980) ICHA,IY1CH(ICHA),IY2CH(ICHA)
          CALL CSEXIT(1)
C
        ENDIF
        IF ( IY1CH(ICHA).LT. 0 .OR.  IY1CH(ICHA).LE.0 .OR.
     &       IY2CH(ICHA).LT. 0 .OR.  IY2CH(ICHA).LE.0     ) THEN
          WRITE(NFECRA,9981) ICHA,IY1CH(ICHA),IY1CH(ICHA)
          CALL CSEXIT(1)
        ENDIF
      ENDDO
C
C ---- Matieres volatiles legeres : [CH(CHX1); CO]
C       CH(ALPHA)O(BETA) --> A1 CH(CHX1) + B1 CO
C                            + (1-A1-B1) CH(GAMMA)O(DELTA)
C       Si IY1CH = 0 CHX1 fixe , Y1CH calcule
C       Si IY1CH = 1 Y1CH fixe , CHX1 calcule
C       Si IY1CH = 2 Y1CH fixe, CHX1 fixe, on ajoute de l'eau
C       CH(ALPHA)O(BETA) --> A1 CH(CHX1) + B1 CO + C1 H2O
C                            + (1-A1-B1) CH(GAMMA)O(DELTA)
C
      DO ICHA = 1, NCHARB
        IF ( IY1CH(ICHA).EQ.0 ) THEN
          CHX1(ICHA)   = 4.D0
          DEN1 = ( (1.D0-DELTA(ICHA))*CHX1(ICHA)-GAMMA(ICHA) )
     &        * WMOLS(ICH(ICHA))
          Y1CH(ICHA) = ( ( GAMMA(ICHA)*(BETA(ICHA)-1.D0) +
     &                     ALPHA(ICHA)*(1.D0-DELTA(ICHA)) )
     &                     *(WMOLAT(IATC)+CHX1(ICHA)*WMOLAT(IATH))
     &                 + ( BETA(ICHA)*(CHX1(ICHA)-GAMMA(ICHA)) +
     &                     DELTA(ICHA)*(ALPHA(ICHA)-CHX1(ICHA)) )
     &                     *WMOLE(ICO) )
     &               / DEN1
        ELSEIF(IY1CH(ICHA) .EQ. 1) THEN
          DEN1 = Y1CH(ICHA)*(1.D0-DELTA(ICHA))*WMOLS(ICH(ICHA))
     &         - WMOLAT(IATH)*(GAMMA(ICHA)*(BETA(ICHA)-1.D0) +
     &                         ALPHA(ICHA)*(1.D0-DELTA(ICHA)))
     &         + WMOLE(ICO)*(DELTA(ICHA)-BETA(ICHA))
          CHX1(ICHA) = ( WMOLAT(IATC)*(ALPHA(ICHA)-GAMMA(ICHA))
     &                 + WMOLAT(IATO)*(DELTA(ICHA)*ALPHA(ICHA)-
     &                                 BETA(ICHA) *GAMMA(ICHA)  )
     &                 + Y1CH(ICHA)*GAMMA(ICHA)*WMOLS(ICH(ICHA)) )
     &               / DEN1
        ELSEIF(IY1CH(ICHA) .EQ. 2) THEN
          CHX1(ICHA)   = 4.D0
        ENDIF
C
        IF ( IY1CH(ICHA).LT.2) THEN
         DEN1     =  (1.D0-DELTA(ICHA))*CHX1(ICHA)-GAMMA(ICHA)
         A1(ICHA) = ( GAMMA(ICHA)*(BETA(ICHA)-1.D0) +
     &                ALPHA(ICHA)*(1.D0-DELTA(ICHA)) )       / DEN1
         B1(ICHA) = ( BETA(ICHA)*(CHX1(ICHA)-GAMMA(ICHA)) +
     &                DELTA(ICHA)*(ALPHA(ICHA)-CHX1(ICHA)) ) / DEN1
         C1(ICHA) = ZERO
        ELSE
         WMV1 = Y1CH(ICHA)*( WMOLAT(IATC)
     &                      +ALPHA(ICHA)*WMOLAT(IATH)
     &                      +BETA(ICHA)*WMOLAT(IATO) )
         WMVCH1 = WMOLAT(IATC) + CHX1(ICHA)*WMOLAT(IATH)
C
         A11 = 4.D0-GAMMA(ICHA)
         A12 = -GAMMA(ICHA)
         A13 = 2.D0
         A21 = -DELTA(ICHA)
         A22 = 1.D0-DELTA(ICHA)
         A23 = 1.D0
         A31 = WMVCH1
         A32 = WMOLE(ICO)
         A33 = WMOLE(IH2O)
         DET = MATDET(A11, A12, A13, A21, A22, A23, A31, A32, A33)
C
         DEN1 = 1.D0 / DET
C
         A11 =  (ALPHA(ICHA)-GAMMA(ICHA))
         A12 = -GAMMA(ICHA)
         A13 =  2.D0
         A21 =  (BETA(ICHA)-DELTA(ICHA) )
         A22 =  (1.D0-DELTA(ICHA))
         A23 =  1.D0
         A31 =  WMV1
         A32 =  WMOLE(ICO)
         A33 =  WMOLE(IH2O)
         DET = MATDET(A11, A12, A13, A21, A22, A23, A31, A32, A33)
         A1(ICHA) = DEN1 * DET
C
         A11 = 4.D0-GAMMA(ICHA)
         A12 = ALPHA(ICHA)-GAMMA(ICHA)
         A13 = 2.D0
         A21 = -DELTA(ICHA)
         A22 = BETA(ICHA)-DELTA(ICHA)
         A23 = 1.D0
         A31 = WMVCH1
         A32 = WMV1
         A33 = WMOLE(IH2O)
         DET = MATDET(A11, A12, A13, A21, A22, A23, A31, A32, A33)
         B1(ICHA) = DEN1 * DET
C
         A11 = 4.D0-GAMMA(ICHA)
         A12 = -GAMMA(ICHA)
         A13 = ALPHA(ICHA)-GAMMA(ICHA)
         A21 = -DELTA(ICHA)
         A22 = 1.D0-DELTA(ICHA)
         A23 = BETA(ICHA)-DELTA(ICHA)
         A31 = WMVCH1
         A32 = WMOLE(ICO)
         A33 = WMV1
         DET = MATDET(A11, A12, A13, A21, A22, A23, A31, A32, A33)
         C1(ICHA) = DEN1 * DET
        ENDIF
       ENDDO
C
C ---- Matieres volatiles lourdes : [CH(CHX2); CO]
C       CH(ALPHA)O(BETA) --> A2 CH(CHX2) + B2 CO
C                            + (1-A2-B2) CH(GAMMA)O(DELTA)
C       Si IY2CH = 0 CHX2 fixe, Y2CH calcule
C       Si IY2CH = 1 Y2CH fixe, CHX2 calcule
C       Si IY2CH = 2 Y2CH fixe, CHX2 fixe, on ajoute de l'eau
C       CH(ALPHA)O(BETA) --> A2 CH(CHX2) + B2 CO + C2 H2O
C                            + (1-A2-B2) CH(GAMMA)O(DELTA)
C
      DO ICHA = 1, NCHARB
        IF ( IY2CH(ICHA).EQ.0 ) THEN
          CHX2(ICHA)   = 2.D0
          DEN2 = ( (1.D0-DELTA(ICHA))*CHX2(ICHA)-GAMMA(ICHA) )
     &        * WMOLS(ICH(ICHA))
          Y2CH(ICHA) = ( ( GAMMA(ICHA)*(BETA(ICHA)-1.D0) +
     &                     ALPHA(ICHA)*(1.D0-DELTA(ICHA)) )
     &                     *(WMOLAT(IATC)+CHX2(ICHA)*WMOLAT(IATH))
     &                 + ( BETA(ICHA)*(CHX2(ICHA)-GAMMA(ICHA)) +
     &                     DELTA(ICHA)*(ALPHA(ICHA)-CHX2(ICHA)) )
     &                     *WMOLE(ICO) )
     &               / DEN2
        ELSEIF(IY2CH(ICHA).EQ.1) THEN
          DEN2 = Y2CH(ICHA)*(1.D0-DELTA(ICHA))*WMOLS(ICH(ICHA))
     &         - WMOLAT(IATH)*(GAMMA(ICHA)*(BETA(ICHA)-1.D0) +
     &                         ALPHA(ICHA)*(1.D0-DELTA(ICHA)))
     &         + WMOLE(ICO)*(DELTA(ICHA)-BETA(ICHA))
          CHX2(ICHA) = ( WMOLAT(IATC)*(ALPHA(ICHA)-GAMMA(ICHA))
     &                 + WMOLAT(IATO)*(DELTA(ICHA)*ALPHA(ICHA)-
     &                                 BETA(ICHA) *GAMMA(ICHA)  )
     &                 + Y2CH(ICHA)*GAMMA(ICHA)*WMOLS(ICH(ICHA)) )
     &               / DEN2
        ELSE
          Y2CH(ICHA) = MIN(2.D0*Y1CH(ICHA),(1.D0+Y1CH(ICHA))/2.D0)
          CHX2(ICHA) = 2.D0
        ENDIF
C
        IF( IY2CH(ICHA).LT.2) THEN
         DEN2     = (1.D0-DELTA(ICHA))*CHX2(ICHA)-GAMMA(ICHA)
         A2(ICHA) = ( GAMMA(ICHA)*(BETA(ICHA)-1.D0) +
     &                ALPHA(ICHA)*(1.D0-DELTA(ICHA)) )       / DEN2
         B2(ICHA) = ( BETA(ICHA)*(CHX2(ICHA)-GAMMA(ICHA)) +
     &                DELTA(ICHA)*(ALPHA(ICHA)-CHX2(ICHA)) ) / DEN2
         C2(ICHA) = ZERO
        ELSE
         WMV2 = Y2CH(ICHA)*( WMOLAT(IATC)
     &                      +ALPHA(ICHA)*WMOLAT(IATH)
     &                      +BETA(ICHA)*WMOLAT(IATO) )
         WMVCH2 = WMOLAT(IATC) + CHX2(ICHA)*WMOLAT(IATH)
C
         A11 = 2.D0-GAMMA(ICHA)
         A12 = - GAMMA(ICHA)
         A13 = 2.D0
         A21 = -DELTA(ICHA)
         A22 = 1.D0-DELTA(ICHA)
         A23 = 1.D0
         A31 = WMVCH2
         A32 = WMOLE(ICO)
         A33 = WMOLE(IH2O)
         DET = MATDET(A11, A12, A13, A21, A22, A23, A31, A32, A33)
         DEN2 = 1.D0 / DET
C
         A11 = (ALPHA(ICHA)-GAMMA(ICHA))
         A12 = - GAMMA(ICHA)
         A13 = 2.D0
         A21 = (BETA(ICHA)-DELTA(ICHA) )
         A22 = (1.D0-DELTA(ICHA))
         A23 = 1.D0
         A31 = WMV2
         A32 = WMOLE(ICO)
         A33 = WMOLE(IH2O)
         DET = MATDET(A11, A12, A13, A21, A22, A23, A31, A32, A33)
         A2(ICHA) = DEN2 * DET
C
         A11 = (2.D0-GAMMA(ICHA))
         A12 = (ALPHA(ICHA)-GAMMA(ICHA))
         A13 = 2.D0
         A21 = -DELTA(ICHA)
         A22 = (BETA(ICHA)-DELTA(ICHA))
         A23 = 1.D0
         A31 = WMVCH2
         A32 = WMV2
         A33 = WMOLE(IH2O)
         DET = MATDET(A11, A12, A13, A21, A22, A23, A31, A32, A33)
         B2(ICHA) = DEN2 * DET
C
         A11 = 2.D0-GAMMA(ICHA)
         A12 = - GAMMA(ICHA)
         A13 = ALPHA(ICHA)-GAMMA(ICHA)
         A21 = -DELTA(ICHA)
         A22 = 1.D0-DELTA(ICHA)
         A23 = (BETA(ICHA)-DELTA(ICHA))
         A31 = WMVCH2
         A32 = WMOLE(ICO)
         A33 = WMV2
         DET = MATDET(A11, A12, A13, A21, A22, A23, A31, A32, A33)
         C2(ICHA) = DEN2 * DET
        ENDIF
      ENDDO
C
C --> Calcul de EHGAZE et de WMOLE
C     pour les especes CH(CHX1) et CH(CHX2)
C
C ---- Especes CH(CHX1)
C
      DO ICHA = 1, NCHARB
        WMOLE(ICHX1C(ICHA)) = WMOLAT(IATC)+CHX1(ICHA)*WMOLAT(IATH)
        IF ( IY1CH(ICHA).EQ.0 .OR. IY1CH(ICHA).EQ.2 ) THEN
          DO IT = 1, NPO
            EHGAZE(ICHX1C(ICHA),IT) = EHGAZE(ICHX1,IT)
          ENDDO
        ELSE
C        On a suppose D(HDEV,1,ICHA) = 0
          DO IT = 1, NPO
            DEN1 = A1(ICHA)*WMOLE(ICHX1C(ICHA))
            EHGAZE(ICHX1C(ICHA),IT) =
     &      ( ( EHSOLI(ICH(ICHA),IT) -
     &          (1.D0-Y1CH(ICHA))*EHSOLI(ICK(ICHA),IT) )
     &          * WMOLS(ICH(ICHA))
     &        - B1(ICHA)*WMOLE(ICO)*EHGAZE(ICO,IT) )
     &      / DEN1
          ENDDO
        ENDIF
      ENDDO
C
C ---- Especes CH(CHX2)
C
      DO ICHA = 1, NCHARB
        WMOLE(ICHX2C(ICHA)) = WMOLAT(IATC)+CHX2(ICHA)*WMOLAT(IATH)
        IF ( IY2CH(ICHA).EQ.0 .OR. IY2CH(ICHA).EQ.2 ) THEN
          DO IT = 1, NPO
            EHGAZE(ICHX2C(ICHA),IT) = EHGAZE(ICHX2,IT)
          ENDDO
        ELSE
C        On a suppose D(HDEV,2,ICHA) = 0
          DO IT = 1, NPO
            DEN2 = A2(ICHA)*WMOLE(ICHX2C(ICHA))
            EHGAZE(ICHX2C(ICHA),IT) =
     &      ( ( EHSOLI(ICH(ICHA),IT) -
     &           (1.D0-Y2CH(ICHA))*EHSOLI(ICK(ICHA),IT) )
     &           * WMOLS(ICH(ICHA))
     &        - B2(ICHA)*WMOLE(ICO)*EHGAZE(ICO,IT) )
     &      / DEN2
          ENDDO
        ENDIF
      ENDDO
C
C
C --> Calcul pour les differentes classes
C
      DO ICLA = 1, NCLACP
C
C ---- Diametre min (m)
C
        DIA2MN(ICLA) = ZERO
C
C ---- Masse volumique (kg/m3)
C
        RHO20(ICLA)  = RHO0CH(ICHCOR(ICLA))
        RHO2MN(ICLA) = RHO20(ICLA)*XASHCH(ICHCOR(ICLA))
C
C ---- Masse initiale de la particule (m)
C
       XMP0(ICLA) = RHO20(ICLA)*PI*(DIAM20(ICLA)**3)/6.D0
C
C ---- Masse de cendres de la particule (m)
C
       XMASH(ICLA) = XMP0(ICLA)*XASHCH(ICHCOR(ICLA))
C
C
      ENDDO
C
C AFFICHAGE RECAPITULATIF
C =======================
C
      WRITE(NFECRA,8000)
C
C  Constante du modele de devolatilisation
C
      WRITE(NFECRA,8100)
      DO ICHA=1,NCHARB
        WRITE(NFECRA,8101) ICHA
        IF ( IY1CH(ICHA) .EQ.0 ) THEN
          WRITE(NFECRA,8102)
        ELSE IF ( IY1CH(ICHA) .EQ.1 ) THEN
          WRITE(NFECRA,8103)
        ELSE IF  ( IY1CH(ICHA) .EQ.2 ) THEN
          WRITE(NFECRA,8104)
        ENDIF
C
        WRITE(NFECRA,8105) Y1CH(ICHA),Y2CH(ICHA)
      ENDDO
C
C
C  Composition matieres volatiles
C
      WRITE(NFECRA,8010)
C
      WRITE(NFECRA,8051)
      IERROR = 0
      DO ICHA=1,NCHARB
        WRITE(NFECRA,8052) ICHA,CHX1(ICHA),CHX2(ICHA)
        IF ( CHX1(ICHA) .LE. 0.D0 .OR. CHX2(ICHA) .LE. 0 ) THEN
          IERROR = 1
        ENDIF
      ENDDO
      IF ( IERROR .EQ. 1 ) THEN
        WRITE(NFECRA,9970)
        CALL CSEXIT(1)
      ENDIF
C
      DO ICHA=1,NCHARB
        WRITE(NFECRA,8011)  ICHA
        WRITE(NFECRA,8012)
        WRITE(NFECRA,8013) A1(ICHA),B1(ICHA),C1(ICHA)
        WRITE(NFECRA,8014) A2(ICHA),B2(ICHA),C2(ICHA)
      ENDDO
C
C  Calcul de DHdev
C
      write(NFECRA,8001)
      DO ICHA=1,NCHARB
C
C    Enthalie du Charbon
C
       HCHAR = EHSOLI(ICH(ICHA),1)
       HCOKE = EHSOLI(ICK(ICHA),1)
C
C    Enthalpie des matieres volatiles
C
       EHVOL1 = ( A1(ICHA)*EHGAZE(ICHX1C(ICHA),1)*WMOLE(ICHX1C(ICHA))
     &           +B1(ICHA)*EHGAZE(ICO ,1)*WMOLE(ICO)
     &           +C1(ICHA)*EHGAZE(IH2O,1)*WMOLE(IH2O) )
     &         /( A1(ICHA)*WMOLE(ICHX1C(ICHA))
     &           +B1(ICHA)*WMOLE(ICO)
     &           +C1(ICHA)*WMOLE(IH2O) )
C
       EHVOL2 = ( A2(ICHA)*EHGAZE(ICHX2C(ICHA),1)*WMOLE(ICHX2C(ICHA))
     &           +B2(ICHA)*EHGAZE(ICO ,1)*WMOLE(ICO)
     &           +C2(ICHA)*EHGAZE(IH2O,1)*WMOLE(IH2O) )
     &         /( A2(ICHA)*WMOLE(ICHX2C(ICHA))
     &           +B2(ICHA)*WMOLE(ICO)
     &           +C2(ICHA)*WMOLE(IH2O) )
C
       DHVOL1 = HCHAR-Y1CH(ICHA)*EHVOL1-(1.D0-Y1CH(ICHA))*HCOKE
       DHVOL2 = HCHAR-Y2CH(ICHA)*EHVOL2-(1.D0-Y2CH(ICHA))*HCOKE
C
       WRITE(NFECRA,8002) ICHA,DHVOL1,DHVOL2,PCICH(ICHA)
C
      ENDDO
C
C  Loi enthalpie/temperature
C
      WRITE(NFECRA,8020)
      DO ICHA = 1, NCHARB
        WRITE(NFECRA,8011) ICHA
        WRITE(NFECRA,8021)
        DO II = 1, NPOC
          write(NFECRA,8022)THC(II),EHSOLI(ICH (ICHA),II)
     &                             ,EHSOLI(ICK (ICHA),II)
     &                             ,EHSOLI(IASH(ICHA),II)
     &                             ,EHSOLI(IWAT(ICHA),II)
        ENDDO
      ENDDO
C
C  Calcul des AiFj : nbre de mole de i par kg de j a l'origine
C
      WMCO   = WMOLE(ICO)
      WMO2   = WMOLE(IO2)
      WMCO2  = WMOLE(ICO2)
      WMH2O  = WMOLE(IH2O)
      WMN2   = WMOLE(IN2)
      WMC    = WMOLAT(IATC)
C
      ACOF3 =  1.D0/WMC
      AO2F3 = -.5D0/WMC
C
      DMF4  = ( OXYO2 (1)*WMO2 +OXYN2 (1)*WMN2
     &         +OXYH2O(1)*WMH2O+OXYCO2(1)*WMCO2 )
C
      IF ( DMF4 .LE. 0.D0 ) THEN
        WRITE(NFECRA,9896) OXYO2(1) ,OXYN2(1) ,
     &                     OXYH2O(1),OXYCO2(1)
        CALL CSEXIT(1)
      ENDIF
C
      AO2F4  = OXYO2(1)  / DMF4
      AN2F4  = OXYN2(1)  / DMF4
      AH2OF4 = OXYH2O(1) / DMF4
      ACO2F4 = OXYCO2(1) / DMF4
C
      IF ( NOXYD .GE. 2.D0 ) THEN
        DMF6  = ( OXYO2 (2)*WMO2 +OXYN2 (2)*WMN2
     &           +OXYH2O(2)*WMH2O+OXYCO2(2)*WMCO2 )
        IF ( DMF6 .LE. 0.D0 ) THEN
          WRITE(NFECRA,9897) OXYO2(2) ,OXYN2(2) ,
     &                       OXYH2O(2),OXYCO2(2)
          CALL CSEXIT(1)
        ENDIF
C
        AO2F6  = OXYO2(2)  / DMF6
        AN2F6  = OXYN2(2)  / DMF6
        AH2OF6 = OXYH2O(2) / DMF6
        ACO2F6 = OXYCO2(2) / DMF6
C
      ELSE
        AO2F6  = 0.D0
        AN2F6  = 0.D0
        AH2OF6 = 0.D0
        ACO2F6 = 0.D0
      ENDIF
C
      IF ( NOXYD .EQ. 3.D0 ) THEN
        DMF7  = ( OXYO2 (3)*WMO2 +OXYN2 (3)*WMN2
     &           +OXYH2O(3)*WMH2O+OXYCO2(3)*WMCO2 )
        IF ( DMF7 .LE. 0.D0 ) THEN
          WRITE(NFECRA,9898) OXYO2(3) ,OXYN2(3) ,
     &                       OXYH2O(3),OXYCO2(3)
          CALL CSEXIT(1)
        ENDIF
C
        AO2F7  = OXYO2(3)  / DMF7
        AN2F7  = OXYN2(3)  / DMF7
        AH2OF7 = OXYH2O(3) / DMF7
        ACO2F7 = OXYCO2(3) / DMF7
C
      ELSE
        AO2F7  = 0.D0
        AN2F7  = 0.D0
        AH2OF7 = 0.D0
        ACO2F7 = 0.D0
      ENDIF
C
      AH2OF5  = 1.D0/WMH2O
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
C
C
C--------
C FORMATS
C--------
C
C
 8000 FORMAT(1X,' RECAPITULATIF SUR LES CHARBONS :'/,
     &       1X,' ==============================  '  )
 8001 FORMAT(/,3X,'Delta Hdev  ',/,
     &       5X,' Numero Charbon',6X,'Legeres',7X,'Lourdes',9X,'PCI')
 8002 FORMAT(10X,I3,8X,G15.7,1X,G15.7,G15.7)
C
 8010 FORMAT(/,3X,'COMPOSITION MATIERES VOLATILES ')
 8011 FORMAT(/,5X,' Numero Charbon',I6)
 8012 FORMAT(18X,'CHX1 ',12X,'CO',12X,'H2O')
 8013 FORMAT( 8X,' MV1 ',G15.7,G15.7,G15.7)
 8014 FORMAT( 8X,' MV2 ',G15.7,G15.7,G15.7)
 8020 FORMAT(/,3X,'LOI ENTHALPIE/TEMERATURE ')
 8021 FORMAT( 10X,'TEMPERATURE',9X,'Hch ',10X,'Hcoke',10X,
     &                             'Hash',10X,'Hwat')
 8022 FORMAT( 8X,5(G15.7))
C
 8050 FORMAT(/,3X,'COMPOSITION DU CHARBON ')
 8051 FORMAT(/,5X,'Numero Charbon',7X,'CHX1',12X,'CHX2')
 8052 FORMAT(10X,I3,8X,G15.7,1X,G15.7)
C
 8100 FORMAT(/,3X,'CONSTANTE DU MODELE DE DEVOLATILISATION')
 8101 FORMAT(/,8X,'CHARBON ',I2)
 8102 FORMAT(/,12X,'OPTION 0 : Y1 et Y2 CALCULES     ')
 8103 FORMAT(/,12X,'OPTION 1 : Y1 et Y2 FIXES        ')
 8104 FORMAT(/,12X,'OPTION 2 : CHX1 et CHX2 FIXES    ')

 8105 FORMAT(12X,'Y1CH = ',G15.7,3X,'Y2CH = ',G15.7)
C
 9970 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  CERTAINES VALEURS DE CHX1 ET CHX2 SONT NEGATIVES          ',/,
     &'@  OU NULLES                                                 ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9980 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  Vous utilisez une option differente pour Y1 et Y2         ',/,
     &'@  pour le charbon ',I2,'                                    ',/,
     &'@         IY1CH = ',I2,'                                     ',/,
     &'@         IY2CH = ',I2,'                                     ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9981 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  Vous utilisez une option  pour Y1 et Y2 non disponible    ',/,
     &'@  pour le charbon ',I2,'                                    ',/,
     &'@         IY1CH = ',I2,'                                     ',/,
     &'@         IY2CH = ',I2,'                                     ',/,
     &'@                                                            ',/,
     &'@  seules les valeur : 0 , 1 ou 2 sont admissibles.          ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9991 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
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
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
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
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
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
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
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
 9995 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  Le nombre de charbons est limite a      ',I10              ,/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9996 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  Le nombre de classes par charbon est limite a ',I10        ,/,
     &'@   Il vaut ',I10   ,' pour le charbon ',I10                  ,/,
     &'@                      dans le fichier parametrique          ',/,
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
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
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
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
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
 9895 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  Le nombre d''Oxydants doit etre compris entre 1 et 3      ',/,
     &'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 9896 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  LA COMPOSITION DE L''OXYDANT 1 EST ERRONEE                ',/,
     &'@     O2  :  ',G15.7,'                                       ',/,
     &'@     N2  :  ',G15.7,'                                       ',/,
     &'@     H2O :  ',G15.7,'                                       ',/,
     &'@     CO2 :  ',G15.7,'                                       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9897 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  LA COMPOSITION DE L''OXYDANT 2 EST ERRONEE                ',/,
     &'@     O2  :  ',G15.7,'                                       ',/,
     &'@     N2  :  ',G15.7,'                                       ',/,
     &'@     H2O :  ',G15.7,'                                       ',/,
     &'@     CO2 :  ',G15.7,'                                       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9898 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  LA COMPOSITION DE L''OXYDANT 3 EST ERRONEE                ',/,
     &'@     O2  :  ',G15.7,'                                       ',/,
     &'@     N2  :  ',G15.7,'                                       ',/,
     &'@     H2O :  ',G15.7,'                                       ',/,
     &'@     CO2 :  ',G15.7,'                                       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9899 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (CPLECD)      ',/,
     &'@    *********                                               ',/,
     &'@      PHYSIQUE PARTICULIERE (CHARBON PULVERISE)             ',/,
     &'@                                                            ',/,
     &'@  LA COMPOSITION DE L''OXYDANT 1 EST ERRONEE                ',/,
     &'@     O2  :  ',G15.7,'                                       ',/,
     &'@     N2  :  ',G15.7,'                                       ',/,
     &'@     H2O :  ',G15.7,'                                       ',/,
     &'@     CO2 :  ',G15.7,'                                       ',/,
     &'@                                                            ',/,
     &'@  ALORS QUE L''OXYDANT 1 DOIT OBLIGATOIREMENT CONTENIR      ',/,
     &'@  DE L''OXYGENE                                             ',/,
     &'@                                                            ',/,
     &'@  Verifier le fichier parametrique.                         ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
      END
C
C
c@z
                            FUNCTION MATDET
C                           ***************
C     -------------------------------------------------------------
     & ( A11 , A12 , A13 , A21 , A22 , A23 , A31 , A32 , A33 )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  CALCUL DU DETERMINANT DE LA MATRICE 3x3
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! AA1 a A33    ! R  !  -> ! COEFFICIENT DE LA MATRICE            !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION A11, A12, A13, A21, A22, A23, A31, A32, A33
      DOUBLE PRECISION MATDET
C
C***********************************************************************
C
C=======================================================================
C 1. CALCUL DU DETERMINANT
C=======================================================================
C
      MATDET = A11*A22*A33 + A21*A32*A13 + A31*A12*A23
     &       - A11*A32*A23 - A21*A12*A33 - A31*A22*A13
C
      END
c@z
