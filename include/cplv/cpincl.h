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
C                             cpincl.h
C
C***********************************************************************
C
C            INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF
C                A LA COMBUSTION DU CHARBON PULVERISE
C
C Necessite ppppar.h et ppthch
C-----------------------------------------------------------------------
C
C--> POINTEURS VARIABLES COMBUSTION CHARBON PULVERISE cpincl, ppincl
C
C -------voir ppppar
CC       NCHARM        --> Nombre maximal de charbons
CC       NCPCMX        --> Nombre maximal de classes par charbon
CC       NCLCPM        --> Nombre total de classes
CC
C      INTEGER    NCHARM  , NCPCMX   , NCLCPM
C      PARAMETER (NCHARM=3, NCPCMX=10, NCLCPM=NCHARM*NCPCMX)
C
CC ------
C      EPSICP : Precision pour les tests
C
        DOUBLE PRECISION EPSICP
        PARAMETER ( EPSICP = 1.D-8 )
C
C
C--> DONNEES RELATIVES AU CHARBON
C
C       NCHARB        --> Nombre de charbons
      INTEGER          NCHARB
C
C ---- PAR CHARBON (grandeurs fournies)
C
C      - Distribution granulometrique
C        NCLPCH(CH)   --> Nombre de classes par charbon
C
      INTEGER          NCLPCH(NCHARM)
C
C      - Proprietes sur charbon sec
C        CCH(CH)      --> Composition elementaire en C, H, O sur sec (%)
C        HCH(CH)          du charbon
C        OCH(CH)
C        ALPHA(CH)    --> Composition du charbon reactif
C        BETA(CH)         sous la forme CH(ALPHA)O(BETA)
C                         ALPHA(CH) = HCH(CH)/CCH(CH)
C                         BETA(CH)  = OCH(CH)/CCH(CH)
C        ALPHAM       --> ALPHA moyen
C        BETAM        --> BETA moyen
C        PCICH(CH)    --> PCI (J/kg) charbon
C        RHO0CH(CH)   --> Masse volumique initiale (kg/m3)
C      - Proprietes sur charbon sec du coke
C        CCK(CH)      --> Composition elementaire en C, H, O sur sec (%)
C        HCK(CH)          du coke
C        OCK(CH)
C        GAMMA(CH)    --> Composition du coke
C        DELTA(CH)        sous la forme CH(GAMMA)O(DELTA)
C                         GAMMA(CH) = HCK(CH)/CCK(CH)
C                         DELTA(CH) = OCK(CH)/CCK(CH)
C        PCICK(CH)    --> PCI (J/kg) coke
C        RHOCK(CH)    --> Masse volumique coke
C      - Proprietes sur charbon sec des cendres (ou humide)
C         XASHCH(CH)   --> Taux de cendre (kg/kg)
C        CPASHC(CH)   --> CP des cendres (J/kg/K)
C        H0ASHC(CH)   --> Enthalpie de formation des cendres (J/kg)
C        H02CH        --> H0 du Charbon
C        CPCH         --> CP du Charbon
C        XWATCH(CH)   --> Taux d'humidite (kg/kg)
C
      DOUBLE PRECISION CCH   (NCHARM), HCH   (NCHARM), OCH   (NCHARM),
     &                 ALPHA (NCHARM), BETA  (NCHARM), ALPHAM, BETAM ,
     &                 PCICH (NCHARM), RHO0CH(NCHARM),
     &                 CCK   (NCHARM), HCK   (NCHARM), OCK   (NCHARM),
     &                 GAMMA (NCHARM), DELTA (NCHARM),
     &                 RHOCK (NCHARM), PCICK (NCHARM),
     &                 XASHCH(NCHARM), CPASHC(NCHARM),
     &                 H0ASHC(NCHARM),
     &                 H02CH (NCHARM), CP2CH (NCHARM),
     &                 XWATCH(NCHARM), CP2WAT(NCHARM)
C
C      - Parametres cinetiques pour la devolatilisation
C         (Modele de Kobayashi)
C        IY1CH(CH)    --> Indicateur : 0 si MVl = {CH4;CO}
C                                      1 si MVl = {CHz;CO}
C        Y1CH(CH)     --> Coefficient stoechiometrique (adim)
C                         calcule si IY1CH = 0 ; donne si IY1CH = 1
C        A1CH(CH)     --> Facteur pre-exponetielle (1/s)
C         E1CH(CH)     --> Energie d'activation (J/mol)
C        IY2CH(CH)    --> Indicateur : 0 si MVL = {C2H4;CO}
C                                      1 si MVL = {CxHy;CO}
C        Y2CH(CH)     --> Coefficient stoechiometrique (adim)
C                         calcule si IY2CH = 0 ; donne si IY2CH = 1
C        A2CH(CH)     --> Constante preexponetielle (1/s)
C        E2CH(CH)     --> Energie d'activation (J/mol)
C         - Parametres cinetiques pour la combustion heterogene du coke
C           (Modele a sphere retrecissante)
C        AHETCH(CH)   --> Constante pre-exponentielle (kg/m2/s/atm)
C        EHETCH(CH)   --> Energie d'activation (kcal/mol)
C        IOCHET(CH)   --> Ordre de la reaction 0.5 si = 0 1 si = 1
C
      INTEGER          IY1CH (NCHARM), IY2CH (NCHARM)
      INTEGER          IOCHET (NCHARM)
      DOUBLE PRECISION Y1CH  (NCHARM), A1CH  (NCHARM), E1CH  (NCHARM),
     &                 Y2CH  (NCHARM), A2CH  (NCHARM), E2CH  (NCHARM),
     &                 AHETCH(NCHARM), EHETCH(NCHARM)
C
C      - Enthalpie du charbon reactif, coke et cendres
C     ICH(CH)      --> Pointeur dans le tableau EHSOLI pour
C                         le Charbon Reactif
C     ICK(CH)      --> Pointeur dans le tableau EHSOLI pour
C                         le Coke
C     IASH(CH)     --> Pointeur dans le tableau EHSOLI pour
C                         les cendres
C     IWAT(CH)     --> Pointeur dans le tableau EHSOLI pour
C                         l'humidite
C     NSOLID       --> Nb constituants solides (Ch.Reactif, Coke, Ash)
C     NSOLIM       --> Nb maximal de constituants solides
C     EHSOLI(S,IT) --> Enthalpie massique (J/kg) du constituant solide
C                         no S a la temperature T(IT)
C     WMOLS(S)     --> Masse molaire du constituant solide
C     EH0SOL(S)    --- Enthalpie de formation (J/kg) du constituant solide
C                      no S
C
      INTEGER    NSOLIM
      PARAMETER( NSOLIM = 4*NCHARM )
C
      INTEGER          NSOLID, ICH(NCHARM), ICK(NCHARM), IASH(NCHARM)
      INTEGER          IWAT(NCHARM)
      DOUBLE PRECISION EHSOLI(NSOLIM,NPOT), WMOLS(NSOLIM)
      DOUBLE PRECISION EH0SOL(NSOLIM)
C
C ---- PAR CLASSES (grandeurs deduites)
C
C        NCLACP     --> Nb de classes
C
      INTEGER          NCLACP
C
C      - Proprietes
C        ICHCOR(CL)  --> = ICH si la classe consideree appartient
C                        au charbon ICH (1, 2, ...)
C        DIAM20(CL)  --> Diametre initial (m)
C        DIA2MN(CL)  --> Diametre minimum (m)
C        RHO20(CL)   --> Masse volumique initiale (kg/m3)
C        RHO2MN(CL)  --> Masse volumique minimale (kg/m3)
C        XMP0(CL)    --> Masse initiale de la particule (m)
C        XMASH(CL)   --> Masse de cendres de la particule (m)
C
      INTEGER          ICHCOR(NCLCPM)
      DOUBLE PRECISION DIAM20(NCLCPM), DIA2MN(NCLCPM),
     &                 RHO20 (NCLCPM), RHO2MN(NCLCPM),
     &                 XMP0  (NCLCPM), XMASH (NCLCPM)
C
C
C--> DONNEES RELATIVES A LA COMBUSTION DES ESPECES GAZEUSES
C
C        ICHX1C(CH)  --> Pointeur CHx1  pour EHGAZE et WMOLE
C        ICHX2C(CH)  --> Pointeur CHx2  pour EHGAZE et WMOLE
C        ICHX1       --> Pointeur CHx1m pour EHGAZE et WMOLE
C        ICHX2       --> Pointeur CHx2m pour EHGAZE et WMOLE
C        ICO         --> Pointeur CO    pour EHGAZE et WMOLE
C        IO2         --> Pointeur O2    pour EHGAZE et WMOLE
C        ICO2        --> Pointeur CO2   pour EHGAZE et WMOLE
C        IH2O        --> Pointeur H2O   pour EHGAZE et WMOLE
C        IN2         --> Pointeur N2    pour EHGAZE et WMOLE
C        XSI         --> XSI = 3,76 pour de l'air
C        F3MAX       --> Maximum pour le traceur F3
C        CHX1(CH)    --> Composition de l'hydrocarbure relatif
C                        au MVl : CH(X1)
C        CHX2(CH)    --> Composition de l'hydrocarbure relatif
C                        au MVL : CH(X2)
C        A1(CH),     --> Coefficients stoechiometriques molaires pour
C        B1(CH)          la reaction de devolatilisation a basses T
C        A2(CH),     --> Coefficients stoechiometriques molaires pour
C        B2(CH)          la reaction de devolatilisation a basses T
C
      INTEGER          ICHX1C(NCHARM), ICHX2C(NCHARM),
     &                 ICHX1, ICHX2, ICO, IO2, ICO2, IH2O, IN2
      DOUBLE PRECISION XSI, F3MAX,
     &                 CHX1(NCHARM), CHX2(NCHARM),
     &                 A1(NCHARM), B1(NCHARM),
     &                 A2(NCHARM), B2(NCHARM)
C
C--> DONNEES COMPLEMENTAIRES RELATIVES AU CALCUL DE RHO
C    SUR LES FACETTES DE BORD
C
C       IENTAT(IENT) --> Indicateur air par type de facette d'entree
C       IENTCP(IENT) --> Indicateur CP  par type de facette d'entree
C       TIMPAT(IENT) --> Temperature en K pour l'air relative
C                         a l'entree IENT
C       X20(IENT,    --> Fraction massique dans le melange de charbon
C           ICLA   )     de la classe ICLA relative a l'entree IENT
C
      INTEGER          IENTAT(NOZPPM), IENTCP(NOZPPM)
      DOUBLE PRECISION TIMPAT(NOZPPM), X20(NOZPPM,NCLCPM)
C
C--> POINTEURS DANS LE TABLEAU TBMCR
C
      INTEGER IF1MC(NCHARM) , IF2MC(NCHARM)
      INTEGER IX1MC ,IX2MC, ICHX1F1, ICHX2F2
      INTEGER ICOF1, ICOF2
C
C--> DEFINITION DES COMMONS
C
      COMMON / ICPCOM / NCHARB, NCLPCH, NCLACP, IY1CH , IY2CH ,
     &                  ICHCOR,
     &                  NSOLID, ICH   , ICK   , IASH  , IWAT  ,
     &                  ICHX1C, ICHX2C,
     &                  ICHX1 , ICHX2 ,
     &                  ICO   , IO2   , ICO2  , IH2O  , IN2   ,
     &                  IENTAT, IENTCP,
     &                  IF1MC , IF2MC ,
     &                  IX1MC , IX2MC , ICHX1F1, ICHX2F2,
     &                  ICOF1 , ICOF2 , IOCHET
C
      COMMON / RCPCOM / CCH   , HCH   , OCH   ,
     &                  ALPHA , BETA  , ALPHAM, BETAM ,
     &                  PCICH , RHO0CH,
     &                  CCK   , HCK   , OCK   ,
     &                  GAMMA , DELTA ,
     &                  RHOCK , PCICK ,
     &                  XASHCH, CPASHC, H0ASHC, XWATCH ,
     &                  H02CH , CP2CH ,
     &                  Y1CH  , A1CH  , E1CH  , Y2CH  , A2CH  , E2CH  ,
     &                  AHETCH, EHETCH,
     &                  EHSOLI, WMOLS ,
     &                  DIAM20, DIA2MN, RHO20 , RHO2MN,
     &                  XMP0  , XMASH ,
     &                  XSI   , F3MAX ,
     &                  CHX1  , CHX2  ,
     &                  A1    , B1    , A2    , B2    ,
     &                  TIMPAT, X20   ,
     &                  EH0SOL, CP2WAT
C
C--> GRANDEURS FOURNIES PAR L'UTILISATEUR EN CONDITIONS AUX LIMITES
C      PERMETTANT DE CALCULER AUTOMATIQUEMENT LA VITESSE, LA TURBULENCE,
C      L'ENTHALPIE D'ENTREE.
C
C    POUR LES ENTREES UNIQUEMENT , IENT ETANT LE NUMERO DE ZONE FRONTIERE
C
C       QIMPAT(IENT)           --> Debit       Air          en kg/s
C       TIMPAT(IENT)           --> Temperature Air          en K
C       QIMPCP(IENT,ICHA)      --> Debit       Charbon ICHA en kg/s
C       TIMPCP(IENT,ICHA)      --> Temperature Charbon ICHA en K
C       DISTCH(IENT,ICHA,ICLA) --> Distribution en %masse de la classe ICLA
C                                  pour le charbon ICHA
C
      DOUBLE PRECISION  QIMPAT(NOZPPM)
      DOUBLE PRECISION  QIMPCP(NOZPPM,NCHARM), TIMPCP(NOZPPM,NCHARM)
      DOUBLE PRECISION  DISTCH(NOZPPM,NCHARM,NCPCMX)
c     DOUBLE PRECISION  COEFE(NGAZEM), XSOLID(NSOLIM)
c     DOUBLE PRECISION  F1MC(NCHARM) , F2MC(NCHARM)
C
      COMMON / RCPCLI / QIMPAT       ,
     &                  QIMPCP               , TIMPCP               ,
     &                  DISTCH
c    &                  COEFE        , XSOLID        ,
c    &                  F1MC         , F2MC
C
C Complement Table
C
       DOUBLE PRECISION THC(NPOT)
       COMMON/TABLEC/   THC
C
       INTEGER          NPOC
       COMMON/ITABLC/   NPOC
C
C Equation sur YCO2
C
       INTEGER         IEQCO2 , IYCO2
       COMMON/EQUCO2 / IEQCO2 , IYCO2





C FIN
c@z
