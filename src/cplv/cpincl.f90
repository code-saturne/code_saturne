!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

!                             cpincl.h

!===============================================================================

!            INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF
!                A LA COMBUSTION DU CHARBON PULVERISE

! Necessite ppppar.h et ppthch
!-------------------------------------------------------------------------------

!--> POINTEURS VARIABLES COMBUSTION CHARBON PULVERISE cpincl, ppincl

! -------voir ppppar
!C       NCHARM        --> Nombre maximal de charbons
!C       NCPCMX        --> Nombre maximal de classes par charbon
!C       NCLCPM        --> Nombre total de classes
!C
!      INTEGER    NCHARM  , NCPCMX   , NCLCPM
!      PARAMETER (NCHARM=3, NCPCMX=10, NCLCPM=NCHARM*NCPCMX)

!C ------
!      EPSICP : Precision pour les tests

  double precision epsicp
  parameter ( epsicp = 1.d-8 )


!--> DONNEES RELATIVES AU CHARBON

!       NCHARB        --> Nombre de charbons
integer          ncharb

! ---- PAR CHARBON (grandeurs fournies)

!      - Distribution granulometrique
!        NCLPCH(CH)   --> Nombre de classes par charbon

integer          nclpch(ncharm)

!      - Proprietes sur charbon sec
!        CCH(CH)      --> Composition elementaire en C, H, O sur sec (%)
!        HCH(CH)          du charbon
!        OCH(CH)
!        ALPHA(CH)    --> Composition du charbon reactif
!        BETA(CH)         sous la forme CH(ALPHA)O(BETA)
!                         ALPHA(CH) = HCH(CH)/CCH(CH)
!                         BETA(CH)  = OCH(CH)/CCH(CH)
!        ALPHAM       --> ALPHA moyen
!        BETAM        --> BETA moyen
!        PCICH(CH)    --> PCI (J/kg) charbon
!        RHO0CH(CH)   --> Masse volumique initiale (kg/m3)
!      - Proprietes sur charbon sec du coke
!        CCK(CH)      --> Composition elementaire en C, H, O sur sec (%)
!        HCK(CH)          du coke
!        OCK(CH)
!        GAMMA(CH)    --> Composition du coke
!        DELTA(CH)        sous la forme CH(GAMMA)O(DELTA)
!                         GAMMA(CH) = HCK(CH)/CCK(CH)
!                         DELTA(CH) = OCK(CH)/CCK(CH)
!        PCICK(CH)    --> PCI (J/kg) coke
!        RHOCK(CH)    --> Masse volumique coke
!      - Proprietes sur charbon sec des cendres (ou humide)
!         XASHCH(CH)   --> Taux de cendre (kg/kg)
!        CPASHC(CH)   --> CP des cendres (J/kg/K)
!        H0ASHC(CH)   --> Enthalpie de formation des cendres (J/kg)
!        H02CH        --> H0 du Charbon
!        CPCH         --> CP du Charbon
!        XWATCH(CH)   --> Taux d'humidite (kg/kg)

double precision cch   (ncharm), hch   (ncharm), och   (ncharm),  &
                 alpha (ncharm), beta  (ncharm), alpham, betam ,  &
                 pcich (ncharm), rho0ch(ncharm),                  &
                 cck   (ncharm), hck   (ncharm), ock   (ncharm),  &
                 gamma (ncharm), delta (ncharm),                  &
                 rhock (ncharm), pcick (ncharm),                  &
                 xashch(ncharm), cpashc(ncharm),                  &
                 h0ashc(ncharm),                                  &
                 h02ch (ncharm), cp2ch (ncharm),                  &
                 xwatch(ncharm), cp2wat(ncharm)

!      - Parametres cinetiques pour la devolatilisation
!         (Modele de Kobayashi)
!        IY1CH(CH)    --> Indicateur : 0 si MVl = {CH4;CO}
!                                      1 si MVl = {CHz;CO}
!        Y1CH(CH)     --> Coefficient stoechiometrique (adim)
!                         calcule si IY1CH = 0 ; donne si IY1CH = 1
!        A1CH(CH)     --> Facteur pre-exponetielle (1/s)
!         E1CH(CH)     --> Energie d'activation (J/mol)
!        IY2CH(CH)    --> Indicateur : 0 si MVL = {C2H4;CO}
!                                      1 si MVL = {CxHy;CO}
!        Y2CH(CH)     --> Coefficient stoechiometrique (adim)
!                         calcule si IY2CH = 0 ; donne si IY2CH = 1
!        A2CH(CH)     --> Constante preexponetielle (1/s)
!        E2CH(CH)     --> Energie d'activation (J/mol)
!         - Parametres cinetiques pour la combustion heterogene du coke avec O2
!           (Modele a sphere retrecissante)
!        AHETCH(CH)   --> Constante pre-exponentielle (kg/m2/s/atm)
!        EHETCH(CH)   --> Energie d'activation (kcal/mol)
!        IOCHET(CH)   --> Ordre de la reaction 0.5 si = 0 1 si = 1

!         - Parametres cinetiques pour la combustion heterogene du coke avec CO2
!           (Modele a sphere retrecissante)
!        AHETC2(CH)   --> Constante pre-exponentielle (kg/m2/s/atm)
!        EHETC2(CH)   --> Energie d'activation (kcal/mol)
!        IOETC2(CH)   --> Ordre de la reaction 0.5 si = 0 1 si = 1

integer          iy1ch (ncharm), iy2ch (ncharm)
integer          iochet (ncharm) , ioetc2(ncharm)
double precision y1ch  (ncharm), a1ch  (ncharm), e1ch  (ncharm),  &
                 y2ch  (ncharm), a2ch  (ncharm), e2ch  (ncharm),  &
                 ahetch(ncharm), ehetch(ncharm),                  &
                 ahetc2(ncharm), ehetc2(ncharm)

!      - Enthalpie du charbon reactif, coke et cendres
!     ICH(CH)      --> Pointeur dans le tableau EHSOLI pour
!                         le Charbon Reactif
!     ICK(CH)      --> Pointeur dans le tableau EHSOLI pour
!                         le Coke
!     IASH(CH)     --> Pointeur dans le tableau EHSOLI pour
!                         les cendres
!     IWAT(CH)     --> Pointeur dans le tableau EHSOLI pour
!                         l'humidite
!     NSOLID       --> Nb constituants solides (Ch.Reactif, Coke, Ash)
!     NSOLIM       --> Nb maximal de constituants solides
!     EHSOLI(S,IT) --> Enthalpie massique (J/kg) du constituant solide
!                         no S a la temperature T(IT)
!     WMOLS(S)     --> Masse molaire du constituant solide
!     EH0SOL(S)    --- Enthalpie de formation (J/kg) du constituant solide
!                      no S

integer    nsolim
parameter( nsolim = 4*ncharm )

integer          nsolid, ich(ncharm), ick(ncharm), iash(ncharm)
integer          iwat(ncharm)
double precision ehsoli(nsolim,npot), wmols(nsolim)
double precision eh0sol(nsolim)

! ---- PAR CLASSES (grandeurs deduites)

!        NCLACP     --> Nb de classes

integer          nclacp

!      - Proprietes
!        ICHCOR(CL)  --> = ICH si la classe consideree appartient
!                        au charbon ICH (1, 2, ...)
!        DIAM20(CL)  --> Diametre initial (m)
!        DIA2MN(CL)  --> Diametre minimum (m)
!        RHO20(CL)   --> Masse volumique initiale (kg/m3)
!        RHO2MN(CL)  --> Masse volumique minimale (kg/m3)
!        XMP0(CL)    --> Masse initiale de la particule (m)
!        XMASH(CL)   --> Masse de cendres de la particule (m)

integer          ichcor(nclcpm)
double precision diam20(nclcpm), dia2mn(nclcpm),                  &
                 rho20 (nclcpm), rho2mn(nclcpm),                  &
                 xmp0  (nclcpm), xmash (nclcpm)

!--> DONNEES RELATIVES AUX OXYDANTS


!       NOXYD        --> Nombre d'Oxydants (Maximum 3)

integer          noxyd

!       OXYO2       --> composition des oxydants en O2
!       OXYN2       --> composition des oxydants en N2
!       OXYH2O      --> composition des oxydants en H2O
!       OXYCO2      --> composition des oxydants en CO2

double precision  oxyo2(3),oxyn2(3),oxyh2o(3),oxyco2(3)

!--> DONNEES RELATIVES A LA COMBUSTION DES ESPECES GAZEUSES

!        ICHX1C(CH)  --> Pointeur CHx1  pour EHGAZE et WMOLE
!        ICHX2C(CH)  --> Pointeur CHx2  pour EHGAZE et WMOLE
!        ICHX1       --> Pointeur CHx1m pour EHGAZE et WMOLE
!        ICHX2       --> Pointeur CHx2m pour EHGAZE et WMOLE
!        ICO         --> Pointeur CO    pour EHGAZE et WMOLE
!        IO2         --> Pointeur O2    pour EHGAZE et WMOLE
!        ICO2        --> Pointeur CO2   pour EHGAZE et WMOLE
!        IH2O        --> Pointeur H2O   pour EHGAZE et WMOLE
!        IN2         --> Pointeur N2    pour EHGAZE et WMOLE
!        CHX1(CH)    --> Composition de l'hydrocarbure relatif
!                        au MVl : CH(X1)
!        CHX2(CH)    --> Composition de l'hydrocarbure relatif
!                        au MVL : CH(X2)
!        A1(CH),     --> Coefficients stoechiometriques molaires pour
!        B1(CH)          la reaction de devolatilisation a basses T
!        C1(CH)
!        A2(CH),     --> Coefficients stoechiometriques molaires pour
!        B2(CH)          la reaction de devolatilisation a basses T
!        C2(CH)

integer          ichx1c(ncharm), ichx2c(ncharm),                  &
                 ichx1, ichx2, ico, io2, ico2, ih2o, in2
double precision chx1(ncharm), chx2(ncharm),                      &
                 a1(ncharm), b1(ncharm),c1(ncharm),               &
                 a2(ncharm), b2(ncharm),c2(ncharm)

!--> DONNEES COMPLEMENTAIRES RELATIVES AU CALCUL DE RHO
!    SUR LES FACETTES DE BORD

!       IENTAT(IENT) --> Indicateur air par type de facette d'entree
!       IENTCP(IENT) --> Indicateur CP  par type de facette d'entree
!       INMOXY(IENT) --> Numero de l'oxydant
!       TIMPAT(IENT) --> Temperature en K pour l'air relative
!                         a l'entree IENT
!       X20(IENT,    --> Fraction massique dans le melange de charbon
!           ICLA   )     de la classe ICLA relative a l'entree IENT

integer          ientat(nozppm), ientcp(nozppm)
integer          inmoxy(nozppm)
double precision timpat(nozppm), x20(nozppm,nclcpm)

!--> POINTEURS DANS LE TABLEAU TBMCR

integer if1mc(ncharm) , if2mc(ncharm)
integer ix1mc ,ix2mc, ichx1f1, ichx2f2
integer icof1, icof2, ih2of1 , ih2of2

!--> DEFINITION DES COMMONS

common / icpcom / ncharb, nclpch, nclacp, iy1ch , iy2ch ,         &
                  ichcor,                                         &
                  nsolid, ich   , ick   , iash  , iwat  ,         &
                  ichx1c, ichx2c,                                 &
                  ichx1 , ichx2 ,                                 &
                  ico   , io2   , ico2  , ih2o  , in2   ,         &
                  ientat, ientcp, inmoxy,                         &
                  if1mc , if2mc ,                                 &
                  ix1mc , ix2mc , ichx1f1, ichx2f2,               &
                  icof1 , icof2 , iochet , ioetc2 ,               &
                  ih2of1 , ih2of2 , noxyd

common / rcpcom / cch   , hch   , och   ,                         &
                  alpha , beta  , alpham, betam ,                 &
                  pcich , rho0ch,                                 &
                  cck   , hck   , ock   ,                         &
                  gamma , delta ,                                 &
                  rhock , pcick ,                                 &
                  xashch, cpashc, h0ashc, xwatch ,                &
                  h02ch , cp2ch ,                                 &
                  y1ch  , a1ch  , e1ch  , y2ch  , a2ch  , e2ch  , &
                  ahetch, ehetch, ahetc2 , ehetc2 ,               &
                  ehsoli, wmols ,                                 &
                  diam20, dia2mn, rho20 , rho2mn,                 &
                  xmp0  , xmash ,                                 &
                  chx1  , chx2  ,                                 &
                  a1    , b1    , c1    , a2    , b2    , c2    , &
                  timpat, x20   ,                                 &
                  eh0sol, cp2wat, oxyo2 , oxyn2 , oxyh2o, oxyco2

!--> GRANDEURS FOURNIES PAR L'UTILISATEUR EN CONDITIONS AUX LIMITES
!      PERMETTANT DE CALCULER AUTOMATIQUEMENT LA VITESSE, LA TURBULENCE,
!      L'ENTHALPIE D'ENTREE.

!    POUR LES ENTREES UNIQUEMENT , IENT ETANT LE NUMERO DE ZONE FRONTIERE

!       QIMPAT(IENT)           --> Debit       Air          en kg/s
!       TIMPAT(IENT)           --> Temperature Air          en K
!       QIMPCP(IENT,ICHA)      --> Debit       Charbon ICHA en kg/s
!       TIMPCP(IENT,ICHA)      --> Temperature Charbon ICHA en K
!       DISTCH(IENT,ICHA,ICLA) --> Distribution en %masse de la classe ICLA
!                                  pour le charbon ICHA

double precision  qimpat(nozppm)
double precision  qimpcp(nozppm,ncharm), timpcp(nozppm,ncharm)
double precision  distch(nozppm,ncharm,ncpcmx)
!     DOUBLE PRECISION  COEFE(NGAZEM), XSOLID(NSOLIM)
!     DOUBLE PRECISION  F1MC(NCHARM) , F2MC(NCHARM)

common / rcpcli / qimpat       ,                                  &
                  qimpcp               , timpcp               ,   &
                  distch
!    &                  COEFE        , XSOLID        ,
!    &                  F1MC         , F2MC

! Complement Table

 double precision thc(npot)
 common/tablec/   thc

 integer          npoc
 common/itablc/   npoc

! FIN
