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
C                             fuincl.h
C
C***********************************************************************
C
C            INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF
C                A LA COMBUSTION DU FUEL
C
C Necessite ppppar.h et ppthch
C-----------------------------------------------------------------------
C      EPSIFL : Precision pour les tests
C
        DOUBLE PRECISION EPSIFL
        PARAMETER ( EPSIFL = 1.D-8 )
C
C--> DONNEES RELATIVES AU FUEL
C
C      - Proprietes du fuel
C        CFOL      --> fractions massiques elementaires en C, H, O, S, In (%)
C        HFOL          du fuel oil liquid
C        OFOL
C        SFOL
C        XInFOL
C        PCIFOL    --> PCI (J/kg) fuel oil liquid
C        RHO0FL   --> Masse volumique initiale (kg/m3)
C      - Proprietes du coke
C        CKF      --> Fractions massiques elementaires en C, H, O, S, In (%)
C        HKF          du coke
C        OKF
C        SKF
C        XInKF
C        GAMMA    --> Composition du coke
C        DELTA        sous la forme CH(GAMMA)O(DELTA)
C                         GAMMA = HCK/CCK
C                         DELTA = OCK/CCK
C        PCIKF     --> PCI (J/kg) coke
C        RHOKF     --> Masse volumique coke
C        FKC       --> Fraction massique initiale de coke dans le fuel
C        H02FOL    --> H0 du fuel oil liquid
C        CPFOL     --> CP du fuel oil liquid
C        HRFVAP    --> H formation vapeur a Tebu
C        Fractions massiques dans les vapeurs
C        HSFOV     --> H2S
C        COFOV     --> CO
C        CHFOV     --> CHn
C        nHCFOV    --> n dans la formule CHn (un réel, car formule molaire moyenne)
C
C        DFOL      --> densite du fuel liquide
C
      DOUBLE PRECISION CFOL , HFOL , OFOL , SFOL, XInFOL,
     &                 PCIFOL , RHO0FL , RHOKF,
     &                 H02FOL , CP2FOL , HRFVAP, DFOL,
     &                 CKF , HKF , OKF , SKF, XInKF, PCIKF, FKC,
     &                 HSFOV, COFOV, CHFOV, nHCFOV
C
C      - Parametres pour l'evaporation
C      TEVAP1      --> temperature de debut d'evaporation
C      TEVAP2      --> temperature de fin d'evaporation
c
c
c
C        - Parametres cinetiques pour la combustion heterogene du coke
C          (Modele a sphere retrecissante)
C        AHETFL   --> Constante pre-exponentielle (kg/m2/s/atm)
C        EHETFL   --> Energie d'activation (kcal/mol)
C        IOFHET   --> Ordre de la reaction 0.5 si = 0 1 si = 1
C
      DOUBLE PRECISION YFOL , AFOL  , EFOL  ,
     &                 AHETFL , EHETFL, TEVAP1, TEVAP2
      INTEGER          IOFHET
C
C      - Enthalpie du fuel et coke
C     IFOL         --> Pointeur dans le tableau EHSOLI pour
C                         le fuel oil liquid
C     ICK          --> Pointeur dans le tableau EHSOLI pour
C                         le Coke
C     EHSOLI(S,IT) --> Enthalpie massique (J/kg) du constituant solide
C                         no S a la temperature T(IT)
C
      INTEGER          IFOL, IKF
C
C--> DONNEES RELATIVES A LA COMBUSTION DES ESPECES GAZEUSES
C
C        IIFOV        --> Pointeur FOV   pour EHGAZE et WMOLE
C        IICO         --> Pointeur CO    pour EHGAZE et WMOLE
C        IIO2         --> Pointeur O2    pour EHGAZE et WMOLE
C        IICO2        --> Pointeur CO2   pour EHGAZE et WMOLE
C        IIH2O        --> Pointeur H2O   pour EHGAZE et WMOLE
C        IIN2         --> Pointeur N2    pour EHGAZE et WMOLE
C        IIH2S        --> Pointeur H2S   pour EHGAZE et WMOLE
C        IISO2        --> Pointeur SO2   pour EHGAZE et WMOLE
C
C        XSI         --> XSI = 3,76 pour de l'air
C        FVAPMX     --> Maximum pour le traceur F3
C        FOV         --> Composition de l'hydrocarbure relatif
C                        aux matieres volatiles
C        A,     --> Coefficients stoechiometriques molaires pour
C        B          la reaction d'evaporation
C
C        Concentrations dans les espèces globales
C        AFOVF1         nb de moles de vapeur associées à un kg de traceur 1
C        ACOF1                          CO
C        AH2SF1                         H2S
C
C        ACOF3                          CO                                 3
C        AO2F3                          O2
C        AH2SF3                         H2S
C        AH2OF3                         H2O
C        AN2F3                          N2
C        FF3MAX fraction massqique maximale du traceur F3
C               (correspondant à la masse libérée par combustion hétérogène il
C                ne peut exister pur)
C
C        AO2F4         nb de moles de O2 dans l'oxydant associé au traceur 4
C        AN2F4                        N2
C
      DOUBLE PRECISION FVAPMX,
     &                 FOV,
     &                 A, B
      INTEGER IFOV,IH2S,ISO2
C
C--> DONNEES COMPLEMENTAIRES RELATIVES AU CALCUL DE RHO
C    SUR LES FACETTES DE BORD
C
C       IENTAT(IENT) --> Indicateur air par type de facette d'entree
C       IENTFL(IENT) --> Indicateur CFOL  par type de facette d'entree
C       TIMPAT(IENT) --> Temperature en K pour l'air relative
C                         a l'entree IENT
C      X30(IENT)    --> Fraction massique  relative a l'entree IENT
C      XMG0(IENT)
C
      INTEGER          IENTFL(NOZPPM)
      DOUBLE PRECISION X30(NOZPPM), XMG0(NOZPPM)
C
C
C--> POINTEURS DANS LE TABLEAU TBMCR
C
      DOUBLE PRECISION AFOVF1,ACOF1,AH2SF1,ACOF3,AO2F3,AH2SF3,AH2OF3
      DOUBLE PRECISION AN2F3,FF3MAX,AO2F4,AN2F4
C
C--> DEFINITION DES COMMONS
C
      COMMON / IFUCOM / IFOL ,
     &                  IENTFL,
     &                  IOFHET
C
      COMMON / IFUESP / IFOV,IH2S,ISO2
      COMMON / RFUESP / AFOVF1,ACOF1,AH2SF1,ACOF3,AO2F3,AH2SF3,AH2OF3,
     &                  AN2F3,FF3MAX,AO2F4,AN2F4
C
      COMMON / RFUCOM / CFOL   , HFOL  , OFOL , SFOL, XInFOL,
     &                  PCIFOL , RHO0FL, RHOKF ,
     &                  CKF    , HKF   , OKF  , SKF , XInKF,
     &                  PCIKF,FKC,
     &                  HSFOV,COFOV,CHFOV,nHCFOV,
     &                  H02FOL , CP2FOL , HRFVAP,
     &                  YFOL  , AFOL  , EFOL  , DFOL,
     &                  AHETFL, EHETFL,TEVAP1, TEVAP2 ,
     &                  FVAPMX, FOV, A   , B ,
     &                  X30 , XMG0
C
C--> GRANDEURS FOURNIES PAR L'UTILISATEUR EN CONDITIONS AUX LIMITES
C      PERMETTANT DE CALCULER AUTOMATIQUEMENT LA VITESSE, LA TURBULENCE,
C      L'ENTHALPIE D'ENTREE.
C
C    POUR LES ENTREES UNIQUEMENT , IENT ETANT LE NUMERO DE ZONE FRONTIERE
C
C       QIMPAT(IENT)      --> Debit       Air                 en kg/s
C       TIMPAT(IENT)      --> Temperature Air                 en K
C       QIMPFL(IENT)      --> Debit      Fuel Oil Liquid     en kg/s
C       TIMPFL(IENT)      --> Temperature  FOL               en K
C       DINIFL(IENT)      --> Diametre initial des gouttes de FOL  en microns
C
      DOUBLE PRECISION  DINIFL,DINIKF,DINIIN
      DOUBLE PRECISION  QIMPFL(NOZPPM), TIMPFL(NOZPPM)
      DOUBLE PRECISION  HLFM
C
      COMMON / RCPCLI / QIMPFL  , TIMPFL ,
     &                  HLFM     , DINIFL ,DINIKF,DINIIN
C
C FIN
c@z
