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
C                              ppincl.h
C
C***********************************************************************
C
C            INCLUDE GENERAL PROPRE A LA PHYSIQUE PARTICULIERE
C
C necessite ppthch et ppppar
C-----------------------------------------------------------------------
C
C
C--> TABLEAU INDICATEURS DU CHOIX DE LA PHYSIQUE PARTICULIERE CHOISIE
C
      INTEGER   NMODMX
      PARAMETER(NMODMX = 50)
      INTEGER           IPPMOD(NMODMX)
      COMMON / IIPPMD / IPPMOD
C
C ---- Indicateur global de physique particuliere
C        IPPMOD(IPHPAR) = 0 : pas de physique particuliere
C                         1 : physique particuliere enclenchee
C                         2 : physique particuliere avec pilotage du
C                             rayonnement par fichier parametrique
      INTEGER IPHPAR
C
C ---- Modeles propres a la combustion gaz ICO...
      INTEGER ICOD3P, ICODEQ, ICOEBU, ICOBML, ICOLWC
C
C ---- Modeles propres a la combustion charbon pulverise ICP...
      INTEGER ICP3PL
C
C ---- Modeles propres aux versions effet Joule et conduction ionique
      INTEGER  IELJOU, IELARC, IELION
C
C ---- Modeles propres a la combustion charbon pulverise couplee Lagrangien
      INTEGER ICPL3C
C ---- Modeles propres a la combustion fuel
      INTEGER ICFUEL
C
C ---- Modele compressible
      INTEGER  ICOMPF
C
      PARAMETER       (IPHPAR = 1 , ICOD3P = 2 , ICODEQ = 3 ,
     &                 ICOEBU = 4 , ICOBML = 5 , ICOLWC = 6 ,
     &                 ICP3PL = 7 , ICPL3C = 8 , ICFUEL = 9 ,
     &                 IELJOU = 10, IELARC = 11, IELION = 12,
     &                 ICOMPF = 13)
C
C--> MODULE RAYONNEMEMT
C    IRAYPP =  0 : pas de rayonnement
C           =  1 : DOM + calcul a partir des CABS des especes
C           =  2 : DOM + CABS par Modak
C           =  3 : P-1 + calcul a partir des CABS des especes
C           =  4 : P-1 + CABS par Modak
      INTEGER IRAYPP
C
C
C--> NOMBRE DE VARIABLES ALGEBRIQUES OU D'ETAT
C    pour la physique particuliere NSALPP
C    total NSALTO
      INTEGER NSALPP, NSALTO
C
C
C--> POINTEURS VARIABLES COMBUSTION GAZ
C
C ---- Variables transportees
      INTEGER IFM, IFP2M, IYGFM, ICM, ICP2M, IFPCPM
      INTEGER IYFM, IYFP2M, ICOYFP
C
C ---- Variables d'etat
      INTEGER IYM(NGAZGM), ITEMP, IFMIN, IFMAX
      INTEGER ICKABS, IT4M, IT3M
C
C --- Pointeurs proprietes (PROPCE)
      INTEGER ITSC
C
C--> POINTEURS VARIABLES COMBUSTION CHARBON PULVERISE
C
C ---- Variables transportees
C        Phase continue (melange gazeux)
      INTEGER IF1M(NCHARM), IF2M(NCHARM)
      INTEGER IF3M, IF4M, IF4P2M, IF5M, IF6M, IF7M
      INTEGER IF3MC2
C        Phase dispersee (classe de particules)
      INTEGER IXCK(NCLCPM), IXCH(NCLCPM), INP(NCLCPM)
      INTEGER IH2(NCLCPM) , IXWT(NCLCPM)
C
C ---- Variables d'etat
C        Phase continue (melange gazeux)
      INTEGER IYM1(NGAZEM), ITEMP1 , IROM1 , IMMEL
C        Phase dispersee (classes de particules)
      INTEGER ITEMP2(NCLCPM), IROM2(NCLCPM), IDIAM2(NCLCPM)
      INTEGER IX2(NCLCPM)
      INTEGER IGMDCH(NCLCPM), IGMHET(NCLCPM) , IGHCO2(NCLCPM)
      INTEGER IGMDV1(NCLCPM), IGMDV2(NCLCPM)
      INTEGER IGMSEC(NCLCPM)
C        Rayonnement : integrale de la luminance sur 4*PI
      INTEGER ILUMI
C
C--> POINTEURS VARIABLES COMBUSTION FUEL
C
C ---- Variables transportees
C        Phase continue
       INTEGER IFVAP, IFHTF
C        Phase dispersee
       INTEGER IHLF(NCLCPM)
       INTEGER IXKF(NCLCPM), IXFOL(NCLCPM), ING(NCLCPM)
       INTEGER ITEMP3(NCLCPM), IROM3(NCLCPM), IDIAM3(NCLCPM)
       INTEGER IX3(NCLCPM)
C
C ---- Variables d'etat
C        Phase continue
       INTEGER IYFOL(NCLCPM)
C        Phase dispersee
       INTEGER IH1HLF(NCLCPM), IGMHTF(NCLCPM), IGMEVA(NCLCPM)
C
C--> POINTEURS VARIABLES VERSION ELECTRIQUES
C
C        Dimension des 'vecteurs electriques'      = NDIMVE
      INTEGER    NDIMVE
      PARAMETER (NDIMVE = 3)
C
C ---- Variables transportees
C        Potentiel reel       = IPOTR
C        Potentiel imaginaire = IPOTI
C        Composantes du potentiel vecteur magnetique = IPOTVA()
C        Fraction massique des constituants = IYCOEL()
C
      INTEGER IPOTR, IPOTI, IPOTVA(NDIMVE), IYCOEL(NGAZGM)
C
C ---- Variables d'etat dans PROPCE
C        Puissance volumique dissipee par effet Joule W/m3 = IEFJOU
C        Forces electromagnetiques de Laplace en N/m3      = ILAPLA()
C        Charge electrique volumique C/m3                  = IQELEC
C        Densite de courant electrique reelle A/m2         = IDJR()
C        Puissance volumique rayonnee W/m3
C        Densite de courant electrique imaginaire en A/m2  = IDJI()
C         ou coeff d'absorption en m-1                     = IDRAD
C
C      Les autres variables deduites seront locales pour
C        economiser de la place memoire.
C        Cela n'empeche pas de les sortir en post-traitement.
C        Gradient du potentiel reel en V/m                 = IGPOTR()
C        Gradient du potentiel imaginaire en V/m           = IGPOTI()
C        Composantes du champ magnetique en Tesla          = IBMG()
C
C
      INTEGER IEFJOU , IQELEC , ILAPLA(NDIMVE)
      INTEGER IDJR(NDIMVE)    , IDJI(NDIMVE)   , IDRAD
C
C--> POINTEURS VARIABLES CONDUCTION IONIQUE
C
      INTEGER   NESIOM
      PARAMETER (NESIOM = 10)
      INTEGER   NESPIO
C
C ---- Variables transportees
C        par espece
      INTEGER IYMION(NESIOM)
C
C ---- Variables d'etat
C
C
C--> POINTEURS COMPRESSIBLE
C
C ---- Variables transportees par phase
      INTEGER IRHO(NPHSMX), IENERG(NPHSMX), ITEMPK(NPHSMX)
C ---- Proprietes supplementaires par phase
      INTEGER ICV(NPHSMX), IVISCV(NPHSMX), IEOS(NPHSMX)
C
C     COMMON complete plus bas
C
C ---- Aliases pour les conditions aux limites
      INTEGER IRUN(NPHSMX), IRUNH(NPHSMX)
C
      COMMON / IACLCF / IRUN , IRUNH
C
C ---- Proprietes supplementaires par phase
      DOUBLE PRECISION CV0(NPHSMX), VISCV0(NPHSMX)
C
      COMMON / RPOCFP / CV0 , VISCV0
C
C ---- Prediction de pression par une equation d'evolution
      INTEGER IPPRED(NPHSMX)
C ---- Flux de masse specifique pour la vitesse
      INTEGER IFLMAU(NPHSMX)
C ---- Utilisation de la pression predite pour resoudre Navier-Stokes
      INTEGER IGRDPP(NPHSMX)
C --- Conditions aux limites prenant en compte l'equilibre hydrostatique
      INTEGER ICFGRP(NPHSMX)
C
      COMMON / IPOCFO / IPPRED , IFLMAU , IGRDPP , ICFGRP
C
C ---- Flux de bord convectifs QDM et energie (numero de PROPFB)
      INTEGER           IFBRHU(NPHSMX) , IFBRHV(NPHSMX) ,
     &                  IFBRHW(NPHSMX) , IFBENE(NPHSMX)
C
      COMMON / IPOBCF/  IFBRHU         , IFBRHV         ,
     &                  IFBRHW         , IFBENE
C
C--> POINTEUR RELATIF A LA VARIABLE ENTHALPIE
C
      INTEGER IHM
C
C--> REMPLISSAGE COMMON RAYONNEMENT
C
      COMMON / IPORAY / IRAYPP
C
C--> REMPLISSAGE COMMON POINTEURS VARIABLES TRANSPORTEES
C                                 VARIABLES D'ETAT
C
      COMMON / IPOVST /
C
C ---- Combustion gaz
     &                  IFM, IFP2M, IYGFM, ICM, ICP2M, IFPCPM,
     &                  IYFM, IYFP2M, ICOYFP,
C ---- Combustion charbon pulverise
     &                  IF1M, IF2M, IF3M, IF4M, IF4P2M, IF5M,
     &                  IF6M, IF7M, IF3MC2 ,
     &                  IXCK, IXCH, INP , IH2 , IXWT  ,
C ---- Combustion fuel
     &                  IHLF, IFVAP, IFHTF,
     &                  IXKF, IXFOL, ING, IX3 ,
C
C ---- Versions electriques
     &                  IPOTI, IPOTR, IPOTVA, IYCOEL,
C ---- Conduction ionique
     &                  NESPIO, IYMION,
C
C ---- Compressible
     &                  IRHO   , IENERG , ITEMPK ,
     &                  ICV    , IVISCV , IEOS   ,
C ---- Enthalpie
     &                  IHM
C
      COMMON / IPOVSA /
C
C ---- Nb de variables d'etat ou algebriques
     &                  NSALPP, NSALTO,
C ---- Combustion gaz
     &                  IYM, ITEMP, IFMIN, IFMAX, ICKABS, IT4M,
     &                                                    IT3M,
C ---- Combustion charbon pulverise
     &                  IYM1, ITEMP1, IROM1 ,IMMEL,
     &                  ITEMP2, IROM2, IDIAM2, IX2,
     &                  IGMDCH, IGMDV1, IGMDV2, IGMHET, IGHCO2 ,
     &                  IGMSEC,
C ---- Combustion fuel
     &                  IYFOL, ITEMP3, IROM3 , IDIAM3,
     &                  IH1HLF, IGMEVA, IGMHTF,
C ---- Rayonnement
     &                  ILUMI,
C ---- Versions electriques
     &                  IEFJOU, ILAPLA , IQELEC ,
     &                  IDJR  , IDJI   , IDRAD
C
C--> Modele de flamme de premelange LWC
C
      COMMON / ILWCPP / ITSC
C
C--> OPTIONS NUMERIQUES
C
C ---- Coefficient de relaxation de la masse volumique
C      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
      DOUBLE PRECISION SRROM
C
      COMMON / ROPTCP / SRROM

C
C--> GRANDEURS FOURNIES PAR L'UTILISATEUR EN CONDITIONS AUX LIMITES
C      PERMETTANT DE CALCULER AUTOMATIQUEMENT LA VITESSE, LA TURBULENCE,
C      L'ENTHALPIE D'ENTREE.
C    LES GRANDEURS CI-DESSOUS SONT COMMUNES A LA COMBUSTION GAZ ET AU
C      CHARBON.
C
C    POUR LES ENTREES UNIQUEMENT , IENT ETANT LE NUMERO DE ZONE FRONTIERE
C
C       IQIMP (IENT) --> Indicateur zone a debit impose
C       ICALKE(IENT) --> Indicateur type de condition sur la turbulence
C         0 : Utilisateur donne les valeurs
C         1 : Automatique a partir de DH
C                                         et de la vitesse d'entree
C         2 : Automatique a partir de l'intensite turbulente
C                                         et de la vitesse d'entree
C       XINTUR(IENT) --> Intensite turbulente (k=1.5(UREF*XINTUR)**2)
C       DH    (IENT) --> Diametre hydraulique
C       QCALC (IENT) --> Debit calcule  : raf la ; direc ds le sspgm
C
      INTEGER          IQIMP(NOZPPM)  , ICALKE(NOZPPM)
      DOUBLE PRECISION XINTUR(NOZPPM) , DH(NOZPPM)
c, QCALC(NOZPPM)
C
      COMMON / IPPCLI / IQIMP         , ICALKE
      COMMON / RPPCLI / XINTUR        , DH
c, QCALC
C
C
C Pointeur dans IA sur IZFPPP pour reperage des zones frontieres associees
C   aux faces de bord
C Peut etre serait il plus approprie de le verser dans pointe
C
C
      INTEGER           IIZFPP
      COMMON / IFROPP / IIZFPP
C
C NZFPPP Nombre de zones de bord (sur le proc courant)
C ILZPPP Liste des numeros de zone de bord (du proc courant)
C NOZAPM Numero de zone de bord atteint max
C   exemple zones 1 4 2 : NZFPPP=3,NOZAPM=4
C
      INTEGER           NOZAPM, NZFPPP, ILZPPP(NBZPPM)
      COMMON / IZONPP / NOZAPM, NZFPPP, ILZPPP
C
C FIN
c@z
