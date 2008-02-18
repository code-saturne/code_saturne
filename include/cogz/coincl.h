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
C                             coincl.h
C
C***********************************************************************
C
C            INCLUDE POUR LA PHYSIQUE PARTICULIERE RELATIF A
C                           LA COMBUSTION GAZ
C
C Necessite ppppar.h
C
C-----------------------------------------------------------------------
C
C--> MODELE FLAMME DE DIFFUSION (CHIMIE 3 POINTS)
C
C ---- Grandeurs fournies par l'utilisateur dans usd3pc.F
C
C       TINOXY       --> Temperature d'entree pour l'oxydant en K
C       TINFUE       --> Temperature d'entree pour le fuel en K
C       IENTOX       --> indicateur oxydant par type de facette d'entree
C       IENTFU       --> indicateur fuel    par type de facette d'entree
C
C ---- Grandeurs deduiites
C
C       HINOXY       --> Enthalpie massique d'entree pour l'oxydant
C       HINFUE       --> Enthalpie massique d'entree pour le fuel
C       HSTOEA       --> Temperature a la stoechiometrie adiabatique
C       NMAXF        --> Nb de points de tabulation en F
C       NMAXFM       --> Nb maximal de points de tabulation en F
C       NMAXH        --> Nb de points de tabulation en H
C       NMAXHM       --> Nb maximal de points de tabulation en H
C       HH           --> Enthalpie stoechiometrique tabulee
C       FF           --> Richesse tabulee
C       TFH(IF,IH)   --> Tabulation richesse - enthalpie stoechiometrique
C
      INTEGER    NMAXF, NMAXFM, NMAXH, NMAXHM
      PARAMETER( NMAXFM = 15 , NMAXHM = 15)
      INTEGER    IENTOX(NOZPPM), IENTFU(NOZPPM)
C
      DOUBLE PRECISION TINOXY, TINFUE, HINFUE, HINOXY, HSTOEA
      DOUBLE PRECISION HH(NMAXHM), FF(NMAXFM), TFH(NMAXFM,NMAXHM)
C
C
C--> MODELE FLAMME DE PREMELANGE (MODELE EBU)
C
C ---- Grandeurs fournies par l'utilisateur dans usebuc.F
C
C       IENTGF       --> indicateur gaz frais  par type de facette d'entree
C       IENTGB       --> indicateur gaz brules par type de facette d'entree
C       QIMP         --> Debit impose en kg/s
C       FMENT        --> Taux de melange par type de facette d'entree
C       TKENT        --> Temperature en K par type de facette d'entree
C       FRMEL        --> Taux de melange constant pour modeles 0 et 1
C       TGF          --> Temperature gaz frais en K identique
C                        pour premelange frais et dilution
C       CEBU         --> Constante Eddy break-Up
C
C ---- Grandeurs deduites
C
C       HGF          --> Enthalpie massique gaz frais identique
C                        pour premelange frais et dilution
C       TGBAD        --> Temperature adiabatique gaz brules en K
C
C
      INTEGER          IENTGF(NOZPPM), IENTGB(NOZPPM)
      DOUBLE PRECISION FMENT(NOZPPM), TKENT(NOZPPM), QIMP(NOZPPM)
      DOUBLE PRECISION FRMEL, TGF, CEBU, HGF, TGBAD
C
C
C--> DEFINITION DES COMMONS
C
      COMMON / ICOCOM / NMAXF , NMAXH , IENTGF, IENTGB, IENTOX, IENTFU
      COMMON / RCOCOM / TINOXY, TINFUE, HINFUE, HINOXY, HSTOEA,
     &                  HH    , FF    , TFH   ,
     &                  FMENT , TKENT , QIMP  ,
     &                  FRMEL , TGF   , CEBU  ,
     &                  HGF   , TGBAD
C
C
C--> MODELE DE FLAMME DE PREMELANGE LWC
C
C       NDRACM : nombre de pics de Dirac maximum
C       NDIRAC : nombre de Dirac (en fonction du modele)
        INTEGER NDRACM
        PARAMETER (NDRACM = 5)
C
        INTEGER NDIRAC
        COMMON / ILWCDI / NDIRAC
C
C --- Grandeurs fournies par l'utilisateur dans uslwc1.F
C
C       VREF : Vitesse de reference
C       LREF : Longueur de reference
C         TA : Temperature d'activation
C      TSTAR : Temperature de cross-over
C
      INTEGER IRHOL(NDRACM), ITEML(NDRACM), IFMEL(NDRACM)
      INTEGER IFMAL(NDRACM), IAMPL(NDRACM), ITSCL(NDRACM)
      INTEGER IMAML(NDRACM), IHHHH(NDRACM), IMAM
C
      COMMON / RLWCET / IRHOL, ITEML, IFMEL, IFMAL, IAMPL,
     &                  ITSCL, IMAML, IHHHH, IMAM
C
      DOUBLE PRECISION VREF, LREF, TA, TSTAR
      DOUBLE PRECISION FMIN, FMAX, HMIN, HMAX
      DOUBLE PRECISION COEFF1, COEFF2, COEFF3
C
      COMMON / RLWCST / VREF, LREF, TA, TSTAR,
     &                  FMIN, FMAX, HMIN, HMAX,
     &                  COEFF1, COEFF2, COEFF3
C
C
c@z
