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
C                             ppthch.h
C
C***********************************************************************
C
C            INCLUDE THERMOCHIMIE POUR LA PHYSIQUE PARTICULIERE
C
C-----------------------------------------------------------------------
C
C
C
C--> CONSTANTES THERMOCHIMIE
C
C       RR           --> Constante des gaz parfaits en J/mol/K
C       TREFTH       --> Temperature de reference (K)
C       VOLMOL       --> Volume molaire dans les conditions NTP
C                        T = 0 C et P = 1 atm
C
      DOUBLE PRECISION RR
      DOUBLE PRECISION TREFTH, PREFTH, VOLMOL
      PARAMETER ( RR     = 8.31434D0      ,
     &            TREFTH = 25.D0 + TKELVI ,
     &            PREFTH = 1.01325D5      ,
     &            VOLMOL = 22.41D-3       )
C
C--> DONNEES
C
C       NRGAZ        --> Nb de reactions globales en phase gaz
C       NRGAZM       --> Nb maximal de reactions globales en phase gaz
C       NATO         --> Nb d especes atomiques (C,H,..)
C       NATOM        --> Nb maximal d especes atomiques (C,H,..)
C       NGAZE        --> Nb de constituants gazeux elementaires
C       NGAZEM       --> Nb maximal de constituants gazeux elementaires
C       NGAZG        --> Nb d especes globales (ex:Fuel,Oxyd,Prod1,Prod2)
C       NGAZGM       --> Nb maximal d especes globales
C       NPO          --> Nb de points de tabulation
C       NPOT         --> Nb maximal de points de tabulation
C       TH           --> Temperature en Kelvin
C       EHGAZG(G,IT) --> Enthalpie massique (J/kg) de l espece globale
C                        no G a la temperature T(IT)
C       WMOLG(G)     --> Masse molaire de l espece globale
C       EHGAZE(G)    --> Enthalpie massique (J/kg) constituant gazeux
C                        elementaire no E a la temperature T(IT)
C       WMOLE(G)     --> Masse molaire du constituant gazeux elementaire
C       WMOLAT(E)    --> Masse molaire des atomes (C,H,..)
C       IATC, IATH   --> Pointeur dans WMOLEL pour les ecpeces
C       IATO, IATN, IATS       elementaires (C,H,..)
C       FS(R)        --> Taux de melange pour la reaction gloable R
C       STOEG(G,R)   --> Stoechio en especes globales des reactions
C                        pour l espece no G et pour la reaction no R
C       CKABSG(G)    --> Coefficient d'absorption des especes globales
C       CKABS1       --> Coefficient d'absorption du melange gazeux
C                        (en CP)
C       DIFTL0       --> Diffusivite dynamique en kg/(m s)
C
      INTEGER    NGAZGM, NGAZEM, NPOT, NATOM, NRGAZM
      PARAMETER( NGAZGM = 25 , NGAZEM = 20 ,
     &           NPOT  = 500 , NATOM  = 5   , NRGAZM = 1 )
      INTEGER    IATC, IATH, IATO, IATN , IATS
      PARAMETER( IATC = 1, IATH = 2, IATO = 3, IATN = 4 , IATS = 5 )
C
      INTEGER           NPO, NGAZE, NGAZG, NATO, NRGAZ
      COMMON / TCHPPI / NPO, NGAZE, NGAZG, NATO, NRGAZ
C
      DOUBLE PRECISION  TH(NPOT),
     &                  EHGAZE(NGAZEM,NPOT), EHGAZG(NGAZGM,NPOT),
     &                  WMOLE(NGAZEM), WMOLG(NGAZGM), WMOLAT(NATOM),
     &                  STOEG(NGAZGM,NRGAZM), FS(NRGAZM),
     &                  CKABSG(NGAZGM), CKABS1,
     &                  DIFTL0, XCO2, XH2O
C ..v.7..1....v....2....v....3....v....4....v....5....v....6....v....7.I
      COMMON / TCHPPR / TH, EHGAZE, EHGAZG, WMOLE, WMOLG, WMOLAT,
     &                  STOEG, FS, CKABSG, CKABS1, DIFTL0, XCO2, XH2O
C
C
c@z
