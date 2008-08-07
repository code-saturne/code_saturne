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
C            INCLUDE POUR LA PHYSIQUE PARTICULIERE
C                        VARIABLE COMMUNE ENTRE
C                    COMBUSTION DU CHARBON PULVERISE
C                    COMBUSTION DU FIOUL LOURD
C
C        XSI         --> XSI = 3,76 pour de l'air
C
      DOUBLE PRECISION  XSI
      COMMON / RCPFU1 / XSI
C
C   nb de moles de I dans J
C
      DOUBLE PRECISION AO2F3,ACOF3,AN2F3,AH2OF3
      DOUBLE PRECISION AO2F4,AN2F4,AH2OF4,ACO2F4
      DOUBLE PRECISION AH2OF5
      DOUBLE PRECISION AO2F6,AN2F6,AH2OF6,ACO2F6
      DOUBLE PRECISION AO2F7,AN2F7,AH2OF7,ACO2F7
C
      COMMON / RCPFU2 / AO2F3,ACOF3,AN2F3,AH2OF3,
     &                  AO2F4,AN2F4,AH2OF4,ACO2F4,
     &                  AH2OF5,
     &                  AO2F6,AN2F6,AH2OF6,ACO2F6,
     &                  AO2F7,AN2F7,AH2OF7,ACO2F7
C
C Equation sur YCO2
C
       INTEGER         IEQCO2 , IYCO2
       COMMON/EQUCO2 / IEQCO2 , IYCO2
C
C Combustion heterogene avec le  CO2
C
       INTEGER         IHTCO2
       COMMON/EHTCO2 / IHTCO2
C
C Equation sur NOX :
C ================
C
C   IEQNOX = 0 pas de NOx
C          = 1 calcul du NOx
C
       INTEGER         IEQNOX
       COMMON/EQUNOX / IEQNOX
C
C   Scalaires supplementaires : fraction massique de HCN et NO
C                               temperature air
C
       INTEGER         IYHCN , IYNO , ITAIRE
       COMMON/EQUNOX / IYHCN , IYNO , ITAIRE
C
C   Propce supplementaires :
C
C         Conversion HCN en NO       : EXP(-E1/RT)
C         Conversion HCN en NO       : EXP(-E2/RT)
C         NO thermique (Zel'dovitch) : EXP(-E3/RT)
C
C
       INTEGER         IGHCN1 , IGHCN2 , IGNOTH
       COMMON/PRONOX / IGHCN1 , IGHCN2 , IGNOTH
C
C   Temperature moyenne d'entree
C   Taux de vapeur moyen
C
       DOUBLE PRECISION TAIRE
       COMMON /NOXDBL/  TAIRE

C
C FIN
c@z
