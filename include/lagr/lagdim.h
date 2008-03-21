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
C                              lagdim.h
C***********************************************************************
C
C=======================================================================
C
C     Include pour le module Lagrangien : dimensions
C
C         Trois fichiers complementaires
C                            lagran.h qui porte les non dimensions
C                            lagdim.h qui porte les dimensions variables
C                            lagpar.h qui porte les parametres
C
C=======================================================================
C 1. Connectivite
C
C     LONGUEUR DU TABLEAU DU CONNECTIVITE CELLULES -> FACES
C     (calcule dans le sous-programme LAGINI)
C
      INTEGER           LNDNOD
      COMMON / ILAGD1 / LNDNOD
C
C=======================================================================
C 2. Classes et particules
C
C     NBPMAX : NOMBRE MAXIMAL DE PARTICULES AUTORISE DANS LE DOMAINE
C              AU COUR DU CALCUL (UTILE SI INJECTION INSTATIONNAIRE)
C
      INTEGER           NBPMAX
      COMMON / ILAGD2 / NBPMAX
C
C=======================================================================
C 3. Dimensions des tableaux particulaires
C
C     NVP          : NOMBRE DE VARIABLES SUR LES PARTICULES
C
C     NVP1         : NOMBRE DE VARIABLES SUR LES PARTICULES
C                     EN ENLEVANT POSITION, VITESSE PARTICULE
C                     ET VITESSE FLUIDE
C
C     NVEP         : NOMBRE D'INFO SUR LES PARTICULES (REELS)
C
C     NIVEP        : NOMBRE D'INFO SUR LES PARTICULES (ENTIERS)
C
C     NTERSL       : NOMBRE DE TERMES SOURCES POUR COUPLAGE RETOUR
C
C     NVLSTA       : NOMBRE DE VARIABLES STATISTIQUES
C
C     NVISBR       : NOMBRE DE VARIABLES A ENREGISTRER SUR LES FRONTIERES
C
C
      INTEGER           NVP    , NVP1   , NVEP   , NIVEP  ,
     &                  NTERSL , NVLSTA , NVISBR
      COMMON / ILAGD3 / NVP    , NVP1   , NVEP   , NIVEP  ,
     &                  NTERSL , NVLSTA , NVISBR
C
C FIN
c@z

