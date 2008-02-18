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
C                              albase.h
C***********************************************************************
C
C  METHODE ALE
C  IALE   : UTILISATION DE LA METHODE ALE
C         = 0 SANS METHODE ALE
C         = 1 AVEC METHODE ALE
C  IIMPAL : POINTEUR SUR IMPALE, INDICATEUR DE DEPLACEMENT IMPOSE
C  IXYZN0 : POINTEUR SUR XYZNO0, POSITION INITIALE DU MAILLAGE
C  IDEPAL : POINTEUR SUR DEPALE, DEPLACEMENT DU MAILLAGE
C  IIALTY : POINTEUR SUR IALTYB, TYPE DE BORD
C  NALINF : NOMBRE D'ITERATIONS D'INITIALISATION DU FLUIDE
C  NALIMX : NOMBRE MAXIMAL D'ITERATIONS D'IMPLICITATION DU DEPLACEMENT
C           DES STRUCTURES
C  IORTVM : TYPE DE VISCOSITE DE MAILLAGE
C         = 0 ISOTROPE
C         = 1 ORTHOTROPE
C  EPALIM : PRECISION RELATIVE D'IMPLICITATION DU DEPLACEMENT DES
C           STRUCTURES
C  ITALIN : ITERATION D'INITIALISATION DE l'ALE
C         = 0 NON
C         = 1 OUI
C
      INTEGER           IALE  , IIMPAL, IXYZN0, IDEPAL, IIALTY, NALINF
      INTEGER           NALIMX, IORTVM, ITALIN
      COMMON / ICMALE / IALE  , IIMPAL, IXYZN0, IDEPAL, IIALTY, NALINF,
     &                  NALIMX, IORTVM, ITALIN
C
      DOUBLE PRECISION EPALIM
      COMMON / RCMALE / EPALIM
C FIN
C
c@z
