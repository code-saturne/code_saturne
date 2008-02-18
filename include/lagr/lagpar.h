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
C                              lagpar.h
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
C 1. Classes et particules
C
C     NCLAGM : NOMBRE MAXIMAL DE CLASSES DE PARTICULES
C
      INTEGER         NCLAGM
      PARAMETER      (NCLAGM = 20)
C
C     NCHARM2 : NOMBRE MAXIMAL DE CLASSES DE CHARBON (voir cpincl.h)
C
      INTEGER         NCHARM2
      PARAMETER      (NCHARM2 = 3)
C
C     NCLSTM : NOMBRE MAXIMUM DE STATISTIQUES VOLUMIQUE PAR GROUPE
C
      INTEGER         NCLSTM
      PARAMETER      (NCLSTM = 100)
C
C=======================================================================
C 2. Conditions aux limites
C
C     NFLAGM : NOMBRE MAXIMAL DE ZONE FRONTIERES
C
      INTEGER         NFLAGM
      PARAMETER      (NFLAGM = 100)
C
C=======================================================================
C 3. Conditions aux limites
C
C     NDLAGM : NOMBRE MAXIMAL DE DONNEES SUR LES PARTICULES (REELS)
C
      INTEGER         NDLAGM
      PARAMETER      (NDLAGM = 50)
C
C     NDLAIM : NOMBRE MAXIMAL DE DONNEES SUR LES PARTICULES (ENTIERS)
C
      INTEGER         NDLAIM
      PARAMETER      (NDLAIM = 10)
C
C=======================================================================
C 4. Schema en temps
C
C     NVGAUS : NOMBRE DE VARIABLES ALEATOIRES GAUSSIENNES PAR PARTICULES
C
      INTEGER         NVGAUS
      PARAMETER      (NVGAUS = 9)
C
C
C=======================================================================
C 5. Mouvement Brownien
C
C     NVGAUS : NOMBRE DE VARIABLES ALEATOIRES GAUSSIENNES PAR PARTICULES
C
      INTEGER         NBRGAU
      PARAMETER      (NBRGAU = 6)
C
C
C=======================================================================
C 6. Variables utilisateurs supplementaires
C
C     NUSVAR : Nombre maximum de variables utilisateur supplementaires
C
      INTEGER         NUSVAR
      PARAMETER      (NUSVAR = 10)
C
C     NUSSTA : Nombre maximum de stats utilisateur supplementaires
C
      INTEGER         NUSSTA
      PARAMETER      (NUSSTA = 20)
C
C     NUSBRD : Nombre maximum interactions particules/frontieres
C              utilisateur supplementaires
C
      INTEGER         NUSBRD
      PARAMETER      (NUSBRD = 10)
C
C=======================================================================
C 7. Affichages et fichiers suite
C
C     NVPLMX : Nombre maximum de variables
C
      INTEGER         NVPLMX
      PARAMETER      (NVPLMX = 50)
C
C=======================================================================
C 8. Visualisation particulaires
C
C     NLISTE : Nombre maximum de particules visualisable
C
      INTEGER         NLISTE
      PARAMETER      (NLISTE = 500)
C
C=======================================================================
C 9. Types d'interaction au bord
C
C     NLISTE : Nombre maximum de particules visualisable
C
      INTEGER         IENTRL     , ISORTL     , IREBOL
      INTEGER         IDEPO1     , IDEPO2     , IDEPO3
      INTEGER         IENCRL     , JBORD1     , JBORD2
      INTEGER         JBORD3     , JBORD4     , JBORD5
      INTEGER         IDEPFA
C
      PARAMETER      (IENTRL =  1, ISORTL =  2, IREBOL =  3)
      PARAMETER      (IDEPO1 =  4, IDEPO2 =  5, IDEPO3 =  6)
      PARAMETER      (IENCRL =  7, JBORD1 =  8, JBORD2 =  9)
      PARAMETER      (JBORD3 = 10, JBORD4 = 11, JBORD5 = 12)
      PARAMETER      (IDEPFA = 13)
C
C=======================================================================
C
C FIN
c@z

