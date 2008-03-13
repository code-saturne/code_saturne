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
C                              alstru.h
C***********************************************************************
C
C  METHODE ALE - MOUVEMENT DE STRUCTURES EN COUPLAGE INTERNE
C
C NBSTRU : NOMBRE DE STRUCTURES MOBILES
C
C XMSTRU : MATRICE DE MASSE DES STRUCTURES (kg)
C XCSTRU : MATRICE DE FRICTION DES STRUCTURES (kg/s)
C XKSTRU : MATRICE DE RAIDEUR DES STRUCTURES (kg/s2 = N/m)
C XSTREQ : VECTEUR ECART DE LA POSITION DES STRUCTURE DANS LE MAILLAGE
C          INITIAL PAR RAPPORT A LEUR POSITION D'EQUILIBRE (m)
C XSTR   : VECTEUR DEPLACEMENT DES STRUCTURES PAR RAPPORT A LEUR POSITION
C          DANS LE MAILLAGE INITIAL (m)
C XPSTR  : VECTEUR VITESSE DES STRUCTURES (m/s)
C XPPSTR : VECTEUR ACCELERATION DES STRUCTURES (m/s2)
C FORSTR : VECTEUR FORCE EXERCE SUR LES STRUCTURES (N)
C XSTA   : VALEUR DE XSTR AU PAS DE TEMPS PRECEDENT
C XPSTA  : VALEUR DE XPSTR AU PAS DE TEMPS PRECEDENT
C XPPSTA : VALEUR DE XPPSTR AU PAS DE TEMPS PRECEDENT
C FORSTA : VALEUR DE FORSTR AU PAS DE TEMPS PRECEDENT
C XSTP   : VALEUR PREDITE DE XSTR
C FORSTP : VALEUR PREDITE DE FORSTR
C DTSTR  : PAS DE TEMPS ASSOCIE AU MOUVEMENT DES STRUCTURES
C AEXXST : COEFFICIENT DE PREDICTION DU DEPLACEMENT (SUR XPSTR)
C BEXXST : COEFFICIENT DE PREDICTION DU DEPLACEMENT (SUR XPSTR-XPSTA)
C CFOPRE : COEFFICIENT DE PREDICTION DES EFFORTS
C
C
      INTEGER          NBSTRU
C
      COMMON / ISTRUC / NBSTRU
C
      DOUBLE PRECISION XMSTRU(3,3,NSTRMX)
      DOUBLE PRECISION XCSTRU(3,3,NSTRMX)
      DOUBLE PRECISION XKSTRU(3,3,NSTRMX)
      DOUBLE PRECISION XSTR(3,NSTRMX)  ,XSTA(3,NSTRMX)
      DOUBLE PRECISION XSTP(3,NSTRMX)  ,XSTREQ(3,NSTRMX)
      DOUBLE PRECISION XPSTR(3,NSTRMX) ,XPSTA(3,NSTRMX)
      DOUBLE PRECISION XPPSTR(3,NSTRMX),XPPSTA(3,NSTRMX)
      DOUBLE PRECISION FORSTR(3,NSTRMX),FORSTA(3,NSTRMX)
      DOUBLE PRECISION FORSTP(3,NSTRMX)
      DOUBLE PRECISION DTSTR(NSTRMX)
      DOUBLE PRECISION AEXXST, BEXXST, CFOPRE
C
      COMMON / RSTRUC / XMSTRU, XCSTRU, XKSTRU, XSTR  , XSTA  , XSTP,
     &                  XSTREQ, XPSTR , XPSTA , XPPSTR, XPPSTA,
     &                  FORSTR, FORSTA, FORSTP, DTSTR,
     &                  AEXXST, BEXXST, CFOPRE
C
C  METHODE DE NEWMARK HHT
C
C  ALPNMK : COEFFICIENT ALPHA DE LA METHODE DE NEWMARK HHT
C  BETNMK : COEFFICIENT BETA  DE LA METHODE DE NEWMARK HHT
C  GAMNMK : COEFFICIENT GAMMA DE LA METHODE DE NEWMARK HHT
      DOUBLE PRECISION ALPNMK, GAMNMK, BETNMK
      COMMON / RNEMRK / ALPNMK, GAMNMK, BETNMK
C
C FIN
C
c@z
