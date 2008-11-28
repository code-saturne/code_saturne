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
C                              alaste.h
C***********************************************************************
C
C  METHODE ALE - MOUVEMENT DE STRUCTURES EN COUPLAGE AVEC CODE_ASTER
C
C NTCAST : NUMERO D'ITERATION DE COUPLAGE AVEC CODE_ASTER
C NBASTE : NOMBRE DE STRUCTURES MOBILES
C NBFAST : NOMBRE DE FACES COUPLEES
C NBNAST : NOMBRE DE NOEUDS COUPLES
C ISYNCP : INDICATEUR D'IMPRESSION DES RESULTATS DES DEUX CODES
C          AUX MEMES INSTANTS (SORTIE ENSIGHT POUR ASTER)
C ASDDLF : BLOCAGE DES DDL DE FORCE
C ASDDLC : BLOCAGE DES DDL CINEMATIQUES
C
      INTEGER           NTCAST
      INTEGER           NBASTE, NBFAST, NBNAST
      INTEGER           IFORAS
      INTEGER           ISYNCP
      INTEGER           ASDDLF(3,NASTMX), ASDDLC(3,NASTMX)
C
      COMMON / IASTER / NTCAST, NBASTE, NBFAST, NBNAST,
     &                  IFORAS,
     &                  ISYNCP,
     &                  ASDDLF, ASDDLC
C
C FIN
C
c@z
