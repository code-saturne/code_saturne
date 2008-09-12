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
C                              mltgrd.h
C***********************************************************************
C
C MULTIGRILLE
C -----------
C   NCEGRM : NOMBRE MAX DE CELLULES SUR MAILLAGE LE PLUS GROSSIER
C   NGRMAX : NOMBRE MAX DE NIVEAUX DE MAILLAGES
C   NAGMX0 : PARAMETRE CONSTRUCTION DE  MAILLAGE AUTOMATIQUE
C   IAGMX0 : PARAMETRE CONSTRUCTION DE  MAILLAGE AUTOMATIQUE
C   NCPMGR : Si > 0, active le post traitement de l'agglomeration, en
C            projetant les numeros de cellules grossieres sur le
C            maillage fin (modulo NCPMGR(IVAR))
C
      INTEGER           NCEGRM, NGRMAX,
     &                  NAGMX0(NVARMX), IAGMX0(NVARMX), NCPMGR(NVARMX)
      COMMON / IMULTG / NCEGRM, NGRMAX,
     &                  NAGMX0, IAGMX0, NCPMGR
C
C FIN
c@z
