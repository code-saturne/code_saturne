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
C     Copyright (C) 1998-2007 EDF S.A., France
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
C                              vector.h
C***********************************************************************
C
C Longueur de registre LREGIS
C           VR number    VR lenght
C               4       4096
C               8       2048
C              16       1024
C              32        512
C              64        256
C             128        128
C             256         64
C
C IVECTI,IVECTB             INDICATEUR (1/0) DE VECTORISATION
C                           FORCEE OU NON (FACES INTERN ET DE BRD
C
      INTEGER   LREGIS
      PARAMETER(LREGIS=1024)
      INTEGER           IVECTI , IVECTB
      COMMON / IVECTO / IVECTI , IVECTB
c@z
