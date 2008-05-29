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
C                              dimens.h
C***********************************************************************
C
C... DONNEES DU MAILLAGE
C
      INTEGER           NCELET,NCEL,NFAC,NFABOR,NCELBR,
     &                  NPRFML,NFML,NNOD,LNDFAC,LNDFBR
      COMMON / IGEOME / NCELET,NCEL,NFAC,NFABOR,NCELBR,
     &                  NPRFML,NFML,NNOD,LNDFAC,LNDFBR
C
      INTEGER           NDIM
      COMMON / INDIM3 / NDIM
C
      INTEGER           NVAR, NSCAL, NVISLS, NPHAS
      COMMON / IVARIA / NVAR, NSCAL, NVISLS, NPHAS
C
      INTEGER           NCOFAB
      COMMON / IPCLIM / NCOFAB
C
      INTEGER           NPROCE, NPROFA, NPROFB, NFLUMA
      COMMON / IPRPPH / NPROCE, NPROFA, NPROFB, NFLUMA
C
C FIN
c@z
