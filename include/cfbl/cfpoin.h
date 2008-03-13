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
C                              cfpoin.h
C***********************************************************************
C
C
C      ACCESSIBLES DIRECTEMENT DANS IA, RA
C
C Pointeur Dimension       Description
C IIFBET ! NFABOR*NPHAS  ! INDICATEUR FLUX THERMIQUE AU BORD IMPOSE
C                          (il faut annuler des contributions de bord
C                           de l'eq de E)
C IIFBRU ! NFABOR*NPHAS  ! INDICATEUR FLUX DE BORD CALCULE PAR RUSANOV
C                          (il faut annuler des contributions de bord
C                           de l'eq de Qdm et de l'eq de E)
C
      INTEGER           IIFBET , IIFBRU
      COMMON / IPFBCF/  IIFBET , IIFBRU
C
C FIN
c@z
