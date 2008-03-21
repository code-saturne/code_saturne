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
C NCEGRM :
C   NOMBRE MAX DE CELLULES SUR MAILLAGE LE + GROSSIER
C NGRMAX :
C   NOMBRE MAX DE NIVEAUX DE MAILLAGES
C JNCEL,JNFAC,JIFACE,JIRESP,JDA,JXA :
C   POINTEURS PARTIELS SUR LES DESCRIPTEURS DES MAILLAGES GROSSIERS
C IPGREN :
C   POINTEURS SUR DESCRIPTEURS ENTIERS DES MAILLAGES GROSSIERS
C IPGRRE :
C   POINTEURS SUR DESCRIPTEURS REELS DES MAILLAGES GROSSIERS
C
      INTEGER           NCEGRM , NGRMAX ,
     &                  JNCEL  , JNFAC  , JIFACE , JIRESP ,
     &                  JDA    , JXA    ,
     &                  IPGREN(NTGREN,NGRMMX),IPGRRE(NTGRRE,NGRMMX)
      COMMON / IPAMGR / NCEGRM , NGRMAX ,
     &                  JNCEL  , JNFAC  , JIFACE , JIRESP ,
     &                  JDA    , JXA    ,
     &                  IPGREN               ,IPGRRE
C
C FIN
c@z
