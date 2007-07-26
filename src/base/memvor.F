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
                        SUBROUTINE MEMVOR
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , IAPPEL , NFABOR , IFINIA , IFINRA )

C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE POUR LA METHODE DES VORTEX
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! IAPPEL       ! E  !  -> ! INDIQUE LES DONNES A RENVOYER        !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE DAN!
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU !______________.____._____.______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "vortex.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0, IDBRA0, IFINIA, IFINRA
      INTEGER          IAPPEL, NFABOR
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0

      IF(IAPPEL.EQ.1) THEN
C
        IIREPV = IDEBIA
        IFINIA = IIREPV + NFABOR
C
        IFINRA = IDEBRA
C
        CALL IASIZE('MEMVOR',IFINIA)
C
      ELSEIF(IAPPEL.EQ.2) THEN
C
        IIFAGL = IDEBIA
        IFINIA = IIFAGL + NNENT*ICVMAX
C
        IXYZV  = IDEBRA
        IVISV  = IXYZV  + NNENT*ICVMAX*3
        IW1X   = IVISV  + NNENT*ICVMAX
        IW1Y   = IW1X   + NNENT*ICVMAX
        IW1Z   = IW1Y   + NNENT*ICVMAX
        IW1V   = IW1Z   + NNENT*ICVMAX
        IW2X   = IW1V   + NNENT*ICVMAX
        IW2Y   = IW2X   + NNENT*ICVMAX
        IW2Z   = IW2Y   + NNENT*ICVMAX
        IW2V   = IW2Z   + NNENT*ICVMAX
        IFINRA = IW2V   + NNENT*ICVMAX
C
        CALL IASIZE('MEMVOR',IFINIA)
        CALL RASIZE('MEMVOR',IFINRA)
C
      ELSEIF(IAPPEL.EQ.3) THEN
C
        IIVRCE = IDEBIA
        IFINIA = IIVRCE + NNENT*NVOMAX
C
        IYZCEL = IDEBRA
        IUVORT = IYZCEL + NNENT*ICVMAX*2
        IVVORT = IUVORT + NNENT*ICVMAX
        IWVORT = IVVORT + NNENT*ICVMAX
        IYZVOR = IWVORT + NNENT*ICVMAX
        IYZVOA = IYZVOR + NNENT*NVOMAX*2
        ISIGNV = IYZVOA + NNENT*NVOMAX*2
        IXSIGM = ISIGNV + NNENT*NVOMAX
        IXGAMM = IXSIGM + NNENT*NVOMAX
        IXTMP  = IXGAMM + NNENT*NVOMAX*2
        IXTMPL = IXTMP  + NNENT*NVOMAX
        IFINRA = IXTMPL + NNENT*NVOMAX
C
        CALL IASIZE('MEMVOR',IFINIA)
        CALL RASIZE('MEMVOR',IFINRA)
      ENDIF
C
      RETURN
      END
c@z
