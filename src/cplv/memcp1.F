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
                        SUBROUTINE MEMCP1
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NVAR   , NCELET , NCEL   , NFAC   , NFABOR ,
     &   NTBCPI , ICPWI  ,
     &   NTBCPR , ICPWR  ,
     &   NTBMCI , IMCWI  ,
     &   NTBMCR , IMCWR  ,
     &   NTBWOI , IWORI  ,
     &   NTBWOR , IWORR  ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   PHYSIQUE PARTICULIERE : FLAMME CHARBON PULVERISE
CFONC     GESTION MEMOIRE CALCUL DES PROPRIETES PHYSIQUES PHASE GAZ
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0/IDBRA0! E  !  -> ! POINTEUR DE LA PREMIERE CAS LIBRE DES!
CARGU !              !    !     !  TABLEAUX IA/RA                      !
CARGU ! NVAR         ! E  !  -> ! NOMBRE DE VARIABLES                  !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NTBCPI       ! E  !  -> ! TAILLE DU MACRO TABLEAU CP ENTIERS   !
CARGU ! ICPWI        ! E  ! <-  ! POINTEUR  MACRO TABLEAU CP ENTIERS   !
CARGU ! NTBCPR       ! E  !  -> ! TAILLE DU MACRO TABLEAU CP REELS     !
CARGU ! ICPWR        ! E  ! <-  ! POINTEUR  MACRO TABLEAU CP REELS     !
CARGU ! NTBMCI       ! E  !  -> ! TAILLE DU MACRO TABLEAU MC ENTIERS   !
CARGU ! IMCWI        ! E  ! <-  ! POINTEUR  MACRO TABLEAU MC ENTIERS   !
CARGU ! NTBMCR       ! E  !  -> ! TAILLE DU MACRO TABLEAU MC REELS     !
CARGU ! IMCWR        ! E  ! <-  ! POINTEUR  MACRO TABLEAU MC REELS     !
CARGU ! NTBWOI       ! E  !  -> ! TAILLE DU MACRO TABLEAU WORK ENTIERS !
CARGU ! IWORI        ! E  ! <-  ! POINTEUR  MACRO TABLEAU WORK ENTIERS !
CARGU ! NTBWOR       ! E  !  -> ! TAILLE DU MACRO TABLEAU WORK REELS   !
CARGU ! IWORR        ! E  ! <-  ! POINTEUR  MACRO TABLEAU WORK REELS   !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS RA EN SORTIE                   !
CARGU .______________.____._____.______________________________________.
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
      INCLUDE "optcal.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 ,IDBRA0
      INTEGER          NVAR
      INTEGER          NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NTBCPI , ICPWI  , NTBCPR , ICPWR
      INTEGER          NTBMCI , IMCWI  , NTBMCR , IMCWR
      INTEGER          NTBWOI , IWORI  , NTBWOR , IWORR
      INTEGER          IFINIA , IFINRA
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA , IDEBRA
C
C***********************************************************************
C
C---> INITIALISATION
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
C
      ICPWI   =       IDEBIA
      IMCWI   =       ICPWI   + NTBCPI*NCELET
      IWORI   =       IMCWI   + NTBMCI*NCELET
      IFINIA =        IWORI   + NTBWOI*NCELET
C
      ICPWR   =       IDEBRA
      IMCWR   =       ICPWR   + NTBCPR*NCELET
      IWORR   =       IMCWR   + NTBMCR*NCELET
      IFINRA =        IWORR   + NTBWOR*NCELET
C
C---> VERIFICATION
C
      CALL IASIZE('MEMPH1',IFINIA)
C     ===========
C
      CALL RASIZE('MEMPH1',IFINRA)
C     ===========
C
      RETURN
      END
c@z
