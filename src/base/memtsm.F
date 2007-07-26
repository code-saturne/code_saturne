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
                        SUBROUTINE MEMTSM
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCELET , NCEL   , NVAR   , NPHAS  ,
     &   IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE POUR LES TERMES SOURCES DE MASSE
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
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! IFINIA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS IA EN SORTIE                   !
CARGU ! IFINRA       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS IA EN SORTIE                   !
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
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C
      INTEGER          IDBIA0 ,IDBRA0
      INTEGER          NCELET , NCEL  ,NVAR   , NPHAS
      INTEGER          IFINIA , IFINRA
C
      INTEGER          IDEBIA, IDEBRA, IOK1, IPHAS
      INTEGER          IPASS
      DATA             IPASS /0/
      SAVE             IPASS
C
C***********************************************************************
C
C---> INITIALISATION
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      IPASS = IPASS + 1
C
C---> VERIFICATION DES DIMENSIONS
C
      IOK1 = 0
      DO IPHAS = 1, NPHAS
        IF(NCETSM(IPHAS).GT.NCELET .OR. NCETSM(IPHAS).LT.0) THEN
          WRITE(NFECRA,1000) IPHAS, NCETSM(IPHAS)
          IOK1 = 1
        ENDIF
      ENDDO
      IF(IOK1.NE.0) THEN
        CALL CSEXIT (1)
      ENDIF
C
C---> CALCUL DU NOMBRE DE CELLULES AVEC TSM TOTAL
C
      DO IPHAS = 1, NPHAS
        NCTSMT(IPHAS) = NCETSM(IPHAS)
      ENDDO
      IF (IRANGP.GE.0) THEN
        CALL PARISM(NPHAS,NCTSMT)
      ENDIF
C
C---> QUELQUES MESSAGES
C
      DO IPHAS = 1, NPHAS
        IF(NCTSMT(IPHAS).EQ.0) THEN
          WRITE(NFECRA,2000) IPHAS, NCTSMT(IPHAS)
          WRITE(NFECRA,3000)
        ELSE
          WRITE(NFECRA,2001) IPHAS, NCTSMT(IPHAS)
          WRITE(NFECRA,3000)
        ENDIF
      ENDDO
C
C---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
C
      IFINIA = IDEBIA
      IFINRA = IDEBRA
C
      DO IPHAS = 1, NPHAS
C
        IICESM(IPHAS) = IFINIA
        IITPSM(IPHAS) = IICESM(IPHAS) + NCETSM(IPHAS)
        IFINIA        = IITPSM(IPHAS) + NCETSM(IPHAS)*NVAR
C
        ISMACE(IPHAS) = IFINRA
        IFINRA        = ISMACE(IPHAS) + NCETSM(IPHAS)*NVAR
C
      ENDDO
C
C---> VERIFICATION
C
      CALL IASIZE('MEMTSM',IFINIA)
C     ==========
C
      CALL RASIZE('MEMTSM',IFINRA)
C     ==========
C
C---> FORMATS
C
 1000 FORMAT(/,' SORTIE DANS MEMTSM CAR LA DIMENSIONS DU TABLEAU ',/,
     &         '   RELATIF AUX SOURCES DE MASSE EST INCORRECTE ',/,
     &         '   PHASE ',I10,/,
     &         '   NCETSM = ',I10)
C
 2000 FORMAT(
     & /,'PHASE ',I6,' : TRAITEMENT DES SOURCES DE MASSE NON ACTIVE ',/,
     &   '                 NCETSM = ',I10,/)
 2001 FORMAT(
     & /,'PHASE ',I6,' : TRAITEMENT DES SOURCES DE MASSE ACTIVE ',/,
     &   '                 SUR  UN TOTAL DE ',I10,' CELLULES')
C
 3000 FORMAT(
     &'-------------------------------------------------------------',/)
C
      RETURN
      END
c@z
