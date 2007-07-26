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
                        SUBROUTINE MEMPDC
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NCELET , NCEL   , NPHAS  , NDIM   , IFINIA , IFINRA )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE TABLEAU PDC : CKUPDC
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
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE (3)            !
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
C***********************************************************************
C
C
      INTEGER          IDBIA0 ,IDBRA0
      INTEGER          NCELET , NCEL , NPHAS, NDIM
      INTEGER          IFINIA , IFINRA
C
      INTEGER          IDEBIA, IDEBRA, IOK, IOK1, IPHAS
C
C***********************************************************************
C
C---> INITIALISATION
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C
C---> VERIFICATION DES DIMENSIONS
C
      IOK1 = 0
      DO IPHAS = 1, NPHAS
C
        IOK = 0
        IF(NCKPDC(IPHAS).NE.3 .AND. NCKPDC(IPHAS).NE.6) THEN
          IOK = 1
        ENDIF
        IF(NCEPDC(IPHAS).GT.NCELET .OR. NCEPDC(IPHAS).LT.0) THEN
          IOK = 1
        ENDIF
        IF(IOK.NE.0) THEN
          WRITE(NFECRA,1000) IPHAS, NCEPDC(IPHAS), NCKPDC(IPHAS)
          IOK1 = 1
        ENDIF
      ENDDO
C
      IF(IOK1.NE.0) THEN
        CALL CSEXIT (1)
      ENDIF
C
C---> CALCUL DU NOMBRE DE CELLULES AVEC PDC TOTAL
C
      DO IPHAS = 1, NPHAS
        NCPDCT(IPHAS) = NCEPDC(IPHAS)
      ENDDO
      IF (IRANGP.GE.0) THEN
        CALL PARISM(NPHAS,NCPDCT)
      ENDIF
C
C---> QUELQUES MESSAGES
C
      DO IPHAS = 1, NPHAS
        IF(NCPDCT(IPHAS).EQ.0) THEN
          WRITE(NFECRA,2000) IPHAS, NCPDCT(IPHAS)
          WRITE(NFECRA,3000)
        ELSE
          IF(NCKPDC(IPHAS).EQ.3) THEN
            WRITE(NFECRA,2001) IPHAS, NCPDCT(IPHAS), NCKPDC(IPHAS)
            WRITE(NFECRA,3000)
          ELSEIF(NCKPDC(IPHAS).EQ.6) THEN
            WRITE(NFECRA,2002) IPHAS, NCPDCT(IPHAS), NCKPDC(IPHAS)
            WRITE(NFECRA,3000)
          ENDIF
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
        IICEPD(IPHAS) = IFINIA
        IFINIA        = IICEPD(IPHAS) + NCEPDC(IPHAS)
C
        ICKUPD(IPHAS) = IFINRA
        IFINRA        = ICKUPD(IPHAS) + NCEPDC(IPHAS)*NCKPDC(IPHAS)
C
      ENDDO
C
C     Si pour une des phases on a des pertes de charge
C       sur un des processeurs
C     et que les tableaux TPUCOU n'ont pas encore ete definis (IPUCOU=0)
C     il faut les dimensionner
      IF (IPUCOU.EQ.0) THEN
        IOK=0
        DO IPHAS = 1, NPHAS
          IF (NCPDCT(IPHAS).GT.0) IOK = 1
        ENDDO
        IF (IOK.EQ.1) THEN
          ITPUCO = IFINRA
          IFINRA = ITPUCO + NCELET *NDIM
        ENDIF
      ENDIF


C
C---> VERIFICATION
C
      CALL IASIZE('MEMPDC',IFINIA)
C     ==========
C
      CALL RASIZE('MEMPDC',IFINRA)
C     ==========
C
C---> FORMATS
C
 1000 FORMAT(/,' SORTIE DANS MEMPDC CAR LES DIMENSIONS DES TABLEAUX ',/,
     &         '   RELATIFS AUX PERTES DE CHARGES SONT INCORRECTES ',/,
     &         '   PHASE ',I10,/,
     &         '     NCEPDC = ',I10,/,
     &         '     NCKPDC = ',I10)
C
 2000 FORMAT(/,
     & 'PHASE ',I6,' : TRAITEMENT DES PERTES DE CHARGES NON ACTIVE ',/,
     & '                 NCEPDC = ',I10,/)
 2001 FORMAT(
     & /,/,
     & 'PHASE ',I6,' : TRAITEMENT DES PERTES DE CHARGES ACTIVE ',/,
     & '                 SUR  UN TOTAL DE NCEPDC = ',I10,' CELLULES',/,
     & '                 AVEC UN TENSEUR DIAGONAL (NCKPDC = ',I1,')',/)
 2002 FORMAT(
     & /,/,
     & 'PHASE ',I6,' : TRAITEMENT DES PERTES DE CHARGES ACTIVE ',/,
     & '                 SUR  UN TOTAL DE NCEPDC = ',I10,' CELLULES',/,
     & '                 AVEC UN TENSEUR COMPLET SYMETRIQUE (NCKPDC = ',
     &  I1,')',/)
C
 3000 FORMAT(
     &'-------------------------------------------------------------',/)
C
      RETURN
      END
c@z
