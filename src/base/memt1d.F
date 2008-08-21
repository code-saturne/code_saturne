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
                        SUBROUTINE MEMT1D
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NFABOR , IFNIA1 , IFNRA1 , IFNIA2 , IFNRA2 ,
     &   IFINIA , IFINRA , IA     , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C
C  FONCTION
C  --------
c@foncb
CFONC
CFONC   GESTION MEMOIRE POUR LE MODULE THERMIQUE 1D EN PAROI
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
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! IFNIA1       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS IA APRES LIBERATION TOUT SAUF  !
CARGU !              !    !     !                               IFPT1D !
CARGU ! IFNRA1         E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS RA APRES LIBERATION TOUT SAUF  !
CARGU !              !    !     !                               TPPT1D !
CARGU ! IFNIA2       ! E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS IA APRES LIBERATION DE NPPT1D  !
CARGU ! IFNRA2         E  ! <-  ! POINTEUR DE LA PREMIERE CAS LIBRE    !
CARGU !              !    !     !  DANS RA APRES LIBERATION DE RGPT1D  !
CARGU !              !    !     !  ET EPPT1D                           !
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
      INCLUDE "cstnum.h"
      INCLUDE "optcal.h"
      INCLUDE "numvar.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "parall.h"
C***********************************************************************
C
C     ARGUMENTS
      INTEGER          IDBIA0 ,IDBRA0
      INTEGER          NFABOR
      INTEGER          IFINIA, IFINRA, IFNIA1, IFNRA1, IFNIA2, IFNRA2
      INTEGER          IA(*)
      DOUBLE PRECISION RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IOK1, IFAC

C
C***********************************************************************

C---> INITIALISATION
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0

C---> VERIFICATION DES DIMENSIONS
C
      IOK1 = 0
      IF(NFPT1D.GT.NFABOR .OR. NFPT1D.LT.0) THEN
        WRITE(NFECRA,1000) NFPT1D
        IOK1 = 1
      ENDIF
      IF(IOK1.NE.0) THEN
        CALL CSEXIT (1)
      ENDIF
C
C---> CALCUL DU NOMBRE DE FACES DE PAROI AVEC MODULE 1D TOTAL
C
      NFPT1T = NFPT1D
C
      IF (IRANGP.GE.0) THEN
        CALL PARCPT(NFPT1T)
      ENDIF
C
C---> QUELQUES MESSAGES
C
      IF(NFPT1T.EQ.0) THEN
        WRITE(NFECRA,2000) NFPT1T
      ELSE
        WRITE(NFECRA,2001) NFPT1T, NFPT1D
      ENDIF
      WRITE(NFECRA,3000)
C
C---> PLACE MEMOIRE RESERVEE AVEC DEFINITION DE IFINIA IFINRA
C
      IFINIA = IDEBIA
      IFINRA = IDEBRA

C
      IIFPT1 = IFINIA
      IFNIA1 = IIFPT1 + NFPT1D
      IICLT1 = IFNIA1
      IFNIA2 = IICLT1 + NFPT1D
      INPPT1 = IFNIA2
      IFINIA = INPPT1 + NFPT1D
C
      ITPPT1 = IFINRA
      IFNRA1 = ITPPT1 + NFPT1D
      ITEPT1 = IFNRA1
      IHEPT1 = ITEPT1 + NFPT1D
      IFEPT1 = IHEPT1 + NFPT1D
      IXLMT1 = IFEPT1 + NFPT1D
      IRCPT1 = IXLMT1 + NFPT1D
      IDTPT1 = IRCPT1 + NFPT1D
      IFNRA2 = IDTPT1 + NFPT1D
      IEPPT1 = IFNRA2
      IRGPT1 = IEPPT1 + NFPT1D
      IFINRA = IRGPT1 + NFPT1D
C
C---> VERIFICATION
C
      CALL IASIZE('MEMT1D',IFINIA)
C     ==========
C
      CALL RASIZE('MEMT1D',IFINRA)
C     ==========
C
C---> INITIALISATION DES TABLEAUX
C     a des valeurs sortant en erreur dans vert1d
C     sauf pour les variables de conditions aux limites
C     qui sont initialisees a des valeurs standard
C     (flux nul, Timpose=0, coef d'echange infini)
      DO IFAC = 1, NFPT1D
        IA(IIFPT1+IFAC-1) = -999
        IA(INPPT1+IFAC-1) = -999
        IA(IICLT1+IFAC-1) = 3
        RA(IEPPT1+IFAC-1) = -999.D0
        RA(IRGPT1+IFAC-1) = -999.D0
        RA(ITPPT1+IFAC-1) = 0.D0
        RA(ITEPT1+IFAC-1) = 0.D0
        RA(IHEPT1+IFAC-1) = RINFIN
        RA(IFEPT1+IFAC-1) = 0.D0
        RA(IXLMT1+IFAC-1) = -999.D0
        RA(IRCPT1+IFAC-1) = -999.D0
        RA(IDTPT1+IFAC-1) = -999.D0
      ENDDO

C---> FORMATS
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT(/,' SORTIE DANS MEMT1D CAR LA DIMENSIONS DU TABLEAU ',/,
     &         '   RELATIF AUX FACES DE BORD COUPLEES AU MODULE ',/,
     &         '   THERMIQUE 1D DE PAROI EST INCORRECTE ',/,
     &         '   NFPT1D = ',I10)
C
 2000 FORMAT(
     &    /,'TTES PHASES  : MODULE THERMIQUE 1D DE PAROI NON ACTIVE ',/,
     &      '                 NFPT1D = ',I10,/)
 2001 FORMAT(
     &    /,'TTES PHASES  : MODULE THERMIQUE 1D EN PAROI ACTIVE     ',/,
     &      '   SUR UN TOTAL DE ',I10,' FACES DE BORD',/,
     &      '   (',I10,' FACES DE BORD EN LOCAL)',/)
C
 3000 FORMAT(
     &'-------------------------------------------------------------',/)
C
#else
C
 1000 FORMAT(/,' ABORT IN MEMT1D BECAUSE THE DIMENSION OF THE ARRAY ',/,
     &         '   RELATIVE TO THE COUPLED FACES OF THE 1D-WALL ',/,
     &         '   THERMAL MODULE IS INCORRECT ',/,
     &         '   NFPT1D = ',I10)
C
 2000 FORMAT(
     & /,'ALL PHASES  : 1D-WALL THERMAL MODULE NOT ACTIVATED ',/,
     &   '                 NFPT1D = ',I10,/)
 2001 FORMAT(
     & /,'ALL PHASES  : 1D-WALL THERMAL MODULE ACTIVATED ',/,
     &   '   ON A TOTAL OF ',I10,' BOUNDARY FACES',/,
     &   '   (',I10,' LOCAL BOUNDARY FACES)',/)
C
 3000 FORMAT(
     &'-------------------------------------------------------------',/)
C
#endif
C
      RETURN
      END
c@z
