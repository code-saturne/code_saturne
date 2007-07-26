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
                        SUBROUTINE LAGTRI
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC  La liste des particules a visualiser est donnee par l'utilisateur
CFONC  dans USLAG1. Cette subroutine verifie cette liste :
CFONC  1. elle suprime les trous dans le tableau,
CFONC  2. suprime les doublons,
CFONC  3. range les particules a visualiser par ordre croissant.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU !              !    !     !                                      !
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
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER I , I1 , I2 , II , IND , J , TLIST(NLISTE)
C
C***********************************************************************
C
C... VERIFICATION DES NUMEROS DES PARTICULES DESIREES
C
      IND = 0
      DO I = 1,NLISTE
        TLIST (I) = -1
        IF ( LISTE (I).GT.0 ) THEN
          IND = IND + 1
          TLIST (IND) = LISTE (I)
        ENDIF
      ENDDO
C
C... RANGEMENT DU TABLEAU LISTE SANS TROU
C
      DO I = 1,NLISTE
        LISTE (I) = TLIST (I)
      ENDDO
C
C... RANGEMENT DU TABLEAU LISTE SANS REPETITION
C
      I1 = 0
      DO I = 1,IND
        II = LISTE (I)
        IF ( II .GT. 0 ) THEN
          DO J = I+1, IND
            IF ( II.EQ.LISTE (J) ) THEN
              LISTE (J) = -1
              I1 = I1 + 1
            ENDIF
          ENDDO
        ENDIF
      ENDDO
C
      IF ( I1.GT.0 ) THEN
        I2 = 0
        DO I = 1,IND
          TLIST (I) = -1
          IF ( LISTE (I).GT.0 ) THEN
            I2 = I2 + 1
            TLIST (I2) = LISTE (I)
          ENDIF
        ENDDO
        IND = I2
        DO I = 1,IND
          LISTE (I) = TLIST (I)
        ENDDO
        DO I = IND+1,NLISTE
          LISTE (I) = -1
        ENDDO
      ENDIF
C
C... RANGEMENT DU TABLEAU LISTE PAR ORDRE CROISSANT
C
      I = 1
  50  CONTINUE
      I2 = I + 1
      IF (  LISTE (I2).GT.0  .AND.
     &      LISTE (I2).LT.LISTE (I) ) THEN
        I1 = LISTE (I)
        LISTE (I) = LISTE (I2)
        LISTE (I2) = I1
        I = 1
        GOTO 50
      ELSE
        I = I2
        IF ( I.LT.IND ) GOTO 50
      ENDIF
C
      IF ( IND.GT.NBVIS ) THEN
        DO I = IND+1,NLISTE
          LISTE (I) = -1
        ENDDO
      ENDIF
C
C----
C FIN
C----
C
      END
c@z
