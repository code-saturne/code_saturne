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
                      INTEGER FUNCTION IDRBLA
C                       *****************
C     -------------------------------------------------------------
     & ( CHAINE, LCH )
C     -------------------------------------------------------------
C**********************************************************************
C
C  FONCTION :
C  --------
c@foncb
CFONC
CFONC  DETERMINER LA POSITION DU DERNIER CARACTERE NON BLANC DANS
CFONC          CHAINE DE LONGUEUR LCH AVEC LA CONVENTION DE ZERO SI
CFONC          LA CHAINE EST BLANCHE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! CHAINE       ! E  !  -> ! CHAINE A VERIFIER                    !
CARGU ! LCH          ! E  !  -> ! LONGUEUR DE LA CHAINE                !
CARGU !______________!____!_____!______________________________________!
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
C ARGUMENTS
C
      CHARACTER CHAINE*(*)
      INTEGER   LCH
C
C VARIABLES LOCALES
C
      INTEGER   II
C
C***********************************************************************
C
C---------------
C POSITIONNEMENT
C---------------
C
      DO 10 II = LCH, 1, -1
         IF ( CHAINE (II:II) .NE. ' ' ) THEN
            IDRBLA = II
            GOTO 20
         ENDIF
   10 CONTINUE
      IDRBLA = 0
C
   20 CONTINUE
C
C-----
C FIN
C-----
C
      END
c@z
