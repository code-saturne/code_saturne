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
                        SUBROUTINE VERLON
C                       *****************
C     -------------------------------------------------------------
     & ( CHAINE, II1, II2, LPOS )
C     -------------------------------------------------------------
C**********************************************************************
C  FONCTION :
C  --------
c@foncb
CFONC
CFONC   VERIFICATION DE LA LONGUEUR D'UNE CHAINE DE CARACTERES
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! CHAINE       ! A  !  -> ! CHAINE DE CARACTERES A VERIFIER      !
CARGU ! II1          ! E  ! <-  ! POSITION  PREMIER CARACTERE NON BLANC!
CARGU ! II2          ! E  ! <-  ! POSITION  DERNIER CARACTERE NON BLANC!
CARGU ! LPOS         ! E  ! <-  ! LONGUEUR EFFECTIVE DE LA CHAINE      !
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
C
C***********************************************************************
C
      IMPLICIT NONE
C
C ARGUMENTS
C
      CHARACTER CHAINE*(*)
      INTEGER   II1,II2,LPOS
C
C VARIABLES LOCALES
C
      INTEGER   N1,IPRBLA,IDRBLA
C
C***********************************************************************
C
      II1  = 0
      II2  = 0
      LPOS = 0
      N1   = LEN ( CHAINE )
      IF ( N1 .LE. 0 ) RETURN
C
      II1  = IPRBLA ( CHAINE, N1 )
      II2  = IDRBLA ( CHAINE, N1 )
      LPOS = II2 - II1 + 1
C
      END
c@z
