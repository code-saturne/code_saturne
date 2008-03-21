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
                        SUBROUTINE GETFBR
C                       ******************
C     -------------------------------------------------------------
     & ( FSTR , FACNB, FACES)
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CONSTRUCTION DE LA LISTE DES FACES DE BORD CORRESPONDANT A
CFONC A L'INTERPRETATION DE LA CHAINE DE CARACTERE FSTR

c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! FSTR         ! A  !  -> ! CHAINE A INTERPRETER                 !
CARGU ! FACES        ! E  !  <- ! FACES SELECTIONNEES                  !
CARGU ! FACNB        ! E  !  <- ! NOMBRE DE FACES SELECTIONNEES        !
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
C***********************************************************************
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*(*)    FSTR
      INTEGER      FACES(*), FACNB
C
C VARIABLES LOCALES
C
      INTEGER          LENSTR
C
C***********************************************************************
C
      LENSTR=LEN(FSTR)
      CALL CSGFBR(FSTR, LENSTR, FACNB, FACES)
C
      RETURN
C
      END
c@z
c@a
                        SUBROUTINE GETCEL
C                       ******************
C     -------------------------------------------------------------
     & ( FSTR , CELLNB, CELLS)
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CONSTRUCTION DE LA LISTE DES CELLULES CORRESPONDANT A
CFONC A L'INTERPRETATION DE LA CHAINE DE CARACTERE FSTR

c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! FSTR         ! A  !  -> ! CHAINE A INTERPRETER                 !
CARGU ! CELLS        ! E  !  <- ! FACES SELECTIONNEES                  !
CARGU ! CELLNB       ! E  !  <- ! NOMBRE DE FACES SELECTIONNEES        !
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
C***********************************************************************
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*(*) FSTR
      INTEGER       CELLS(*), CELLNB
C
C VARIABLES LOCALES
C
      INTEGER          LENSTR
C
C***********************************************************************
C
      LENSTR=LEN(FSTR)
      CALL CSGCEL(FSTR, LENSTR, CELLNB, CELLS)
C
      RETURN
C
      END
c@z
c@a
                        SUBROUTINE GETFAC
C                       ******************
C     -------------------------------------------------------------
     & ( FSTR , FACNB, FACES)
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CONSTRUCTION DE LA LISTE DES FACES INTERNES CORRESPONDANT A
CFONC A L'INTERPRETATION DE LA CHAINE DE CARACTERE FSTR

c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! FSTR         ! A  !  -> ! CHAINE A INTERPRETER                 !
CARGU ! CELLS        ! E  !  <- ! FACES SELECTIONNEES                  !
CARGU ! CELLNB       ! E  !  <- ! NOMBRE DE FACES SELECTIONNEES        !
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
C***********************************************************************
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C***********************************************************************
C
C ARGUMENTS
C
      CHARACTER*(*) FSTR
      INTEGER       FACES(*), FACNB
C
C VARIABLES LOCALES
C
      INTEGER          LENSTR
C
C***********************************************************************
C
      LENSTR=LEN(FSTR)
      CALL CSGFAC(FSTR, LENSTR, FACNB, FACES)
C
      RETURN
C
      END
