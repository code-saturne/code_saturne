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
                        SUBROUTINE DEFSYR
C                       *****************
C     -------------------------------------------------------------
     & ( NUMSYR, NOMSYR , CPROJ, CRITBO, CRITVL, IWARNI )
C
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC      DEFINITION DE COUPLAGE(S) AVEC LE CODE SYRTHES
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NUMSYR       ! A  !  -> ! NUMERO DU CAS SYRTHES ASSOCIE        !
CARGU ! NOMSYR       ! A  !  -> ! NOM DU CAS SYRTHES ASSOCIE           !
CARGU ! CPROJ        ! C  !  -> ! DIRECTION DE PROJECTION (' ' POUR    !
CARGU !              !    !     ! 3D, 'x', 'y', ou 'z' POUR 2D)        !
CARGU ! CRITBO       ! A  !  -> ! CRITERE DE SELECTION DES FACES DE    !
CARGU !              !    !     ! BORD COUPLEES (OU VIDE)              !
CARGU ! CRITVL       ! A  !  -> ! CRITERE DE SELECTION DES CELLULES    !
CARGU !              !    !     ! VOLUMIQUES COUPLEES (OU VIDE)        !
CARGU ! IWARNI       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
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
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
C***********************************************************************
C
C Arguments
C
      CHARACTER*(*) NOMSYR, CRITBO, CRITVL
      CHARACTER     CPROJ
      INTEGER       NUMSYR, IWARNI
C
C Variables locales
C
      INTEGER       LNOMSY, LCRITB, LCRITV
C
C***********************************************************************
C
      LNOMSY = LEN(NOMSYR)
      LCRITB = LEN(CRITBO)
      LCRITV = LEN(CRITVL)
C
      CALL DEFSY1(NUMSYR, NOMSYR, CPROJ,  CRITBO, CRITVL, IWARNI,
C     ===========
     &            LNOMSY, LCRITB, LCRITV)
C
C
      RETURN
C
      END
c@z

