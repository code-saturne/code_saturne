c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE USSYRC
C                       *****************
C     -------------------------------------------------------------
     & ( )
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
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C Arguments
C
C Variables locales
C
      CHARACTER        CPRJSY
      INTEGER          NUMSYR, NBCSYR, II
      INTEGER          IBOSYR, IVLSYR, IWARNS
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C
      IF(1.EQ.1) RETURN
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
      IWARNS = 1
C
      NUMSYR = -1
      NBCSYR = 3
C
      DO II = 1, NBCSYR
C
C       IPRJSY : ' ' : couplage 3D standard
C                'x', 'y', ou 'z' : axe de projection pour couplage avec
C                                   SYRTHES 2D.
C       IWARNS : niveau d'impression associe
C
C       Exemple : couplage surfacique 3D aux les faces de couleur 3
C                 avec une instance nommee SYRTHES_01
C
        IF (II .EQ. 1) THEN
C
          CPRJSY= ' '
          IBOSYR = 1
          IVLSYR = 0
C
          CALL DEFSYR(NUMSYR, 'SYRTHES_01', CPRJSY, '3', ' ', IWARNS)
C         ===========

C       Exemple : couplage surfacique 2D aux les faces du groupe 'wall'
C                 avec une instance nommee SYRTHES_02
C
        ELSE IF (II .EQ. 2) THEN
C
          CPRJSY= 'z'
          IBOSYR = 1
          IVLSYR = 0
C
          CALL DEFSYR(NUMSYR, 'SYRTHES_02', CPRJSY,
C         ===========
     &                'wall', ' ', IWARNS)
C
C       Exemple : couplage volumique 3D au niveau de la boite
C                 de coins (0, 0, 0) et (1, 1, 1), avec une instance
C                 nommee SOLID

        ELSE IF (II .EQ. 3) THEN
C
          CPRJSY= ' '
          IBOSYR = 1
          IVLSYR = 0
C
          CALL DEFSYR(NUMSYR, 'SOLID', CPRJSY,
C         ===========
     &                ' ', 'box[0., 0., 0., 1., 1., 1.]', IWARNS)
C
C
        ENDIF
C
      ENDDO
C
      RETURN
C
      END
c@z

