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
                        SUBROUTINE ORDITA
C                       *****************
C     -------------------------------------------------------------
     & ( NFABOR , IFABOR, ICLASS)
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC ORDONNE LES NFABOR PAR IFABOR DECROISSANT
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                 ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NFABOR       ! E  !  -> ! NB D'EL A ORDONNER                   !
CARGU ! IFABOR(NFABOR! TE !  -> ! CRITERE D'ORDONNANCEMENT             !
CARGU ! ICLASS(NFABOR! TE ! <-  ! TABLE DE RENUM                       !
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
C     DONNEES EN COMMON
C***********************************************************************
C
C
      INTEGER          NFABOR
      INTEGER          IFABOR(NFABOR)
      INTEGER          ICLASS(NFABOR)
C
      INTEGER          JJ, ITMP, II
C
C***********************************************************************
C
      IF (NFABOR.EQ.0) RETURN
C
      DO JJ = 1, NFABOR
        ICLASS(JJ) = JJ
      ENDDO
      DO JJ = NFABOR/2, 1, -1
        CALL TDESI1(JJ,NFABOR,NFABOR,IFABOR,ICLASS)
      ENDDO
      DO II = NFABOR, 3, -1
        ITMP       = ICLASS(1)
        ICLASS(1)  = ICLASS(II)
        ICLASS(II) = ITMP
        CALL TDESI1(1,NFABOR,II-1,IFABOR,ICLASS)
      ENDDO
      ITMP      = ICLASS(1)
      ICLASS(1) = ICLASS(2)
      ICLASS(2) = ITMP
C
      RETURN
      END
c@z
