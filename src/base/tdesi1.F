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
                        SUBROUTINE TDESI1
C                       *****************
C     -------------------------------------------------------------
     & ( JJ     , NFABOR , NN     , IFABOR , ICLASS )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION
C  --------
c@foncb
CFONC
CFONC               SOUS-PROGRAMME DE DESCENTE D'UN ARBRE BINAIRE POUR
CFONC                L'ALGORITHME DE TRI PAR ARBRE (HEAPSORT)
CFONC
CFONC               VERSION 2D
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! JJ           ! E  !  -> ! NIVEAU DE L'ARBRE BINAIRE TRAITE     !
CARGU ! NN           ! E  !  -> ! NOMBRE DE FACES    ACTIVES           !
CARGU ! NFABOR       ! E  !  -> ! TAILLE DU TABLEAU A DESCENDRE        !
CARGU ! ICLASS       ! E  ! <-> ! POINTEURS SUR LES NUMEROS D'ELEMENTS !
CARGU !              !    !     ! DANS L'ARBRE BINAIRE                 !
CARGU ! IFABOR       ! TR !  -> ! VOISINS DES FACES                    !
CARGU !  (  NFABOR)  !    !     !                                      !
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
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          JJ, NN, NFABOR
C
      INTEGER          ICLASS(NFABOR)
      INTEGER          IFABOR(NFABOR)
C
C VARAIBLES LOCALES
C
      INTEGER      II, LL, ISVNUM
      LOGICAL      DESC, ISUP
C
C=======================================================================
C
      II = JJ
      LL = 2*JJ
C
      DESC = .TRUE.
C
      IF ((LL.GT.NN).OR.(.NOT. DESC)) GOTO 20
C
 10   CONTINUE
C
         IF (LL .LT. NN) THEN
C
C--->       ISUP = .TRUE. SI ICLASS (LL) > ICLASS (LL+1)
C
            ISUP = .TRUE.
            IF ( IFABOR (ICLASS (LL+1)) .GT.
     &                                     IFABOR (ICLASS (LL)) ) THEN
               ISUP = .FALSE.
            ENDIF
C
            IF (ISUP) LL = LL + 1
C
         ENDIF
C
C--->    ISUP = .TRUE. SI ICLASS (II) > ICLASS (LL)
C
         ISUP = .TRUE.
         IF ( IFABOR (ICLASS (LL)) .GT.
     &                                    IFABOR (ICLASS (II)) ) THEN
            ISUP = .FALSE.
         ENDIF
C
         IF (ISUP) THEN
            ISVNUM = ICLASS (II)
            ICLASS (II) = ICLASS (LL)
            ICLASS (LL) = ISVNUM
            II = LL
            LL = 2*LL
            DESC = .TRUE.
         ELSE
            DESC = .FALSE.
         ENDIF
C
         IF ((LL .LE. NN) .AND. DESC) GOTO 10
C
 20   CONTINUE
C
      RETURN
      END
c@z
