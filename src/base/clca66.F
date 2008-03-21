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
                        SUBROUTINE CLCA66
C                       *****************
C     -------------------------------------------------------------
     & ( CLSYME , ELOGLO , ALPHA )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC       CALCUL DE LA MATRICE ALPHA (VOIR CNDRIJ)
CFONC POUR LES CONDITIONS AUX LIMITES EN PAROI ET SYMETRIE EN RIJ-EPS
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! CLSYME       !  R !  -> ! INDICATEUR SYMETRIE = 1., PAROI = 0. !
CARGU ! ELOGLO       ! TR !  -> ! MATRICE DE CHANGEMENT DE BASE        !
CARGU !  (6,6)       !    !     !  LOCAL -> GLOBAL                     !
CARGU ! ALPHA        ! TR ! <-  ! VOIR CNDRIJ                          !
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
C ARGUMENTS
C
      DOUBLE PRECISION CLSYME
      DOUBLE PRECISION ELOGLO(3,3),ALPHA(6,6)
C
C VARIABLES LOCALES
C
      INTEGER II,JJ,KK,PP,JJ1,JJ2
C
C***********************************************************************
C
C=======================================================================
C CALCUL DE ALPHA(I,J) POUR I DANS [1,3] ET J DANS [1,3]  : 9 TERMES
C=======================================================================
C
       DO II = 1,3
         DO JJ = 1,3
C
           ALPHA(II,JJ) =   ELOGLO(II,1)**2*ELOGLO(JJ,1)**2 +
     &                      ELOGLO(II,2)**2*ELOGLO(JJ,2)**2 +
     &                      ELOGLO(II,3)**2*ELOGLO(JJ,3)**2 +
     &   2.D0*CLSYME*ELOGLO(II,1)*ELOGLO(II,3)*ELOGLO(JJ,1)*ELOGLO(JJ,3)
C
         ENDDO
       ENDDO
C
C=======================================================================
C CALCUL DE ALPHA(I,J) POUR I DANS [1,3] ET J DANS [4,6]  : 9 TERMES
C=======================================================================
C
       DO II=1,3
         DO JJ=1,3
C
           IF (JJ.EQ.1) THEN
             KK=1
             PP=2
           ELSE IF (JJ.EQ.2) THEN
             KK=1
             PP=3
           ELSE IF (JJ.EQ.3) THEN
             KK=2
             PP=3
           ENDIF
C
           ALPHA(II,JJ+3) = 2.D0*(
     &       ELOGLO(II,1)**2*ELOGLO(KK,1)*ELOGLO(PP,1) +
     &       ELOGLO(II,2)**2*ELOGLO(KK,2)*ELOGLO(PP,2) +
     &       ELOGLO(II,3)**2*ELOGLO(KK,3)*ELOGLO(PP,3) +
     &       CLSYME*ELOGLO(II,3)*ELOGLO(II,1)*
     &        (ELOGLO(KK,1)*ELOGLO(PP,3)+ELOGLO(KK,3)*ELOGLO(PP,1)))
C
         ENDDO
       ENDDO
C
C=======================================================================
C CALCUL DE ALPHA(I,J) POUR I DANS [4,6] ET J DANS [1,3]  : 9 TERMES
C=======================================================================
C
       DO II=1,3
         DO JJ=1,3
C
           IF (II.EQ.1) THEN
             KK=1
             PP=2
           ELSE IF (II.EQ.2) THEN
             KK=1
             PP=3
           ELSE IF (II.EQ.3) THEN
             KK=2
             PP=3
           ENDIF
           ALPHA(II+3,JJ) =
     &       ELOGLO(KK,1)*ELOGLO(PP,1)*ELOGLO(JJ,1)**2 +
     &       ELOGLO(KK,2)*ELOGLO(PP,2)*ELOGLO(JJ,2)**2 +
     &       ELOGLO(KK,3)*ELOGLO(PP,3)*ELOGLO(JJ,3)**2 +
     &       CLSYME*ELOGLO(JJ,3)*ELOGLO(JJ,1)*
     &        (ELOGLO(KK,1)*ELOGLO(PP,3)+ELOGLO(KK,3)*ELOGLO(PP,1))
C
         ENDDO
       ENDDO
C
C
C=======================================================================
C CALCUL DE ALPHA(I,J) POUR I DANS [4,6] ET J DANS [4,6]  : 9 TERMES
C=======================================================================
C
       DO II=1,3
         DO JJ=1,3
C
           IF (II.EQ.1) THEN
             KK=1
             PP=2
           ELSE IF (II.EQ.2) THEN
             KK=1
             PP=3
           ELSE IF (II.EQ.3) THEN
             KK=2
             PP=3
           ENDIF
           IF (JJ.EQ.1) THEN
             JJ1=1
             JJ2=2
           ELSE IF (JJ.EQ.2) THEN
             JJ1=1
             JJ2=3
           ELSE IF (JJ.EQ.3) THEN
             JJ1=2
             JJ2=3
           ENDIF
           ALPHA(II+3,JJ+3)=
     &     2.D0*(ELOGLO(KK,1)*ELOGLO(PP,1)*ELOGLO(JJ1,1)*ELOGLO(JJ2,1)+
     &           ELOGLO(KK,2)*ELOGLO(PP,2)*ELOGLO(JJ1,2)*ELOGLO(JJ2,2)+
     &           ELOGLO(KK,3)*ELOGLO(PP,3)*ELOGLO(JJ1,3)*ELOGLO(JJ2,3))+
     &     CLSYME*(ELOGLO(KK,1)*ELOGLO(PP,3)+ELOGLO(KK,3)*ELOGLO(PP,1))*
     &         (ELOGLO(JJ1,3)*ELOGLO(JJ2,1)+ELOGLO(JJ1,1)*ELOGLO(JJ2,3))
C
         ENDDO
       ENDDO
C
       RETURN
       END
c@z
