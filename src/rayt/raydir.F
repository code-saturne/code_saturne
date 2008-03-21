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
                        SUBROUTINE RAYDIR
C                       *****************
C     -------------------------------------------------------------
     &  ( SX,SY,SZ,NDIRS )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
C
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE RAYONNEMENT :
CFONC    --------------------------------------
CFONC
CFONC   CALCUL DES COSINUS DIRECTEURS DE LA DIRECTION
CFONC   DE PROPAGATION DU RAYONNEMENT
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU !SX,SY,SZ      ! R  ! <-  ! COSINUS DIRECTEURS DU RAYONNEMENT    !
CARGU !  NDIRS       ! E  ! ->  ! NOMBRE DE DIRECTIONS PAR 1/8 DE      !
CARGU !              !    !     ! SPHERE                               !
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
      INTEGER          NDIRS
      DOUBLE PRECISION SX(NDIRS),SY(NDIRS),SZ(NDIRS)
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION  TETA(6,3)
      DOUBLE PRECISION  PHI1,TETA1
      INTEGER           II,KK,NRAY,JRAY,IIA,IIB
C
C***********************************************************************
C
      DATA TETA
     & / .785398D0,  .785398D0,  .785398D0, 1.256637D0,  .300162D0,
     &                                                       .448799D0,
     &   .785398D0, 1.256637D0, 1.427996D0,  .999598D0, 1.374447D0,
     &                                                      1.121997D0,
     &   .785398D0,  .314159D0,  .142800D0,  .205855D0,  .571199D0,
     &                                                       .785398D0 /
C
C***********************************************************************
C
C    0 - INITIALISATION
C        --------------
C  FORMULE : NDIRS = 3*NRAY -2
C
      IF (NDIRS.EQ.4)  NRAY = 2
      IF (NDIRS.EQ.16) NRAY = 6
C
C    1 - CALCUL  DES COSINUS DIRECTEURS
C        ------------------------------
C
      SZ(1 ) = COS(ATAN(TAN(TETA(1,3))/COS(TETA(1,1))))
      SY(1 ) = SZ(1 )
      SX(1 ) = SZ(1 )
C
      JRAY = 1
C
      DO II = 0,2
C
        IIA = 3 + II
        IIB = 1 + II
        IF (IIA.GT.3) IIA = IIA-3
C
        DO KK = 2,NRAY
          JRAY = JRAY+1
          TETA1 = TETA(KK,IIB)
          PHI1 = ATAN(TAN(TETA(KK,IIA))/COS(TETA(KK,IIB)))
          SX(JRAY) = SIN(PHI1)*COS(TETA1)
          SY(JRAY) = SIN(PHI1)*SIN(TETA1)
          SZ(JRAY) = COS(PHI1)
        ENDDO
C
      ENDDO
C
      RETURN
C
      END
c@z
