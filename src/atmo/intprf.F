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
                        SUBROUTINE INTPRF
C                       *****************
C     ------------------------------------------------------
     & ( NPROFZ , NPROFT ,
     &   PROFZ  , PROFT  , PROFV  , XZ     , TEMPS  , VAR    )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC       INTERPOLATION VERTICALE-TEMPORELLE OPTIMISEE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! VI           ! R  ! <-  ! RESULTAT DE L'INTERPOLATION          !
CARGU ! ZI           ! R  !  -> ! ALTITUDE A LAQUELLE ON INTERPOLE     !
CARGU ! TI           ! R  !  -> ! TEMPS PHYSIQUE AUQUEL ON INTERPOLE   !
CARGU ! APROM        ! TR !  -> ! VARIABLES DES PROFILS METEO A        !
CARGU !              !    !     !  INTERPOLER                          !
CARGU ! TMPROM       ! TR !  -> ! INSTANTS DES PROFILS METEO A         !
CARGU !              !    !     !  INTERPOLER                          !
CARGU ! ZDPROM       ! E  !  -> ! ALTITUDES DES PROFILS METEO A        !
CARGU !              !    !     !  INTERPOLER                          !
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
      INTEGER          NPROFZ, NPROFV, NPROFT
      DOUBLE PRECISION PROFZ(NPROFZ), PROFT(NPROFT)
      DOUBLE PRECISION PROFV(NPROFZ,NPROFT)
      DOUBLE PRECISION XZ, TEMPS, VAR
C
C VARIABLES LOCALES
C
      INTEGER          IT, IT1, IT2
      INTEGER          IZ, IZ1, IZ2
      DOUBLE PRECISION ALPHAZ, ALPHAT, VAR1, VAR2
C
C***********************************************************************
C
C=======================================================================
C 1. INTERPOLATION EN TEMPS
C=======================================================================
C
      IF (TEMPS.LE.PROFT(1)) THEN
        IT1 = 1
        IT2 = 1
        ALPHAT = 1.D0
      ELSE IF (TEMPS.GE.PROFT(NPROFT)) THEN
        IT1 = NPROFT
        IT2 = NPROFT
        ALPHAT = 1.D0
C     sinon on a forcement NPROFT>1
      ELSE
        IT = 1
 102    CONTINUE
        IF (TEMPS.GT.PROFT(IT+1)) THEN
          IT = IT + 1
          GOTO 102
        ELSE
          IT1 = IT
          IT2 = IT+1
          ALPHAT = (PROFT(IT2)-TEMPS)/(PROFT(IT2)-PROFT(IT1))
        ENDIF
      ENDIF
C
C=======================================================================
C 2. INTERPOLATION VERTICALE EN ESPACE
C=======================================================================
C
      IF (XZ.LE.PROFZ(1)) THEN
        IZ1 = 1
        IZ2 = 1
        ALPHAZ = 1.D0
      ELSE IF (XZ.GE.PROFZ(NPROFZ)) THEN
        IZ1 = NPROFZ
        IZ2 = NPROFZ
        ALPHAZ = 1.D0
C     sinon on a forcement NPROFZ>1
      ELSE
        IZ = 1
 103    CONTINUE
        IF (XZ.GT.PROFZ(IZ+1)) THEN
          IZ = IZ + 1
          GOTO 103
        ELSE
          IZ1 = IZ
          IZ2 = IZ+1
          ALPHAZ = (PROFZ(IZ2)-XZ)/(PROFZ(IZ2)-PROFZ(IZ1))
        ENDIF
      ENDIF
C
C=======================================================================
C 3. INTERPOLATION
C=======================================================================
C
      VAR1 = ALPHAZ*PROFV(IZ1,IT1) + (1.D0-ALPHAZ)*PROFV(IZ2,IT1)
      VAR2 = ALPHAZ*PROFV(IZ1,IT2) + (1.D0-ALPHAZ)*PROFV(IZ2,IT2)
C
      VAR = ALPHAT*VAR1 + (1.D0-ALPHAT)*VAR2
C
      END

C


