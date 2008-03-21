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
                SUBROUTINE LWCGFU
C               *****************
C     --------------------------------------------------------
     & ( GFUNC , F    , FM    , YFP2M  , FP2M )
C     --------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL DES VALEURS DE LA FONCTION G
CFONC  SUIVANT LES PARAMETRES F, FM, FP2M. YP2M
CFONC
CFONC  LE RESULTAT EST :
CFONC  ---------------
CFONC     CALCUL DE LA VALEUR DE G AU POINT F
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU !              !    !     !                                      !
CARGU ! GFUNC        ! R  ! <-  ! Valeur de G au Point F               !
CARGU ! F            ! R  !  -> ! VALEUR DE LA FRACTION DE MELANGE     !
CARGU ! FM           ! R  !  -> ! MOYENNE DE LA FRACTION DE MELANGE    !
CARGU ! FP2M         ! R  !  -> ! VARIANCE DE LA FRACTION DE MELANGE   !
CARGU ! YP2M         ! R  !  -> ! VARIANCE DE LA FRACTION MASSIQUE     !
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
C ARGUMENTS
C
      DOUBLE PRECISION GFUNC , F, FM, YFP2M, FP2M

C
C VARIABLES LOCALES
C
      DOUBLE PRECISION EPSI
C
C
C***********************************************************************
C
C=======================================================================
C 1.  CALCULS PRELIMINAIRES
C=======================================================================
C
C --> Initialisation
C
      GFUNC = 0.D0
      EPSI = 1.D-09
C
C=======================================================================
C 2.  CALCUL DES VALEURS DE LA FONCTION G
C=======================================================================
C
      IF (FP2M .LE. EPSI) THEN
        GFUNC = 1.D0
      ELSE
        GFUNC = (F-FM) * SQRT( 1.D0 + YFP2M/FP2M )
      ENDIF
C
C----
C FIN
C----
C
      END
c@z
