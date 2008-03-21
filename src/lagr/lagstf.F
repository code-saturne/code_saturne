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
                        SUBROUTINE LAGSTF
C                       *****************
C     -----------------------------------------------------------------
     &  ( NCELET , NFABOR , NVISBR ,
     &    IVAR   ,
     &    GMIN   , GMAX   , GMOY   ,
     &    PARBOR , UNSNBR )
C     -----------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC      CALCUL DES VALEURS MIN MAX ET MOYENNE POUR LES
CFONC      STATISTIQUES AUX FRONTIERES.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NVISBR       ! E  !  -> ! NOMBRE DE STATISTIQUES AUX FRONTIERES!
CARGU ! IVAR         ! E  !  -> ! NUMERO DE LA VARIABLE A INTEGRER     !
CARGU ! GMIN...GMOY  ! E  ! <-  ! VARIATIONS MIN MAX ET MOYENNE        !
CARGU ! PARBOR       ! TR !  -> ! INFOS SUR INTERACTION DES PARTICULES !
CARGU !(NFABOR,NVISBR!    !     !   AUX FACES DE BORD                  !
CARGU ! UNSNBR(NFABOR! E  !  -> ! INVERSE DU NOMBRE INTERACTION        !
CARGU !              !    !     !   FRONTIERE POUR MOYENNE             !
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
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NFABOR , NVISBR
      INTEGER          IVAR
      DOUBLE PRECISION GMAX   , GMIN   , GMOY
      DOUBLE PRECISION PARBOR(NFABOR, NVISBR), UNSNBR(NFABOR)
C
C  VARIABLES LOCALES
C
      INTEGER          NBRFAC , IFAC
C
C**********************************************************************
C
C
      NBRFAC = 0
C
      GMAX = -GRAND
      GMIN =  GRAND
      GMOY =  0.D0
C
      IF (IMOYBR(IVAR).EQ.2) THEN
C
        DO IFAC = 1,NFABOR
          IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
            NBRFAC = NBRFAC + 1
            GMAX = MAX (GMAX, PARBOR(IFAC,IVAR)*UNSNBR(IFAC))
            GMIN = MIN (GMIN, PARBOR(IFAC,IVAR)*UNSNBR(IFAC))
            GMOY = GMOY + (PARBOR(IFAC,IVAR)*UNSNBR(IFAC))
          ENDIF
        ENDDO
C
      ELSE IF (IMOYBR(IVAR).EQ.1) THEN
C
        DO IFAC = 1,NFABOR
          IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
            NBRFAC = NBRFAC + 1
            GMAX = MAX (GMAX, PARBOR(IFAC,IVAR)/TSTATP)
            GMIN = MIN (GMIN, PARBOR(IFAC,IVAR)/TSTATP)
            GMOY = GMOY + (PARBOR(IFAC,IVAR)/TSTATP)
          ENDIF
        ENDDO
C
      ELSE IF (IMOYBR(IVAR).EQ.0) THEN
C
        DO IFAC = 1,NFABOR
          IF (PARBOR(IFAC,INBR).GT.SEUILF) THEN
            NBRFAC = NBRFAC + 1
            GMAX = MAX (GMAX, PARBOR(IFAC,IVAR))
            GMIN = MIN (GMIN, PARBOR(IFAC,IVAR))
            GMOY = GMOY + PARBOR(IFAC,IVAR)
          ENDIF
        ENDDO
C
      ENDIF
C
      IF (NBRFAC.GT.0) THEN
        GMOY = GMOY /DBLE(NBRFAC)
      ELSE
        GMAX =  0.D0
        GMIN =  0.D0
        GMOY =  0.D0
      ENDIF
C
C-------
C FORMAT
C-------
C
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
