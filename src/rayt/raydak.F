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
                        SUBROUTINE RAYDAK
C                       *****************
C     -------------------------------------------------------------
     & ( NCEL   , NCELET ,
     &   CK     , PCO2   , PH2O   , FV     , TEMP   )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE RAYONNEMENT :
CFONC    --------------------------------------
CFONC
CFONC   CALCUL DES PROPRIETES RADIATIVES D'UN GAZ EN FONCTION DE LA
CFONC   TEMPERATURE,DE LA COMPOSITION DES PRODUITS EN CO2, H2O ET SUIES
CFONC   EN UTILISANT LES REGRESSIONS LINEAIRES ETABLIES PAR MODAK.
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! CK (NCELET)  ! TR ! <-  ! COEFFICIENT D'ABSORPTION DU MILIEU   !
CARGU !              !    !     ! (NUL SI TRANSPARENT)                 !
CARGU ! PCO2(NCELET) ! TR !  -> ! PRESSION PARTIELLE DE CO2            !
CARGU ! PCO2(NCELET) ! TR !  -> ! PRESSION PARTIELLE DE H2O            !
CARGU ! FV  (NCELET) ! TR !  -> ! FRACTION VOLUMIQUE DE SUIES          !
CARGU ! TEMP(NCELET) ! TR !  -> ! TEMPERATURE                          !
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
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET, NCEL
      DOUBLE PRECISION CK(NCELET), TEMP(NCELET), FV(NCELET)
      DOUBLE PRECISION PCO2(NCELET), PH2O(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL
      DOUBLE PRECISION ALPHA, PATH, TE, TS, SOOTK, TMIN, TMAX
C
C***********************************************************************
C=======================================================================
C CALCULS
C=======================================================================
C
C --- Longueur moyenne de penetration du rayonnement
C
      PATH = 15.D0
C
      TMAX = 2000.D0
      TMIN = 300.D0
C
C ATTENTION : LES TEMPERATURES UTILISEES DANS MODAK SONT MISES EN KELVIN
C =========
C
      DO IEL = 1, NCEL
C
C --- Temperature du melange gazeux
C
        TE = TEMP(IEL)
C
C --- Temperature du corps noir
C
        TS = TEMP(IEL)
C
C --- Limitation A TMAX = 2000 K et TMIN = 300 K
C
        IF ( TEMP(IEL).GT.TMAX ) THEN
          TS = TMAX
          TE = TMAX
        ENDIF
        IF ( TEMP(IEL).LT.TMIN ) THEN
          TS = TMIN
          TE = TMIN
        ENDIF
C
C --- Fraction volumique de suies
C
        SOOTK = 7.D0*FV(IEL)/0.95D-6
C
C --- Calcul de l'absorptivite du fluide
C
        CALL ABSORB
C       ===========
     &  ( TS , TE , PATH , SOOTK , PCO2(IEL) , PH2O(IEL) , ALPHA )
C
C --- Test d'erreur
C
        IF ( (1.D0-ALPHA).LE.EPZERO ) THEN
          WRITE(NFECRA,1000) IEL, ALPHA, PCO2(IEL), PH2O(IEL),
     &                       SOOTK, TE, PATH, FV(IEL)
          CALL CSEXIT(1)
        ENDIF
C
C --- Calcul du coeffcient d'absorption
C
        CK(IEL) = - LOG(1.D0-ALPHA)/PATH
C
      ENDDO
C
C========
C FORMATS
C========
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ERREUR RAYDAK : CALCUL DE ABSORPTIVITE                  ',/,
     &'@    *************                                           ',/,
     &'@ IEL   = ', I10                                              ,/,
     &'@ ALPHA = ', G15.7                                            ,/,
     &'@ PCO2  = ', G15.7                                            ,/,
     &'@ PH2O  = ', G15.7                                            ,/,
     &'@ SOOTK = ', G15.7                                            ,/,
     &'@ TE    = ', G15.7                                            ,/,
     &'@ PATH  = ', G15.7                                            ,/,
     &'@ FV    = ', G15.7                                            ,/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
      RETURN
C
      END
c@z
