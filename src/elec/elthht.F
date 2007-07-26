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
                        SUBROUTINE ELTHHT
C                       *****************
C     -------------------------------------------------------------
     & ( MODE   , NESP , YESP , ENTHAL , TEMPER )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  --------
c@foncb
CFONC
CFONC CALCULE L'ENTHALPIE OU LA TEMPERATURE  A PARTIR DE LA
CFONC   COMPOSITION ET DE LA VALEUR DE LA TEMPERATURE OU DE
CFONC                     ENTHALPIE
CFONC         SPECIFIQUE AU MODULE ELECTRIQUE
c@fonce
C
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! MODE         ! E  !  -> !  -1 : T -> H  ;   1 : H -> T         !
CARGU ! NESP         ! E  !  -> ! NB DE CONSTITUANTS                   !
CARGU ! YESP         ! TR !  -> ! FRACTION MASSIQUE DES CONSTITUANTS   !
CARGU ! ENTHAL       ! R  ! <-> ! ENTHALPIE                            !
CARGU ! TEMPER       ! R  ! <-> ! TEMPERATURE                          !
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "ppthch.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C
C ARGUMENTS
C
      INTEGER          MODE , NESP
C
      DOUBLE PRECISION ENTHAL , TEMPER , YESP(NESP)
C
C
C VARIABLES LOCALES
C
      INTEGER          IT , IESP
C
      DOUBLE PRECISION EH1 , EH0
C
C
C***********************************************************************
C
C=======================================================================
C 1. CALCUL DE L'ENTHALPIE A PARTIR DE LA TEMPERATURE
C=======================================================================
C
      IF ( MODE.EQ.-1 ) THEN
C
        IT=NPO
        IF (TEMPER.GE.TH(IT)) THEN
          ENTHAL = 0
          DO IESP=1,NESP
            ENTHAL = ENTHAL + YESP(IESP)*EHGAZG(IESP,IT)
          ENDDO
          RETURN
        ENDIF
C
        IT=1
        IF (TEMPER.LE.TH(IT)) THEN
          ENTHAL = 0
          DO IESP=1,NESP
            ENTHAL = ENTHAL + YESP(IESP)*EHGAZG(IESP,IT)
          ENDDO
          RETURN
        ENDIF
C
 10     CONTINUE
          IT=IT+1
          IF (TEMPER.LE.TH(IT)) THEN
            EH0 = 0
            EH1 = 0
            DO IESP=1,NESP
              EH0 = EH0 + YESP(IESP)*EHGAZG(IESP,IT-1)
              EH1 = EH1 + YESP(IESP)*EHGAZG(IESP,IT)
            ENDDO
            ENTHAL=EH0+(EH1-EH0)*(TEMPER-TH(IT-1))
     &                          /(TH(IT)-TH(IT-1))
            RETURN
          ENDIF
        GOTO 10
C
C
C=======================================================================
C 2. CALCUL DE LA TEMPERATURE A PARTIR DE l'ENTHALPIE
C=======================================================================
C
      ELSE IF ( MODE.EQ.1 ) THEN
C
        IT  = NPO
        EH1 = 0
        DO IESP=1,NESP
          EH1 = EH1+YESP(IESP)*EHGAZG(IESP,IT)
        ENDDO
        IF ( ENTHAL .GE. EH1 ) THEN
          TEMPER = TH (IT)
          RETURN
        ENDIF
C
        IT  = 1
        EH1 = 0
        DO IESP=1,NESP
          EH1 = EH1+YESP(IESP)*EHGAZG(IESP,IT)
        ENDDO
        IF ( ENTHAL .LE. EH1 ) THEN
           TEMPER = TH (IT)
           RETURN
        ENDIF
C
 20     CONTINUE
          IT  = IT+1
          EH0 = EH1
          EH1 = 0
          DO IESP=1,NESP
            EH1 = EH1+YESP(IESP)*EHGAZG(IESP,IT)
          ENDDO
          IF ( ENTHAL .LE. EH1 ) THEN
            TEMPER = TH(IT-1)+ (ENTHAL-EH0)
     &                        *(TH(IT)-TH(IT-1))/(EH1-EH0)
            RETURN
          ENDIF
        GOTO 20
C
C
      ELSE
C
        WRITE(NFECRA,1000) MODE
        CALL CSEXIT (1)
C       ===========
C
C
      ENDIF
C
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR DANS ELTHHT                          ',/,
     &'@    *********                                               ',/,
     &'@    VALEUR INCORRECTE DE L''ARGUMENT MODE                   ',/,
     &'@    CE DOIT ETRE UN ENTIER EGAL A 1 OU -1                   ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
