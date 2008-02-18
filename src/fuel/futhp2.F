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
                        SUBROUTINE FUTHP2
C                       *****************
C     --------------------------------------------------------------
     & ( MODE   , ENTHAL , XSOLID , TEMPER )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  --------
c@foncb
CFONC CALCUL DE LA TEMPERATURE DES PARTICULES
CFONC   EN FONCTION DE L'ENTHALPIE ET DES CONCENTRATIONS
CFONC   SI IMODE = 1
CFONC CALCUL DE L'ENTHALPIE DES PARTICULES
CFONC   EN FONCTION DE LA TEMPERATURE ET DES CONCENTRATIONS
CFONC   SI IMODE = -1
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! MODE         ! E  !  -> !  -1 : T -> H  ;   1 : H -> T         !
CARGU ! ENTHAL       ! R  ! <-> ! ENTHALPIE MASSIQUE J/kg              !
CARGU ! XSOLID       ! TR !  -> ! FRACTION MASSIQUE DES CONSTITUANTS   !
CARGU ! TEMPER       ! R  ! <-> ! TEMPERATURE EN KELVIN                !
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
      INCLUDE "pointe.h"
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C
C ARGUMENTS
C
      INTEGER          MODE
C
      DOUBLE PRECISION XSOLID(2)
      DOUBLE PRECISION TEMPER , ENTHAL
C
C VARIABLES LOCALES
C
      INTEGER          IT , ISOL , IHFLT2
C
      DOUBLE PRECISION EH1 , EH0
C
C***********************************************************************
C=======================================================================
C 1. RQ IMPORTANTE : On suppose pour l'instant que H2 = H02 + CP2(T2-TREF)
C=======================================================================
C
      IHFLT2 = 0
C
      IF ( IHFLT2.EQ.0 ) THEN
C
C=======================================================================
C 2. H2 FONCTION LINEAIRE T2
C=======================================================================
C
C
        IF ( MODE.EQ.-1 ) THEN
C
C --> Loi temperature -> enthalpie (MODE = -1)
C
          ENTHAL = H02FOL + CP2FOL*(TEMPER-TREFTH)
C
        ELSEIF ( MODE.EQ.1 ) THEN
C
C --> Loi enthalpie -> temperature (MODE = 1)
C
          TEMPER =  (ENTHAL-H02FOL)/CP2FOL + TREFTH
C
        ELSE
C
          WRITE(NFECRA,1000) MODE
          CALL CSEXIT (1)
C         ===========
C
        ENDIF
C
C
      ELSE
C
C=======================================================================
C 3. H2 TABULE
C=======================================================================
C
        IF ( MODE.EQ.-1 ) THEN
C
C --> Loi temperature -> enthalpie (MODE = -1)
C
          IT = NPO
          IF ( TEMPER.GE.TH(IT) ) THEN
            ENTHAL = ZERO
            DO ISOL = 1, 2
              ENTHAL = ENTHAL + XSOLID(ISOL)*EHSOLI(ISOL,IT)
            ENDDO
            GO TO 11
          ENDIF
C
          IT = 1
          IF ( TEMPER.LE.TH(IT) ) THEN
            ENTHAL = ZERO
            DO ISOL = 1, 2
              ENTHAL = ENTHAL + XSOLID(ISOL)*EHSOLI(ISOL,IT)
            ENDDO
            GO TO 11
          ENDIF
          IT = 1
 10       CONTINUE
C
          IT = IT + 1
          IF ( TEMPER.LE.TH(IT) ) THEN
            EH0 = ZERO
            EH1 = ZERO
            DO ISOL = 1, 2
              EH0 = EH0 + XSOLID(ISOL)*EHSOLI(ISOL,IT-1)
              EH1 = EH1 + XSOLID(ISOL)*EHSOLI(ISOL,IT  )
            ENDDO
            ENTHAL = EH0
     &             + (EH1-EH0)*(TEMPER-TH(IT-1))/(TH(IT)-TH(IT-1))
            GOTO 11
          ENDIF
          GOTO 10
 11       CONTINUE
C
        ELSEIF ( MODE.EQ.1 ) THEN
C
C --> Loi enthalpie -> temperature (MODE = 1)
C
          IT  = NPO-1
          EH1 = ZERO
          DO ISOL = 1, 2
            EH1 = EH1 + XSOLID(ISOL)*EHSOLI(ISOL,IT+1)
          ENDDO
          IF ( ENTHAL.GE.EH1 ) TEMPER = TH(IT+1)
C
          IT  = 1
          EH0 = ZERO
          DO ISOL = 1, 2
            EH0 = EH0 + XSOLID(ISOL)*EHSOLI(ISOL,IT  )
          ENDDO
          IF ( ENTHAL.LE.EH0 ) TEMPER = TH(IT)
C
          DO IT = 1, NPO-1
            EH0 = ZERO
            EH1 = ZERO
            DO ISOL = 1, 2
              EH0 = EH0 + XSOLID(ISOL)*EHSOLI(ISOL,IT  )
              EH1 = EH1 + XSOLID(ISOL)*EHSOLI(ISOL,IT+1)
            ENDDO
            IF ( ENTHAL.GE.EH0 .AND. ENTHAL.LE.EH1 )
     &        TEMPER = TH(IT)
     &               + (ENTHAL-EH0)*(TH(IT+1)-TH(IT))/(EH1-EH0)
          ENDDO
C
        ELSE
C
          WRITE(NFECRA,1000) MODE
          CALL CSEXIT (1)
C         ===========
C
        ENDIF
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
     &'@ @@ ATTENTION : ERREUR DANS FUTHP2                          ',/,
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
