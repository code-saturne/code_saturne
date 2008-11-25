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
                        SUBROUTINE PPCSCA
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC        CALCUL DE NSCAPP SUIVANT LA PHYSIQUE PARTICULIERE
CFONC        DEMANDEE PAR l'UTILISATEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
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
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppcpfu.h"
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER          IOK
C
C
C***********************************************************************
C
C=======================================================================
C 1. MODELES DE COMBUSTION GAZ (ICO...)
C=======================================================================
C
C --> Flamme de diffusion
C
      IF ( IPPMOD(ICOD3P).EQ.0 ) NSCAPP = 2
      IF ( IPPMOD(ICOD3P).EQ.1 ) NSCAPP = 3
C
      IF ( IPPMOD(ICODEQ).EQ.0 ) NSCAPP = 2
      IF ( IPPMOD(ICODEQ).EQ.1 ) NSCAPP = 3
C
C --> Flamme de premelange (modele EBU)
C
      IF ( IPPMOD(ICOEBU).EQ.0 ) NSCAPP = 1
      IF ( IPPMOD(ICOEBU).EQ.1 ) NSCAPP = 2
      IF ( IPPMOD(ICOEBU).EQ.2 ) NSCAPP = 2
      IF ( IPPMOD(ICOEBU).EQ.3 ) NSCAPP = 3
C
C ---> Modele BML
C
      IF ( IPPMOD(ICOBML).EQ.0 ) NSCAPP = 2
      IF ( IPPMOD(ICOBML).EQ.1 ) NSCAPP = 3
C
C ---> Modele LWC
C
      IF ( IPPMOD(ICOLWC).EQ.0 ) NSCAPP = 4
      IF ( IPPMOD(ICOLWC).EQ.1 ) NSCAPP = 5
      IF ( IPPMOD(ICOLWC).EQ.2 ) NSCAPP = 5
      IF ( IPPMOD(ICOLWC).EQ.3 ) NSCAPP = 6
      IF ( IPPMOD(ICOLWC).EQ.4 ) NSCAPP = 5
      IF ( IPPMOD(ICOLWC).EQ.5 ) NSCAPP = 6
C
C ---nombre de Dirac pour le modele LWC
C
      IF ( IPPMOD(ICOLWC).EQ.0 .OR.
     &     IPPMOD(ICOLWC).EQ.1           ) THEN
        NDIRAC = 2
C
      ELSE IF ( IPPMOD(ICOLWC).EQ.2 .OR.
     &          IPPMOD(ICOLWC).EQ.3      ) THEN
        NDIRAC = 3
C
      ELSE IF ( IPPMOD(ICOLWC).EQ.4 .OR.
     &          IPPMOD(ICOLWC).EQ.5      ) THEN
        NDIRAC = 4
      ENDIF
C
C=======================================================================
C 2. MODELES DE COMBUSTION CHARBON PULVERISE (ICP...)
C=======================================================================
C
C --> Flamme charbon pulverise
C
C ------ Transport d'H2
      IF ( IPPMOD(ICP3PL).EQ.0 ) THEN
C
        NSCAPP = 1 + 2*NCHARB + 2 + 4*NCLACP + (NOXYD-1)
C
        IF ( IHTCO2.EQ. 1) NSCAPP = NSCAPP + 1
        IF ( IEQCO2.GE. 1) NSCAPP = NSCAPP + 1
C
      ENDIF
C
C ------ Transport d'H2 + phase de sechage
C
      IF ( IPPMOD(ICP3PL).EQ.1 ) THEN
C
        NSCAPP = 1 + 2*NCHARB + 3 +5*NCLACP + (NOXYD-1)
C
        IF ( IHTCO2.EQ. 1) NSCAPP = NSCAPP + 1
        IF ( IEQCO2.EQ. 1) NSCAPP = NSCAPP + 1
C
      ENDIF
C
C --> Flamme charbon pulverise couple Lagrangien
C
      IF ( IPPMOD(ICPL3C).EQ.0 ) THEN
        NSCAPP = 1 + 2*NCHARB + 2
      ENDIF
C
C=======================================================================
C 3. MODELE COMPRESSIBLE SANS CHOC : rho, e, s
C=======================================================================
C
      IF ( IPPMOD(ICOMPF).GE.0 ) NSCAPP = 3*NPHAS
C
C
C=======================================================================
C 4. MODELES ELECTRIQUES : Effet Joule        (IELJOU)
C                          Arc electrique     (IELARC)
C                          Conduction ionique (IELION)
C=======================================================================
C
C
C --> Effet Joule
C
        IF ( IPPMOD(IELJOU).EQ.1 ) NSCAPP = 2 + NGAZG-1
        IF ( IPPMOD(IELJOU).EQ.2 ) NSCAPP = 3 + NGAZG-1
C
        IF ( IPPMOD(IELJOU).EQ.3 ) NSCAPP = 2 + NGAZG-1
        IF ( IPPMOD(IELJOU).EQ.4 ) NSCAPP = 3 + NGAZG-1
C
C --> Arc electrique
C
        IF ( IPPMOD(IELARC).EQ.1 ) NSCAPP = 2 + NGAZG-1
        IF ( IPPMOD(IELARC).EQ.2 ) NSCAPP = 5 + NGAZG-1
C
C --> Conduction ionique
C
        IF ( IPPMOD(IELION).EQ.1 ) NSCAPP = 2 + NGAZG-1
C
C=======================================================================
C 5. MODELES DE COMBUSTION FUEL (ICFUEL)
C=======================================================================
C
      IF ( IPPMOD(ICFUEL).GE.0 ) THEN
        NSCAPP = 4 + 3*NCLAFU
        IF ( IEQCO2.GE. 1) NSCAPP = NSCAPP + 1
        IF ( IEQNOX.EQ. 1) NSCAPP = NSCAPP + 3
      ENDIF
C
C=======================================================================
C 6. MODELE ATMOSPHERIQUE (IATMOS)
C=======================================================================
C
      IF ( IPPMOD(IATMOS).GE.0 ) NSCAPP = 0
C
C=======================================================================
C 7. VERIFICATION : UNE SEULE PHASE
C=======================================================================
C
      IOK = 0
C
      IF( NSCAPP.GT.0 .AND. IPPMOD(ICOMPF).LT.0 ) THEN
C                     ^^^^^^^^^^^^^^^^^^^^^^^^^
C     On peut penser faire du multiphasique avec la physique particuliere
C     compressible sans choc... (mettons deux pressions indépendantes
C     par exemple) ... a mettre en place si necessaire ...
C
        IF(NPHAS.NE.1) THEN
          WRITE(NFECRA,1000)
          IOK = IOK + 1
        ENDIF
      ENDIF
C
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    SEULE UNE PHASE PERMISE EN PHYSIQUE PARTICULIERE        ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@  Arret dans ppcsca.                                        ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)

C
      RETURN
      END
c@z
