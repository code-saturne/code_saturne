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
                        SUBROUTINE PPLECD
C                       *****************
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  LECTURE DU FICHIER DE DONNEES PHYSIQUE PARTICULIERE
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
      INCLUDE "pointe.h"
      INCLUDE "vector.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
C
C VARIABLES LOCALES
C
C
C***********************************************************************
C
C=======================================================================
C 1. AIGUILLAGE VERS LE MODELE ADEQUAT
C=======================================================================
C
C ---> Flamme de diffusion  - Chimie 3 points
C      Flamme de premelange - Modele EBU
C      Flamme de premelange - Modele LWC
C
      IF ( IPPMOD(ICOD3P).GE.0 .OR. IPPMOD(ICOEBU).GE.0
     &                         .OR. IPPMOD(ICOLWC).GE.0  ) THEN
        CALL COLECD
C       ===========
      ENDIF
C
C
C ---> Flamme charbon pulverise
C
      IF ( IPPMOD(ICP3PL).GE.0 ) THEN
        CALL CPLECD
C       ===========
      ENDIF
C
C ---> Combustion charbon pulverise couple transport Lagrangien
C      des particules de charbon
C
      IF ( IPPMOD(ICPL3C).GE.0 ) THEN
        CALL CPLECD
C       ===========
      ENDIF
C
C ---> Flamme fuel
C
      IF ( IPPMOD(ICFUEL).GE.0 ) THEN
        CALL FULECD
C       ===========
      ENDIF
C
C ---> Version Electrique : Effet Joule, Arc Electrique,
C                           Conduction Ionique
C
      IF ( IPPMOD(IELJOU).GE.1 .OR.
     &     IPPMOD(IELARC).GE.1 .OR.
     &     IPPMOD(IELION).GE.1      ) THEN
        CALL ELLECD
C       ===========
      ENDIF
C
C----
C FIN
C----
C
      RETURN
      END
C
c@z
