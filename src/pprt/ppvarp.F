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
                        SUBROUTINE PPVARP
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INIT DES POSITIONS DES VARIABLES SELON
CFONC    LE TYPE DE PHYSIQUE PARTICULIERE
CFONC  REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP
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
C
C***********************************************************************
C
C
C***********************************************************************
C
C ---> Physique particuliere : Combustion Gaz
C
      IF (  IPPMOD(ICOD3P).GE.0 .OR. IPPMOD(ICOEBU).GE.0
     &                          .OR. IPPMOD(ICOLWC).GE.0 ) THEN
        CALL COVARP
C       ===========
      ENDIF
C
C ---> Physique particuliere :  Combustion Charbon Pulverise
C
      IF ( IPPMOD(ICP3PL).GE.0 ) THEN
        CALL CPVARP
C       ===========
      ENDIF
C
C ---> Physique particuliere :  Combustion Charbon Pulverise
C      Couplee Transport Lagrangien des particules de charbon
C
      IF ( IPPMOD(ICPL3C).GE.0 ) THEN
        CALL CPLVAR
C       ===========
      ENDIF
C
C ---> Physique particuliere :  Combustion Fuel
C
      IF ( IPPMOD(ICFUEL).GE.0 ) THEN
        CALL FUVARP
C       ===========
      ENDIF
C
C ---> Physique particuliere : Compressible
C
      IF ( IPPMOD(ICOMPF).GE.0) THEN
        CALL CFVARP
C       ===========
      ENDIF
C
C ---> Physique particuliere : Versions Electriques
C
      IF ( IPPMOD(IELJOU).GE.1 .OR.
     &     IPPMOD(IELARC).GE.1 .OR.
     &     IPPMOD(IELION).GE.1       ) THEN
        CALL ELVARP
C       ===========
      ENDIF
C
      RETURN
      END
c@z
