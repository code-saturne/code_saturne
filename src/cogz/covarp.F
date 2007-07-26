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
                        SUBROUTINE COVARP
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC       INIT DES POSITIONS DES VARIABLES SELON
CFONC               POUR LA COMBUSTION
CFONC         FLAMME DE DIFFUSION ET DE PREMELANGE
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
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER        ISC, IPHAS
C
C***********************************************************************
C=======================================================================
C 1. DEFINITION DES POINTEURS
C=======================================================================
C
C 1.1 Flamme de diffusion : chimie 3 points
C =========================================
C
      IF ( IPPMOD(ICOD3P).GE.0 ) THEN
C
C ---- Taux de melange
        IFM = ISCAPP(1)
C
C ---- Variance du taux de melange
        IFP2M = ISCAPP(2)
        ISCAVR(IFP2M) = IFM
C
C ---- Enthalpie
        IF ( IPPMOD(ICOD3P).EQ.1 ) IHM = ISCAPP(3)
      ENDIF
C
C
C 1.2 Flamme de premelange : modele EBU
C =====================================
C
      IF ( IPPMOD(ICOEBU).GE.0 ) THEN
C
C ---- Fraction massique des gaz frais
        IYGFM   = ISCAPP(1)
        IF ( IPPMOD(ICOEBU).EQ.2 .OR.
     &       IPPMOD(ICOEBU).EQ.3      ) THEN
C
C ---- Taux de melange
          IFM   = ISCAPP(2)
        ENDIF
        IF ( IPPMOD(ICOEBU).EQ.1 ) IHM = ISCAPP(2)
        IF ( IPPMOD(ICOEBU).EQ.3 ) IHM = ISCAPP(3)
      ENDIF
C
C
C 1.3 Flamme de premelange : modele BML A DEVELOPPER
C ==================================================
C
C 1.4 Flamme de premelange : modele LWC
C =====================================
C
      IF (IPPMOD(ICOLWC).GE.0 ) THEN
C
        IFM           = ISCAPP(1)
        IFP2M         = ISCAPP(2)
        ISCAVR(IFP2M) = IFM
C
        IYFM          = ISCAPP(3)
        IYFP2M        = ISCAPP(4)
        ISCAVR(IYFP2M)= IYFM
C
        IF (IPPMOD(ICOLWC).GE.2 ) THEN
          ICOYFP = ISCAPP(5)
        ENDIF
C
        IF (IPPMOD(ICOLWC).EQ.1 ) IHM = ISCAPP(5)
        IF (IPPMOD(ICOLWC).EQ.3 .OR.
     &      IPPMOD(ICOLWC).EQ.5 ) IHM = ISCAPP(6)
C
      ENDIF
C
C=======================================================================
C 2. PROPRIETES PHYSIQUES
C    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
C      IPHSCA, IVISLS, ICP
C=======================================================================
C
      DO ISC = 1, NSCAPP
C
        IF ( ISCAVR(ISCAPP(ISC)).LE.0 ) THEN
C
C ---- Notre physque particuliere est monophasique
          IPHSCA(ISCAPP(ISC)) = 1
C
C ---- Viscosite dynamique de reference relative au scalaire
C      ISCAPP(ISC)
          IVISLS(ISCAPP(ISC)) = 0
C
        ENDIF
C
      ENDDO
C
      IF ( IPPMOD(ICOD3P).EQ.1 .OR.
     &     IPPMOD(ICOEBU).EQ.1 .OR.
     &     IPPMOD(ICOEBU).EQ.3 .OR.
     &     IPPMOD(ICOLWC).EQ.1 .OR.
     &     IPPMOD(ICOLWC).EQ.3 .OR.
     &     IPPMOD(ICOLWC).EQ.5   ) THEN
C
C
C ---- Bien que l on soit en enthalpie on conserve un CP constant
C
        IPHAS = IPHSCA(IHM)
        ICP(IPHAS)    = 0
C
      ENDIF
C
      RETURN
      END
C
c@z
