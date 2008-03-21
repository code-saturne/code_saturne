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
                        SUBROUTINE ELVARP
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC       INIT DES POSITIONS DES VARIABLES POUR LE MODULE ELECTRIQUE
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
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER        IS, IESP , IDIMVE, ISC, IPHAS
C
C***********************************************************************
C=======================================================================
C 1. DEFINITION DES POINTEURS
C=======================================================================
C
C 1.0 Dans toutes les versions electriques
C ========================================
C
C ---- Enthalpie
      IS  = 1
      IHM = ISCAPP(IS)
C
C ---- Potentiel reel
      IS = IS+1
      IPOTR = ISCAPP(IS)
C
C 1.1 Effet Joule (cas potentiel imaginaire)
C ==========================================
C
      IF(IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4) THEN
C
C ---- Potentiel imaginaire
        IS = IS+1
        IPOTI  = ISCAPP(IS)
C
      ENDIF
C
C 1.2 Arc electrique
C ==================
C
      IF ( IPPMOD(IELARC).GE.2 ) THEN
C
C ---- Potentiel vecteur
        DO IDIMVE = 1, NDIMVE
          IS = IS+1
          IPOTVA(IDIMVE) = ISCAPP(IS)
        ENDDO
      ENDIF
C
C 1.3 Conduction ionique
C ======================
C
C
C 1.4 Dans toutes les versions electriques
C ========================================
C
C ---- Fractions massiques des constituants
C
      IF ( NGAZG .GT. 1 ) THEN
        DO IESP = 1, NGAZG-1
          IS = IS+1
          IYCOEL(IESP)=ISCAPP(IS)
        ENDDO
      ENDIF
C
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
C ---- Notre physique particuliere est monophasique
          IPHSCA(ISCAPP(ISC)) = 1
C
C ---- Viscosite dynamique moleculaire variable pour les
C                                              scalaires ISCAPP(ISC)
C        Pour l'enthalpie en particulier.
C        Pour le potentiel vecteur, voir plus bas
          IVISLS(ISCAPP(ISC)) = 1
C
        ENDIF
C
      ENDDO
C
C ---- "Viscosite dynamique moleculaire" = 1
C                                  pour le potentiel vecteur en Arc
      IF ( IPPMOD(IELARC).GE.2 ) THEN
        DO IDIMVE = 1, NDIMVE
          IVISLS(IPOTVA(IDIMVE)) = 0
        ENDDO
      ENDIF
C
C ---- Cp est variable ; pas sur que ce soit indispensable pour le verre
C                              mais pour le moment c'est comme ca.
      IPHAS      = IPHSCA(IHM)
      ICP(IPHAS) = 1
C
      RETURN
      END
C
c@z
