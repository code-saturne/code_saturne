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
                        SUBROUTINE CPLVAR
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
CFONC    --------------------------------------------------------------
CFONC
CFONC     ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE
CFONC
CFONC       COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
CFONC       TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON
CFONC
CFONC       INIT DES POSITIONS DES VARIABLES TRANSPORTEES
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
      INTEGER        IS, IPHAS, ICHA, ISC
C
C***********************************************************************
C
C=======================================================================
C 1. DEFINITION DES POINTEURS
C=======================================================================
C
C ---> Variables propres a la phase continue
C
      IHM   = ISCAPP(1)
C
C ---> Variables propres a la phase continue
C
      DO ICHA = 1, NCHARB
        IS          = 1+ICHA
        IF1M(ICHA)  = ISCAPP(IS)
        IS          = 1+NCHARB+ICHA
        IF2M(ICHA)  = ISCAPP(IS)
      ENDDO
C
      IS = 1+NCHARB*2
      IS = IS+1
      IF3M  = ISCAPP(IS)
      IS = IS+1
      IF4P2M = ISCAPP(IS)
C
C=======================================================================
C 2. PROPRIETES PHYSIQUES
C    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
C    - PROPRES AUX SCALAIRES   : IPHSCA, IVISLS, ISCAVR
C      Rq : pas de variance associee a un scalaire dans notre cas
C    - PROPRES A LA SUSPENSION : ICP
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
C ---- Bien que l on soit en enthalpie on conserve un CP constant
C
      IPHAS = IPHSCA(IHM)
      ICP(IPHAS)    = 0
C
C
      RETURN
      END
c@z
