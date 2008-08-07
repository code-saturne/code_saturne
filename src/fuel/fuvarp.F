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
                        SUBROUTINE FUVARP
C                       *****************
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC        INIT DES POSITIONS DES VARIABLES TRANSPORTEES POUR
CFONC                 COMBUSTION FUEL
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
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
      INCLUDE "ppcpfu.h"
C
C***********************************************************************
C
      INTEGER        IS, IPHAS, ISC, ICLA
C
C***********************************************************************
C=======================================================================
C 1. DEFINITION DES POINTEURS
C=======================================================================
C
C ---> Variables propres a la suspension gaz - particules
C
      IHM   = ISCAPP(1)
C
C ---> Variables propres a la phase dispersee
C
      DO ICLA = 1, NCLAFU
        IS          = 1+ICLA
        ING(ICLA)   = ISCAPP(IS)
        IS          = 1+1*NCLAFU+ICLA
        IYFOL(ICLA) = ISCAPP(IS)
        IS          = 1+2*NCLAFU+ICLA
        IHLF(ICLA)  = ISCAPP(IS)
      ENDDO
C
C ---> Variables propres a la phase gaz
C
      IS       = IS + 1
      IFVAP    = ISCAPP(IS)
      IS       = IS + 1
      IFHTF    = ISCAPP(IS)
      IS       = IS + 1
      IF4P2M   = ISCAPP(IS)
C
      IF ( IEQCO2 .GE. 1 ) THEN
        IS = IS+1
        IYCO2 = ISCAPP(IS)
      ENDIF
C
      IF ( IEQNOX .EQ. 1 ) THEN
        IS = IS+1
        IYHCN = ISCAPP(IS)
        IS = IS+1
        IYNO  = ISCAPP(IS)
        IS = IS+1
        ITAIRE = ISCAPP(IS)
      ENDIF
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
        IF ( ISCAVR(ISCAPP(ISC)) .LE. 0 ) THEN
C
C ---- Notre physique particuliere est monophasique
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
      RETURN
      END
c@z
