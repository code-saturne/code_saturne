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
                        SUBROUTINE CPVARP
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC        INIT DES POSITIONS DES VARIABLES TRANSPORTEES POUR
CFONC                 COMBUSTION CHARBON PULVERISE
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
      INCLUDE "ppcpfu.h"
      INCLUDE "ihmpre.h"
C
C***********************************************************************
C
      INTEGER        ICLA,  IS, IPHAS, ICHA, ISC , IS1
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
      DO ICLA = 1, NCLACP
        IS = 1+ICLA
        INP(ICLA) = ISCAPP(IS)
        IS = 1+1*NCLACP+ICLA
        IXCH(ICLA)= ISCAPP(IS)
        IS = 1+2*NCLACP+ICLA
        IXCK(ICLA) = ISCAPP(IS)
        IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
          IS = 1+3*NCLACP+ICLA
          IXWT(ICLA) = ISCAPP(IS)
C
          IS = 1+4*NCLACP+ICLA
          IH2(ICLA) = ISCAPP(IS)
C
        ELSE
          IS = 1+3*NCLACP+ICLA
          IH2(ICLA) = ISCAPP(IS)
        ENDIF
C
      ENDDO
C
C ---> Variables propres a la phase continue
C
      IS1 = IS
      DO ICHA = 1, NCHARB
        IS          = IS1+ICHA
        IF1M(ICHA)  = ISCAPP(IS)
        IS          = IS1+NCHARB+ICHA
        IF2M(ICHA)  = ISCAPP(IS)
      ENDDO
C
      IS = IS+1
      IF3M  = ISCAPP(IS)
      IF ( IHTCO2 .EQ. 1) THEN
        IS = IS+1
        IF3MC2 = ISCAPP(IS)
      ENDIF
      IS = IS+1
      IF4P2M = ISCAPP(IS)
      IF ( IPPMOD(ICP3PL) .EQ. 1 ) THEN
        IS = IS+1
        IF5M  = ISCAPP(IS)
      ENDIF
C
      IF ( NOXYD .GE. 2 ) THEN
        IS = IS+1
        IF6M  = ISCAPP(IS)
      ENDIF
      IF ( NOXYD .EQ. 3 ) THEN
        IS = IS+1
        IF7M  = ISCAPP(IS)
      ENDIF
C
      IF ( IEQCO2 .GE. 1 ) THEN
        IS = IS+1
        IYCO2 = ISCAPP(IS)
      ENDIF
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C     Construction de l'indirection entre la numerotation du noyau et XML
      IF(IIHMPR.EQ.1) THEN
        CALL UICPSC (NCHARB, NCLACP, IPPMOD, ICP3PL, IEQCO2,
     &               IHM, INP, IXCH, IXCK, IXWT, IH2,
     &               IF1M, IF2M, IF3M, IF4P2M, IF5M, IYCO2)
      ENDIF

#endif
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
C ---- Bien que l'on soit en enthalpie on conserve un CP constant
C
      IPHAS = IPHSCA(IHM)
      ICP(IPHAS)    = 0
C
      RETURN
      END
c@z
