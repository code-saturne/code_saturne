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
                        SUBROUTINE CLPSCA
C                       *****************
C     --------------------------------------------------------------
     & ( NCELET , NCEL   , NVAR   , NSCAL  , ISCAL  ,
     &   PROPCE , SCANDD , RTP    )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CLIPPING
CFONC    POUR UN SCALAIRE OU VARIANCE
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! ISCAL        ! E  !  -> ! NUMERO DU SCALAIRE                   !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT        )!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! SCANDD       ! TR !  -> ! SCALAIRE AUQUEL EST ASSOCIE LA       !
CARGU ! (NCELET)     !    !     !    VARIANCE TRAITEE (SI C'EN EST UNE)!
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
C-----------------------------------------------------------------------
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
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL
      INTEGER          NVAR   , NSCAL
      INTEGER          ISCAL
C
      DOUBLE PRECISION RTP(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION SCANDD(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IVAR, IEL
      INTEGER          ICLMAX, ICLMIN, IISCAV
      INTEGER          IPPVAR
      DOUBLE PRECISION VMIN, VMAX, VFMIN, VFMAX
      DOUBLE PRECISION SCMAX, SCMIN
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C --- Numero de variable de calcul et de post associe au scalaire traite
      IVAR   = ISCA(ISCAL)
      IPPVAR = IPPRTP(IVAR)
C
C --- Numero du scalaire eventuel associe dans le cas fluctuation
      IISCAV = ISCAVR(ISCAL)
C
C=======================================================================
C 2. IMPRESSIONS ET CLIPPINGS
C=======================================================================
C
C      IL Y A TOUJOURS CLIPPING DES VARIANCES A DES VALEURS POSITIVES
C
C --- Calcul du min et max
      VMIN = RTP(1,IVAR)
      VMAX = RTP(1,IVAR)
      DO IEL = 1, NCEL
        VMIN = MIN(VMIN,RTP(IEL,IVAR))
        VMAX = MAX(VMAX,RTP(IEL,IVAR))
      ENDDO
      IF (IRANGP.GE.0) THEN
        CALL PARMIN(VMIN)
C       ===========
        CALL PARMAX(VMAX)
C       ===========
      ENDIF
      VARMNA(IPPVAR) = VMIN
      VARMXA(IPPVAR) = VMAX
C
      IF(IISCAV.EQ.0) THEN
C
C --- Clipping des scalaires non variances
C
        ICLMAX = 0
        ICLMIN = 0
        IF(SCAMAX(ISCAL).GT.SCAMIN(ISCAL))THEN
          DO IEL = 1, NCEL
            IF(RTP(IEL,IVAR).GT.SCAMAX(ISCAL))THEN
              ICLMAX = ICLMAX + 1
              RTP(IEL,IVAR) = SCAMAX(ISCAL)
            ENDIF
            IF(RTP(IEL,IVAR).LT.SCAMIN(ISCAL))THEN
              ICLMIN = ICLMIN + 1
              RTP(IEL,IVAR) = SCAMIN(ISCAL)
            ENDIF
          ENDDO
        ENDIF
C
        IF (IRANGP.GE.0) THEN
          CALL PARCPT (ICLMIN)
C         ===========
          CALL PARCPT (ICLMAX)
C         ===========
        ENDIF
C
        ICLPMN(IPPVAR) = ICLMIN
        ICLPMX(IPPVAR) = ICLMAX
C
      ELSE
C
C --- Clipping des variances
C
        ICLMAX = 0
        ICLMIN = 0
C
C   -- Clipping minimal au minimum 0.
        IF(ICLVFL(ISCAL).EQ.0) THEN
          DO IEL = 1, NCEL
            IF(RTP(IEL,IVAR).LT.0.D0) THEN
              ICLMIN = ICLMIN + 1
              RTP(IEL,IVAR) = 0.D0
            ENDIF
          ENDDO
C
C   -- Clipping a partir des valeurs du scalaire (ou 0 au min)
        ELSEIF(ICLVFL(ISCAL).EQ.1) THEN
          DO IEL = 1, NCEL
            IF(RTP(IEL,IVAR).LT.0.D0) THEN
              ICLMIN = ICLMIN + 1
              RTP(IEL,IVAR) = 0.D0
            ENDIF
          ENDDO
C
          SCMAX = SCAMAX(ISCAVR(ISCAL))
          SCMIN = SCAMIN(ISCAVR(ISCAL))
          DO IEL = 1, NCEL
            VFMAX = (SCANDD(IEL)-SCMIN)*(SCMAX-SCANDD(IEL))
            IF(RTP(IEL,IVAR).GT.VFMAX) THEN
              ICLMAX = ICLMAX + 1
              RTP(IEL,IVAR) = VFMAX
            ENDIF
          ENDDO
C
C   -- Clipping a partir des valeurs donnees par l'utilisateur
C        (ou 0 au min)
        ELSEIF(ICLVFL(ISCAL).EQ.2) THEN
          VFMIN = 0.D0
          VFMIN = MAX(SCAMIN(ISCAL),VFMIN)
          VFMAX = SCAMAX(ISCAL)
          IF(VFMAX.GT.VFMIN)THEN
            DO IEL = 1, NCEL
              IF(RTP(IEL,IVAR).GT.VFMAX)THEN
                ICLMAX = ICLMAX + 1
                RTP(IEL,IVAR) = VFMAX
              ENDIF
              IF(RTP(IEL,IVAR).LT.VFMIN)THEN
                ICLMIN = ICLMIN + 1
                RTP(IEL,IVAR) = VFMIN
              ENDIF
            ENDDO
          ENDIF
        ENDIF
C
        IF (IRANGP.GE.0) THEN
          CALL PARCPT (ICLMIN)
C         ===========
          CALL PARCPT (ICLMAX)
C         ===========
        ENDIF
C
        ICLPMN(IPPVAR) = ICLMIN
        ICLPMX(IPPVAR) = ICLMAX
C
      ENDIF
C
C
C--------
C FORMATS
C--------
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
