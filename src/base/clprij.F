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
                        SUBROUTINE CLPRIJ
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NVAR   , NPHAS  ,
     &   IPHAS  , ICLIP  ,
     &   PROPCE , RTPA   , RTP    )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CLIPPING DE Rij ET EPSILON
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE DE CELLULES                   !
CARGU ! NVAR         ! E  !  -> ! NOMBRE DE VARIABLES                  !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! IPHAS        ! E  !  -> ! NUMERO DE LA PHASE A TRAITER         !
CARGU ! ICLIP        ! E  !  -> ! INDICATEUR = 1 ON N'UTILISE PAS RTPA !
CARGU !              !    !     !  (INIVAR)                            !
CARGU !              !    !     !            SINON ON PEUT (TURRIJ)    !
CARGU ! PROPCE       ! TR ! <-> ! TABLEAUX DES VARIABLES AU PDT COURANT!
CARGU !(NCELET,*     !    !     !                                      !
CARGU ! RTPA         ! TR !  -> ! TABLEAUX DES VARIABLES AU PDT PRECEDT!
CARGU ! (NCELET,NVAR)!    !     !                                      !
CARGU ! RTP          ! TR ! <-> ! TABLEAUX DES VARIABLES AU PDT COURANT!
CARGU ! (NCELET,NVAR)!    !     !                                      !
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
C
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
      INCLUDE "entsor.h"
      INCLUDE "numvar.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NVAR, NCELET, NCEL, NPHAS
      INTEGER          IPHAS, ICLIP
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION RTPA(NCELET,NVAR)
      DOUBLE PRECISION RTP(NCELET,NVAR)
C
C VARIABLES LOCALES
C
      INTEGER          ICLEPS, IEL, IVAR, IVAR1, IVAR2, ISOU, IPP
      INTEGER          IR11IP, IR22IP, IR33IP, IR12IP, IR13IP, IR23IP
      INTEGER          IEPIP
      INTEGER          ICLRIJ(6)
      DOUBLE PRECISION VMIN, VMAX, VAR, RIJMIN, VARREL, UND0, EPZ2
C
C***********************************************************************
C
C
      IR11IP = IR11(IPHAS)
      IR22IP = IR22(IPHAS)
      IR33IP = IR33(IPHAS)
      IR12IP = IR12(IPHAS)
      IR13IP = IR13(IPHAS)
      IR23IP = IR23(IPHAS)
      IEPIP  = IEP(IPHAS)
C
C Une petite valeur pour eviter des valeurs exactement nulles.
C
      EPZ2 = EPZERO**2
C
C=======================================================================
C  ---> Stockage Min et Max pour listing
C=======================================================================
C
      DO ISOU = 1, 7
        IF    (ISOU.EQ.1) THEN
          IVAR = IR11IP
        ELSEIF(ISOU.EQ.2) THEN
          IVAR = IR22IP
        ELSEIF(ISOU.EQ.3) THEN
          IVAR = IR33IP
        ELSEIF(ISOU.EQ.4) THEN
          IVAR = IR12IP
        ELSEIF(ISOU.EQ.5) THEN
          IVAR = IR13IP
        ELSEIF(ISOU.EQ.6) THEN
          IVAR = IR23IP
        ELSEIF(ISOU.EQ.7) THEN
          IVAR = IEPIP
        ENDIF
        IPP = IPPRTP(IVAR)
C
        VMIN =  GRAND
        VMAX = -GRAND
        DO IEL = 1, NCEL
          VAR = RTP(IEL,IVAR)
          VMIN = MIN(VMIN,VAR)
          VMAX = MAX(VMAX,VAR)
        ENDDO
        IF (IRANGP.GE.0) THEN
          CALL PARMIN(VMIN)
C         ===========
          CALL PARMAX(VMAX)
C         ===========
        ENDIF
        VARMNA(IPP) = VMIN
        VARMXA(IPP) = VMAX
C
      ENDDO
C
C ---> Clipping (modif pour eviter les valeurs exactement nulles)
C
      IF(ICLIP.EQ.1) THEN
C
        DO ISOU = 1, 3
C
          ICLRIJ(ISOU) = 0
C
          IF(ISOU.EQ.1) IVAR=IR11IP
          IF(ISOU.EQ.2) IVAR=IR22IP
          IF(ISOU.EQ.3) IVAR=IR33IP
C
          DO IEL = 1, NCEL
            IF (ABS(RTP(IEL,IVAR)).LE.EPZ2) THEN
              ICLRIJ(ISOU) = ICLRIJ(ISOU) + 1
              RTP(IEL,IVAR) = MAX(RTP(IEL,IVAR),EPZ2)
            ELSEIF(RTP(IEL,IVAR).LE.0.D0) THEN
              ICLRIJ(ISOU) = ICLRIJ(ISOU) + 1
              RTP(IEL,IVAR) = ABS(RTP(IEL,IVAR))
            ENDIF
          ENDDO
C
        ENDDO
C
        ICLEPS = 0
        DO IEL = 1, NCEL
          IF (ABS(RTP(IEL,IEPIP)).LE.EPZ2) THEN
            ICLEPS = ICLEPS + 1
            RTP(IEL,IEPIP) = MAX(RTP(IEL,IEPIP),EPZ2)
          ELSEIF(RTP(IEL,IEPIP).LE.0.D0) THEN
            ICLEPS = ICLEPS + 1
            RTP(IEL,IEPIP) = ABS(RTP(IEL,IEPIP))
          ENDIF
        ENDDO
C
      ELSE
C
        VARREL = 1.1D0
C
        DO ISOU = 1, 3
C
          ICLRIJ(ISOU) = 0
C
          IF(ISOU.EQ.1) IVAR=IR11IP
          IF(ISOU.EQ.2) IVAR=IR22IP
          IF(ISOU.EQ.3) IVAR=IR33IP
C
          DO IEL = 1, NCEL
            IF (ABS(RTP(IEL,IVAR)).LE.EPZ2) THEN
              ICLRIJ(ISOU) = ICLRIJ(ISOU) + 1
              RTP(IEL,IVAR) = MAX(RTP(IEL,IVAR),EPZ2)
            ELSEIF(RTP(IEL,IVAR).LE.0.D0) THEN
              ICLRIJ(ISOU) = ICLRIJ(ISOU) + 1
              RTP(IEL,IVAR) =
     &          MIN(ABS(RTP(IEL,IVAR)), VARREL*ABS(RTPA(IEL,IVAR)))
            ENDIF
          ENDDO
C
        ENDDO
C
        ICLEPS = 0
        DO IEL = 1, NCEL
          IF (ABS(RTP(IEL,IEPIP)).LT.EPZ2) THEN
            ICLEPS = ICLEPS + 1
            RTP(IEL,IEPIP) = MAX(RTP(IEL,IEPIP),EPZ2)
          ELSEIF(RTP(IEL,IEPIP).LE.0.D0) THEN
            ICLEPS = ICLEPS + 1
            RTP(IEL,IEPIP) =
     &          MIN(ABS(RTP(IEL,IEPIP )), VARREL*ABS(RTPA(IEL,IEPIP )))
          ENDIF
        ENDDO
C
      ENDIF
C
C On force l'inegalite de Cauchy Schwarz
C
      DO ISOU = 4, 6
C
        ICLRIJ(ISOU) = 0
C
        IF(ISOU.EQ.4) THEN
          IVAR  = IR12IP
          IVAR1 = IR11IP
          IVAR2 = IR22IP
        ELSEIF(ISOU.EQ.5) THEN
          IVAR  = IR13IP
          IVAR1 = IR11IP
          IVAR2 = IR33IP
        ELSEIF(ISOU.EQ.6) THEN
          IVAR  = IR23IP
          IVAR1 = IR22IP
          IVAR2 = IR33IP
        ENDIF
        UND0 = 1.D0
        DO IEL = 1 , NCEL
          RIJMIN = SQRT(RTP(IEL,IVAR1)*RTP(IEL,IVAR2))
          IF(RIJMIN.LT.ABS(RTP(IEL,IVAR))) THEN
            RTP(IEL,IVAR) = SIGN(UND0,RTP(IEL,IVAR)) * RIJMIN
            ICLRIJ(ISOU) = ICLRIJ(ISOU) + 1
          ENDIF
        ENDDO
C
      ENDDO

C
C ---> Stockage nb de clippings pour listing
C
      IF (IRANGP.GE.0) THEN
        CALL PARCPT (ICLRIJ(1))
C       ===========
        CALL PARCPT (ICLRIJ(2))
C       ===========
        CALL PARCPT (ICLRIJ(3))
C       ===========
        CALL PARCPT (ICLRIJ(4))
C       ===========
        CALL PARCPT (ICLRIJ(5))
C       ===========
        CALL PARCPT (ICLRIJ(6))
C       ===========
        CALL PARCPT (ICLEPS)
C       ===========
      ENDIF
C
      ICLPMN(IPPRTP(IR11IP)) = ICLRIJ(1)
      ICLPMN(IPPRTP(IR22IP)) = ICLRIJ(2)
      ICLPMN(IPPRTP(IR33IP)) = ICLRIJ(3)
      ICLPMN(IPPRTP(IR12IP)) = ICLRIJ(4)
      ICLPMN(IPPRTP(IR13IP)) = ICLRIJ(5)
      ICLPMN(IPPRTP(IR23IP)) = ICLRIJ(6)
      ICLPMN(IPPRTP(IEPIP )) = ICLEPS
C
      RETURN
C
      END
c@z
