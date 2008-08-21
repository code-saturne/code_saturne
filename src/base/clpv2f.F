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
                        SUBROUTINE CLPV2F
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NVAR   , NPHAS  ,
     &   IPHAS  , IWAPHI ,
     &   PROPCE , RTP    )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CLIPPING DE PHI EN V2F (PAS DE CLIPPING SUR F_BARRE)
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
CARGU ! IWAPHI       ! E  !  -> ! NIVEAU D'IMPRESSION                  !
CARGU ! PROPCE       ! TR ! <-> ! TABLEAUX DES VARIABLES AU PDT COURANT!
CARGU !(NCELET,*     !    !     !                                      !
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
      INTEGER          IPHAS, IWAPHI
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION RTP(NCELET,NVAR)
C
C VARIABLES LOCALES
C
      INTEGER          IEL, IPP
      INTEGER          IPHIPH
      INTEGER          NCLPMX, NCLPMN
      DOUBLE PRECISION XPHI, VMIN, VMAX, VAR
C
C***********************************************************************
C
C
      IPHIPH = IPHI(IPHAS)
C
C=======================================================================
C  ---> Stockage Min et Max pour listing
C=======================================================================
C
      IPP = IPPRTP(IPHIPH)
C
      VMIN =  GRAND
      VMAX = -GRAND
      DO IEL = 1, NCEL
        VAR = RTP(IEL,IPHIPH)
        VMIN = MIN(VMIN,VAR)
        VMAX = MAX(VMAX,VAR)
      ENDDO
      IF (IRANGP.GE.0) THEN
        CALL PARMIN(VMIN)
C       ===========
        CALL PARMAX(VMAX)
C       ===========
      ENDIF
      VARMNA(IPP) = VMIN
      VARMXA(IPP) = VMAX
C
C=====================================================================
C  ---> Reperage des valeurs superieures a 2, pour affichage seulement
C=====================================================================
C
      IF (IWAPHI.GE.2) THEN
        NCLPMX = 0
        DO IEL = 1, NCEL
          IF (RTP(IEL,IPHIPH).GT.2.D0) NCLPMX = NCLPMX+1
        ENDDO
        IF(IRANGP.GE.0) CALL PARCPT(NCLPMX)
C                       ===========
        IF (NCLPMX.GT.0) WRITE(NFECRA,1000) IPHAS,NCLPMX
      ENDIF

C=====================================================================
C  ---> Clipping en valeur absolue pour les valeurs negatives
C=====================================================================
C
      NCLPMN = 0
      DO IEL = 1, NCEL
        XPHI = RTP(IEL,IPHIPH)
        IF (XPHI.LT.0.D0) THEN
          RTP(IEL,IPHIPH) = -XPHI
          NCLPMN = NCLPMN + 1
        ENDIF
      ENDDO
      IF(IRANGP.GE.0) CALL PARCPT(NCLPMN)
C                     ===========
      ICLPMN(IPP) = NCLPMN
C
C=======================================================================
C ---> Formats
C=======================================================================
C
#if defined(_CS_LANG_FR)
C
 1000 FORMAT('ATTENTION VARIABLE PHI, PHASE ',I3,/,
     &     'VALEUR MAXIMALE PHYSIQUE DE 2 DEPASSEE SUR ',I10,
     &     ' CELLULES')
C
#else
C
 1000 FORMAT('WARNING VARIABLE PHI, PHASE ',I3,/,
     &     'MAXIMUM PHYSICAL VALUE OF 2 EXCEEDED FOR ',I10,
     &     ' CELLS')
C
#endif
C
      RETURN
C
      END
c@z
