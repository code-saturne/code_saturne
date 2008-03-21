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
                        SUBROUTINE FUPHY2
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   RTP    , PROPCE )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC  CALCUL DES PROPRIETES PHYSIQUES DE LA PHASE DISPERSEE
CFONC  VALEURS CELLULES
CFONC  ----------------
CFONC
CFONC    FRACTION MASSIQUE DE LIQUIDE
CFONC      ET CLIPPING EVENTUELS
CFONC    DIAMETRE
CFONC    MASSE VOLUMIQUE
CFONC      ET CLIPPING EVENTUELS
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! RTP          ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT)        !
CARGU ! PROPCE       ! TR ! <-> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
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
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL
C
      DOUBLE PRECISION RTP(NCELET,*) , PROPCE(NCELET,*)
C
C VARIABLES LOCALES
C
      INTEGER          IEL
      INTEGER          N1     , N2     , IPCDIA
      DOUBLE PRECISION XNG    ,  D1S3
      DOUBLE PRECISION RHOFOL , DIAM3M , DIAM3X
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
      D1S3 = 1.D0/3.D0
C
C=======================================================================
C 2. CALCUL POUR CHAQUE CLASSE
C    DE LA MASSE VOLUMIQUE DU FOL
C    DE LA FRACTION MASSIQUE DE FOL
C    DU DIAMETRE DU COKE
C=======================================================================
C
      N1 = 0
      N2 = 0
      DIAM3M =  1.D0
      DIAM3X =  0.D0
C
      IPCDIA = IPPROC(IDIAM3)
C
      DO IEL = 1, NCEL
C
C      Masse Volumique
C
        PROPCE(IEL,IPPROC(IROM3)) = RHO0FL
C
        YFOL   = RTP(IEL,ISCA(IYFOL))
        XNG    = RTP(IEL,ISCA(ING))
        RHOFOL = PROPCE(IEL,IPPROC(IROM3))
C
C --- Calcul du diametre de la particule
C
C     Attention !
C     Diam en mm
C     Xng en Ggouttes (en milliards de gouttes/kg de mélange)
C     ça marche mais il faudrait un facteur mis à UN pour s'y retrouver
C
        IF ( YFOL .GT. EPSIFL .AND. (XNG*YFOL) .GT. 0.D0) THEN
          PROPCE(IEL,IPCDIA) = ((YFOL / RHOFOL)
     &                        /(PI/6.D0 * XNG) ) ** D1S3
        ELSE
          PROPCE(IEL,IPCDIA) = DINIFL
        ENDIF
C
        IF ( PROPCE(IEL,IPCDIA) .GT. DINIFL ) THEN
          N1 = N1+1
          DIAM3X = MAX(DIAM3X,PROPCE(IEL,IPCDIA))
          PROPCE(IEL,IPCDIA) = DINIFL
        ENDIF
C
        IF ( PROPCE(IEL,IPCDIA) .LT. DINIIN ) THEN
          N2 = N2+1
          DIAM3M = MIN(DIAM3M,PROPCE(IEL,IPCDIA))
          PROPCE(IEL,IPCDIA) = DINIIN
        ENDIF
C
      ENDDO
C
      IF (IRANGP.GE.0) THEN
C
        CALL PARCPT (N1)
C       ===========
        CALL PARCPT (N2)
C       ===========
C
        CALL PARMAX (DIAM3X)
C       ===========
        CALL PARMIN (DIAM3M)
C       ===========
      ENDIF
C
      IF ( N1 .GT. 0 ) THEN
         WRITE(NFECRA,1001) N1, DIAM3X
      ENDIF
      IF ( N2 .GT. 0 ) THEN
         WRITE(NFECRA,1002)  N2, DIAM3M
      ENDIF
C
C----
C FORMATS
C----
C
 1001 FORMAT(/,1X,' CLIPPING EN MAX DU DIAMETRE :',
     &       /,10X,' Nombre de points : ',I8,
     &       /,10X,' Valeur Max       : ',G15.7)
 1002 FORMAT(/,1X,' CLIPPING EN MIN DU DIAMETRE :',
     &       /,10X,' Nombre de points : ',I8,
     &       /,10X,' Valeur Min       : ',G15.7)
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
