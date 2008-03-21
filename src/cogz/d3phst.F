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
                        SUBROUTINE D3PHST
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , INDPDF ,
     &   DIRMIN , DIRMAX , FDEB   , FFIN   , HREC   ,
     &   FM     , HM     ,
     &   HSTOE  )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC ROUTINE PHYSIQUE PARTICULIERE : FLAMME DE DIFFUSION
CFONC  CALCUL DE L'ENTHALPIE STOECHIOMETRIQUE LOCALE
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! INDPDF       ! TE !  -> ! INDICATEUR PASSAGE OU NON PAR LES PDF!
CARGU ! DIRMIN       ! TR !  -> ! PDF : DIRAC EN FMIN                  !
CARGU ! DIRMAX       ! TR !  -> ! PDF : DIRAC EN FMAX                  !
CARGU ! FDEB         ! TR !  -> ! PDF : ABSCISSE DEBUT RECTANGLE       !
CARGU ! FFIN         ! TR !  -> ! PDF : ABSCISSE FIN RECTANGLE         !
CARGU ! HREC         ! TR !  -> ! PDF : HAUTEUR RECTANGLE              !
CARGU ! FM           ! TR !  -> ! FRACTION DE MELANGE MOYENNE          !
CARGU ! HM           ! TR !  -> ! ENTHALPIE MASSIQUE MOYENNE           !
CARGU !              !    !     !  SI ECOULEMENT PERMEATIQUE           !
CARGU ! HSTOE        ! TR ! <-> ! ENTHALPIE STOECHIOMETRIQUE           !
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
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C ARGUMENTS
C
      INTEGER          NCELET, NCEL
      INTEGER          INDPDF(NCELET)
      DOUBLE PRECISION DIRMIN(NCELET), DIRMAX(NCELET)
      DOUBLE PRECISION FDEB(NCELET), FFIN(NCELET), HREC(NCELET)
      DOUBLE PRECISION FM(NCELET), HM(NCELET), HSTOE(NCELET)
C
C
C VARIABLES LOCALES
C
      INTEGER          ICEL
      DOUBLE PRECISION FSIR, HHH, HCT, F1, F2
      DOUBLE PRECISION EPSI
C
      INTEGER          N1, N2
      DOUBLE PRECISION HSMAX, HSMIN
C
C
C***********************************************************************
C
C=======================================================================
C 0. INITIALISATION
C=======================================================================
C
      EPSI  = 1.D-6
      FSIR =  FS(1)
C
      N1 = 0
      N2 = 0
      HSMIN = GRAND
      HSMAX =-GRAND
C
C
      DO ICEL = 1, NCEL
C
        IF ( INDPDF(ICEL) .EQ. 0 ) THEN
C
C=======================================================================
C 1. DETERMINATION DE HSTOE SANS INTEGRATION
C=======================================================================
C
          IF ( FM(ICEL).LE.FSIR .AND. FM(ICEL).GT.EPSI ) THEN
            HSTOE(ICEL) = ( FSIR*HM(ICEL)+(FM(ICEL)-FSIR)*HINOXY )
     &                  / FM(ICEL)
          ELSEIF( FM(ICEL).LT.(1.D0-EPSI) ) THEN
            HSTOE(ICEL) = ((1.D0-FSIR)*HM(ICEL)+(FSIR-FM(ICEL))*HINFUE)
     &                 / (1.D0-FM(ICEL))
          ENDIF
C
        ELSE
C
C=======================================================================
C 2. DETERMINATION DE HSTOE AVEC INTEGRATION
C=======================================================================
C
          HCT = DIRMIN(ICEL)*HINOXY+DIRMAX(ICEL)*HINFUE
          HHH = 0.D0
          IF ( HREC(ICEL).GT.EPSI ) THEN
C
            IF (FDEB(ICEL).LE.FSIR) THEN
              F1 = FDEB(ICEL)
              F2 = MIN(FSIR,FFIN(ICEL))
              HCT = HCT + HREC(ICEL)*
     &              (F2-F1)*HINOXY*(2.D0*FSIR-F1-F2)/(2.D0*FSIR)
              HHH = HHH + HREC(ICEL)*(F2**2-F1**2)/(2.D0*FSIR)
            ENDIF
            IF (FFIN(ICEL).GT.FSIR) THEN
              F1 = MAX( FSIR,FDEB(ICEL))
              F2 = FFIN(ICEL)
              HCT = HCT + HREC(ICEL) *
     &             (F2-F1)*HINFUE*(F2+F1-2.D0*FSIR)/(2.D0*(1.D0-FSIR))
              HHH = HHH +
     &              HREC(ICEL)*(F2-F1)*(2.D0-F1-F2)/(2.D0*(1.D0-FSIR))
            ENDIF
            HSTOE(ICEL) = (HM(ICEL)-HCT)/ HHH
C
C Clipping a HSTOEA = HH(1)     en max
C Clipping a HSTOEA = HH(NMAXH) em min
C
            IF ( HSTOE(ICEL) .GT. HH(1) ) THEN
              N1 = N1 + 1
              HSMAX = MAX(HSTOE(ICEL),HSMAX)
              HSTOE(ICEL) = HH(1)
            ENDIF
C
            IF ( HSTOE(ICEL) .LT. HH(NMAXH) ) THEN
              N2 = N2 + 1
              HSMIN = MIN(HSTOE(ICEL),HSMIN)
              HSTOE(ICEL) = HH(NMAXH)
            ENDIF
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
      IF (IRANGP.GE.0) THEN
        CALL PARCPT (N1)
C       ===========
        CALL PARCPT (N2)
C       ===========
        CALL PARMAX (HSMAX)
C       ===========
        CALL PARMIN (HSMIN)
C       ===========
      ENDIF

      IF ( N1.GT.0 ) THEN
        WRITE(NFECRA,1000) N1,HSMAX,HH(1)
      ENDIF
      IF ( N2.GT.0 ) THEN
        WRITE(NFECRA,1001) N2,HSMIN,HH(NMAXH)
      ENDIF
C
C----
C FORMATS
C----
C
 1000   FORMAT(1X,' Clipping de HSTOE EN MAX EN ',I8,' POINTS',/,
     &         1X,'     Valeur Max : ',G15.7,/,
     &         1X,'     Valeur De Clipping : ',G15.7,/)
 1001   FORMAT(1X,' Clipping de HSTOE EN MIN EN ',I8,' POINTS',/,
     &         1X,'     Valeur Max : ',G15.7,/,
     &         1X,'     Valeur De Clipping : ',G15.7,/)
C
      RETURN
      END
C
c@z
