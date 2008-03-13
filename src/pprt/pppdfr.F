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
                        SUBROUTINE PPPDFR
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , INDPDF ,
     &   FM     , FP2M   ,
     &   FMINI  , FMAXI  ,
     &   DIRMIN , DIRMAX , FDEB   , FFIN   , HREC )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL DES PARAMETRES DE LA PDF
CFONC  PDF RECTANGLE - PICS DE DIRAC "GENERALISEE" PPl - AE
CFONC
CFONC  LE RESULTAT EST :
CFONC  ---------------
CFONC     CALCUL DES PARAMETRES ASSOCIES AUX FONCTIONS RECTANGLE - DIRAC
CFONC
CFONC          INDPDF indique le passage ou non par la pdf
CFONC          DIRMIN contient le Dirac en FMINI
CFONC          DIRMAX - - - -  le Dirac en FMAXI
CFONC          FDEB - - - -  l'abcisse de debut du rectangle
CFONC          FFIN - - - - - - - - - - - fin   - - - - - -
CFONC          HREC - - - -  la hauteur du rectangle
CFONC
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
CARGU ! INDPDF       ! TE ! <-  ! INDICATEUR PASSAGE OU NON PAR PDF    !
CARGU ! FM           ! TR !  -> ! MOYENNE DE LA FRACTION DE MELANGE    !
CARGU ! FP2M         ! TR !  -> ! VARIANCE DE LA FRACTION DE MELANGE   !
CARGU ! FMINI        ! TR !  -> ! BORNE MIN DE LA FRACTION DE MELANGE !
CARGU ! FMAXI        ! TR !  -> ! BORNE MAX DE LA FRACTION DE MELANGE !
CARGU ! DIRMIN       ! TR !  <- ! DIRAC EN FMINI                       !
CARGU ! DIRMAX       ! TR !  <- ! DIRAC EN FMAXI                       !
CARGU ! FDEB         ! TR !  <- ! ABSCISSE DEBUT RECTANGLE             !
CARGU ! FFIN         ! TR !  <- ! ABSCISSE FIN RECTANGLE               !
CARGU ! HREC         ! TR !  <- ! HAUTEUR RECTANGLE                    !
CARGU ! XINPDF       ! TR !  <- ! INDICATEUR PASSAGE OU NON PAR PDF    !
CARGU !              !    !     ! POURRAIT ETRE ENTIER MAIS PUISQU'ON A!
CARGU !              !    !     ! UN TABLEAU DE TRAVAIL REEL ON S"EN   !
CARGU !              !    !     ! SERT                                 !
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
C**********************************************************************
C     DONNEES EN COMMON
C**********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "pointe.h"
      INCLUDE "parall.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET, NCEL
      INTEGER          INDPDF(NCELET)
C
      DOUBLE PRECISION FM(NCELET), FP2M(NCELET)
      DOUBLE PRECISION FMINI(NCELET), FMAXI(NCELET)
      DOUBLE PRECISION DIRMIN(NCELET), DIRMAX(NCELET)
      DOUBLE PRECISION FDEB(NCELET), FFIN(NCELET)
      DOUBLE PRECISION HREC(NCELET)

C
C VARIABLES LOCALES
C
      INTEGER          IEL, N1, N2, N3, N4, N5
      DOUBLE PRECISION T1, T2, T3, T1MOD, T2MOD
C
C
C***********************************************************************
C
C=======================================================================
C 1.  CALCULS PRELIMINAIRES
C=======================================================================
C
C Parametre relatif a la variance
      T1 = 1.D-08
C Parametre relatif a la moyenne
      T2 = 5.D-07
C
      DO IEL = 1, NCEL
        INDPDF(IEL) = 0
C Modifs des parametres T1 et T2 afin de tenir compte du fait que
C   FMINI < FM < FMAXI
        T1MOD = T1*(FMAXI(IEL)-FMINI(IEL))**2
        T2MOD = T2*(FMAXI(IEL)-FMINI(IEL))
        IF ( (FP2M(IEL).GT.T1MOD)
     &       .AND.(FM(IEL) .GE. (FMINI(IEL) + T2MOD))
     &       .AND.(FM(IEL) .LE. (FMAXI(IEL) - T2MOD)) ) THEN
          INDPDF(IEL) = 1
        ENDIF
      ENDDO
C
C
C=======================================================================
C 2.  CALCUL DES PARAMETRES DE LA FONCTION DENSITE DE PROBABILITE
C=======================================================================
C
      DO IEL = 1, NCEL
C
        IF ( INDPDF(IEL).EQ.1 ) THEN
          IF (    (     (FM(IEL) .LE.(FMINI(IEL) + FMAXI(IEL))*0.5D0)
     &            .AND.(FP2M(IEL).LE.(FM(IEL) - FMINI(IEL))**2/3.D0))
     &       .OR. (     (FM(IEL) .GT.(FMINI(IEL) + FMAXI(IEL))*0.5D0)
     &            .AND.(FP2M(IEL).LE.(FMAXI(IEL) -FM(IEL))**2/3.D0)) )
     &      THEN
C
C --> Rectangle seul
C
            HREC(IEL)   = SQRT(3.D0*FP2M(IEL))
            DIRMIN(IEL) = 0.D0
            DIRMAX(IEL) = 0.D0
            FDEB(IEL)   = FM(IEL) - HREC(IEL)
            FFIN(IEL)   = FM(IEL) + HREC(IEL)
C
          ELSEIF(      (FM(IEL)  .LE.(FMINI(IEL) + FMAXI(IEL))*0.5D0)
     &           .AND. (FP2M(IEL).LE.((FM(IEL) - FMINI(IEL))
     &              *(2.D0*FMAXI(IEL) +FMINI(IEL)-3.D0*FM(IEL))/3.D0)) )
     &      THEN
C
C --> Rectangle et un Dirac en FMINI
C
            FDEB(IEL)   = FMINI(IEL)
            DIRMAX(IEL) = 0.D0
            FFIN(IEL)   = FMINI(IEL) +1.5D0*( (FM(IEL) - FMINI(IEL))**2
     &                                       + FP2M(IEL) )
     &                                    /(FM(IEL) - FMINI(IEL))
            DIRMIN(IEL) = (3.D0*FP2M(IEL) -(FM(IEL) - FMINI(IEL))**2)
     &                  / (3.D0*((FM(IEL) - FMINI(IEL))**2 +FP2M(IEL)))
C
          ELSEIF(      (FM(IEL)  .GT.(FMINI(IEL) + FMAXI(IEL))*0.5D0)
     &           .AND. (FP2M(IEL).LE.((FMAXI(IEL) - FM(IEL))
     &               *(3.D0*FM(IEL)-FMAXI(IEL)-2.D0*FMINI(IEL))/3.D0)) )
     &      THEN
C
C --> Rectangle et un Dirac en FMAXI (c'est juste ;
C                          le HI/81/02/03/A contient une erreur  p 12)
C
            FFIN(IEL)   = FMAXI(IEL)
            DIRMIN(IEL) = 0.D0
            FDEB(IEL)   = FMINI(IEL)
     &                  + ( 3.D0*( (FM(IEL)-FMINI(IEL))**2+FP2M(IEL) )
     &                     + (FMAXI(IEL) - FMINI(IEL))**2
     &                     - 4.D0*(FM(IEL) - FMINI(IEL))
     &                           *(FMAXI(IEL) - FMINI(IEL))
     &                    )
     &                    / (2.D0*(FM(IEL) - FMAXI(IEL)))
            DIRMAX(IEL) = ( 3.D0*FP2M(IEL) -(FM(IEL) - FMAXI(IEL))**2 )
     &                 / ( 3.D0*((FM(IEL) - FMAXI(IEL))**2 +FP2M(IEL)) )
C
          ELSE
C
C --> Rectangle et deux Diracs
C
            FDEB(IEL)   = FMINI(IEL)
            FFIN(IEL)   = FMAXI(IEL)
            DIRMAX(IEL) = 3.D0*((FM(IEL) - FMINI(IEL))**2 +FP2M(IEL))
     &                    /(FMAXI(IEL) - FMINI(IEL))**2
     &                   -2.D0*(FM(IEL) - FMINI(IEL))
     &                    /(FMAXI(IEL) - FMINI(IEL))
            DIRMIN(IEL) = DIRMAX(IEL) + 1.D0 - 2.D0*(FM(IEL)-FMINI(IEL))
     &                                          /(FMAXI(IEL)-FMINI(IEL))
C
          ENDIF
C
          IF ( ABS(FFIN(IEL) - FDEB(IEL)).GT.EPZERO ) THEN
            HREC(IEL) = ( 1.D0-DIRMIN(IEL)-DIRMAX(IEL) )
     &                / ( FFIN(IEL)-FDEB(IEL) )
          ELSE
            T3 = SQRT(3.D0*T1*(FMAXI(IEL)-FMINI(IEL))**2)
            FDEB(IEL) = MIN(FMAXI(IEL),MAX(FMINI(IEL),FM(IEL) - T3))
            FFIN(IEL) = MIN(FMAXI(IEL),MAX(FMINI(IEL),FM(IEL) + T3))
            HREC(IEL) = ( 1.D0-DIRMIN(IEL)-DIRMAX(IEL) )
     &                 /( FFIN(IEL) - FDEB(IEL) )
          ENDIF
C
        ELSE
          DIRMIN(IEL) = 0.D0
          DIRMAX(IEL) = 0.D0
          FDEB(IEL)   = 0.D0
          FFIN(IEL)   = 0.D0
          HREC(IEL)   = 0.D0
        ENDIF
C
      ENDDO
C
C
C=======================================================================
C 4.  IMPRESSION
C=======================================================================
C
      N1 = 0
      N2 = 0
      N3 = 0
      N4 = 0
      N5 = 0
      DO IEL = 1, NCEL
        IF ( INDPDF(IEL).EQ.1 ) THEN
          N1 = N1+1
          IF ( DIRMIN(IEL).GT.EPZERO
     &         .AND. DIRMAX(IEL).LT.EPZERO ) THEN
            N2 = N2+1
          ELSE IF ( DIRMIN(IEL).LT.EPZERO
     &         .AND. DIRMAX(IEL).GT.EPZERO ) THEN
            N3 = N3+1
          ELSE IF ( DIRMIN(IEL).GT.EPZERO
     &         .AND. DIRMAX(IEL).GT.EPZERO ) THEN
            N4 = N4+1
          ELSE IF ( DIRMIN(IEL).LT.EPZERO
     &         .AND. DIRMAX(IEL).LT.EPZERO ) THEN
            N5 = N5+1
          ENDIF
        ENDIF
      ENDDO
C
      IF ( IRANGP.GE.0 ) THEN
C
        CALL PARCPT (N1)
C       ===========
        CALL PARCPT (N2)
C       ===========
        CALL PARCPT (N3)
C       ===========
        CALL PARCPT (N4)
C       ===========
        CALL PARCPT (N5)
C       ===========
      ENDIF
C
      WRITE(NFECRA,1000) N1, NCEL
      WRITE(NFECRA,2000) N5, N2, N3, N4
C
C
C----
C FORMATS
C----
C
 1000 FORMAT ( /,
     &'PDF RECTANGLE - PICS DE DIRAC COPDFR',/,
     &'MOYENNE, VARIANCE DU TRACEUR TRANPORTES',/,
     &'NOMBRE DE POINTS TURBULENTS (PASSAGE PAR LES PDF)   = ',I6,/,
     &'NOMBRE DE POINTS DE CALCULS                         = ',I6 )
 2000 FORMAT(
     &' Nb points avec PDF rectangle sans Dirac                = ',I6,/,
     &' - - - - - - - - - -- - - -  et Dirac en FMINI          = ',I6,/,
     &' - - - - - - - - - -- - - - - - - - - -  FMAXI          = ',I6,/,
     &' - - - - - - - - - - - - - - - Diracs en FMINI et FMAXI = ',I6,/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
