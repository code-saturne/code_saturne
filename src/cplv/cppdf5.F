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
                        SUBROUTINE CPPDF5
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   F1M    , F2M    , F3M    , F4M    , F3P2M  , F4P2M  ,
     &   INDPDF ,
     &   SI7    , SI8    , SP2M   , F4I7   )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CALCUL DU TYPE DE PDF
CFONC   PDF CONJOINTE DEGENERRE EN UNE PDF 1D DE TYPE RECTANGLE - DIRAC
CFONC   SUPPORT (I7,I8) passant par M
CFONC
CFONC   --> RECONSTITUTION DE 5 MOMENTS
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! F1M          ! TR !  -> ! MOYENNE DU TRACEUR 1 (MV CHZ +CO)    !
CARGU ! F2M          ! TR !  -> ! MOYENNE DU TRACEUR 2 (MV CXHY +CO)   !
CARGU ! F3M          ! TR !  -> ! MOYENNE DU TRACEUR 3 (CO C.HET)      !
CARGU ! F4M          ! TR !  -> ! MOYENNE DU TRACEUR 4 (AIR)           !
CARGU ! F3P2M        ! TR !  -> ! VARIANCE DU TRACEUR 3 (CO C.HET)     !
CARGU ! F4P2M        ! TR !  -> ! VARIANCE DU TRACEUR 4 (AIR)          !
CARGU ! INDPDF       ! TE !  <- ! PASSAGE PAR LES PDF                  !
CARGU ! SI7          ! TR !  <- ! ABSCISSE CURVILIGNE AU POINT I7      !
CARGU ! SI8          ! TR !  <- ! ABSCISSE CURVILIGNE AU POINT I8      !
CARGU ! SP2M         ! TR !  <- ! VARIANCE DROITE SUPPORT              !
CARGU ! F4I7         ! TR !  <- ! F4 AU POINT I7                       !
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
C
      INCLUDE "paramx.h"
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
C ARGUMENTS
C
      INTEGER          NCELET, NCEL
      INTEGER          INDPDF(NCELET)
C
      DOUBLE PRECISION F1M(NCELET), F2M(NCELET), F3M(NCELET)
      DOUBLE PRECISION F4M(NCELET), F3P2M(NCELET), F4P2M(NCELET)
      DOUBLE PRECISION SI7(NCELET), SI8(NCELET), SP2M(NCELET)
      DOUBLE PRECISION F4I7(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL, II
      INTEGER          ITYP1, ITYPA, ITYP(2), ICPTS1
C
      DOUBLE PRECISION S2MAX, T1, T2
      DOUBLE PRECISION F4L5, RAV3S4
      DOUBLE PRECISION F1I7T(2), F2I7T(2), F3I7T(2), F4I7T(2), F4I8T(2)
      DOUBLE PRECISION SI8T(2), SI7T(2), SP2MT(2)
      DOUBLE PRECISION A, B, C, ANGALP, ANGGAM
      DOUBLE PRECISION F4I7TS, TETA(2)
      DOUBLE PRECISION XR1, F4I7T1, SI7T1, SI8T1, SP2MT1
      DOUBLE PRECISION XRA, F4I7TA, SI7TA, SI8TA, SP2MTA
      DOUBLE PRECISION F1I3, F2I3, F3I3, F1I5, F2I5, F3I5
      DOUBLE PRECISION XI4I5, XI4I3
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C --- Parametres relatifs a la variance T1 et a la moyenne T2
C
      T1 = 1.D-04
      T2 = 5.D-03
C
C --- Initialisation des tableaux
C
      DO IEL = 1, NCEL
        F4I7(IEL)   = 0.D0
        SI7(IEL)    = 0.D0
        SI8(IEL)    = 0.D0
        SP2M(IEL)   = 0.D0
        INDPDF(IEL) = 0.D0
      ENDDO
C
C
C=======================================================================
C 2. TYPE DE PDF
C    INDPDF = 0  : PAS DE PDF
C    INDPDF = 10 : SUPPORT (I7,I8) FORME AVEC TETA VARIABLE
C                  RECONSTITUTION DE 5 MOMENTS
C                    INDPDF = 11 : I7 APP. [L5,I3max]
C                    INDPDF = 12 : I7 APP. [I4,L3]
C                    INDPDF = 13 : I7 APP. [I4,L5]
C=======================================================================
C
      DO IEL = 1, NCEL
        IF ( F4P2M(IEL).GT.T1 .AND.
     &       F4M(IEL)  .GE.T2 .AND. F4M(IEL).LE.(1.D0-T2) ) THEN
          INDPDF(IEL) = 10
        ELSE
          INDPDF(IEL) = 0
        ENDIF
      ENDDO
C
C
C=======================================================================
C 3. CALCULS DES PARAMETRES DE LA PDF : SI7, SI8 et SP2M
C    CALCUL DE F4I7
C=======================================================================
C
!OCL SCALAR
C
      DO IEL = 1, NCEL
C
        IF ( INDPDF(IEL).EQ.10 ) THEN
C
          DO II = 1, 2
            F4I7T(II) = 0.D0
            SI7T(II)  = 0.D0
            SI8T(II)  = 0.D0
            SP2MT(II) = 0.D0
          ENDDO
C
C --> Calculs preliminaires
C
          RAV3S4 = F3P2M(IEL)/F4P2M(IEL)
C     Droite limite = (I5,M)
          F4L5 = F4M(IEL)/(F4M(IEL)+F3M(IEL))
C     Calcul de XI4I5
          F1I5 = F1M(IEL)/(F1M(IEL)+F2M(IEL))
          F2I5 = F2M(IEL)/(F1M(IEL)+F2M(IEL))
          F3I5 = 0.D0
          XI4I5 = SQRT (
     &            ( SQRT(6.D0)/2.D0*(0.D0-F1I5) +
     &              SQRT(6.D0)/4.D0*(0.D0-F2I5-F3I5)   )**2 +
     &            ( 3.D0*SQRT(2.D0)/4.D0*(0.D0-F2I5) +
     &                   SQRT(2.D0)/4.D0*(0.D0-F3I5)   )**2 +
     &            ( 0.D0-F3I5                          )**2 )
          F1I3 = 0.D0
          F2I3 = 0.D0
          F3I3 = 1.D0
C      Calcul de XI4I3
          XI4I3 = SQRT (
     &            ( SQRT(6.D0)/2.D0*(0.D0-F1I3) +
     &              SQRT(6.D0)/4.D0*(0.D0-F2I3-F3I3)  )**2 +
     &            ( 3.D0*SQRT(2.D0)/4.D0*(0.D0-F2I3) +
     &                   SQRT(2.D0)/4.D0*(0.D0-F3I3)  )**2 +
     &            ( 0.D0-F3I3                         )**2 )
          ANGALP = ACOS(XI4I3/(2.D0*XI4I5))
          ANGGAM = PI - 2.D0*ANGALP
C
C --> Determination des deux angles TETA
C
          A = 1.D0/(SIN(ANGALP)*XI4I3)
          B = 1.D0/XI4I5
          C = 1.D0/(XI4I5*TAN(ANGGAM))
          TETA(1) = ATAN(-SQRT(RAV3S4)*B/(A+SQRT(RAV3S4)*C))
          TETA(2) = ATAN( SQRT(RAV3S4)*B/(A-SQRT(RAV3S4)*C))
C
C --> Determination du type de PDF
C       - TYPE DE PDF : 13, I7 appartient au segment [I4,L5]
C       - TYPE DE PDF : 12, I7 appartient au segment [I4,L3]
C       - TYPE DE PDF : 11, I7 appartient au segment [L5,I3max]
C
          DO II = 1, 2
C
            A = SIN(TETA(II))/SIN(ANGALP)/XI4I3
            B = COS(TETA(II))/XI4I5+SIN(TETA(II))/XI4I5/TAN(ANGGAM)
            F4I7TS = (-A/B*F4M(IEL)+(1.D0-F3M(IEL)))/(1.D0-A/B)
C ---- Cette formulation permet de degenerer de facon naturelle
C      si RAV3S4 = 0.
C
            IF ( F4I7TS.LE.1.D0 .AND. F4I7TS.GT.F4L5 ) THEN
C
C ---- TYPE DE PDF : 13, I7 appartient au segment [I4,L5]
C                        I8 appartient au segment [I5,L4]
              ITYP(II) = 3
              F4I7T(II) = (-A/B*F4M(IEL)+(1.-F3M(IEL)))/(1.-A/B)
              F1I7T(II) = 0.D0
              F2I7T(II) = 0.D0
              F3I7T(II) = 1.D0 -F4I7T(II)
C
            ELSEIF ( F4I7TS.LE.F4L5 .AND.
     &               F4I7TS.GT.(1.D0-F3MAX) ) THEN
C
C ---- TYPE DE PDF : 11, I7 appartient au segment [L5,I3max]
C                        I8 appartient au segment [L3,I5]
              ITYP(II) = 1
              F4I7T(II) = (F4M(IEL)-(1.D0-F3M(IEL))*B/A)/(1.D0-B/A)
              F1I7T(II) = 0.D0
              F2I7T(II) = 0.D0
              F3I7T(II) = 1.D0 - F4I7T(II)
C
            ELSE
C
C TYPE DE PDF : 12, C appartient au segment [I4,L3]
C                   D appartient au segment [L4,I3max]
              ITYP(II) = 2
              F4I7T(II) = F3M(IEL)*B/A+F4M(IEL)
              F1I7T(II) = (1.D0-F4I7T(II))*F1M(IEL)/(F1M(IEL)+F2M(IEL))
              F2I7T(II) = (1.D0-F4I7T(II))*F2M(IEL)/(F1M(IEL)+F2M(IEL))
              F3I7T(II) = 0.D0
C
            ENDIF
C
            SI7T(II) = - SQRT (
     &  ( SQRT(6.D0)/2.D0*(F1M(IEL)-F1I7T(II)) +
     &    SQRT(6.D0)/4.D0*(F2M(IEL)+F3M(IEL)-F2I7T(II)-F3I7T(II)) )**2 +
     &  ( 3.D0*SQRT(2.D0)/4.D0*(F2M(IEL)-F2I7T(II)) +
     &         SQRT(2.D0)/4.D0*(F3M(IEL)-F3I7T(II))               )**2 +
     &  ( F3M(IEL)-F3I7T(II)                                      )**2 )
C
            F4I8T(II) =
     &  (1.D0-F3MAX)*(F4I7T(II)*(-F1M(IEL)-F2M(IEL)+1.D0)-F4M(IEL))
     &  / (F4I7T(II)-1.D0+F3M(IEL)+F3MAX*(F1M(IEL)+F2M(IEL)))
C
            SI8T(II) = SI7T(II)*
     &  (F4M(IEL)-F4I8T(II))/(F4M(IEL)-F4I7T(II))
C
            SP2MT(II) = F4P2M(IEL)/
     &  (F4M(IEL)-F4I7T(II))**2*(SI7T(II))**2
            S2MAX = -SI7T(II)*SI8T(II)
            IF (SP2MT(II).GT.S2MAX) SP2MT(II) = 0.D0
C
          ENDDO
C
C
C --> Vu que l'on a deux angles, on a deux supports possibles
C     On priviligie alors la forme 11 avec le support le plus long
C
          XRA     = ZERO
          SI7TA   = 0.D0
          SI8TA   = 0.D0
          SP2MTA  = 0.D0
          F4I7TA  = 0.D0
          ITYPA   = -10
C
          XR1     = ZERO
          SI7T1   = 0.D0
          SI8T1   = 0.D0
          SP2MT1  = 0.D0
          F4I7T1  = 0.D0
          ITYP1   = -10
          ICPTS1  = 0
C
          DO II = 1, 2
C ---- On priviligie la forme 11 avec le support le plus long
            IF ( (ITYP(II).EQ.1) ) THEN
              IF ( (SI8T(II)-SI7T(II)).GT.XR1 .AND.
     &             (SP2MT(II).GT.EPZERO)         ) THEN
                XR1     = (SI8T(II)-SI7T(II))
                F4I7T1  = F4I7T(II)
                ITYP1   = ITYP(II)
                SI7T1   = SI7T(II)
                SI8T1   = SI8T(II)
                SP2MT1  = SP2MT(II)
                ICPTS1  = ICPTS1 + 1
              ENDIF
            ELSE
              IF ( (SI8T(II)-SI7T(II)).GT.XRA .AND.
     &             (SP2MT(II).GT.EPZERO)         ) THEN
                XRA     = (SI8T(II)-SI7T(II))
                F4I7TA  = F4I7T(II)
                ITYPA   = ITYP(II)
                SI7TA   = SI7T(II)
                SI8TA   = SI8T(II)
                SP2MTA  = SP2MT(II)
              ENDIF
            ENDIF
          ENDDO
C
          IF ( ICPTS1.GE.1 ) THEN
            F4I7(IEL)   = F4I7T1
            SI7(IEL)    = SI7T1
            SI8(IEL)    = SI8T1
            SP2M(IEL)   = SP2MT1
            INDPDF(IEL) = INDPDF(IEL) + ITYP1
          ELSE
            F4I7(IEL)   = F4I7TA
            SI7(IEL)    = SI7TA
            SI8(IEL)    = SI8TA
            SP2M(IEL)   = SP2MTA
            INDPDF(IEL) = INDPDF(IEL) + ITYPA
          ENDIF
C
C
        ENDIF
C
      ENDDO
C
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
