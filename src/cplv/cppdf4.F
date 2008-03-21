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
                        SUBROUTINE CPPDF4
C                       *****************
C     -------------------------------------------------------------
     & (  NCELET , NCEL   ,
     &    F1M    , F2M    , F3M    , F4M    , F4P2M  ,
     &    INDPDF ,
     &    SI7    , SI8    , SP2M   , F4I7   )
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
CFONC   --> RECONSTITUTION DE 4 MOMENTS : I7 = I4 et I8 = I6
CFONC   -->
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
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
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
      DOUBLE PRECISION F4M(NCELET), F4P2M(NCELET)
      DOUBLE PRECISION SI7(NCELET), SI8(NCELET), SP2M(NCELET)
      DOUBLE PRECISION F4I7(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL
C
      DOUBLE PRECISION F1I7, F2I7, F3I7, F4I8, S2MAX, T1, T2
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
C --- Parametres relatifs a la variance T1 et a la moyenne T2
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
C    INDPDF = 3  : SUPPORT (I7,I8)= (I4,M)
C=======================================================================
C
      DO IEL = 1, NCEL
        IF ( F4P2M(IEL).GT.T1 .AND.
     &       F4M(IEL)  .GE.T2 .AND. F4M(IEL).LE.(1.D0-T2) ) THEN
          INDPDF(IEL) = 3
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
      DO IEL = 1, NCEL
C
        IF ( INDPDF(IEL).EQ.3 ) THEN
          F4I7(IEL) = 1.D0
          F1I7       = 0.D0
          F2I7       = 0.D0
          F3I7       = 0.D0
          SI7(IEL) = - SQRT (
     &      ( SQRT(6.D0)/2.D0*(F1M(IEL)-F1I7) +
     &        SQRT(6.D0)/4.D0*(F2M(IEL)+F3M(IEL)-F2I7-F3I7) )**2 +
     &      ( 3.D0*SQRT(2.D0)/4.D0*(F2M(IEL)-F2I7) +
     &        SQRT(2.D0)/4.D0*(F3M(IEL)-F3I7) )**2 +
     &       ( F3M(IEL)-F3I7 )**2 )
          F4I8 = (1.D0-F3MAX)*F3M(IEL)
     &        / (F3M(IEL)+F3MAX*(1.D0-F3M(IEL)-F4M(IEL)))
          SI8(IEL)  = SI7(IEL)*(F4M(IEL)-F4I8)/(F4M(IEL)-F4I7(IEL))
          SP2M(IEL) = F4P2M(IEL)/(F4M(IEL)-F4I7(IEL))**2*(SI7(IEL))**2
          S2MAX     = -SI7(IEL)*SI8(IEL)
          IF (SP2M(IEL).GT.S2MAX) INDPDF(IEL) = 0
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
