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
                        SUBROUTINE FUPDFR
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   FVAP   , FHTF   , F4P2M  ,
     &   INDPDF ,
     &   F4M1   , F4M2   , D4CL   , D4F4     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC CALCUL DU TYPE DE PDF
CFONC
c@fonce
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! FVAP         ! TR !  -> ! MOYENNE DU TRACEUR 1 MVl [FOVm+CO]   !
CARGU ! FHET         ! TR !  -> ! MOYENNE DU TRACEUR 3 C héterogène    !
CARGU ! F4P2M        ! TR !  -> ! VARIANCE DU TRACEUR 4 (air)          !
CARGU ! INDPDF       ! TE !  -> ! PASSAGE PAR LES PDF                  !
CARGU ! F4M1         ! TR ! <-> ! BORNE MINIMUM                        !
CARGU ! F4M2         ! TR ! <-> ! BORNE MAX                            !
CARGU ! D4CL         ! TR ! <-> ! AMPLITUDE DU PIC DE DIRAC EN F4CL    !
CARGU ! D4F4         ! TR ! <-> ! AMPLITUDE DU PIC DE DIRAC EN 1       !
CARGU !              !    !     ! (AIR PUR)                            !
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
C     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
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
      INCLUDE "parall.h"
      INCLUDE "pointe.h"
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
      INTEGER          INDPDF(NCELET)
C
      DOUBLE PRECISION FVAP(NCELET) , FHTF(NCELET) , F4P2M(NCELET)
      DOUBLE PRECISION F4M1(NCELET) , F4M2(NCELET) , D4CL(NCELET)
      DOUBLE PRECISION D4F4(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL
C
C   Variables pour le support de la Pdf et sa description
C
      DOUBLE PRECISION T1,T2
      DOUBLE PRECISION F4M,F42M,F1CL,F3CL,F4CL,F1M,F3M
C
      DOUBLE PRECISION D2PD1,D2MD1
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATIONS
C=======================================================================
C
C     Bornes pour le passage par les pdf T1 sur var, T2 sur moy
C
      T1 = 1.D-4
      T2 = 5.D-3
C
C --- Initialisation des tableaux
C
      DO IEL = 1, NCEL
        F4M1(IEL)   = 0.D0
        F4M2(IEL)   = 0.D0
        D4CL(IEL)   = 0.D0
        D4F4(IEL)  = 0.D0
C
        INDPDF(IEL) = 0
      ENDDO
C
C=======================================================================
C 2. TYPE DE PDF
C    INDPDF = 0  : PAS DE PDF
C    INDPDF = 1  : Passage par les PDF
C
C    CALCULS DES PARAMETRES DE LA PDF : F4M1,F4M2,D4CL,D4F4
C    CALCUL DE F4I7
C
C
C=======================================================================
C
      DO IEL = 1, NCEL
C
C       Traceur virtuel au point moyen
        F1M =  FVAP(IEL)
        F3M =  FHTF(IEL)/FF3MAX
        F4M = 1.D0 - F1M - F3M
C
C       Calcul des caractéristiques du point correspondant au
C       combustible moyen
C       F3cl : fraction de masse provenant de F3max et non de F3
C
        F1CL = F1M*FF3MAX/(F3M+F1M*FF3MAX)
        F3CL = 1.D0-F1CL
        F4CL = (1.D0-FF3MAX)*F3CL
C
C       L'intervalle permis pour F4 est [F4cl , 1 ]
C       La valeur maximale de la variance est donc
C       F4p2max = (F4m-F4cl)*(1-F4m)
C       On ne passe par les PDF que si la variance représente
C       une fraction significative de son maximum
C
        IF (  F4M .GT. ( F4CL + T2 ) .AND.
     &       F4M .LT. ( 1.D0 - T2 ) .AND.
     &       F4P2M(IEL) .GT. (T1*(F4M-F4CL)*(1.D0-F4M)) ) THEN
          INDPDF(IEL) = 1
        ELSE
          INDPDF(IEL) = 0
        ENDIF
C
        IF ( INDPDF(IEL) .EQ.1  )THEN
C
C        Calcul préliminaire pour une pdf rectangle et pics de Dirac
C        sur une droite passant par l'entrée d'air (F4) et le
C        combsutible moyen local Fcl
C        Soit F4m1 la borne Min
C             F4m2 la borne Max
C             D4cl amplitude du pic de Dirac en F4cl
C             D4f4 amplitude du pic de Dirac en 1 (air pur)
C
          F42M = F4M**2 + F4P2M(IEL)
C
C         rectangle seul
C
          F4M1(IEL) = F4M - SQRT(3.D0*F4P2M(IEL))
          F4M2(IEL) = F4M + SQRT(3.D0*F4P2M(IEL))
          D4CL(IEL) = ZERO
          D4F4(IEL) = ZERO
C
          IF ( F4M1(IEL).LE.F4cl .OR. F4M2(IEL).GE. 1.D0) THEN
C
C         pic en Fcl
C
            F4M1(IEL) = F4cl
            F4M2(IEL) = (3.D0*F42M + F4CL**2 - 4.D0*F4CL*F4M)
     &                 /(2.D0*(F4M-F4CL))
            D4CL(IEL) = (F4CL+F4M2(IEL)-2.D0*F4M)
     &                 /(F4M2(IEL)-F4CL)
C
            IF (F4M2(IEL).GE.1.D0 .OR. D4CL(IEL).LE.ZERO) THEN
C
C           pic en F4
C
              F4M2(IEL) = 1.D0
              F4M1(IEL) = (3.D0*F42M + 1.D0 - 4.D0*F4M)
     &                   /(2.D0*(F4M-1.D0))
              D4F4(IEL) = (2.D0*F4M-1.D0-F4M1(IEL))
     &                   /(1.D0-F4M1(IEL))
              D4CL(IEL) = 0.D0
C
              IF (F4M1(IEL).LE.F4CL .OR. D4F4(IEL).LE.ZERO) THEN
C
                F4M1(IEL) = F4CL
                F4M2(IEL) = 1.D0
                D2PD1 = (2.D0*F4M-1.D0-F4CL)/(1.D0-F4CL)
                D2MD1 = ( 6.D0*(F42M-F4M*(1.D0+F4CL))
     &                  + 1.d0+F4CL**2+4.D0*F4CL)
     &                 / (1.d0-F4CL)**2
                D4F4(IEL) = 0.5D0*(D2PD1+D2MD1)
                D4CL(IEL) = 0.5D0*(D2PD1-D2MD1)
              ENDIF
C
            ENDIF
C
          ENDIF
C
        ENDIF
C
      ENDDO
C
C----
C FIN
C----
C
      RETURN
      END
c@z
