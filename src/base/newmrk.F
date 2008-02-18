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
                        SUBROUTINE NEWMRK
C                       *****************
C     -------------------------------------------------------------
     & (ISTR  , ALPNMK, BETNMK, GAMNMK,
     &  XM    , XC    , XK    , XN0   , XN    , XPN   , XPPN  , XNM1  ,
     &  XPNM1 , XPPNM1, XFN   , XFNM1 , DT    )

C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  RESOLUTION PAR LA METHODE DE NEWMARK HHT D'UNE EQUATION
CFONC    DIFFERENTIELLE LINEAIRE DE SECOND ORDRE DU TYPE
CFONC
CFONC     M.X'' + C.X' + K.(X+X0) = F
CFONC
CFONC  X EST UN CHAMP DE VECTEUR TRIDIMENSIONNEL ET M, C ET K
CFONC    SONT DES MATRICES 3x3 QUELCONQUES
CFONC
CFONC  LA RESOLUTION EST FAITE EN DEPLACEMENT
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! ISTR         ! E  !  -> ! NUMERO DE STRUCTURE                  !
CARGU ! XM(3,3)      ! TR !  -> ! MATRICE DE MASSE DU SYSTEME          !
CARGU ! XC(3,3)      ! TR !  -> ! MATRICE DE FRICTION DU SYSTEME       !
CARGU ! XK(3,3)      ! TR !  -> ! MATRICE DE RAIDEUR DU SYSTEME        !
CARGU ! XN0(3)       ! TR !  -> ! CHAMP DE DEPLACEMENT INITIAL         !
CARGU ! XN(3)        ! TR ! <-  ! CHAMP DE DEPLACEMENT AU TEMPS N      !
CARGU ! XNM1(3)      ! TR !  -> ! CHAMP DE DEPLACEMENT AU TEMPS N-1    !
CARGU ! XPN(3)       ! TR ! <-  ! CHAMP DE VITESSE AU TEMPS N          !
CARGU ! XPNM1(3)     ! TR !  -> ! CHAMP DE VITESSE AU TEMPS N-1        !
CARGU ! XPPN(3)      ! TR ! <-  ! CHAMP D'ACCELERATION AU TEMPS N      !
CARGU ! XPPNM1(3)    ! TR !  -> ! CHAMP D'ACCELERATION AU TEMPS N-1    !
CARGU ! XFTLD(3)     ! TR !  -> ! CHAMP DE FORCE AU TEMPS INTERMEDIAIRE!
CARGU ! DT           ! R  !  -> ! PAS DE TEMPS                         !
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
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "optcal.h"
      INCLUDE "entsor.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          ISTR
C
      DOUBLE PRECISION ALPNMK, BETNMK, GAMNMK
      DOUBLE PRECISION XM(3,3),XC(3,3),XK(3,3)
      DOUBLE PRECISION XN0(3),XN(3),XPN(3),XPPN(3)
      DOUBLE PRECISION XNM1(3),XPNM1(3),XPPNM1(3)
      DOUBLE PRECISION XFN(3),XFNM1(3)
      DOUBLE PRECISION DT
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ
C
      DOUBLE PRECISION A(3,3), B1(3), B2(3), B(3)
      DOUBLE PRECISION A0, A1, A2, A3, A4, A5, A6, A7
      DOUBLE PRECISION DET, DET1, DET2, DET3, EPSDET
C
C***********************************************************************
C
C     Critere de nullite du determinant
      EPSDET = 1.D-12
C
C     Ceofficients des equations
      A0 = 1.D0/BETNMK/DT**2
      A1 = (1.D0+ALPNMK)*GAMNMK/BETNMK/DT
      A2 = 1.D0/BETNMK/DT
      A3 = 1.D0/2.D0/BETNMK - 1.D0
      A4 = (1.D0+ALPNMK)*GAMNMK/BETNMK - 1.D0
      A5 = (1.D0+ALPNMK)*DT*(GAMNMK/2.D0/BETNMK - 1.D0)
      A6 = DT*(1.D0 - GAMNMK)
      A7 = GAMNMK*DT
C
      DO II = 1, 3
        DO JJ = 1, 3
          A(II,JJ) = (1.D0+ALPNMK)*XK(II,JJ)
     &         + A1*XC(II,JJ) + A0*XM(II,JJ)
        ENDDO
        B(II)  = (1.D0+ALPNMK)*XFN(II) - ALPNMK*XFNM1(II)
        B1(II) = A0*XNM1(II) + A2*XPNM1(II) + A3*XPPNM1(II)
        B2(II) = A1*XNM1(II) + A4*XPNM1(II) + A5*XPPNM1(II)
      ENDDO
C
      DO II = 1, 3
        DO JJ = 1, 3
          B(II) = B(II)
     &         + XM(II,JJ)*B1(JJ) + XC(II,JJ)*B2(JJ)
     &         + XK(II,JJ)*( ALPNMK*XNM1(JJ) + XN0(JJ) )
        ENDDO
      ENDDO
C
      DET = A(1,1)*A(2,2)*A(3,3)
     &    + A(2,1)*A(3,2)*A(1,3)
     &    + A(3,1)*A(1,2)*A(2,3)
     &    - A(3,1)*A(2,2)*A(1,3)
     &    - A(2,1)*A(1,2)*A(3,3)
     &    - A(1,1)*A(3,2)*A(2,3)
C
      IF (DABS(DET).LT.EPSDET) THEN
        WRITE(NFECRA,1000) ISTR,DABS(DET),EPSDET
        NTMABS = NTCABS
      ENDIF
C
      DET1 = B(1)*A(2,2)*A(3,3)
     &     + B(2)*A(3,2)*A(1,3)
     &     + B(3)*A(1,2)*A(2,3)
     &     - B(3)*A(2,2)*A(1,3)
     &     - B(2)*A(1,2)*A(3,3)
     &     - B(1)*A(3,2)*A(2,3)
C
      DET2 = A(1,1)*B(2)*A(3,3)
     &     + A(2,1)*B(3)*A(1,3)
     &     + A(3,1)*B(1)*A(2,3)
     &     - A(3,1)*B(2)*A(1,3)
     &     - A(2,1)*B(1)*A(3,3)
     &     - A(1,1)*B(3)*A(2,3)
C
      DET3 = A(1,1)*A(2,2)*B(3)
     &     + A(2,1)*A(3,2)*B(1)
     &     + A(3,1)*A(1,2)*B(2)
     &     - A(3,1)*A(2,2)*B(1)
     &     - A(2,1)*A(1,2)*B(3)
     &     - A(1,1)*A(3,2)*B(2)
C
      XN(1) = DET1/DET
      XN(2) = DET2/DET
      XN(3) = DET3/DET
C
      DO II = 1, 3
        XPPN(II) = A0*(XN(II)-XNM1(II)) - A2*XPNM1(II)  - A3*XPPNM1(II)
        XPN(II)  = XPNM1(II)            + A6*XPPNM1(II) + A7*XPPN(II)
      ENDDO
C
C----
C FORMATS
C----
C
 1000 FORMAT (
     &'@                                                            ',/,
     &'@ @@ ATTENTION : DEPLACEMENT DE STRUCTURES INTERNES ALE      ',/,
     &'@    *********                                               ',/,
     &'@  Structure : ',I10                                          ,/,
     &'@  La valeur absolue du determinant de la matrice de         ',/,
     &'@    deplacement vaut : ',E14.5                               ,/,
     &'@  La matrice est consideree comme non inversible            ',/,
     &'@    (valeur limite fixee a ',E14.5     ,')                  ',/,
     &'@                                                            ',/,
     &'@  Arret du calcul                                           ',/,
     &'@                                                            '  )
C
C----
C FIN
C----
C
      END
c@z
