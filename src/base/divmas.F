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
                        SUBROUTINE DIVMAS
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NFAC   , NFABOR ,
     &   INIT   , NFECRA ,
     &   IFACEL , IFABOR ,
     &   FLUMAS , FLUMAB ,
     &   DIVERG )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC INTEGRATION DU FLUX DE MASSE SUR LES CELLULES
CFONC
CFONC   .    .       --  .
CFONC   m =  m     + \   m
CFONC    i    i      /__  ij
CFONC                   j
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
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! INIT         ! E  !  -> ! INDICATEUR > 0 REMISE A 0 DE DIVERG  !
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! IFABOR(NFABOR! TE !  -> ! No DE L'ELT VOISIN D'UNE FACE DE BORD!
CARGU ! FLUMAS(NFAC) ! TR !  -> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR !  -> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! DIVERG(NCELET! TR ! <-> ! DIVERGENCE DE FLUMAS FLUMAB          !
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
      INCLUDE "vector.h"
C
C***********************************************************************
C
C
      INTEGER          NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          INIT   , NFECRA
C
      INTEGER          IFACEL(2,NFAC), IFABOR(NFABOR)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION DIVERG(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          IEL, IFAC, II, JJ
C
C***********************************************************************
C
C=======================================================================
C 1.  INITIALISATION
C=======================================================================
C
      IF( INIT.GE.1 ) THEN
        DO IEL = 1, NCELET
          DIVERG(IEL) = 0.D0
        ENDDO
      ELSEIF( INIT.EQ.0.AND.NCELET.GT.NCEL ) THEN
        DO IEL = NCEL+1, NCELET
          DIVERG(IEL) = 0.D0
        ENDDO
      ELSEIF(INIT.NE.0) THEN
        WRITE(NFECRA,1000) INIT
        CALL CSEXIT (1)
      ENDIF
C
C
C=======================================================================
C 2.  INTEGRATION SUR LES FACETTES INTERNES
C=======================================================================
C
      IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          DIVERG(II) = DIVERG(II) +FLUMAS(IFAC)
          DIVERG(JJ) = DIVERG(JJ) -FLUMAS(IFAC)
        ENDDO
C
      ELSE
C
C VECTORISATION NON FORCEE
        DO IFAC = 1, NFAC
          II = IFACEL(1,IFAC)
          JJ = IFACEL(2,IFAC)
          DIVERG(II) = DIVERG(II) +FLUMAS(IFAC)
          DIVERG(JJ) = DIVERG(JJ) -FLUMAS(IFAC)
        ENDDO
C
      ENDIF
C
C
C=======================================================================
C 3.  INTEGRATION SUR LES FACETTES DE BORD
C=======================================================================
C
      IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          DIVERG(II) = DIVERG(II) +FLUMAB(IFAC)
        ENDDO
C
      ELSE
C
        DO IFAC = 1, NFABOR
          II = IFABOR(IFAC)
          DIVERG(II) = DIVERG(II) +FLUMAB(IFAC)
        ENDDO
C
      ENDIF
C
 1000 FORMAT('DIVMAS APPELE AVEC INIT = ',I10)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
