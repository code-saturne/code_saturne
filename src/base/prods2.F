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
C     Copyright (C) 1998-2007 EDF S.A., France
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
                        SUBROUTINE PRODS2
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , ISQRT  ,
     &   VA1    , VB1    , VA2    , VB2    ,
     &   VAVB1  , VAVB2  )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CALCUL SIMULTANE DE 2 PRODUITS SCALAIRES
CFONC                    ______
CFONC VAPVB = VA.VB OU \/ VA.VB  SI ISQRT=1
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! ISQRT        ! E  !  -> ! INDICATEUR = 1 POUR PRENDRE LA RACINE!
CARGU ! VA1(), VB1() ! TR !  -> ! PREMIERS VECTEURS A MULTIPLIER       !
CARGU ! VA2(), VB2() ! TR !  -> ! SECONDS  VECTEURS A MULTIPLIER       !
CARGU ! VAVB1        ! R  ! <-  ! PREMIER PRODUIT SCALAIRE             !
CARGU ! VAVB2        ! R  ! <-  ! SECOND  PRODUIT SCALAIRE             !
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
      INCLUDE "parall.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET,NCEL,ISQRT
      DOUBLE PRECISION VAVB1, VAVB2
      DOUBLE PRECISION VA1(NCELET),VB1(NCELET)
      DOUBLE PRECISION VA2(NCELET),VB2(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER IEL, NVAVB
      DOUBLE PRECISION VAVB(2)
C
#if defined(_CS_HAVE_BLAS)
      INTEGER INCX, INCY
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
#endif
C
C***********************************************************************
C
C
#if defined(_CS_HAVE_BLAS)
      INCX = 1
      INCY = 1
      VAVB(1) = DDOT(NCEL, VA1, INCX, VB1, INCY)
      VAVB(2) = DDOT(NCEL, VA2, INCX, VB2, INCY)
#else
      VAVB(1) = 0.D0
      VAVB(2) = 0.D0
      DO IEL=1,NCEL
        VAVB(1) = VAVB(1) + VA1(IEL)*VB1(IEL)
        VAVB(2) = VAVB(2) + VA2(IEL)*VB2(IEL)
      ENDDO
#endif
C
      IF (IRANGP.GE.0) THEN
        NVAVB = 2
        CALL PARRSM (NVAVB, VAVB)
C       ===========
      ENDIF
C
      VAVB1 = VAVB(1)
      VAVB2 = VAVB(2)
C
      IF (ISQRT.EQ.1) THEN
        VAVB1 = SQRT(VAVB1)
        VAVB2 = SQRT(VAVB2)
      ENDIF
C
C----
C FIN
C----
C
      RETURN
C
      END
C
c@z
