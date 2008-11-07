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
                        SUBROUTINE PRODSC
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , ISQRT  , VA     , VB     , VAVB   )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC                                     ______
CFONC PRODUIT SCALAIRE VAPVB = VA.VB OU \/ VA.VB  SI ISQRT=1
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
CARGU ! VA, VB(NCELET! TR !  -> ! VECTEURS A MULTIPLIER                !
CARGU ! VAVB         ! R  ! <-  ! PRODUIT SCALAIRE                     !
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
      DOUBLE PRECISION VAVB
      DOUBLE PRECISION VA(NCELET),VB(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER INCX, INCY
      DOUBLE PRECISION DDOT
      EXTERNAL         DDOT
C
C***********************************************************************
C
      INCX = 1
      INCY = 1
      VAVB = DDOT(NCEL, VA, INCX, VB, INCY)
C
      IF (IRANGP.GE.0) CALL PARSOM (VAVB)
C                      ===========
      IF( ISQRT.EQ.1 ) VAVB= SQRT(VAVB)
C----
C FIN
C----
C
      RETURN
C
      END
C
c@z
