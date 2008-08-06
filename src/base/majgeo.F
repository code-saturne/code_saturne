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
                        SUBROUTINE MAJGEO
C                       *****************
C     -------------------------------------------------------------
     & ( NCEL2  , NCELE2 , NFAC2  , NFABO2 , NSOM2 ,
     &   LNDFA2 , LNDFB2 , NCELG2 , NFACG2 , NFBRG2 , NSOMG2 )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC PASSAGE DES DIMENSIONS DU MAILLAGE DU C AU FORTRAN.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                 ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCEL2        ! E  !  -> ! NOMBRE DE CELLULES                   !
CARGU ! NCELE2       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NFAC2        ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABO2       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NSOM2        ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFA2       ! E  !  -> ! TAILLE DE LNDFAC                     !
CARGU ! LNDFB2       ! E  !  -> ! TAILLE DE LNDFBR                     !
CARGU ! NCELG2       ! E  !  -> ! NOMBRE GLOBAL DE CELLULES            !
CARGU ! NFACG2       ! E  !  -> ! NOMBRE GLOBAL DE FACES INTERNES      !
CARGU ! NFBRG2       ! E  !  -> ! NOMBRE GLOBAL DE FACES DE BORD       !
CARGU ! NSOMG2       ! E  !  -> ! NOMBRE GLOBAL DE SOMMETS             !
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
      INCLUDE "dimens.h"
      INCLUDE "dimfbr.h"
      INCLUDE "paramx.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
C
C ARGUMENTS
C
      INTEGER          NCEL2, NCELE2, NFAC2, NFABO2, NSOM2
      INTEGER          LNDFA2, LNDFB2
      INTEGER          NCELG2, NFACG2 , NFBRG2, NSOMG2
C
C
C***********************************************************************
C
C=======================================================================
C 1. MISE A JOUR DU NOMBRE DE CELLULES
C=======================================================================
C
      NCEL = NCEL2
      NCELET = NCELE2
C
C=======================================================================
C 2. MISE A JOUR DU NOMBRE DES FACES
C=======================================================================
C
      NFAC = NFAC2
      NFABOR = NFABO2
C
      LNDFAC = LNDFA2
      LNDFBR = LNDFB2
C
C     On remplit maintenant NDIMFB
      IF (NFABOR.EQ.0) THEN
        NDIMFB = 1
      ELSE
        NDIMFB = NFABOR
      ENDIF
C
C=======================================================================
C 3. MISE A JOUR DU NOMBRE DES SOMMETS
C=======================================================================
C
      NNOD = NSOM2
C
C=======================================================================
C 4. MISE A JOUR DES TAILLES GLOBALES
C=======================================================================
C
      NCELGB = NCELG2
      NFACGB = NFACG2
      NFBRGB = NFBRG2
      NSOMGB = NSOMG2
C
      RETURN
      END
c@z
