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
                        SUBROUTINE MAJGEO
C                       *****************
C     -------------------------------------------------------------
     & ( NCELE2 , NFAC2  , NFABO2 , NFACG2 , NFBRG2 )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ---------
c@foncb
CFONC
CFONC PASSAGE DE NCELET DU C AU FORTRAN.
CFONC
CFONC PASSAGE DE NFAC ET NFABOR DU C AU FORTRAN AU CAS OU CES
CFONC VALEURS SONT MODIFIEES PAR LE DECOUPAGE DES FACES GAUCHES.
CFONC
CFONC PASSAGE DES DIMENSIONS GLOBALES DES FACES DE BORD ET INTERNES
CFONC POUR ASSURER LA COHERENCE DES VALEURS EN MONOPROC CAR EN CAS
CFONC DE PARALLELISME CETTE TACHE EST ASSUREE PAR PARGEO.
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                 ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELE2       ! E  !  -> ! NOUVEAU NOMBRE D'ELEMENTS            !
CARGU !              !    !     ! HALO COMPRIS                         !
CARGU ! NFAC2        ! E  !  -> ! NOUVEAU NOMBRE DE FACES INTERNES     !
CARGU ! NFABO2       ! E  !  -> ! NOUVEAU NOMBRE DE FACES DE BORD      !
CARGU ! NFACG2       ! E  !  -> ! NOUVEAU NOMBRE DE FACES INTERNES GLB.!
CARGU ! NFBRG2       ! E  !  -> ! NOUVEAU NOMBRE DE FACES DE BORD GLOB.!
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
      INCLUDE "parall.h"
C
C ARGUMENTS
C
      INTEGER          NCELE2, NFAC2  , NFABO2 , NFACG2 , NFBRG2
C
C
C***********************************************************************
C
C=======================================================================
C 1. MISE A JOUR DU NOMBRE DE CELLULES ETENDUES
C=======================================================================
C
      NCELET = NCELE2
C
C=======================================================================
C 2. MISE A JOUR DU NOMBRE DES FACES
C=======================================================================
C
      NFAC = NFAC2
      NFABOR = NFABO2
C
C     On remplit maintenant NDIMFB
      IF (NFABOR.EQ.0) THEN
        NDIMFB = 1
      ELSE
        NDIMFB = NFABOR
      ENDIF
C
C=======================================================================
C 3. MISE A JOUR DU NOMBRE DES FACES GLOBALES
C=======================================================================
C
      NFACGB = NFACG2
      NFBRGB = NFBRG2
C
      RETURN
      END
c@z
