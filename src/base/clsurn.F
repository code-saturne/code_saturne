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
                        SUBROUTINE CLSURN
C                       *****************
C     -------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NFAC   , NFABOR ,
     &   SURFAC , SURFBO ,
     &   SURFAN , SURFBN ,
     &   IA     , RA     )
C     -------------------------------------------------------------
C***********************************************************************
C
C  FONCTION  :
C  ----------
c@foncb
CFONC
CFONC             CALCUL DES SURFACES
CFONC
c@fonce
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0/IDBRA0! E  !  -> ! POINTEUR DE LA PREMIERE CAS LIBRE DES!
CARGU !              !    !     !  TABLEAUX IA/RA                      !
CARGU ! NFAC  /NFABOR! E  !  -> ! NOMBRE TOTAL DE FACES INTERNES/DE BRD!
CARGU ! SURFAC       ! TR !  -> ! COORDS DU VECTEUR SURFACE DES NFAGLO !
CARGU ! (3,NFAC  )   !    !     !  FACES INTERNES ; DIRIGE DU VOIS 1   !
CARGU !              !    !     !  VERS LE VOISIN 2 (IFACEL)           !
CARGU !              !    !     !  NON UNITAIRE                        !
CARGU ! SURFBO       ! TR !  -> ! COORDS DU VECTEUR SURFACE DES NFAGBR !
CARGU ! (3,NFABOR)   !    !     !  FACES DE BORD ; DIRIGE VERS L'      !
CARGU !              !    !     !  EXTERIEUR DU DOMAINE ; NON UNITAIRE !
CARGU ! SURFAN       ! TR !  -> ! SURFACE DES FACES INTERNES           !
CARGU ! (NFAC    )   !    !     !     (NORME DE SURFAC)                !
CARGU ! SURFBN       ! TR !  -> ! SURFACE DES FACES DE BORD            !
CARGU ! (NFABOR  )   !    !     !     (NORME DE SURFBO)                !
CARGU ! IA           ! TE !  -  ! TABLEAU DE TRAVAIL ENTIER            !
CARGU ! RA           ! TR !  -  ! TABLEAU DE TRAVAIL REEL              !
CARGU .______________. __ ! ____.______________________________________.
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
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER IDBIA0,IDBRA0
      INTEGER NFAC,NFABOR
      INTEGER IA(*)
      DOUBLE PRECISION SURFAC(3,NFAC  ),SURFBO(3,NFABOR)
      DOUBLE PRECISION SURFAN(NFAC    ),SURFBN(NFABOR  )
      DOUBLE PRECISION RA(*)
C
      INTEGER IFAC
      DOUBLE PRECISION SURFX,SURFY,SURFZ
C
C-----------------------------------------------------------------------
C 1. FACES INTERNES
C-----------------------------------------------------------------------
C
      DO IFAC = 1, NFAC
C
         SURFX = SURFAC(1,IFAC)
         SURFY = SURFAC(2,IFAC)
         SURFZ = SURFAC(3,IFAC)
         SURFAN(IFAC) = SQRT(SURFX**2 + SURFY**2 + SURFZ**2)
C
      ENDDO
C
C-----------------------------------------------------------------------
C 2. FACES DE BORD
C-----------------------------------------------------------------------
C
      DO IFAC = 1, NFABOR
C
         SURFX = SURFBO(1,IFAC)
         SURFY = SURFBO(2,IFAC)
         SURFZ = SURFBO(3,IFAC)
         SURFBN(IFAC) = SQRT(SURFX**2 + SURFY**2 + SURFZ**2)
C
      ENDDO
C
C-----------------------------------------------------------------------
C 3. FIN
C-----------------------------------------------------------------------
C
      RETURN
      END
c@z
