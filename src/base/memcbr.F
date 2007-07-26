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
                        SUBROUTINE MEMCBR
C                       *****************
C     -------------------------------------------------------------
     & ( IICELB , NCELET , NCEL   , NFABOR ,
     &   NCELBR , IFINIA ,
     &   IFABOR ,
     &   IA     )
C     -------------------------------------------------------------
C***********************************************************************
C
C  FONCTION  :
C  --------
c@foncb
CFONC             REPERAGE DES ELEMENTS AYANT AU MOINS UNE FACE DE BORD
CFONC             POUR LA RECONSTRUCTION DES GRADIENTS (RECALCUL DE
CFONC             COCG UNIQUEMENT AUX BORDS)
CFONC
CFONC           ON CALCULE NCELBR ET ON REMPLIT IA(IICELB)
CFONC
CFONC
c@fonce
C
C-----------------------------------------------------------------------
C                                  ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IICELB       ! E  !  -> ! POINTEUR DE LA PREMIERE CAS LIBRE DU !
CARGU !              !    !     !  TABLEAU  IA (ET POINTEUR SUR ICELBR)!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NCELBR       ! E  ! <-  ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! IFINIA       ! E  ! <-  ! IICELB+NCELBR DEBUT DE ZONE LIBRE    !
CARGU !              !    !     ! DANS IA EN SORTIE                    !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IA           ! TE !  -  ! TABLEAU DE TRAVAIL ENTIER            !
CARGU .______________.____._____.______________________________________.
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
      INTEGER          IICELB , NCELET , NCEL   , NFABOR
      INTEGER          NCELBR , IFINIA
      INTEGER          IFABOR(NFABOR)
      INTEGER          IA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IFAC, II, IEL , IIASSE
C
C***********************************************************************
C
C
C=======================================================================
C 1. MEMOIRE
C=======================================================================
C
C On surdimensionne ICELBR (NFABOR)
C     et un tableau d'assemblage NCELET
C
      IIASSE = IICELB + NFABOR
      IFINIA = IIASSE + NCELET
C
      CALL IASIZE('MEMCBR',IFINIA)
C     ===========
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      DO IEL = 1, NCELET
        IA(IIASSE+IEL-1) = 0
      ENDDO
C
C=======================================================================
C 2. REPERAGE DES ELEMENTS DE BORD
C=======================================================================
C
      DO IFAC = 1, NFABOR
        IEL = IFABOR(IFAC)
        IA(IIASSE+IEL-1) = IA(IIASSE+IEL-1) + 1
      ENDDO
C
C=======================================================================
C 3. REMPLISSAGE DE ICELBR
C=======================================================================
C
      II = 0
      DO IEL = 1, NCEL
        IF (IA(IIASSE+IEL-1).GT.0) THEN
          II = II + 1
          IA(IICELB+II-1) = IEL
        ENDIF
      ENDDO
C
      NCELBR = II
C
C=======================================================================
C 4. DESALLOCATION DE LA MEMOIRE EN TROP ET DU TABLEAU DE TRAVAIL IASSE
C=======================================================================
C
      IFINIA = IICELB + NCELBR
C
C=======================================================================
C 5. FIN
C=======================================================================
C
      RETURN
      END
c@z
