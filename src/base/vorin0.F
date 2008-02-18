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
                        SUBROUTINE VORIN0
C                       *****************
C     -------------------------------------------------------------
     & ( NFABOR , IREPVO )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  INITIALISATION PAR DEFAUT DES PARAMETRES POUR LA METHODE
CFONC    DES VORTEX AVANT DE PASSER LA MAIN A L'UTILISATEUR
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! IREPVO       ! TE ! <-  ! NUMERO DE L'ENTREE ASSOCIE A CHAQUE  !
CARGU !     (NFABOR) !    !     ! FACE DE BORD (=0 SI PAS DE VORTEX)   !
CARGU !______________!____!_____!______________________________________!
c@argue

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
      INCLUDE "vortex.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NFABOR, IREPVO(NFABOR)
C
C VARIABLES LOCALES
C
      INTEGER          II, JJ, IFAC
C
C***********************************************************************
C
C
C=======================================================================
C 1. PARAMETRES DE LA METHODE
C=======================================================================
C
C ---> Nombre d'entrees
C
      NNENT = -999
C
C ---> Nombre de vortex
C
      DO II = 1, NENTMX
        NVORT(II) = -999
      ENDDO
C
C ---> Defintion du cas
C
      DO II = 1, NENTMX
        ICAS(II)   = - 999
      ENDDO
C
C ---> Tableau de reperage et repere local
C
      DO IFAC = 1, NFABOR
        IREPVO(IFAC) = 0
      ENDDO
C
      DO II = 1, NENTMX
        DO JJ = 1, 3
          DIR1(JJ,II) = 0.D0
          DIR2(JJ,II) = 0.D0
          CEN(JJ,II)  = 0.D0
        ENDDO
      ENDDO
C
C ---> Conditions aux limites et grandeurs caracteristiques
C
      DO II = 1, NENTMX
        DO JJ = 1, 4
          ICLVOR(JJ,II) = -999
        ENDDO
        LLY(II) = -999.D0
        LLZ(II) = -999.D0
        LLD(II) = -999.D0
      ENDDO
C
C ---> Parametres physiques
C
      DO II = 1, NENTMX
        ITLIVO(II) = -999
        TLIMVO(II) = -999.D0
        ISGMVO(II) = -999
        XSGMVO(II) = -999.D0
        IDEPVO(II) = -999
        UD(II)     = 0.D0
      ENDDO
C
C ---> Donnees utilisateurs
C
      DO II = 1, NENTMX
        WRITE(FICVOR(II),'(1A6,I2.2)') 'vordat',II
        UDEBIT(II) = 0.D0
        KDEBIT(II) = -999.D0
        EDEBIT(II) = -999.D0
      ENDDO
C
      RETURN
      END
c@z
