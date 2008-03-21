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
                        SUBROUTINE ELPROP
C                       *****************
C     -------------------------------------------------------------
     & ( IPROPP , IPPPST )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC      INIT DES POSITIONS DES VARIABLES D'ETAT POUR
CFONC                 LE MODULE ELECTRIQUE
CFONC          (DANS VECTEURS PROPCE, PROPFA, PROPFB)
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IPROPP       ! E  ! ->  ! NUMERO DE LA DERNIERE CASE UTLISEE   !
CARGU !              !    !     ! DANS IPPROC, IPPROB, IPPROF          !
CARGU ! IPPPST       ! E  ! <-> ! POINTEUR INDIQUANT LE RANG DE LA     !
CARGU !              !    !     !  DERNIERE GRANDEUR DEFINIE AUX       !
CARGU !              !    !     !  CELLULES (RTP,PROPCE...) POUR LE    !
CARGU !              !    !     !  POST TRAITEMENT                     !
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
      INCLUDE "paramx.h"
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER       IPROPP , IPPPST
C
C VARIABLES LOCALES
C
      INTEGER       IPROP, IDIMVE
C
C***********************************************************************
C=======================================================================
C 1. DEFINITION DES POINTEURS
C=======================================================================
C
C     Pointeurs dans propce (ca n'implique pas qu'on ne calcule pas
C     les variables non definies ici)
C
      IPROP = IPROPP
C
C ---> Temperature en K
C
      IPROP  = IPROP + 1
      ITEMP  = IPROP
C
C ---> Puissance volumique dissipee par effet Joule W/m3
C
      IPROP  = IPROP + 1
      IEFJOU = IPROP
C
C ---> Densite de courant electrique reelle A/m2
C
      DO IDIMVE = 1, NDIMVE
        IPROP        = IPROP + 1
        IDJR(IDIMVE) = IPROP
      ENDDO
C
C Variables specifiques Effet Joule
C =================================
C
      IF ( IPPMOD(IELJOU).EQ.4 ) THEN
C
C ---> Densite de courant electrique imaginaire A/m2
C
        DO IDIMVE = 1, NDIMVE
          IPROP        = IPROP + 1
          IDJI(IDIMVE) = IPROP
        ENDDO
C
      ENDIF

C
C Variables specifiques Arc Electrique
C ====================================
C
      IF ( IPPMOD(IELARC).GE.1 ) THEN
C
C ---> Forces electromagnetiques de Laplace en N/m3
C
        DO IDIMVE = 1, NDIMVE
          IPROP          = IPROP + 1
          ILAPLA(IDIMVE) = IPROP
        ENDDO
C
C ---> Puissance volumique rayonnee W/m3
C      ou coefficient d'absorption
C
        IF ( IXKABE .GT.0 ) THEN
          IPROP = IPROP + 1
          IDRAD = IPROP
        ENDIF
      ENDIF
C
C Variables specifiques Conduction Ionique
C ========================================
C
      IF ( IPPMOD(IELION).GE.1 ) THEN
C
C ---> Charge electrique volumique C/m3
C
        IPROP  = IPROP + 1
        IQELEC = IPROP
C
      ENDIF
C
C ----  Nb de variables algebriques (ou d'etat)
C         propre a la physique particuliere NSALPP
C         total NSALTO
C
      NSALPP = IPROP - IPROPP
      NSALTO = IPROP
C
C ----  On renvoie IPROPP au cas ou d'autres proprietes devraient
C         etre numerotees ensuite
C
      IPROPP = IPROP
C
C=======================================================================
C 2. POSITIONNEMENT DES PROPRIETES : PROPCE, PROPFB, PROPFA
C=======================================================================
C
C ---> Positionnement dans le tableau PROPCE
C
      IPROP         = NPROCE
C
      IPROP         = IPROP + 1
      IPPROC(ITEMP) = IPROP
      IPPPST        = IPPPST + 1
      IPPPRO(IPROP) = IPPPST
C
      IPROP          = IPROP + 1
      IPPROC(IEFJOU) = IPROP
      IPPPST         = IPPPST + 1
      IPPPRO(IPROP)  = IPPPST
C
      DO IDIMVE = 1, NDIMVE
        IPROP                = IPROP + 1
        IPPROC(IDJR(IDIMVE)) = IPROP
        IPPPST               = IPPPST + 1
        IPPPRO(IPROP)        = IPPPST
      ENDDO
C
      IF ( IPPMOD(IELJOU).EQ.4 ) THEN
C
C ---> Densite de courant electrique imaginaire A/m2
C
        DO IDIMVE = 1, NDIMVE
          IPROP                = IPROP + 1
          IPPROC(IDJI(IDIMVE)) = IPROP
          IPPPST               = IPPPST + 1
          IPPPRO(IPROP)        = IPPPST
        ENDDO
C
      ENDIF
C
      IF ( IPPMOD(IELARC).GE.1 ) THEN
C
        DO IDIMVE = 1, NDIMVE
          IPROP                  = IPROP + 1
          IPPROC(ILAPLA(IDIMVE)) = IPROP
          IPPPST                 = IPPPST + 1
          IPPPRO(IPROP)          = IPPPST
        ENDDO
C
        IF ( IXKABE .GT. 0 ) THEN
          IPROP          = IPROP + 1
          IPPROC(IDRAD)  = IPROP
          IPPPST         = IPPPST + 1
          IPPPRO(IPROP)  = IPPPST
        ENDIF
C
      ENDIF
C
      IF ( IPPMOD(IELION).GE.1 ) THEN
C
C ---> Charge electrique volumique C/m3
C
        IPROP          = IPROP + 1
        IPPROC(IQELEC) = IPROP
        IPPPST         = IPPPST + 1
        IPPPRO(IPROP)  = IPPPST
C
      ENDIF
C
      NPROCE = IPROP
C
C
C=======================================================================
C 3. POSITIONNEMENT DES PROPRIETES : PROPFB, PROPFA
C=======================================================================
C
C ---> Positionnement dans le tableau PROPFB
C      Au centre des faces de bord
C
        IPROP = NPROFB
        NPROFB = IPROP
C
C
C ---> Positionnement dans le tableau PROPFA
C      Au centre des faces internes (flux de masse)
C
       IPROP = NPROFA
       NPROFA = IPROP
C
      RETURN
      END
c@z
