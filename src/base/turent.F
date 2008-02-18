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
                        SUBROUTINE KEENDB
C                       *****************
C     --------------------------------------------------------------
     & ( UREF2, DH, XRHO, XMU , CMU, XKAPPA, USTAR2, XK, XEPS )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     CALCUL DE U*, K ET EPSILON A PARTIR D'UN DIAMETRE ET D'UNE
CFONC       VITESSE DEBITANTE POUR DES ECOULEMENTS EN CONDUITE
CFONC       CIRCULAIRE A PAROI LISSE
CFONC     -> UTILISE POUR LES CONDITIONS AUX LIMITES D'ENTREE
CFONC
CFONC     ON SORT A LA FOIS U* ET (XK,XEPS) POUR PERMETTRE
CFONC       A L'UTILISATEUR DE CALCULER DES XK ET XEPS DIFFERENTS
CFONC       AVEC LE U*, S'IL LE SOUHAITE
CFONC
CFONC
CFONC     ON UTILISE DES LOIS TIREES DE IDEL'CIK
CFONC     LE COEFFICIENT DE PERTE DE CHARGE XLMBDA EST DEFINI PAR
CFONC      |dP/dx| = XLMBDA/DH * 1/2*XRHO*UREF**2
CFONC
CFONC      PUIS U*=UREF*SQRT(XLMBDA/8)
CFONC
CFONC     POUR Re < 2000
CFONC       XLMBDA = 64/Re
CFONC
CFONC     POUR Re > 4000
CFONC       XLMBDA = 1/( 1.8*LOG10(Re)-1.64 )**2
CFONC
CFONC     POUR 2000 < Re < 4000, ON COMPLETE PAR UNE DROITE
CFONC       XLMBDA = 0.021377 + 5.3115D-6*Re
CFONC
CFONC
CFONC     A PARTIR DE U*, ON ESTIME XK ET XEPS A PARTIR DE FORMULES
CFONC       CLASSIQUES DE TURBULENCE DEVELOPPEE
CFONC
CFONC     XK   = USTAR2/SQRT(CMU)
CFONC     XEPS = USTAR2**1.5D0/(XKAPPA*DH*0.1D0)
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! UREF2        ! R  !  -> ! CARRE DE LA VITESSE DEBITANTE DE     !
CARGU !              !    !     !                            REFERENCE !
CARGU ! DH           ! R  !  -> ! DIAMETRE HYDRAULIQUE                 !
CARGU ! XRHO         ! R  !  -> ! MASSE VOLUMIQUE                      !
CARGU ! XMU          ! R  !  -> ! VISCOSITE DYNAMIQUE                  !
CARGU ! CMU          ! R  !  -> ! CONSTANTE CMU                        !
CARGU ! XKAPPA       ! R  !  -> ! CONSTANTE KAPPA                      !
CARGU ! USTAR2       ! R  ! <-  ! CARRE DE LA VITESSE DE FROTTEMENT    !
CARGU ! XK           ! R  ! <-  ! INTENSITE TURBULENTE CALCULEE        !
CARGU ! XEPS         ! R  ! <-  ! DISSIPATION TURBULENTE CALCULEE      !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION UREF2, DH, XRHO, XMU , USTAR2, XK, XEPS
      DOUBLE PRECISION CMU, XKAPPA
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION RE, XLMBDA
C
C***********************************************************************
C
      RE = SQRT(UREF2)*DH*XRHO/XMU
C
      IF (RE.LT.2000) THEN
C     dans ce cas on calcule directement u*^2 pour eviter un probleme
C      sur xlmbda=64/Re quand Re->0
C
        USTAR2 = 8.D0*XMU*SQRT(UREF2)/XRHO/DH
C
      ELSE IF (RE.LT.4000) THEN
C
        XLMBDA = 0.021377D0 + 5.3115D-6*RE
        USTAR2 = UREF2*XLMBDA/8.D0
C
      ELSE
C
        XLMBDA = 1/( 1.8D0*LOG(RE)/LOG(10.D0)-1.64D0)**2
        USTAR2 = UREF2*XLMBDA/8.D0
C
      ENDIF
C
      XK   = USTAR2/SQRT(CMU)
      XEPS = USTAR2**1.5D0/(XKAPPA*DH*0.1D0)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
c@a
                        SUBROUTINE KEENIN
C                       *****************
C     --------------------------------------------------------------
     & ( UREF2, XINTUR, DH, CMU, XKAPPA, XK, XEPS )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     CALCUL DE U*, K ET EPSILON A PARTIR D'UN DIAMETRE ET D'UNE
CFONC       INTENSITE TURBULENTE
CFONC     -> UTILISE POUR LES CONDITIONS AUX LIMITES D'ENTREE
CFONC
CFONC
CFONC     ON ESTIME XK ET XEPS A PARTIR DE FORMULES
CFONC       CLASSIQUES DE TURBULENCE DEVELOPPEE
CFONC
CFONC     XK   = USTAR2/SQRT(CMU)
CFONC     XEPS = USTAR2**1.5D0/(XKAPPA*DH*0.1D0)
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! UREF2        ! R  !  -> ! CARRE DE LA VITESSE DEBITANTE DE     !
CARGU !              !    !     !                            REFERENCE !
CARGU ! XINTUR       ! R  !  -> ! INTENSITE TURBULENTE                 !
CARGU ! DH           ! R  !  -> ! DIAMETRE HYDRAULIQUE                 !
CARGU ! CMU          ! R  !  -> ! CONSTANTE CMU                        !
CARGU ! XKAPPA       ! R  !  -> ! CONSTANTE KAPPA                      !
CARGU ! XK           ! R  ! <-  ! INTENSITE TURBULENTE CALCULEE        !
CARGU ! XEPS         ! R  ! <-  ! DISSIPATION TURBULENTE CALCULEE      !
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
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION UREF2, XINTUR, DH, CMU, XKAPPA, XK, XEPS
C
C VARIABLES LOCALES
C
C***********************************************************************
C
C
      XK   = 1.5D0*UREF2*XINTUR**2
      XEPS =
     &     10.D0*CMU**(0.75D0)*XK**1.5D0/(XKAPPA*DH)
C
C----
C FIN
C----
C
      RETURN
      END
c@z
