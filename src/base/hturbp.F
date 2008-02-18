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
                        SUBROUTINE HTURBP
C                       *****************
C     -------------------------------------------------------------
     & ( PRL    , PRT    , CKARM  , YPLUS  , HTUR   )
C     -------------------------------------------------------------
C***********************************************************************
C
C FONCTION :
C --------
c@foncb
CFONC 1) CALCUL DU COEFFICIENT CORRECTEUR DU COEFFICIENT D'ECHANGE
CFONC    ENTRE LE FLUIDE ET LA PAROI POUR UN ECOULEMENT TURBULENT
CFONC    EN FONCTION DE LA DISTANCE ADIMENSIONELLE YPLUS = USTAR*DP/RNU
CFONC    HTUR = PR*YPLUS/TPLUS
CFONC
CFONC
CFONC CE COEFFICIENT EST CALCULE A L'AIDE D'UN MODELE DE SIMILITUDE
CFONC ENTRE COUCHE LIMITE DYNAMIQUE ET COUCHE LIMITE THERMIQUE
CFONC
CFONC LE  TPLUS EST CALCULE :
CFONC
CFONC  - POUR UN NOMBRE DE PRANDTL  << 0.1 (METAUX LIQUIDES) :
CFONC    PAR LE MODELE STANDARD A DEUX COUCHES (PRANDTL-TAYLOR)
CFONC
CFONC  - POUR UN NOMBRE DE PRANDTL  >> 0.1  (LIQUIDES et GAZ):
CFONC    PAR UN MODELE A TROIS COUCHES (ARPACI-LARSEN)
CFONC
CFONC  -->> LE COEFFICIENT D'ECHANGE FINAL : H = (K/dp)*htur
CFONC
c@fonce
C-----------------------------------------------------------------------
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! YPLUS        ! R  !  -> ! DISTANCE A LA PAROI ADIMENSIONNELLE  !
CARGU ! CKARM        ! R  !  -> ! CONSTANTE DE KARMAN                  !
CARGU ! PRT          ! R  !  -> ! NOMBRE DE PRANDTL TURBULENT          !
CARGU ! PRL          ! R  !  -> ! NOMBRE DE PRANDTL MOLECULAIRE        !
CARGU ! HTUR         ! R  ! <-  ! COEFFICIENT CORRECTEUR D'ECHANGE(ADIM!
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
C-----------------------------------------------------------------------
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION HTUR
      DOUBLE PRECISION PRL,CKARM,PRT,YPLUS
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION TPLUS
      DOUBLE PRECISION BETA2,A2
      DOUBLE PRECISION YP0,YP1,YP2
      DOUBLE PRECISION PRLM1
C
C*********************************************************************
C
C     1)INITIALISATIONS
C     -----------------
C
      HTUR = 1.D0
C
      PRLM1 = 0.1D0
C
      YP0   = PRT/(PRL*CKARM)
      YP2   = CKARM*1000.D0/PRT
      YP2   = SQRT(YP2)
      YP1   = (1000.D0/PRL)**(1.D0/3.D0)
C
C     ====================================================
C     2) CALCUL DU COEFFICIENT CORRECTEUR
C        POUR LES NOMBRES DE PRANDTL TRES PETITS
C     ====================================================
      IF( PRL .LE. PRLM1) THEN
        IF(YPLUS .GT. YP0) THEN
          TPLUS = PRL*YP0 + PRT/CKARM * LOG(YPLUS/YP0)
          HTUR = PRL*YPLUS/TPLUS
        ENDIF
C
      ENDIF
C
C
C     ====================================================
C     3) CALCUL DU COEFFICIENT CORRECTEUR
C        POUR UN MODELE A TROIS COUCHES
C     ====================================================
      IF( PRL .GT. PRLM1) THEN
C
       A2 = 15.D0*(PRL**(2.D0/3.D0))
       BETA2 = A2 - 500.D0/ (YP2**2)
C
       IF( (YPLUS .GE. YP1) .AND. (YPLUS.LT.YP2) )THEN
          TPLUS = A2 - 500.D0/(YPLUS*YPLUS)
          HTUR = PRL*YPLUS/TPLUS
       ENDIF
C
       IF( (YPLUS .GE. YP2) )THEN
          TPLUS = BETA2 + PRT/CKARM*LOG(YPLUS/YP2)
          HTUR = PRL*YPLUS/TPLUS
       ENDIF
C
      ENDIF
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
