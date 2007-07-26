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
                        SUBROUTINE VORLGV
C                       *****************
C     -------------------------------------------------------------
     & ( NCEVOR , IENT   , DTREF  ,
     &   YZCEL  , XU     , XV     , XW     )
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  GENRATION DES FLUCUTUATIONS DE VITESSE DANS LA DIRECTION
CFONC  PRINCIPALE A PARTIR DE L'EQUATION DE LANGEVIN
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCEVOR       ! E  !  -> ! NOMBRE DE FACE A L'ENTREE OU EST     !
CARGU !              !    !     ! UTILISE LA METHODE                   !
CARGU ! IENT         ! E  !  -> ! NUMERO DE L'ENTREE                   !
CARGU ! DTREF        ! R  !  -> ! PAS DE TEMPS                         !
CARGU ! YZCEL        ! TR !  -> ! COORDONNEES DES FACES D'ENTREE DANS  !
CARGU !   (ICVMAX ,2)!    !     ! LE REFERENTIEL LOCAL                 !
CARGU ! XU(ICVMAX)   ! TR ! <-> ! COMPOSANTE DE VITESSE PRINCIPALE     !
CARGU ! XV(ICVMAX)   ! TR !  -> ! COMPOSANTES DE VITESSE TRANSVERSES   !
CARGU ! XW(ICVMAX)   ! TR !  -> !                                      !
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
      INCLUDE "entsor.h"
      INCLUDE "vortex.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCEVOR, IENT
C
      DOUBLE PRECISION DTREF
      DOUBLE PRECISION YZCEL(ICVMAX ,2)
      DOUBLE PRECISION XU(ICVMAX ), XV(ICVMAX ), XW(ICVMAX )
C
C     VARIABLES LOCALES
C
      INTEGER          II, III, IUN
      DOUBLE PRECISION DW1(1)
C
      DOUBLE PRECISION CST1, CST2
      PARAMETER       (CST1 = 1.8D0)
      PARAMETER       (CST2 = 0.6D0)
C
      DOUBLE PRECISION SINTH, COSTH
      DOUBLE PRECISION NORME, UFLUC, EK_VOR, EE_VOR, U_VOR, DU_VOR
      DOUBLE PRECISION PHIDAT, VFLUC, YY, ZZ
C
C=======================================================================
C 1. CALCUL DES FLUCTUATIONS DE VITESSE (SELON LA DIRECTION X)
C=======================================================================
C
      IUN = 1
C
      DO II = 1, NCEVOR
        YY = YZCEL(II,1)
        ZZ = YZCEL(II,2)
C
        III = 0
        U_VOR = PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &          YDAT(1,IENT),ZDAT(1,IENT),UDAT(1,IENT),III)
C
        IF(ICAS(IENT).EQ.2) THEN
          DU_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &              YDAT(1,IENT),ZDAT(1,IENT),DUDAT(1,IENT),III)
          EK_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &              YDAT(1,IENT),ZDAT(1,IENT),KDAT(1,IENT),III)
          EE_VOR =  PHIDAT(NFECRA,ICAS(IENT),NDAT(IENT),YY,ZZ,
     &              YDAT(1,IENT),ZDAT(1,IENT),EPSDAT(1,IENT),III)
C
          UFLUC = XU(II) - U_VOR
C
          NORME = SQRT(YZCEL(II,1)**2+YZCEL(II,2)**2)
          COSTH = YZCEL(II,1) / NORME
          SINTH = YZCEL(II,2) / NORME
          VFLUC = - COSTH*XV(II) - SINTH*XW(II)
C
C le signe - vient du fait que l'on veut
C la fluctuation dans la direction normale
C a la paroi
C
          CALL NORMALEN(IUN,DW1(1))
          UFLUC = ( UFLUC -
     &            (1.D0-2.D0/3.D0*CST2)*DU_VOR*VFLUC*DTREF +
     &             2.D0* SQRT(2.D0/3.D0*(CST1-1.D0)*EE_VOR*DTREF)
     &            *DW1(1))
          UFLUC = UFLUC/(1.D0+(0.5D0*CST1*DTREF*EE_VOR
     &                 /EK_VOR))
          XU(II) = U_VOR + UFLUC
        ELSE
          XU(II) = U_VOR
        ENDIF
      ENDDO
C
C ---
C FIN
C ---
C
      RETURN
C
      END
c@z
