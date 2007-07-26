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
                        SUBROUTINE FUTEH1
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   ,
     &   EH     ,
     &   FUEL1  , FUEL3  , OXYD   , PROD1  , PROD2  ,
     &   XINER  , XH2S   , XSO2   ,
     &   TP     ,
     &   EH0    , EH1    )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC CALCUL DE LA TEMPERATURE DU GAZ
CFONC   EN FONCTION DE L'ENTHALPIE DU GAZ ET DES CONCENTRATIONS
CFONC
c@fonce
C                             ARGUMENTS
c@argub
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! EH           ! TR !  -> ! ENTHALPIE DU GAZ                     !
CARGU !              !    !     ! (J/kg de melange gazeux)             !
CARGU ! FUEL1        ! TR !  -> ! FRACTION MASSIQUE CHx1               !
CARGU ! FUEL2        ! TR !  -> ! FRACTION MASSIQUE CHx2               !
CARGU ! FUEL3        ! TR !  -> ! FRACTION MASSIQUE CO                 !
CARGU ! OXYD         ! TR !  -> ! FRACTION MASSIQUE O2                 !
CARGU ! PROD1        ! TR !  -> ! FRACTION MASSIQUE CO2                !
CARGU ! PROD2        ! TR !  -> ! FRACTION MASSIQUE H2O                !
CARGU ! XINER        ! TR !  -> ! FRACTION MASSIQUE N2                 !
CARGU ! TP           ! TR ! <-  ! TEMPERATURE DU GAZ (Kelvin)          !
CARGU ! EH0          ! TR ! <-> ! TABLEAU REEL DE TRAVAIL              !
CARGU ! EH1          ! TR ! <-> ! TABLEAU REEL DE TRAVAIL              !
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
C     TYPE : E (ENTIER), R (REEL), A (ALPHAMNUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C**********************************************************************
C     DONNEES EN COMMON
C**********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL
C
      DOUBLE PRECISION EH(NCELET)
      DOUBLE PRECISION FUEL1(NCELET), FUEL3(NCELET)
      DOUBLE PRECISION OXYD(NCELET), XINER(NCELET)
      DOUBLE PRECISION PROD1(NCELET),PROD2(NCELET)
      DOUBLE PRECISION XH2S(NCELET) , XSO2(NCELET)
      DOUBLE PRECISION TP(NCELET)
      DOUBLE PRECISION EH0(NCELET), EH1(NCELET)
C
C VARIABLES LOCALES
C
      INTEGER          II, ICEL
C
      DOUBLE PRECISION YFOV0 , EH2OV, EH2OV1
      DOUBLE PRECISION DEN1 , EHCHX1
C
C***********************************************************************
C
      II = NPO-1
      DO ICEL = 1, NCEL
C
C --- Clipping eventuel de TP a TH(NPO) si EH > EH1
C
        EH1(ICEL) = FUEL1(ICEL)*EHGAZE(IFOV,II+1)
     &       + FUEL3(ICEL)*EHGAZE(ICO ,II+1)
     &       + OXYD(ICEL) *EHGAZE(IO2 ,II+1)
     &       + PROD1(ICEL)*EHGAZE(ICO2,II+1)
     &       + PROD2(ICEL)*EHGAZE(IH2O,II+1)
     &       + XINER(ICEL)*EHGAZE(IN2 ,II+1)
     &       + XH2S (ICEL)*EHGAZE(IH2S,II+1)
     &       + XSO2 (ICEL)*EHGAZE(ISO2,II+1)
C
        IF ( EH(ICEL).GE.EH1(ICEL) ) TP(ICEL)= TH(II+1)
      ENDDO
C
      II = 1
      DO ICEL = 1, NCEL
C
C --- Clipping eventuel de TP a TH(1) si EH < EH0
C
        EH0(ICEL) = FUEL1(ICEL)*EHGAZE(IFOV,II)
     &            + FUEL3(ICEL)*EHGAZE(ICO ,II)
     &            + OXYD(ICEL) *EHGAZE(IO2 ,II)
     &            + PROD1(ICEL)*EHGAZE(ICO2,II)
     &            + PROD2(ICEL)*EHGAZE(IH2O,II)
     &            + XINER(ICEL)*EHGAZE(IN2 ,II)
     &            + XH2S (ICEL)*EHGAZE(IH2S,II)
     &            + XSO2 (ICEL)*EHGAZE(ISO2,II)
C
        IF ( EH(ICEL).LE.EH0(ICEL) ) TP(ICEL)= TH(1)
      ENDDO
C
      DO II = 1, NPO-1
        DO ICEL = 1, NCEL
C
          EH0(ICEL) = FUEL1(ICEL)*EHGAZE(IFOV,II)
     &               +FUEL3(ICEL)*EHGAZE(ICO ,II)
     &               +OXYD(ICEL) *EHGAZE(IO2 ,II)
     &               +PROD1(ICEL)*EHGAZE(ICO2,II)
     &               +PROD2(ICEL)*EHGAZE(IH2O,II)
     &               +XINER(ICEL)*EHGAZE(IN2 ,II)
     &               +XH2S (ICEL)*EHGAZE(IH2S,II)
     &               +XSO2 (ICEL)*EHGAZE(ISO2,II)
C
          EH1(ICEL) = FUEL1(ICEL)*EHGAZE(IFOV,II+1)
     &               +FUEL3(ICEL)*EHGAZE(ICO ,II+1)
     &               +OXYD(ICEL) *EHGAZE(IO2 ,II+1)
     &               +PROD1(ICEL)*EHGAZE(ICO2,II+1)
     &               +PROD2(ICEL)*EHGAZE(IH2O,II+1)
     &               +XINER(ICEL)*EHGAZE(IN2 ,II+1)
     &               +XH2S (ICEL)*EHGAZE(IH2S,II+1)
     &               +XSO2 (ICEL)*EHGAZE(ISO2,II+1)
C
          IF ( EH(ICEL).GE.EH0(ICEL) .AND. EH(ICEL).LE.EH1(ICEL) ) THEN
            TP(ICEL)= TH(II) + (EH(ICEL)-EH0(ICEL))
     &                       *(TH(II+1)-TH(II))/(EH1(ICEL)-EH0(ICEL))
          ENDIF
        ENDDO
      ENDDO
C
C----
C FIN
C----
C
      RETURN
      END
c@z
