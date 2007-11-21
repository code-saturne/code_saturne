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
                        SUBROUTINE FUTHP1
C                       *****************
C     -------------------------------------------------------------
     & ( MODE  , EH     , XESP   , TP  )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC CALCUL DE LA TEMPERATURE DU GAZ
CFONC   EN FONCTION DE L'ENTHALPIE DU GAZ ET DES CONCENTRATIONS
CFONC   SI MODE = 1
CFONC CALCUL DE L'ENTHALPIE DU GAZ
CFONC   EN FONCTION DE LA TEMPERATURE DU GAZ ET DES CONCENTRATIONS
CFONC   SI MODE = -1
CFONC
c@fonce
C                             ARGUMENTS
c@argub
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! EH           ! TR !  -> ! ENTHALPIE DU GAZ                     !
CARGU !              !    !     ! (J/kg de melange gazeux)             !
CARGU ! XESP         ! TR !  -> ! FRACTION MASSIQUE DES ESPECES        !
CARGU ! TP           ! TR ! <-  ! TEMPERATURE DU GAZ (Kelvin)          !
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
      INTEGER          MODE
      DOUBLE PRECISION EH,TP
      DOUBLE PRECISION XESP(NGAZEM)
C
C VARIABLES LOCALES
C
      INTEGER          II
C
      DOUBLE PRECISION EH0 , EH1
C
C***********************************************************************
C=======================================================================
C 1. CALCUL DE LA TEMPERATURE A PARTIR DE l'ENTHALPIE
C=======================================================================
C
      IF ( MODE .EQ. 1 ) THEN
C
        II = NPO
C
C
C --- Clipping eventuel de TP a TH(NPO) si EH > EH1
C
        EH1 = XESP(IFOV )*EHGAZE(IFOV,II)
     &      + XESP(ICO  )*EHGAZE(ICO ,II)
     &      + XESP(IO2  )*EHGAZE(IO2 ,II)
     &      + XESP(ICO2 )*EHGAZE(ICO2,II)
     &      + XESP(IH2O )*EHGAZE(IH2O,II)
     &      + XESP(IN2  )*EHGAZE(IN2 ,II)
     &      + XESP(IH2S )*EHGAZE(IH2S,II)
     &      + XESP(ISO2 )*EHGAZE(ISO2,II)
C
        IF ( EH .GE. EH1 ) THEN
          TP = TH(II)
          GOTO 501
        ENDIF
C
        II = 1
C
C --- Clipping eventuel de TP a TH(1) si EH < EH0
C
        EH0 = XESP(IFOV )*EHGAZE(IFOV,II)
     &      + XESP(ICO  )*EHGAZE(ICO ,II)
     &      + XESP(IO2  )*EHGAZE(IO2 ,II)
     &      + XESP(ICO2 )*EHGAZE(ICO2,II)
     &      + XESP(IH2O )*EHGAZE(IH2O,II)
     &      + XESP(IN2  )*EHGAZE(IN2 ,II)
     &      + XESP(IH2S )*EHGAZE(IH2S,II)
     &      + XESP(ISO2 )*EHGAZE(ISO2,II)
C
        IF ( EH .LE. EH0 ) THEN
          TP= TH(II)
          GOTO 501
        ENDIF
C
C
 500    CONTINUE
        II = II + 1
C
C --- Clipping eventuel de TP a TH(II-1) si EH < EH0
C
        EH0 = XESP(IFOV )*EHGAZE(IFOV,II-1)
     &      + XESP(ICO  )*EHGAZE(ICO ,II-1)
     &      + XESP(IO2  )*EHGAZE(IO2 ,II-1)
     &      + XESP(ICO2 )*EHGAZE(ICO2,II-1)
     &      + XESP(IH2O )*EHGAZE(IH2O,II-1)
     &      + XESP(IN2  )*EHGAZE(IN2 ,II-1)
     &      + XESP(IH2S )*EHGAZE(IH2S,II-1)
     &      + XESP(ISO2 )*EHGAZE(ISO2,II-1)
C
        EH1 = XESP(IFOV )*EHGAZE(IFOV,II)
     &      + XESP(ICO  )*EHGAZE(ICO ,II)
     &      + XESP(IO2  )*EHGAZE(IO2 ,II)
     &      + XESP(ICO2 )*EHGAZE(ICO2,II)
     &      + XESP(IH2O )*EHGAZE(IH2O,II)
     &      + XESP(IN2  )*EHGAZE(IN2 ,II)
     &      + XESP(IH2S )*EHGAZE(IH2S,II)
     &      + XESP(ISO2 )*EHGAZE(ISO2,II)
C
        IF ( EH .GE. EH0  .AND. EH .LE. EH1  ) THEN
          TP = TH(II-1) + (EH-EH0) *
     &                 (TH(II)-TH(II-1))/(EH1-EH0)
          GOTO 501
        ENDIF
        GOTO 500
 501    CONTINUE
C
C=======================================================================
C 1. CALCUL DE L'ENTHALPIE A PARTIR DE LA TEMPERATURE
C=======================================================================
C
      ELSE IF ( MODE .EQ. -1 ) THEN
C
        II = NPO
C
C --- Calcul en Max
C
        EH = XESP(IFOV )*EHGAZE(IFOV,II)
     &      +XESP(ICO  )*EHGAZE(ICO ,II)
     &      +XESP(IO2  )*EHGAZE(IO2 ,II)
     &      +XESP(ICO2 )*EHGAZE(ICO2,II)
     &      +XESP(IH2O )*EHGAZE(IH2O,II)
     &      +XESP(IN2  )*EHGAZE(IN2 ,II)
     &      +XESP(IH2S )*EHGAZE(IH2S,II)
     &      +XESP(ISO2 )*EHGAZE(ISO2,II)
        IF (TP.GT.TH(II)) GOTO 601
C
C Clipping en Min
C
        II = 1
C
C
C --- Clipping eventuel de TP a TH(1) si EH < EH0
C
        EH = XESP(IFOV )*EHGAZE(IFOV,II)
     &      +XESP(ICO  )*EHGAZE(ICO ,II)
     &      +XESP(IO2  )*EHGAZE(IO2 ,II)
     &      +XESP(ICO2 )*EHGAZE(ICO2,II)
     &      +XESP(IH2O )*EHGAZE(IH2O,II)
     &      +XESP(IN2  )*EHGAZE(IN2 ,II)
     &      +XESP(IH2S )*EHGAZE(IH2S,II)
     &      +XESP(ISO2 )*EHGAZE(ISO2,II)
        IF (TP.LT.TH(II)) GOTO 601
C
C Interpolation dans la table
C
        II = 1
 600    CONTINUE
C
        II = II + 1
        IF ( TP .LE. TH(II) ) THEN
C
C --- Calcul de l'enthalpie de l'espece gazeuse TH(II-1)
C
          EH0 = XESP(IFOV )*EHGAZE(IFOV,II-1)
     &         +XESP(ICO  )*EHGAZE(ICO ,II-1)
     &         +XESP(IO2  )*EHGAZE(IO2 ,II-1)
     &         +XESP(ICO2 )*EHGAZE(ICO2,II-1)
     &         +XESP(IH2O )*EHGAZE(IH2O,II-1)
     &         +XESP(IN2  )*EHGAZE(IN2 ,II-1)
     &         +XESP(IH2S )*EHGAZE(IH2S,II-1)
     &         +XESP(ISO2 )*EHGAZE(ISO2,II-1)
C
C --- Calcul de l'enthalpie de l'espece gazeuse TH(I)
C
          EH1 = XESP(IFOV )*EHGAZE(IFOV,II)
     &         +XESP(ICO  )*EHGAZE(ICO ,II)
     &         +XESP(IO2  )*EHGAZE(IO2 ,II)
     &         +XESP(ICO2 )*EHGAZE(ICO2,II)
     &         +XESP(IH2O )*EHGAZE(IH2O,II)
     &         +XESP(IN2  )*EHGAZE(IN2 ,II)
     &         +XESP(IH2S )*EHGAZE(IH2S,II)
     &         +XESP(ISO2 )*EHGAZE(ISO2,II)
C
          EH =  EH0
     &         +(EH1-EH0)*(TP-TH(II-1))/(TH(II)-TH(II-1))
          GOTO 601
        ENDIF
        GOTO 600
C
 601    CONTINUE
C
      ELSE
        WRITE(NFECRA,1000) MODE
        CALL CSEXIT(1)
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ERREUR DANS FUTHP1                          ',/,
     &'@    *********                                               ',/,
     &'@    VALEUR INCORRECTE DE L''ARGUMENT MODE                   ',/,
     &'@    CE DOIT ETRE UN ENTIER EGAL A 1 OU -1                   ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut etre execute.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
      END
c@z

