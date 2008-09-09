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
                          SUBROUTINE LWCURL
C                         *****************
C     -------------------------------------------------------------
     & ( AMPEN1 , VALMOY , VALVAR , VALMIN , VALMAX ,
     &   EXIT01 , EXIT02  , AMPL01 , AMPL02  )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC  CALCUL DES PARAMETRES DE LA PDF
CFONC  PDF LIBBY - WILLIAMS 2 POINTS EN UTILISANT LE
CFONC  MOMENT D'ORDRE 3 DEDUIT DES RECURRENCES SUR
CFONC  LES FONCTIONS BETA
CFONC
CFONC  COMMENTAIRES : HYPOTHESE DE CURL MODIFIEE
CFONC  ------------
CFONC                 a partir de la valeur moyenne d'une variable,
CFONC                 des extremas et de la variance de cette variable
CFONC                 on en deduit 2 etat autour de l'etat moyen
CFONC                 et une amplitude pour chaque etat
CFONC
CFONC  LE RESULTAT EST :
CFONC  ---------------
CFONC
CFONC     CALCUL DES PARAMETRES ASSOCIES AUX FONCTIONS DIRAC
CFONC
CFONC       Les Diracs sont en position [F(.,1),Y(.,1)] et [F(.,2),Y(.,2)]
CFONC       Leurs amplitudes respectives sont D(.,1) et D(.,2)
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! AMPEN1       ! R  !  -> !  AMPLITUDE TOTALE DES PICS           !
CARGU ! VALMOY       ! R  !  -> !  VALEUR MOYENNE DE LA VARIABLE       !
CARGU ! VALVAR       ! R  !  -> !  VARIANCE DE LA VARIABLE             !
CARGU ! VALMIN       ! R  !  -> !  MIN DE LA VARIABLE                  !
CARGU ! VALMAX       ! R  !  -> !  MAX DE LA VARIABLE                  !
CARGU ! EXIT01       ! R  ! <-  !  ETAT 1                              !
CARGU ! EXIT02       ! R  ! <-  !  ETAT 2                              !
CARGU ! AMPL01       ! R  ! <-  !  AMPLITUDE 1                         !
CARGU ! AMPL02       ! R  ! <-  !  AMPLITUDE 2                         !
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
C**********************************************************************
C     DONNEES EN COMMON
C**********************************************************************
C
C
      INCLUDE "paramx.h"
      INCLUDE "pointe.h"
      INCLUDE "entsor.h"
      INCLUDE "cstnum.h"
      INCLUDE "cstphy.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      DOUBLE PRECISION VALMOY, VALVAR
      DOUBLE PRECISION VALMIN, VALMAX
      DOUBLE PRECISION EXIT01, EXIT02
      DOUBLE PRECISION AMPL01, AMPL02
      DOUBLE PRECISION AMPEN1
C
C VARIABLES LOCALES
C
      DOUBLE PRECISION TVV, C, D
      DOUBLE PRECISION MOYADM, VARADM, TVVADM
      DOUBLE PRECISION EPSI
C=======================================================================
C 0.  CALCULS PRELIMINAIRES
C=======================================================================
C
C ---> Test sur l'amplitude totale des pics
C       si elle est trop faible on positionne
C        les deux etat sur l'etat moyen
C
      EPSI = 1.D-6
C
      IF (AMPEN1.GT.EPSI) THEN
C
C ---> Test sur la variance,
C       si elle est trop faible on positionne
C        les deux etat sur l'etat moyen
C
        IF ((VALVAR.GT.EPSI)) THEN
C
C ---> on travaille en variable adimentionnelle pour ce calcul
C
          MOYADM = (VALMOY-VALMIN)/(VALMAX-VALMIN)
          VARADM = VALVAR/((VALMAX-VALMIN)**2)
C
C ---> calcul du moment d'ordre 3 adim
C
          TVVADM = 2.D0*VARADM**2
     &         *((1.D0-2.D0*MOYADM)/((MOYADM*(1-MOYADM))+VARADM))
C
C ---> calcul du moment d'ordre 3 non adim
C
          TVV = ((VALMAX-VALMIN)**3)*TVVADM
C
C ---> calcul du termedu polynome du moment 3
C
          C=(4.D0+TVV**2/VALVAR**3)
C
C ---> determination du signe de la racine
C
          IF ((1.D0-MOYADM).GT.MOYADM) THEN
            D = 0.5D0+SQRT((C-4.D0)/(4.D0*C))
          ELSE
            D = 0.5D0-SQRT((C-4.D0)/(4.D0*C))
          ENDIF
C
C ---> calcul des amplitudes des 2 pics
C
          AMPL01 = AMPEN1 * D
          AMPL02 = AMPEN1 - AMPL01
C
C -----> Calcul des positions des deux pics
C
          EXIT01 = VALMOY - SQRT((1.D0-D)/(D)*VALVAR)
          EXIT02 = VALMOY + SQRT((D)/(1.D0-D)*VALVAR)
C
          EXIT01 = MAX(VALMIN,MIN(EXIT01,VALMAX))
          EXIT02 = MAX(VALMIN,MIN(EXIT02,VALMAX))
C
        ELSE
C variance faible
C
          AMPL01 = AMPEN1/2.D0
          AMPL02 = AMPL01
C
          EXIT01 = VALMOY
          EXIT02 = VALMOY
C
        ENDIF
C
      ELSE
C amplitude totale faible
C
        AMPL01 = AMPEN1/2.D0
        AMPL02 = AMPL01
C
        EXIT01 = VALMOY
        EXIT02 = VALMOY
C
      ENDIF
C
      END
