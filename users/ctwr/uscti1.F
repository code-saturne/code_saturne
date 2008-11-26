c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
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
                        SUBROUTINE USCTI1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
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
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "ctincl.h"
C
C***********************************************************************
C
      INTEGER          IPHAS
C
C***********************************************************************
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C=======================================================================
C
      IF(1.EQ.1) THEN
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES DONNEES         ',/,
     &'@    *********                                               ',/,
     &'@                      MODULE AEROREFRIGERANTS               ',/,
     &'@                                                            ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usctin DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@     Ce sous-programme utilisateur permet de definir les    ',/,
     &'@       options generales. Il est indispensable.             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C=======================================================================
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 1.  PARAMETRES POUR L'ECART DE TEMPERATURE IMPOSE
C=======================================================================
C
C     ACTIVATION
      IAEERI = 0
C
C     ECART DE REFRIGERATION A IMPOSER
      VAEERI = 13.D0
C
C     FREQUENCE DE MODIFICATION DE LA TEMPERATURE
      IAEERP = 5
C
C     PAS DE TEMPERATURE POUR LE CALCUL DE LA PENTE DE ECARTREF(TEAU)
      PASERI = 0.015D0
C
C     MAXIMUM DE LA TEMPERATURE D'EAU CHAUDE MOYENNE PONDEREE
      AETEMX = 80.D0
C
C     MINIMUM DE LA TEMPERATURE D'EAU REFROIDIE MOYENNE PONDEREE
      AETEMN = 10.D0
C
C=======================================================================
C 2.  POST-PROCESSING DES ZONES D'ECHANGES
C=======================================================================
C
      ICHRZE = 1
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
