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
                        SUBROUTINE USELI1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC          INIT DES OPTIONS DES VARIABLES POUR
CFONC               LE MODULE ELECTRIQUE
CFONC    EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1
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
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
      INTEGER          IPP, IESP , IDIMVE
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
     &'@                      MODULE ELECTRIQUE                     ',/,
     &'@                                                            ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR useli1 DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@     Ce sous-programme utilisateur permet de definir les    ',/,
     &'@       options generales. Il est indispensable.             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C
C=======================================================================
C 1. VARIABLES TRANSPORTEES
C=======================================================================
C
C  Sortie chrono, suivi listing, sortie histo
C     Si l'on n'affecte pas les tableaux suivants,
C     les valeurs par defaut seront utilisees
C
C       ICHRVR( ) = sortie chono (oui 1/non 0)
C       ILISVR( ) = suivi listing (oui 1/non 0)
C       IHISVR( ) = sortie historique (nombre de sondes et numeros)
C       si IHISVR(.,1)  = -1 sortie sur toutes les sondes definies
C                            dans usini1
C
C
C --> Variables communes aux versions electriques
C
C ---- Enthalpie
      IPP = IPPRTP(ISCA(IHM))
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
c ---- Potentiel reel
      IPP = IPPRTP(ISCA(IPOTR))
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C---- Fractions massiques des constituants
      IF ( NGAZG .GT. 1 ) THEN
        DO IESP = 1, NGAZG-1
          IPP = IPPRTP(ISCA(IYCOEL(IESP)))
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDDO
      ENDIF
C
C --> Version effet Joule
C
      IF ( IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4) THEN
        IPP = IPPRTP(ISCA(IPOTI))
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
C --> Version arc electrique 3D
C
      IF ( IPPMOD(IELARC).GE.2 ) THEN
        DO IDIMVE = 1, NDIMVE
          IPP = IPPRTP(ISCA(IPOTVA(IDIMVE)))
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDDO
      ENDIF
C
C --> Version conduction ionique
C     indisponible dans la version presente
C
C=======================================================================
C 2. VARIABLES ALGEBRIQUES OU D'ETAT
C=======================================================================
C
C ---- Temperature
      IPP = IPPPRO(IPPROC(ITEMP) )
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C ---- Conductivite electrique
      IPP = IPPPRO(IPPROC(IVISLS(IPOTR)))
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C ---- Puissance volumique dissipee par effet Joule
      IPP = IPPPRO(IPPROC(IEFJOU) )
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C ---- densite de courant reelle
      DO IDIMVE = 1, NDIMVE
        IPP = IPPPRO(IPPROC(IDJR(IDIMVE)) )
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDDO
C
C ---- densite de courant imaginaire
      IF ( IPPMOD(IELJOU).EQ.4 ) THEN
        DO IDIMVE = 1, NDIMVE
          IPP = IPPPRO(IPPROC(IDJI(IDIMVE)) )
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDDO
      ENDIF
C
      IF ( IPPMOD(IELARC).GE.1 ) THEN
C
C ---- Forces electromagnetiques de Laplace
        DO IDIMVE = 1, NDIMVE
          IPP = IPPPRO(IPPROC(ILAPLA(IDIMVE)) )
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDDO
C
C ---- Coefficient d'absorption ou TS radiatif
        IF ( IXKABE.GT.0 ) THEN
          IPP = IPPPRO(IPPROC(IDRAD) )
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDIF
      ENDIF
C
C ---- Charge electrique volumique
      IF ( IPPMOD(IELION).GE.1 ) THEN
        IPP = IPPPRO(IPPROC(IQELEC) )
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
C
C=======================================================================
C 3. OPTIONS DE CALCUL
C=======================================================================
C
C --> Coefficient de relaxation de la masse volumique
C      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
      SRROM = 0.D0
C
C --> Recalage des variables electriques (joule ou arc electrique)
C      IELCOR = 0 : pas de correction
C      IELCOR = 1 : correction
      IELCOR = 0
C
C     Intensite de courant imposee (arc electrique) en Ampere
C             et Puissance imposee (effet Joule/verre) en Watt
C       ces valeurs doivent etre positives
C       (et en general elles sont strictement positives)
      COUIMP = 0.D0
      PUISIM = 0.D0
C
C     Differentiel de potentiel initiale
C       la valeur doit etre strictement positive
      DPOT = 0.D0
C
C
C----
C FIN
C----
C
      RETURN
      END
c@z
