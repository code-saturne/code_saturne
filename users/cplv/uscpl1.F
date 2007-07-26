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
                        SUBROUTINE USCPL1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
CFONC    --------------------------------------------------------------
CFONC
CFONC     ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE
CFONC
CFONC       COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
CFONC       TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON
CFONC
CFONC     ROUTINE UTILISATEUR POUR ENTREE DES PARAMETRES DE CALCUL
CFONC       (COMMONS)
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
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
      INTEGER          IPP , ICHA
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
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@     COMBUSTION CHARBON PULVERISE COUPLE AU                 ',/,
     &'@     TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON :       ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR uscpl1 DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C=======================================================================
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
C --> Variables propres a la phase gaz continue
C
C      - Enthalpie de la phase gaz continue
      IPP = IPPRTP(ISCA(IHM))
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C --> Variables propres a la phase continue
C
      DO ICHA = 1, NCHARB
C
C       - Moyenne du traceur 1
C         (representatif des MV legeres du charbon ICHA)
        IPP = IPPRTP(ISCA(IF1M(ICHA)))
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
C       - Moyenne du traceur 2
C         (representatif des MV lourdes du charbon ICHA)
        IPP = IPPRTP(ISCA(IF2M(ICHA)))
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
      ENDDO
C
C     - Moyenne du traceur 3 (representatif du C libere sous forme de CO
C       lors de la combustion heterogene)
      IPP = IPPRTP(ISCA(IF3M))
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C     - Variance associe au traceur 4 (representatif de l'air)
      IPP = IPPRTP(ISCA(IF4P2M))
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C=======================================================================
C 2. VARIABLES ALGEBRIQUES OU D'ETAT
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
C --> Variables algebriques propres a la suspension gaz - particules
C
C     - Masse molaire du melange gazeux
      IPP = IPPPRO(IPPROC(IMMEL))
      ICHRVR(IPP)   = 0
      ILISVR(IPP)   = 0
      IHISVR(IPP,1) = -1
C
C     - Temperature du melange gazeux
      IPP = IPPPRO(IPPROC(ITEMP1))
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     - Fraction massique (dans le melange gazeux) du CHx1m
      IPP = IPPPRO(IPPROC(IYM1(1)))
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     - Fraction massique (dans le melange gazeux) du CHx2m
      IPP = IPPPRO(IPPROC(IYM1(2)))
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     - Fraction massique (dans le melange gazeux) du CO
      IPP = IPPPRO(IPPROC(IYM1(3)))
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     - Fraction massique (dans le melange gazeux) du O2
      IPP = IPPPRO(IPPROC(IYM1(4)))
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     - Fraction massique (dans le melange gazeux) du CO2
      IPP = IPPPRO(IPPROC(IYM1(5)))
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     - Fraction massique (dans le melange gazeux) du H2O
      IPP = IPPPRO(IPPROC(IYM1(6)))
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C     - Fraction massique (dans le melange gazeux) du N2
      IPP = IPPPRO(IPPROC(IYM1(7)))
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C ---> Rayonnement
C
      IPP = IPPPRO(IPPROC(ILUMI))
      NOMVAR(IPP)   = 'IntLuminance_4PI'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C
C=======================================================================
C 3. OPTIONS DE CALCUL
C=======================================================================
C
C --- Coefficient de relaxation de la masse volumique
C      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
C
      SRROM = 0.8D0
C
C
C=======================================================================
C 4. CONSTANTES PHYSIQUES
C=======================================================================
C
C ---> Viscosite laminaire associee au scalaire enthalpie
C       DIFTL0 (diffusivite dynamique en kg/(m s))
      DIFTL0 = 4.25D-5
C
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
