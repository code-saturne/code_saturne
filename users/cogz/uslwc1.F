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
                        SUBROUTINE USLWC1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC          INIT DES OPTIONS DES VARIABLES POUR
CFONC                   POUR LA COMBUSTION
CFONC            FLAMME DE PREMELANGE : MODELE LWC
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
      INTEGER          IPP, IDIRAC
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
     &'@     MODULE COMBUSTION GAZ MODELE LWC :                     ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR uslwc1 DOIT ETRE COMPLETE',/,
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
C ---- Fraction de melange
      IF ( IPPMOD(ICOLWC).GE.0 ) THEN
        IPP = IPPRTP(ISCA(IFM))
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
C ---- Variance de la fraction de melange
        IPP = IPPRTP(ISCA(IFP2M))
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
C ---- Fraction massique de fuel
        IPP = IPPRTP(ISCA(IYFM))
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
C ---- Variance de la fraction massique de fuel
        IPP = IPPRTP(ISCA(IYFP2M))
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
      IF (IPPMOD(ICOLWC).GE.2) THEN
          IPP = IPPRTP(ISCA(ICOYFP))
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
      ENDIF
C
C ---- Enthalpie
      IF ( IPPMOD(ICOLWC).EQ.1 .OR.
     &     IPPMOD(ICOLWC).EQ.3 .OR.
     &     IPPMOD(ICOLWC).EQ.5    ) THEN
        IPP = IPPRTP(ISCA(IHM))
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
C
C=======================================================================
C 2. VARIABLES ALGEBRIQUES OU D'ETAT
C=======================================================================
C
C --- Terme source
        IPP = IPPPRO(IPPROC(ITSC))
        NOMVAR(IPP)   = 'T.SOURCE'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C --- temperature
        IPP = IPPPRO(IPPROC(ITEMP))
        NOMVAR(IPP)   = 'Temperature'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C --- fraction massique Combustible
        IPP = IPPPRO(IPPROC(IYM(1)))
        NOMVAR(IPP)   = 'YM_Fuel'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C --- fraction massique Oxydant
        IPP = IPPPRO(IPPROC(IYM(2)))
        NOMVAR(IPP)   = 'YM_Oxyd'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C --- fraction massique Produit
        IPP = IPPPRO(IPPROC(IYM(3)))
        NOMVAR(IPP)   = 'YM_Prod'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
        DO IDIRAC = 1, NDIRAC
          IPP = IPPPRO(IPPROC(IRHOL(IDIRAC)))
          WRITE(NOMVAR(IPP),'(A4,I1)') 'RHOL',IDIRAC
          ICHRVR(IPP)   = 1
          ILISVR(IPP)   = 1
          IHISVR(IPP,1) = -1
C
          IPP = IPPPRO(IPPROC(ITEML(IDIRAC)))
          WRITE(NOMVAR(IPP),'(A4,I1)') 'TEML',IDIRAC
          ICHRVR(IPP)   = 1
          ILISVR(IPP)   = 1
          IHISVR(IPP,1) = -1
C
          IPP = IPPPRO(IPPROC(IFMEL(IDIRAC)))
          WRITE(NOMVAR(IPP),'(A4,I1)') 'FMEL',IDIRAC
          ICHRVR(IPP)   = 1
          ILISVR(IPP)   = 1
          IHISVR(IPP,1) = -1
C
          IPP = IPPPRO(IPPROC(IFMAL(IDIRAC)))
          WRITE(NOMVAR(IPP),'(A4,I1)') 'FMAL',IDIRAC
          ICHRVR(IPP)   = 1
          ILISVR(IPP)   = 1
          IHISVR(IPP,1) = -1
C
          IPP = IPPPRO(IPPROC(IAMPL(IDIRAC)))
          WRITE(NOMVAR(IPP),'(A4,I1)') 'AMPL',IDIRAC
          ICHRVR(IPP)   = 1
          ILISVR(IPP)   = 1
          IHISVR(IPP,1) = -1
C
          IPP = IPPPRO(IPPROC(ITSCL(IDIRAC)))
          WRITE(NOMVAR(IPP),'(A4,I1)') 'TSCL',IDIRAC
          ICHRVR(IPP)   = 1
          ILISVR(IPP)   = 1
          IHISVR(IPP,1) = -1
C
          IPP = IPPPRO(IPPROC(IMAML(IDIRAC)))
          WRITE(NOMVAR(IPP),'(A4,I1)') 'MAML',IDIRAC
          ICHRVR(IPP)   = 1
          ILISVR(IPP)   = 1
          IHISVR(IPP,1) = -1
        ENDDO
C
C ---- Modele LWC AVEC RAYONNEMENT
C
      IF ( IRAYPP.GT.0 ) THEN
C
C ---- Coeff d'absorption
        IPP = IPPPRO(IPPROC(ICKABS))
        NOMVAR(IPP)   = 'KABS'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C ---- Terme T^4
        IPP = IPPPRO(IPPROC(IT4M))
        NOMVAR(IPP)   = 'TEMP4'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
C
C ---- Terme T^3
        IPP = IPPPRO(IPPROC(IT3M))
        NOMVAR(IPP)   = 'TEMP3'
        ICHRVR(IPP)   = 0
        ILISVR(IPP)   = 0
        IHISVR(IPP,1) = -1
C
      ENDIF
C
C
C=======================================================================
C 3. OPTIONS DE CALCUL
C=======================================================================
C
C --> Coefficient de relaxation de la masse volumique
C      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
C
      SRROM = 0.95D0
C
C
C=======================================================================
C 4. CONSTANTES PHYSIQUES
C=======================================================================
C
C --> Viscosite laminaire associee au scalaire enthalpie
C       DIFTL0 (diffusivite dynamique en kg/(m s))
      DIFTL0 = 4.25D-5
C
C --> Constante du modele LWC
C
C --- Vitesse de reference
       VREF = 60.D0
C --- Longueur de reference
       LREF = 0.1D0
C --- Temperature d'activation
       TA   = 0.2D5
C --- Temperature de cross-over (combustion du propane)
       TSTAR= 0.12D4
C
C----
C FIN
C----
C
      RETURN
      END
c@z
