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
                        SUBROUTINE CFINI1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC          INIT DES OPTIONS DES VARIABLES POUR
CFONC               LE COMPRESSIBLE SANS CHOC
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
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER          IPP , II
      INTEGER          IPHAS, IOK
C
C***********************************************************************
C=======================================================================
C 0. VERIFICATION ISCALT, ISCSTH
C=======================================================================
C     L'utilisateur ne doit pas y avoir touche.
C
      DO IPHAS = 1, NPHAS
        IF(ISCALT(IPHAS).NE.-1) THEN
          WRITE(NFECRA,1000)IPHAS,ISCALT(IPHAS)
          CALL CSEXIT (1)
C         ===========
        ENDIF
      ENDDO
      DO II = 1, NSCAPP
        IF(ISCSTH(ISCAPP(II)).NE.-10) THEN
          WRITE(NFECRA,1001)II,ISCAPP(II),ISCAPP(II),ISCSTH(ISCAPP(II))
          CALL CSEXIT (1)
C         ===========
        ENDIF
      ENDDO
C=======================================================================
C 1. VARIABLES TRANSPORTEES
C=======================================================================
C
C 1.1 Definition des scamin et des scamax des variables transportees
C ==================================================================
C
C
      DO IPHAS = 1, NPHAS
C
        IF(  (ABS(SCAMIN(IRHO  (IPHAS))+GRAND).GT.EPZERO).OR.
     &       (ABS(SCAMIN(IENERG(IPHAS))+GRAND).GT.EPZERO).OR.
     &       (ABS(SCAMIN(ITEMPK(IPHAS))+GRAND).GT.EPZERO).OR.
     &       (ABS(SCAMAX(IRHO  (IPHAS))-GRAND).GT.EPZERO).OR.
     &       (ABS(SCAMAX(IENERG(IPHAS))-GRAND).GT.EPZERO).OR.
     &       (ABS(SCAMAX(ITEMPK(IPHAS))-GRAND).GT.EPZERO) ) THEN
          WRITE(NFECRA,2000)
     &         SCAMIN(IRHO  (IPHAS)),SCAMAX(IRHO  (IPHAS)),
     &         SCAMIN(IENERG(IPHAS)),SCAMAX(IENERG(IPHAS)),
     &         SCAMIN(ITEMPK(IPHAS)),SCAMAX(ITEMPK(IPHAS))
          CALL CSEXIT (1)
        ENDIF
c        SCAMIN(IRHO  (IPHAS))   = -GRAND
c        SCAMAX(IRHO  (IPHAS))   =  GRAND
c        SCAMIN(IENERG(IPHAS))   = -GRAND
c        SCAMAX(IENERG(IPHAS))   =  GRAND
c        SCAMIN(ITEMPK(IPHAS))   = -GRAND
c        SCAMAX(ITEMPK(IPHAS))   =  GRAND
C
      ENDDO
C
C 1.2 Nature des scalaires transportes
C ====================================
C
C ---- Type de scalaire (0 passif, 1 temperature en K
C                                 -1 temperature en C
C                                  2 enthalpie en J
C                                  3 energie totale en J)
C      La distinction -1/1 sert pour le rayonnement
C
      DO IPHAS = 1, NPHAS
C
        ISCSTH(IRHO  (IPHAS)) = 0
        ISCSTH(IENERG(IPHAS)) = 3
        ISCSTH(ITEMPK(IPHAS)) = 0
C
        ISCALT(IPHAS) = IENERG(IPHAS)
C
      ENDDO
C
C
C         - Schema convectif % schema 2ieme ordre
C           = 0 : upwind
C           = 1 : second ordre
      DO II = 1, NVARMX
        BLENCV(II) = 0.D0
      ENDDO
C
C         Upwind necessaire pour le schema utilise
C
C
C 1.3 Variable courante : nom, sortie chrono, suivi listing, sortie hist
C ======================================================================
C
C     Comme pour les autres variables,
C       si l'on n'affecte pas les tableaux suivants,
C       les valeurs par defaut seront utilisees
C
C     NOMVAR( ) = nom de la variable
C     ICHRVR( ) = sortie chono (oui 1/non 0)
C     ILISVR( ) = suivi listing (oui 1/non 0)
C     IHISVR( ) = sortie historique (nombre de sondes et numeros)
C     si IHISVR(.,1)  = -1 sortie sur toutes les sondes
C
C     NB : Les 8 premiers caracteres du noms seront repris dans le
C          listing 'developpeur'
C
C ======================================================================
C
      DO IPHAS = 1, NPHAS
C
        IPP = IPPRTP(ISCA(IRHO  (IPHAS)))
        NOMVAR(IPP)  = 'Rho'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
        IPP = IPPRTP(ISCA(IENERG(IPHAS)))
        NOMVAR(IPP)  = 'EnergieT'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
        IPP = IPPRTP(ISCA(ITEMPK(IPHAS)))
        NOMVAR(IPP)  = 'Temp K'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
      ENDDO
C
C
C=======================================================================
C 2. PARAMETRES GLOBAUX
C=======================================================================
C
C --- Couplage vitesse/pression (0 : algorithme classique,
C                                1 : couplage instationnaire)
C     Uniquement en monophasique et en incompressible
C
      IF( IPUCOU.NE.0 ) THEN
        WRITE(NFECRA,3000) IPUCOU
        CALL CSEXIT (1)
      ENDIF
C
C
C --- Estimateurs pour Navier-Stokes
C
C     Interdits en compressible
C
      DO IPHAS = 1, NPHAS
        IF( (IESCAL(IESPRE,IPHAS).NE.0) .OR.
     &      (IESCAL(IESDER,IPHAS).NE.0) .OR.
     &      (IESCAL(IESCOR,IPHAS).NE.0) .OR.
     &      (IESCAL(IESTOT,IPHAS).NE.0) ) THEN
          WRITE(NFECRA,4000)
          CALL CSEXIT (1)
        ENDIF
c       IESCAL(IESPRE,IPHAS) = 0
c       IESCAL(IESDER,IPHAS) = 0
c       IESCAL(IESCOR,IPHAS) = 0
c       IESCAL(IESTOT,IPHAS) = 0
      ENDDO
C
C
C=======================================================================
C 3. OPTIONS DE CALCUL PAR DEFAUT
C=======================================================================
C
C     Pour chaque phase
C
      DO IPHAS = 1, NPHAS
C
C --> Conditions aux limites prenant en compte l'equilibre hydrostatique
C     (oui = 1 , non = 0)
C
        ICFGRP(IPHAS) = 1
C
C
C ---> Masse volumique variable et viscosite constante (pour les suites)
        IROVAR(IPHAS) = 1
        IVIVAR(IPHAS) = 0
C
      ENDDO
C
C
C=======================================================================
C 4. ON REDONNE LA MAIN A L'UTLISATEUR
C=======================================================================
C
      CALL USCFX1
C     ===========
C
C=======================================================================
C 5. OPTIONS DE CALCUL OBLIGATOIRES
C     qui pourront etre remontees au dessus de uscfx1
C     selon les developpements
C=======================================================================
C
C     Pour chaque phase
C
      DO IPHAS = 1, NPHAS
C
        IDIFF(ISCA(IRHO(IPHAS))) = 1
C
C --> Implicitation du terme de convection de l'equation de masse
C     (oui = 1 , non = 0)
C     On choisit 0 ; c'est la seule option qui a ete testee. Elle
C       facilite le codage pour le respect du flux de masse au bord.
C
        ICONV(ISCA(IRHO(IPHAS))) = 0
C
C --> Prise en compte de la pression predite pour resoudre Navier-Stokes
C     (oui = 1 , non = 0)
C
        IGRDPP(IPHAS) = 0
C
C --> Prediction de pression par une equation d'evolution
C
C     ATTENTION   PAS ENCORE IMPLEMENTE
C     *********   LAISSER IPPRED(IPHAS) = 0
C
        IPPRED(IPHAS) = 0
C
C
      ENDDO
C
C=======================================================================
C 6. VERIFICATIONS
C=======================================================================
C
      IOK = 0
      DO IPHAS = 1, NPHAS
        IF(ICFGRP(IPHAS).NE.0.AND.ICFGRP(IPHAS).NE.1) THEN
          WRITE(NFECRA,5000)IPHAS,'ICFGRP',ICFGRP(IPHAS)
          IOK = 1
        ENDIF
      ENDDO
C
      IF (IOK.NE.0) THEN
        CALL CSEXIT (1)
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
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,
     &'@                                                            ',/,
     &'@  La valeur de ISCALT est renseignee automatiquement.       ',/,
     &'@                                                            ',/,
     &'@  L''utilisateur ne doit pas la renseigner dans usini1, or  ',/,
     &'@    elle a ete affectee comme suit :                        ',/,
     &'@    ISCALT(',I10   ,') = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,
     &'@                                                            ',/,
     &'@  Les valeurs de ISCSTH sont renseignees automatiquement.   ',/,
     &'@                                                            ',/,
     &'@  L''utilisateur ne doit pas les renseigner dans usini1, or ',/,
     &'@    pour le scalaire ',I10   ,' correspondant au scalaire   ',/,
     &'@    physique particuliere ',I10   ,' on a                   ',/,
     &'@    ISCSTH(',I10   ,') = ',I10                               ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,
     &'@                                                            ',/,
     &'@  Les bornes des variables rho, energie ou temperature      ',/,
     &'@    ont ete modifiees dans usini1 :                         ',/,
     &'@                                                            ',/,
     &'@                      SCAMIN        SCAMAX                  ',/,
     &'@  rho         ',2E14.5                                       ,/,
     &'@  energie     ',2E14.5                                       ,/,
     &'@  temperature ',2E14.5                                       ,/,
     &'@                                                            ',/,
     &'@  Les bornes de ces variables ne doivent pas etre modifiees ',/,
     &'@  dans usini1. On peut modifier les bornes des variables    ',/,
     &'@  rho et energie dans uscfx1, mais ce n''est pas conseille. ',/,
     &'@  Il est preferable de gerer les depassements éventuels     ',/,
     &'@  au moyen du sous programme uscfth (arret du calcul en fin ',/,
     &'@  de pas de temps en cas de depassement).                   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier usini1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,
     &'@                                                            ',/,
     &'@  L''option IPUCOU = ',I10                                   ,/,
     &'@    n''est pas compatible avec le module compressible       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Imposer IPUCOU = 0 dans usini1.                           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 4000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,
     &'@                                                            ',/,
     &'@  Les estimateurs ne sont pas compatibles avec le module    ',/,
     &'@    compressible                                            ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Imposer IESCAL(.,.) = 0 dans usini1.                      ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 5000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,
     &'@                                                            ',/,
     &'@    PHASE ',I10                                              ,/,
     &'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,
     &'@    IL VAUT ICI ',I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier uscfx1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
      END
c@z
