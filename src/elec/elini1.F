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
                        SUBROUTINE ELINI1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC    INIT DES OPTIONS DES VARIABLES POUR LE MODULE ELECTRIQUE
CFONC       EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1
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
C VARIABLES LOCALES
C
      INTEGER          IDIMVE , IESP
      INTEGER          IPP , II , IOK
      INTEGER          ISC , IPHAS , IVAR
C
C***********************************************************************
C=======================================================================
C 0. VERIFICATION ISCALT, ISCSTH
C=======================================================================
C     L'utilisateur ne doit pas y avoir touche.
C
      IOK = 0
C
      DO IPHAS = 1, NPHAS
        IF(ISCALT(IPHAS).NE.-1) THEN
          WRITE(NFECRA,1000)IPHAS,ISCALT(IPHAS)
          IOK = IOK + 1
        ENDIF
      ENDDO
      DO II = 1, NSCAPP
        IF(ISCSTH(ISCAPP(II)).NE.-10) THEN
          WRITE(NFECRA,1001)II,ISCAPP(II),ISCAPP(II),ISCSTH(ISCAPP(II))
          IOK = IOK + 1
        ENDIF
      ENDDO
C
      IF(IOK.NE.0) THEN
        CALL CSEXIT (1)
C       ===========
      ENDIF
C
C=======================================================================
C 1. VARIABLES TRANSPORTEES
C=======================================================================
C
C 1.1 Definition des scamin et des scamax des variables transportees
C ==================================================================
C
C --> Dans toutes les versions electriques
C     Enthalpie
      SCAMIN(IHM)   = -GRAND
      SCAMAX(IHM)   = +GRAND
C     Potentiel reel
      SCAMIN(IPOTR) = -GRAND
      SCAMAX(IPOTR) = +GRAND
C     Fractions massiques des constituants
      IF ( NGAZG .GT. 1 ) THEN
        DO IESP = 1, NGAZG-1
          SCAMIN(IYCOEL(IESP)) = 0.D0
          SCAMAX(IYCOEL(IESP)) = 1.D0
        ENDDO
      ENDIF
C
C --> Effet Joule (cas potentiel imaginaire)
C     Potentiel imaginaire
C
      IF ( IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4 ) THEN
        SCAMIN(IPOTI) = -GRAND
        SCAMAX(IPOTI) = +GRAND
      ENDIF
C
C --> Arc electrique
C     Potentiel vecteur
      IF ( IPPMOD(IELARC).GE.2 ) THEN
        DO IDIMVE = 1, NDIMVE
          SCAMIN(IPOTVA(IDIMVE)) = -GRAND
          SCAMAX(IPOTVA(IDIMVE)) = +GRAND
        ENDDO
      ENDIF
C
C --> Conduction ionique (a developper)
C
C 1.2 Nature des scalaires transportes
C ====================================
C
C ---- Type de scalaire (0 passif, 1 temperature en K
C                                 -1 temperature en C
C                                  2 enthalpie)
C      La distinction -1/1 sert pour le rayonnement
C
C     Par defaut, scalaire "passif"
      DO ISC = 1, NSCAPP
        ISCSTH(ISCAPP(ISC)) = 0
      ENDDO
C
C     Pour l'enthalpie
      IPHAS = IPHSCA(IHM)
      ISCALT(IPHAS) = IHM
      ISCSTH(IHM)   = 2
C
C 1.4 Donnees physiques ou numeriques propres aux scalaires ELECTRIQUES
C =====================================================================
C
C
C --> Conditions associees aux potentiels
C     (les autres variables ont des comportements par defaut)
      IVAR = ISCA(IPOTR)
      ICONV (IVAR) = 0
      ISTAT (IVAR) = 0
      IDIFF (IVAR) = 1
      IDIFFT(IVAR) = 0
      IDIRCL(IVAR) = 1
C
      IF(IPPMOD(IELJOU).EQ.2 .OR. IPPMOD(IELJOU).EQ.4) THEN
        IVAR = ISCA(IPOTI)
        ICONV (IVAR) = 0
        ISTAT (IVAR) = 0
        IDIFF (IVAR) = 1
        IDIFFT(IVAR) = 0
        IDIRCL(IVAR) = 1
      ENDIF
C
      IF(IPPMOD(IELARC).GE.2) THEN
        DO IDIMVE = 1, NDIMVE
          IVAR = ISCA(IPOTVA(IDIMVE))
          ICONV (IVAR) = 0
          ISTAT (IVAR) = 0
          IDIFF (IVAR) = 1
          IDIFFT(IVAR) = 0
          IDIRCL(IVAR) = 1
        ENDDO
      ENDIF
C
C --> "Viscosite" associee au potentiel vecteur
C     (c'est la seule qui est constante)
      IF ( IPPMOD(IELARC).GE.2 ) THEN
        VISLS0(IPOTVA(1)) = 1.D0
        VISLS0(IPOTVA(2)) = 1.D0
        VISLS0(IPOTVA(3)) = 1.D0
      ENDIF
C
C --> Schmidt ou Prandtl turbulent
C     (pour les potentiels, c'est inutile puisque IDIFFT=0)
C
      DO ISC = 1, NSCAPP
        SIGMAS(ISCAPP(ISC)) = 0.7D0
      ENDDO
C
C ---> Pour tous les scalaires
C
      DO ISC = 1, NSCAPP
C
C ----- Niveau de detail des impressions pour les variables et
C          donc les scalaires (valeurs 0 ou 1)
C          Si = -10000 non modifie par l'utilisateur -> niveau 1

        IVAR = ISCA(ISCAPP(ISC))
        IF(IWARNI(IVAR).EQ.-10000) THEN
          IWARNI(IVAR) = 1
        ENDIF
C
C ----- Informations relatives a la resolution des scalaires
C
C       - Facteur multiplicatif du pas de temps
C
        CDTVAR(IVAR) = 1.D0
C
C         - Schema convectif % schema 2ieme ordre
C           = 0 : upwind
C           = 1 : second ordre
        BLENCV(IVAR) = 1.D0
C
C         - Type de schema convectif second ordre (utile si BLENCV > 0)
C           = 0 : Second Order Linear Upwind
C           = 1 : Centre
        ISCHCV(IVAR) = 1
C
C         - Test de pente pour basculer d'un schema centre vers l'upwind
C           = 0 : utlisation automatique du test de pente
C           = 1 : calcul sans test de pente
        ISSTPC(IVAR) = 0
C
C         - Reconstruction des flux de convection et de diffusion aux faces
C           = 0 : pas de reconstruction
        IRCFLU(IVAR) = 1
C
      ENDDO
C
C
C 1.5 Variable courante : nom, sortie chrono, suivi listing, sortie histo
C =======================================================================
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
C     NB : Les 8 premiers caracteres du nom seront repris dans le
C          listing 'developpeur'
C
C =======================================================================
C
C --> Variables communes aux versions electriques
C
      IPP = IPPRTP(ISCA(IHM))
      NOMVAR(IPP)  = 'Enthalpy'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
      IPP = IPPRTP(ISCA(IPOTR))
      NOMVAR(IPP)  = 'POT_EL_R'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
      IF ( NGAZG .GT. 1 ) THEN
        DO IESP = 1, NGAZG-1
          IPP = IPPRTP(ISCA(IYCOEL(IESP)))
          WRITE(NOMVAR(IPP),'(A6,I2.2)')'YM_ESL',IESP
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
        NOMVAR(IPP)  = 'POT_EL_I'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
C --> Version arc electrique
C
      IF ( IPPMOD(IELARC).GE.2 ) THEN
        DO IDIMVE = 1, NDIMVE
          IPP = IPPRTP(ISCA(IPOTVA(IDIMVE)))
          WRITE(NOMVAR(IPP),'(A7,I1.1)')'POT_VEC',IDIMVE
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDDO
      ENDIF
C
C --> Version conduction ionique
C
C=======================================================================
C 2. VARIABLES ALGEBRIQUES OU D'ETAT
C=======================================================================
C
      IPP = IPPPRO(IPPROC(ITEMP) )
      NOMVAR(IPP)  = 'Temper'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
      IPP = IPPPRO(IPPROC(IEFJOU) )
      NOMVAR(IPP)  = 'PuisJoul'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
      DO IDIMVE = 1, NDIMVE
        IPP = IPPPRO(IPPROC(IDJR(IDIMVE)) )
        WRITE(NOMVAR(IPP),'(A7,I1.1)')'Cour_re',IDIMVE
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDDO
C
      IF ( IPPMOD(IELJOU).EQ.4 ) THEN
        DO IDIMVE = 1, NDIMVE
          IPP = IPPPRO(IPPROC(IDJI(IDIMVE)) )
          WRITE(NOMVAR(IPP),'(A7,I1.1)')'CouImag',IDIMVE
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDDO
      ENDIF
C
      IF ( IPPMOD(IELARC).GE.1 ) THEN
        DO IDIMVE = 1, NDIMVE
          IPP = IPPPRO(IPPROC(ILAPLA(IDIMVE)) )
          WRITE(NOMVAR(IPP),'(A7,I1.1)')'For_Lap',IDIMVE
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDDO
C
        IF ( IXKABE .EQ.1 ) THEN
          IPP = IPPPRO(IPPROC(IDRAD) )
          NOMVAR(IPP)  = 'Coef_Abso'
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDIF

        IF ( IXKABE .EQ.2 ) THEN
          IPP = IPPPRO(IPPROC(IDRAD) )
          NOMVAR(IPP)  = 'TS_radia'
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDIF
C
      ENDIF
C
      IF ( IPPMOD(IELION).GE.1 ) THEN
        IPP = IPPPRO(IPPROC(IQELEC) )
        NOMVAR(IPP)  = 'Charge'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
C Conductivite Electrique
C
      IPP = IPPPRO(IPPROC(IVISLS(IPOTR)) )
      NOMVAR(IPP)  = 'Sigma'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C     Conductivite electrique imaginaire :
C     La conductivite reelle et imaginaire sont dans le meme tableau.
C     Il convient donc de ne pas renseigner NOMVAR ICHRVR ILISVR IHISVR
C     pour IPP = IPPPRO(IPPROC(IVISLS(IPOTI)) )
C     puisque IPPROC(IVISLS(IPOTI)) = IPPROC(IVISLS(IPOTR))
C
C=======================================================================
C 3. INFORMATIONS COMPLEMENTAIRES
C=======================================================================
C
C --> Coefficient de relaxation de la masse volumique
C      a partir du 2ieme pas de temps, on prend :
C      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
      SRROM = 0.D0
C
C --> Recalage des variables electriques
C      IELCOR = 0 : pas de correction
C      IELCOR = 1 : correction
      IELCOR = 0
C
C     Intensite de courant imposee (arc electrique) ou
C                Puissance imposee (Joule)
      COUIMP = 0.D0
      PUISIM = 0.D0
C
C     Differentiel de potentiel Initial en arc (et Joule)
      DPOT = 0.D0
C
C     Coefficient pour la correction en Joule
      COEJOU = 1.D0
C
C ---> Masse volumique variable et viscosite variable (pour les suites)
      IROVAR(IPHAS) = 1
      IVIVAR(IPHAS) = 1
C
C=======================================================================
C 4. ON REDONNE LA MAIN A L'UTLISATEUR
C=======================================================================
C
      CALL USELI1
C     ===========
C
C=======================================================================
C 5. VERIFICATION DES DONNEES ELECTRIQUES
C=======================================================================
C
      IOK = 0
C
      CALL ELVERI (IOK)
C     ===========
C
      IF(IOK.GT.0) THEN
        WRITE(NFECRA,9999)IOK
        CALL CSEXIT (1)
C       ===========
      ELSE
        WRITE(NFECRA,9998)
      ENDIF
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (JOULE) DEMANDEE                  ',/,
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
     &'@    PHYSIQUE PARTICULIERE (JOULE) DEMANDEE                  ',/,
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
 9998 FORMAT(
     &'                                                             ',/,
     &' Pas d erreur detectee lors de la verification des donnees   ',/,
     &'                                                    (useli1).',/)
 9999 FORMAT(
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,
     &'@                                                            ',/,
     &'@  Se reporter aux impressions precedentes pour plus de      ',/,
     &'@    renseignements.                                         ',/,
     &'@  Verifier useli1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
      END
c@z
