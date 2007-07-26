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
                        SUBROUTINE FUINI1
C                       *****************
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC        INIT DES OPTIONS DES VARIABLES TRANSPORTEES
CFONC                      ET DES VARIABLES ALGEBRIQUES
CFONC                COMBUSTION FUEL
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
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "fuincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
      INTEGER          IPP ,  II , JJ , IOK
      INTEGER          ISC , IS, IPHAS
      DOUBLE PRECISION WMOLME
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
C
C=======================================================================
C 1. VARIABLES TRANSPORTEES
C=======================================================================
C
C --> Definition des scamin et des scamax des variables transportees
C
C ---- Variables propres a la phase dispersee
C
      SCAMIN(IHM) = H02FOL
      SCAMAX(IHM) = +GRAND
      SCAMIN(ING)   = 0.D0
      SCAMAX(ING)   = +RINFIN
      SCAMIN(IYFOL) = 0.D0
      SCAMAX(IYFOL) = 4.D-1
      SCAMIN(IHLF)  = -GRAND
      SCAMAX(IHLF)  = +GRAND
C
C ---- Variables propres a la phase continue
C
      SCAMIN(IFVAP)   = 0.D0
      SCAMAX(IFVAP)   = 1.D0
      SCAMIN(IFHTF)   = 0.D0
      SCAMAX(IFHTF)   = 1.D0
      SCAMIN(IF4P2M)  = 0.D0
      SCAMAX(IF4P2M)  = 0.25D0
C
      IF ( IEQCO2 .EQ. 1 ) THEN
        SCAMIN(IYCO2) = 0.D0
        SCAMAX(IYCO2) = 1.D0
      ENDIF
C
C --> Nature des scalaires transportes
C
      DO ISC = 1, NSCAPP
C
C ---- Type de scalaire (0 passif, 1 temperature en K
C                                 -1 temperature en C
C                                  2 enthalpie)
C      La distinction -1/1 sert pour le rayonnement
        ISCSTH(ISCAPP(ISC)) = 0
C
      ENDDO
C
      IPHAS = 1
C
C ---- On resout en enthalpie avec un CP constant (Cf. cpvarp)
C
      IPHAS = IPHSCA(IHM)
      ISCALT(IPHAS) = IHM
      ISCSTH(IHM)   = 2
C
C --> Donnees physiques ou numeriques propres aux scalaires CP
C
      DO ISC = 1, NSCAPP
C
        JJ = ISCAPP(ISC)
C
        IF ( ISCAVR(JJ) .LE. 0 ) THEN
C
C        En combustion on considere que la viscosite turbulente domine
C        ON S'INTERDIT DONC LE CALCUL DES FLAMMES LAMINAIRES AVEC Le =/= 1
C
          VISLS0(JJ) = VISCL0(IPHSCA(JJ))
C
        ENDIF
C
C ------ Schmidt ou Prandtl turbulent
C
        SIGMAS(JJ) = 0.7D0
C
C ------ Coeff dissipation des fluctuations
C
        RVARFL(JJ) = 0.8D0
C
        II = ISCA(ISCAPP(ISC))
C
C ------ Niveau de detail des impressions pour les variables et
C          donc les scalaires (valeurs 0 ou 1)
C          Si = -10000 non modifie par l'utilisateur -> niveau 1
        IF(IWARNI(II).EQ.-10000) THEN
          IWARNI(II) = 1
        ENDIF
C
C ------ Informations relatives a la resolution des scalaires
C
C         - Facteur multiplicatif du pas de temps
        CDTVAR(II) = 1.D0
C
C         - Schema convectif % schema 2ieme ordre
C           = 0 : upwind
C           = 1 : second ordre
        BLENCV(II) = 0.D0
C
C         - Type de schema convetif second ordre (utile si BLENCV > 0)
C           = 0 : Second Order Linear Upwind
C           = 1 : Centre
        ISCHCV(II) = 1
C
C         - Test de pente pour basculer d'un schema centre vers l'upwind
C           = 0 : utlisation automatique du test de pente
C           = 1 : calcul sans test de pente
        ISSTPC(II) = 0
C
C         - Reconstruction des flux de convetion et de diffusion aux faces
C           = 0 : pas de reconstruction
        IRCFLU(II) = 0
C
      ENDDO
C
C
C---> Variable courante : nom, sortie chrono, suivi listing, sortie histo
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
C
C ---- Variables propres a la suspension gaz - particules
C
      IPP = IPPRTP(ISCA(IHM))
      NOMVAR(IPP)  = 'Enthalpy'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C ---- Variables propres a la phase dispersee
      IPP = IPPRTP(ISCA(IYFOL))
      NOMVAR(IPP) = 'YFOL_FOL'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
      IPP = IPPRTP(ISCA(IHLF))
      NOMVAR(IPP) = 'HLF_FOL'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
      IPP = IPPRTP(ISCA(ING))
      NOMVAR(IPP) = 'NG_FOL'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
C ---- Variables propres a la phase gaz
C
      IPP = IPPRTP(ISCA(IFVAP))
      NOMVAR(IPP)  = 'Fr_VAP'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
      IPP = IPPRTP(ISCA(IFHTF))
      NOMVAR(IPP)  = 'Fr_HET'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
      IPP = IPPRTP(ISCA(IF4P2M))
      NOMVAR(IPP)  = 'Var_AIR'
      ICHRVR(IPP)  = 1
      ILISVR(IPP)  = 1
      IHISVR(IPP,1)= -1
C
      IF ( IEQCO2 .EQ. 1 ) THEN
        IPP = IPPRTP(ISCA(IYCO2))
        NOMVAR(IPP)  = 'FR_CO2'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
C=======================================================================
C 2. VARIABLES ALGEBRIQUES OU D'ETAT
C=======================================================================
C
C
C ---> Variables algebriques propres a la suspension gaz - particules
C
      IPP = IPPPRO(IPPROC(IMMEL))
      NOMVAR(IPP)   = 'XM'
      ICHRVR(IPP)   = 0
      ILISVR(IPP)   = 0
      IHISVR(IPP,1) = -1
C
C ---> Variables algebriques propres a la phase dispersee
C
      IPP = IPPPRO(IPPROC(ITEMP3))
      WRITE(NOMVAR(IPP),'(A6,I2.2)')'Tem_FOL'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IROM3))
      WRITE(NOMVAR(IPP),'(A6,I2.2)')'Rho_FOL'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IDIAM3))
      WRITE(NOMVAR(IPP),'(A6,I2.2)')'Dia_gt'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IH1HLF))
      WRITE(NOMVAR(IPP),'(A6,I2.2)')'H1-Hlf'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IGMEVA))
      WRITE(NOMVAR(IPP),'(A6,I2.2)')'Ga_EVA'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IGMHTF))
      WRITE(NOMVAR(IPP),'(A6,I2.2)')'Ga_HET'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
C
C ---> Variables algebriques propres a la phase continue
C
      IPP = IPPPRO(IPPROC(ITEMP1))
      NOMVAR(IPP)   = 'Temp_GAZ'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IROM1))
      NOMVAR(IPP)   = 'ROM_GAZ'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IYM1(1)))
      NOMVAR(IPP)   = 'YM_FOV'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IYM1(2)))
      NOMVAR(IPP)   = 'YM_CO'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IYM1(3)))
      NOMVAR(IPP)   = 'YM_H2S'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IYM1(4)))
      NOMVAR(IPP)   = 'YM_H2O'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IYM1(5)))
      NOMVAR(IPP)   = 'YM_CO2'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IYM1(6)))
      NOMVAR(IPP)   = 'YM_SO2'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IYM1(7)))
      NOMVAR(IPP)   = 'YM_O2'
      ICHRVR(IPP)   = 1
      ILISVR(IPP)   = 1
      IHISVR(IPP,1) = -1
      IPP = IPPPRO(IPPROC(IYM1(8)))
      NOMVAR(IPP)   = 'YM_N2'
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
C 3. INFORMATIONS COMPLEMENTAIRES
C=======================================================================
C
C ---> Definition des pointeurs du tableau TBMCR utilise dans fuphy1.F
C      et dans les sous-programmes appeles
C
      IS    = 0
      IS    = IS + 1
C     IFOV = IS
C
C ---> Initialisation
      IPHAS = 1
C ---- Calcul de RO0 a partir de T0 et P0
C        (loi des gaz parfaits applliquee a l'air)
C
      WMOLME = (WMOLE(IO2)+XSI*WMOLE(IN2)) / (1.D0+XSI)
      RO0(IPHAS) = P0(IPHAS)*WMOLME / (RR*T0(IPHAS))
C
C ---- Initialisation pour la masse volumique du coke
C    Maintenant c'est fait dans FULECD
C     RHOKF = RHO0FL
C
C ---> Coefficient de relaxation de la masse volumique
C      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
      SRROM = 0.90D0
C
C ---> Viscosite laminaire associee au scalaire enthalpie
C       DIFTL0 (diffusivite dynamique en kg/(m s))
C     C'est cette valeur par defaut qui est TOUJOURS utilisee dans les
C       calculs charbon (un peu etonnant de ne pas prendre en
C       compte les variations de ce parametre physique si on
C       recherche des informations sur les flux thermiques aux parois)
C
      DIFTL0      = 4.25D-5
      VISLS0(IHM) = DIFTL0
C
C ---> Masse volumique variable et viscosite constante (pour les suites)
      IROVAR(IPHAS) = 1
      IVIVAR(IPHAS) = 0
C
C=======================================================================
C 4. ON REDONNE LA MAIN A L'UTLISATEUR
C=======================================================================
C
      CALL USFUI1
C     ===========
C
C=======================================================================
C 5. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
C=======================================================================
C
      IOK = 0
      CALL FUVERI (IOK)
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
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (FUEL) DEMANDEE                   ',/,
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
     &'@    PHYSIQUE PARTICULIERE (FUEL) DEMANDEE                   ',/,
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
     &'                                                    (usfui1).',/)
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
     &'@  Verifier usfui1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
      END
c@z
