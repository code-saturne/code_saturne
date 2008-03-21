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
                        SUBROUTINE COINI1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC          INIT DES OPTIONS DES VARIABLES POUR
CFONC               POUR LA COMBUSTION
CFONC         FLAMME DE DIFFUSION ET DE PREMELANGE
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
      INCLUDE "coincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "ppincl.h"
C
C***********************************************************************
C
C VARIABLES LOCALES
C
      INTEGER          IPP , II , JJ, IOK, IDIRAC
      INTEGER          ISC , IPHAS
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
C=======================================================================
C 1. VARIABLES TRANSPORTEES
C=======================================================================
C
C 1.1 Definition des scamin et des scamax des variables transportees
C ==================================================================
C
C --> Flamme de diffusion : chimie 3 points
C
      IF ( IPPMOD(ICOD3P).GE.0 ) THEN
C
C ---- Taux de melange
        SCAMIN(IFM) = 0.D0
        SCAMAX(IFM) = 1.D0
C
C ---- Variance du taux de melange
C        Type de clipping superieur pour la variance
C        0 pas de clipping, 1 clipping var max de fm, 2 clipping a SCAMAX
       ICLVFL(IFP2M) = 1
C        SCAMIN(IFP2M) = 0.D0
C        SCAMAX(IFP2M) = 0.25D0
C
C ---- Enthalpie
        IF ( IPPMOD(ICOD3P).EQ.1 ) THEN
          SCAMIN(IHM) = -GRAND
          SCAMAX(IHM) = +GRAND
        ENDIF
C
      ENDIF
C
C
C --> Flamme de premelange : modele EBU
C
      IF ( IPPMOD(ICOEBU).GE.0 ) THEN
C
C ---- Fraction massique des gaz frais
        SCAMIN(IYGFM) = 0.D0
        SCAMAX(IYGFM) = 1.D0
C
C ---- Taux de melange
        IF ( IPPMOD(ICOEBU).EQ.2 .OR.
     &       IPPMOD(ICOEBU).EQ.3       ) THEN
          SCAMIN(IFM) = 0.D0
          SCAMAX(IFM) = 1.D0
        ENDIF
C
C ---- Enthalpie
        IF ( IPPMOD(ICOEBU).EQ.1 .OR.
     &       IPPMOD(ICOEBU).EQ.3      ) THEN
          SCAMIN(IHM) = -GRAND
          SCAMAX(IHM) = +GRAND
        ENDIF
C
      ENDIF
C
C
C --> Flamme de premelange : modele LWC
C
      IF ( IPPMOD(ICOLWC).GE.0 ) THEN
        SCAMIN(IFM) = 0.D0
        SCAMAX(IFM) = 1.D0
C
        ICLVFL(IFP2M) = 0
C
        SCAMIN(IYFM) = 0.D0
        SCAMAX(IYFM) = 1.D0
C
        ICLVFL(IYFP2M) = 0
C
        IF ( IPPMOD(ICOLWC).GE.2 ) THEN
          SCAMIN(ICOYFP) =-0.25D0
          SCAMAX(ICOYFP) = 0.25D0
        ENDIF
C
      ENDIF
C
C 1.2 Nature des scalaires transportes
C ====================================
C
      DO ISC = 1, NSCAPP
C
C
C ---- Type de scalaire (0 passif, 1 temperature en K
C                                 -1 temperature en C
C                                  2 enthalpie)
C      La distinction -1/1 sert pour le rayonnement
        ISCSTH(ISCAPP(ISC)) = 0
C
      ENDDO
C
C
C 1.3 L'utilisation de la variable enthalpie necessite un traitement
C     particulier
C ==================================================================
C
C
C ---- On resout en enthalpie avec un CP constant (Cf. cpvarp)
      IF ( IPPMOD(ICOD3P).EQ.1 .OR.
     &     IPPMOD(ICOEBU).EQ.1 .OR.
     &     IPPMOD(ICOEBU).EQ.3 .OR.
     &     IPPMOD(ICOLWC).EQ.1 .OR.
     &     IPPMOD(ICOLWC).EQ.3 .OR.
     &     IPPMOD(ICOLWC).EQ.5     ) THEN
C
        IPHAS = IPHSCA(IHM)
C
        ISCALT(IPHAS) = IHM
        ISCSTH(IHM) = 2
C
      ENDIF
C
C
C 1.4 Donnees physiques ou numeriques propres aux scalaires COMBUSTION
C ====================================================================
C
      DO ISC = 1, NSCAPP
C
        JJ = ISCAPP(ISC)
C
        IF ( ISCAVR(JJ).LE.0 ) THEN
C
C ---- En combustion on considere que la viscosite turbulente domine
C      ON S'INTERDIT DONC LE CALCUL DES FLAMMES LAMINAIRES AVEC Le =/= 1
C
          VISLS0(JJ) = VISCL0(IPHSCA(JJ))
C
        ENDIF
C
C ---- Schmidt ou Prandtl turbulent
C
        SIGMAS(JJ) = 0.7D0
C
C ---- Coeff dissipation des fluctuations
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
C ---- Informations relatives a la resolution des scalaires
C
C       - Facteur multiplicatif du pas de temps
C
        CDTVAR(II) = 1.D0
C
C         - Schema convectif % schema 2ieme ordre
C           = 0 : upwind
C           = 1 : second ordre
        BLENCV(II) = 1.D0
C
C         - Type de schema convetif second ordre (utile si BLENCV > 0)
C           = 0 : Second Order Linear Upwind
C           = 1 : Centre
        ISCHCV(II) = 1
C
C         - Test de pente pour basculer d'un schema centre vers l'upwind
C           = 0 : utilisation automatique du test de pente
C           = 1 : calcul sans test de pente
        ISSTPC(II) = 0
C
C         - Reconstruction des flux de convetion et de diffusion aux faces
C           = 0 : pas de reconstruction
        IRCFLU(II) = 1
C
      ENDDO
C
C
C 1.5 Variable courante : nom, sortie chrono, suivi listing, sortie histo
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
C =======================================================================
C
C --> Flamme de diffusion : chimie 3 points - chmie equilibre
C
      IF ( IPPMOD(ICOD3P).GE.0 ) THEN
C
C ---- Taux de melange
        IPP = IPPRTP(ISCA(IFM))
        NOMVAR(IPP)  = 'Fra_MEL'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
C ---- Variance du taux melange
        IPP = IPPRTP(ISCA(IFP2M))
        NOMVAR(IPP)  = 'Var_FrMe'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C
      ENDIF
C
C
C --> Flamme de premelange : Modele EBU
C
      IF ( IPPMOD(ICOEBU).GE.0 ) THEN
C
C ---- Fraction massique des gaz frais
        IPP = IPPRTP(ISCA(IYGFM))
        NOMVAR(IPP)  = 'Fra_GF'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
C ---- Taux de melange
      IF ( IPPMOD(ICOEBU).GE.2 ) THEN
        IPP = IPPRTP(ISCA(IFM))
        NOMVAR(IPP)  = 'Fra_MEL'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
      ENDIF
C
C
C --> Flamme de premelange : Modeles BML et LWC
C
      IF ( IPPMOD(ICOLWC).NE.-1 ) THEN
C --- Taux de melange
        IPP = IPPRTP(ISCA(IFM))
        NOMVAR(IPP) = 'Fra_Mel'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C --- Variance du taux de melange
        IPP = IPPRTP(ISCA(IFP2M))
        NOMVAR(IPP) = 'Var_FMe'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C --- Fraction massique
        IPP = IPPRTP(ISCA(IYFM))
        NOMVAR(IPP) = 'Fra_Mas'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C --- Variance de la fraction massique
        IPP = IPPRTP(ISCA(IYFP2M))
        NOMVAR(IPP) = 'Var_FMa'
        ICHRVR(IPP)  = 1
        ILISVR(IPP)  = 1
        IHISVR(IPP,1)= -1
C --- Covariance
        IF (IPPMOD(ICOLWC).GE.2) THEN
          IPP = IPPRTP(ISCA(ICOYFP))
          NOMVAR(IPP) = 'COYF_PP4'
          ICHRVR(IPP)  = 1
          ILISVR(IPP)  = 1
          IHISVR(IPP,1)= -1
        ENDIF
C
      ENDIF
C --> Si Transport de l'enthalpie
C
       IF ( IPPMOD(ICOD3P).EQ.1 .OR.
     &      IPPMOD(ICOEBU).EQ.1 .OR.
     &      IPPMOD(ICOEBU).EQ.3 .OR.
     &      IPPMOD(ICOLWC).EQ.1 .OR.
     &      IPPMOD(ICOLWC).EQ.3 .OR.
     &      IPPMOD(ICOLWC).EQ.5   ) THEN
         IPP = IPPRTP(ISCA(IHM))
         NOMVAR(IPP)  = 'Enthalpy'
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
C --> Flamme de diffusion :
C
      IF ( IPPMOD(ICOD3P).GE.0 ) THEN
        IPP = IPPPRO(IPPROC(ITEMP))
        NOMVAR(IPP)   = 'Temperature'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IYM(1)))
        NOMVAR(IPP)   = 'YM_Fuel'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IYM(2)))
        NOMVAR(IPP)   = 'YM_Oxyd'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IYM(3)))
        NOMVAR(IPP)   = 'YM_Prod'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
      ENDIF
C
C ---> Physique particuliere : Flamme de premelange (modele EBU)
C
      IF ( IPPMOD(ICOEBU).GE.0 ) THEN
        IPP = IPPPRO(IPPROC(ITEMP))
        NOMVAR(IPP)   = 'Temperature'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IYM(1)))
        NOMVAR(IPP)   = 'YM_Fuel'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IYM(2)))
        NOMVAR(IPP)   = 'YM_Oxyd'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IYM(3)))
        NOMVAR(IPP)   = 'YM_Prod'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
      ENDIF
C
C ---> Physique particuliere : Flamme de premelange (modele LWC)
C
      IF ( IPPMOD(ICOLWC).GE. 0 ) THEN
        IPP = IPPPRO(IPPROC(ITSC))
        NOMVAR(IPP)   = 'T.SOURCE'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(ITEMP))
        NOMVAR(IPP)   = 'Temperature'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IYM(1)))
        NOMVAR(IPP)   = 'YM_Fuel'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IYM(2)))
        NOMVAR(IPP)   = 'YM_Oxyd'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
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
      ENDIF
C
C ---> Physique particuliere : Flamme de diffusion (chimie 3 points)
C                              Flamme de premelange (modele EBU 1 et 3)
C                              Flamme de premelange (modele LWC 1, 3 et 5)
C                              AVEC RAYONNEMENT
C
      IF ( ( IPPMOD(ICOD3P).EQ.1 .OR.
     &       IPPMOD(ICOEBU).EQ.1 .OR. IPPMOD(ICOEBU).EQ.3 .OR.
     &       IPPMOD(ICOLWC).EQ.1 .OR. IPPMOD(ICOLWC).EQ.3 .OR.
     &       IPPMOD(ICOLWC).EQ.5 )
     &      .AND. (IRAYPP.GE.1) ) THEN
        IPP = IPPPRO(IPPROC(ICKABS))
        NOMVAR(IPP)   = 'KABS'
        ICHRVR(IPP)   = 1
        ILISVR(IPP)   = 1
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IT4M))
        NOMVAR(IPP)   = 'TEMP4'
        ICHRVR(IPP)   = 0
        ILISVR(IPP)   = 0
        IHISVR(IPP,1) = -1
        IPP = IPPPRO(IPPROC(IT3M))
        NOMVAR(IPP)   = 'TEMP3'
        ICHRVR(IPP)   = 0
        ILISVR(IPP)   = 0
        IHISVR(IPP,1) = -1
      ENDIF
C
C
C=======================================================================
C 3. INFORMATIONS COMPLEMENTAIRES
C=======================================================================
C
C --> Calcul de RO0 a partir de T0 et P0
C
       IPHAS = 1
       IF ( IPPMOD(ICOD3P).NE.-1 .OR.
     &      IPPMOD(ICOEBU).NE.-1 .OR.
     &      IPPMOD(ICOLWC).NE.-1     ) THEN
         WMOLME = WMOLG(2)
         RO0(IPHAS) = P0(IPHAS)*WMOLME / (RR*T0(IPHAS))
       ENDIF
C
C On met les constantes a -GRAND pour obliger l'utilisateur a les definir
C  (ou les laisser) dans usebu1, usd3p1 ou uslwc1.
C --> Constante modele EBU par defaut
      CEBU   =-GRAND
C
C --> Constantes modele LWC par defaut
      VREF  =-GRAND
      LREF  =-GRAND
      TA    =-GRAND
      TSTAR =-GRAND
C
C --> Coefficient de relaxation de la masse volumique
C      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
      SRROM =-GRAND
C
C --> Viscosite laminaire associee au scalaire enthalpie
C       DIFTL0 (diffusivite dynamique en kg/(m s))
      DIFTL0      =-GRAND
C
C --> Diffusion 3 points, tableaux HH HF THF
C           (generes dans d3pphy)
C
      NMAXF = 9
      NMAXH = 9
C
C ---> Masse volumique variable et viscosite constante (pour les suites)
      IROVAR(IPHAS) = 1
      IVIVAR(IPHAS) = 0
C
C=======================================================================
C 4. ON REDONNE LA MAIN A L'UTLISATEUR
C=======================================================================
C
      IF ( IPPMOD(ICOEBU).GE.0 ) THEN
        CALL USEBU1
C       ===========
      ELSE IF( IPPMOD(ICOD3P).GE.0 ) THEN
        CALL USD3P1
C       ===========
      ELSE IF( IPPMOD(ICOLWC).GE.0 ) THEN
        CALL USLWC1
C       ==========
      ENDIF
C
C=======================================================================
C 5. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
C=======================================================================
C
      IOK = 0
      IF ( IPPMOD(ICOEBU).GE.0 ) THEN
C
        CALL EBUVER (IOK)
C       ===========
        IF(IOK.GT.0) THEN
          WRITE(NFECRA,9999)IOK
          CALL CSEXIT (1)
C         ===========
        ELSE
          WRITE(NFECRA,9998)
        ENDIF
C
      ELSE IF( IPPMOD(ICOD3P).GE.0 ) THEN
        CALL D3PVER (IOK)
C       ===========
        IF(IOK.GT.0) THEN
          WRITE(NFECRA,9991)IOK
          CALL CSEXIT (1)
        ELSE
          WRITE(NFECRA,9990)
        ENDIF
C
      ELSE IF( IPPMOD(ICOLWC).GE.0 ) THEN
        CALL LWCVER (IOK)
C       ===========
        IF(IOK.GT.0) THEN
          WRITE(NFECRA,9993)IOK
          CALL CSEXIT (1)
        ELSE
          WRITE(NFECRA,9992)
        ENDIF
C
      ENDIF
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
     &'@    *********                                               ',/,
     &'@    PHYSIQUE PARTICULIERE (COMBUSTION) DEMANDEE             ',/,
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
     &'@    PHYSIQUE PARTICULIERE (COMBUSTION) DEMANDEE             ',/,
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
     &'                                                    (usebu1).',/)
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
     &'@  Verifier usebu1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9990 FORMAT(
     &'                                                             ',/,
     &' Pas d erreur detectee lors de la verification des donnees   ',/,
     &'                                                    (usd3p1).',/)
 9991 FORMAT(
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
     &'@  Verifier usd3p1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 9992 FORMAT(
     &'                                                             ',/,
     &' Pas d erreur detectee lors de la verification des donnees   ',/,
     &'                                                    (uslwc1).',/)
 9993 FORMAT(
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
     &'@  Verifier uslwc1.                                          ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
      RETURN
      END
c@z
