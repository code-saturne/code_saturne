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
                        SUBROUTINE INIUSI(IVERIF)
C                       *****************
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC  ROUTINE APPELANT LES ROUTINES UTILISATEUR POUR L'ENTREE DES
CFONC    PARAMETRES DE CALCUL : ON PASSE ICI POUR TOUT CALCUL
CFONC
CFONC  CETTE ROUTINE PERMET DE CACHER A L'UTILISATEUR LES APPELS
CFONC    A VARPOS ET AU LECTEUR XML DE L'IHM
CFONC
CFONC  LE DECOUPAGE DE L'ANCIEN USINI1 PERMET EGALEMENT DE MIEUX
CFONC    CONTROLER LES ZONES OU SONT INITIALISES LES VARIABLES (PAR
CFONC    LE BIAIS DE PARAMETRES PASSES EN ARGUMENT)
CFONC
CFONC  ON INITIALISE EGALEMENT ICI NPHAS ET ISCAPH QUI VALENT 1 EN
CFONC    PRATIQUE DANS TOUS LES CALCULS (C'EST UN PARAMETRE UTILISATEUR,
CFONC    MAIS SA VALEUR EST QUELQUE PEU CONTRAINTE ENCORE...)
CFONC
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IVERIF       ! E  !  -> ! INDICATEUR DES TESTS ELEMENTAIRES    !
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
      INCLUDE "cstnum.h"
      INCLUDE "dimens.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "mltgrd.h"
      INCLUDE "cstphy.h"
      INCLUDE "entsor.h"
      INCLUDE "vector.h"
      INCLUDE "albase.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ihmpre.h"
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
      INTEGER          IVERIF
C
C VARIABLES LOCALES
C
      INTEGER          II, IPHAS , ISCAL , NMODPP, IIIIII
      INTEGER          NPHMAX, NSCMAX, NESMAX, NPHUSI, NSCUSI
      INTEGER          IEEPRE, IEEDER, IEECOR, IEETOT, IIHMPU
      INTEGER          IALGCE
      INTEGER          IAPPEL
      DOUBLE PRECISION RELAXP, EXTRAP
C
C***********************************************************************
C
C=======================================================================
C 0. INITIALISATION DE L'INFORMATION "FICHIER XML (IHM) REQUIS & EXISTE"
C=======================================================================
C
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
C     Avec Xml, on regarde si le fichier a ete ouvert (requis et existe,
C       selon les tests realises dans cs_main)
C     IIHMPR a ete initialise a 0 juste avant (INIINI)
C
      CALL CSIHMP(IIHMPR)
C
#endif
C
C=======================================================================
C 1. INITIALISATION DE PARAMETRES DEPENDANT DU NOMBRE DE PHASES
C=======================================================================
C
C --- Nombre de phases
C
C     Egal a 1 en version 1.2
C     =======================
C
C     Le nombre de phases maximal est donne par NPHSMX dans paramx.h
C     On teste la valeur donnee par l'utilisateur avant de completer
C       les tableaux dimensionnes a NPHSMX (ex. ITURB(NPHSMX)).
C
      NPHAS = 1
C
C
C --- Varpos
C     Verification du nombre de phases
C      1er passage
      CALL VARPOS(IIIIII)
C     ===========
C
C --- Parametres dependant du nombre de phases
C
C     Turbulence
C     Chaleur massique variable ou non
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
      IF(IIHMPR.EQ.1) THEN
C
        CALL CSTURB(ITURB, IDEUCH, IGRAKE, IGRARI, XLOMLG)
C       ===========
C
        CALL CSCPVA(ICP)
C       ===========
C
      ENDIF
C
#endif
C
C
C   - Sous-programme utilisateur
C     ==========================
C
      NPHMAX = NPHSMX
      NPHUSI = NPHAS
      IIHMPU = IIHMPR
      CALL USIPPH
C     ===========
     & (NPHMAX , NPHUSI , IIHMPU , NFECRA , ITURB  , ICP , IVERIF)
C
C=======================================================================
C 2. INITIALISATION DE PARAMETRES DEPENDANT DU NOMBRE DE SCALAIRES
C=======================================================================
C
C --- Nombre de scalaires utilisateurs
C
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
      IF(IIHMPR.EQ.1) THEN
C
        CALL CSNSCA(NSCAUS)
C       ===========
C
      ENDIF
C
#endif
C
C
C   - Sous-programme utilisateur
C     ==========================
C
      IIHMPU = IIHMPR
      CALL USINSC(IIHMPU , NFECRA , NSCAUS , IVERIF)
C     ===========
C
C
C --- Dans le cas de physiques particulieres definies par des modules
C        specifiques du code tels que charbon, combustion, electrique
C        le sous-programme USPPMO doit etre complete imperativement
C        par l'utilisateur
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
      IF (IIHMPR.EQ.1) THEN
C
        CALL UIPPMO
C       ===========
     & ( IPPMOD, ICOD3P, ICODEQ, ICOEBU, ICOBML,
     &   ICOLWC, ICP3PL, ICPL3C, ICFUEL,
     &   IELJOU, IELARC, IELION, ICOMPF, INDJON, IEQCO2)
C
      ENDIF
C
#endif
C
      CALL USPPMO
C     ===========
C
C --- Varpos
C     Verification et construction de ISCAPP
C      2ieme passage
      CALL VARPOS(NMODPP)
C     ===========
C
C --- Parametres dependant du nombre de scalaires utilisateurs
C
C     Moyenne du carre des fluctuations d'un scalaire UTILISATEUR
C     Diffusivite variable ou non
C
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
      IF(IIHMPR.EQ.1) THEN
C
        CALL  CSISCA(ISCAVR)
C       ============
C
        CALL  CSIVIS(ISCAVR, IVISLS, ISCALT, ISCSTH)
C       ============
C
      ENDIF
C
#endif
C
C
C   - Sous-programme utilisateur
C     ==========================
C
      NSCMAX = NSCAMX
      NSCUSI = NSCAUS
      IIHMPU = IIHMPR
      CALL USIPSC
C     ===========
     & ( NSCMAX , NSCUSI , IIHMPU , NFECRA , ISCAVR , IVISLS , IVERIF )
C
C
C --- Parametres dependant du nombre de scalaires : exception
C     IPHSCA indique le numero de la phase porteuse pour chaque
C                                                    scalaire UTILISATEUR
C
C       Dans le cas d'une seule phase, IPHSCA(ISCAL) = 1 : ne rien changer.
C       ==================================================================
C
C       Sinon noter bien que :
C
C       La phase porteuse des scalaires UTILISATEUR ISCAL qui
C         representent la moyenne du carre des fluctuations d'un
C         scalaire utilisateur K sera la meme que celle de ce scalaire
C         utilisateur K.
C         Donc, pour de tels scalaires ISCAL (reperes par ISCAVR(ISCAL)>0),
C                          on ne doit pas renseigner IPHSCA(ISCAL) ici.
C         C'est l'objet du test inclus dans l'exemple ci-dessous.
C
C       Pour les scalaires non utilisateur relatifs a des physiques
C         particulieres, (charbon, combustion, electrique : voir usppmo)
C         implicitement definis selon le modele,
C         les informations sont donnees automatiquement par ailleurs :
C                                         on ne modifie pas IPHSCA ici.
C
      DO ISCAL = 1, NSCAUS
        IF(ISCAVR(ISCAL).LE.0) THEN
          IPHSCA(ISCAL) = 1
        ENDIF
      ENDDO
C
C
C=======================================================================
C 3. INITIALISATION DE PARAMETRES "GLOBAUX"
C=======================================================================
C
C
C --- Parametres globaux
C
C     Pas de temps
C     Couplage vitesse/pression
C     Prise en compte de la pression hydrostatique
C     Estimateurs (pas encore dans l'IHM)
C
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
      IF(IIHMPR.EQ.1) THEN
C
        CALL CSIDTV(IDTVAR)
C       ===========
C
        CALL CSIPHY(IPHYDR)
C       ===========
C
      ENDIF
C
#endif
C
C
C   - Sous-programme utilisateur
C     ==========================
C
      NPHMAX = NPHSMX
      NESMAX = NESTMX
      IEEPRE = IESPRE
      IEEDER = IESDER
      IEECOR = IESCOR
      IEETOT = IESTOT
      NPHUSI = NPHAS
      IIHMPU = IIHMPR
C     IALGCE permet de remplir la variable cs_glob_maillage_grd_cdg_cel dans
C       cs_maillage_grd.c, a travers la routine ALGCEN.
C     cs_glob_maillage_grd_cdg_cel est initialise a 0 dans cs_maillage_grd.c,
C       et on ne change sa valeur ici que si elle a vraiment ete touchee par
C       l'utilisateur (pour garder l'initialisation en un seul endroit).
C     Le blindage en erreur est dans cs_maillage_grd.c (erreur si IALGCE>1,
C       cs_glob_maillage_grd_cdg_cel inchange si IALGCE<0)
      IALGCE = -999
C
      CALL USIPGL
C     ===========
     & ( NPHMAX , NESMAX ,
     &   IEEPRE , IEEDER , IEECOR , IEETOT ,
     &   NPHUSI , IIHMPU , NFECRA ,
     &   IDTVAR , IPUCOU , IPHYDR , IALGCE , IESCAL , IVERIF )
C
      IF (IALGCE.NE.-999) CALL ALGCEN(IALGCE)
C
C --- Parametres de la methode ALE
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
      IF (IIHMPR.EQ.1) THEN
C
        CALL UIALIN (IALE, NALINF, NALIMX, EPALIM, IORTVM)
C       ===========
      ENDIF
C
#endif
C
      CALL USALIN
C     ===========
C
C --- Varpos
C     Positionnement de pointeurs
C     Verifications
C     Determination de IPR, IU ... ISCA, NVAR
C     Determination de IPP...
C
C      3ieme passage
      CALL VARPOS(NMODPP)
C     ===========
C
C
C=======================================================================
C 4. INITIALISATION DE PARAMETRES UTILISATEUR SUPPLEMENTAIRES
C=======================================================================
C
C --- Format des fichiers aval (entsor.h)
C --- Options du calcul (optcal.h)
C --- Constantes physiques (cstphy.h)
C
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
      IF(IIHMPR.EQ.1) THEN
C
        CALL CSVNUM
C       ===========
     &            (NVAR,
     &             IU, IV, IW, IPR,
     &             ITURB, IK, IEP,
     &             IR11, IR22, IR33,
     &             IR12, IR13, IR23,
     &             IOMG, IPHI, IFB,
     &             IALE, IUMA, IVMA, IWMA,
     &             ISCA, ISCAPP)
C
C     Format des fichiers aval
        CALL CSIFOA(IFOAVA, IFOAVX)
C       ===========
C
C     Suite de calcul, relecture fichier auxiliaire, champ de vitesse figé
        CALL CSISUI(ISUITE, ILEAUX, ICCVFG)
C       ===========
C
C     Pas de temps (seulement NTMABS, DTREF, INPDT0)
        CALL CSTIME
C       ===========
     &             (INPDT0, IPTLRO, NTMABS, IDTVAR, DTREF, DTMIN,
     &              DTMAX, COUMAX, FOUMAX, VARRDT, RELXST)
C
C     Temperature ou enthalpie (hors physiques particulieres)
        IF(NMODPP.EQ.0) THEN
          CALL CSSCA1(ISCALT, ISCSTH)
C         ===========
C
        ENDIF
C
C      Options numériques locales
C
        CALL UINUM1
C       ===========
     &        (ISCA, ISCAPP, BLENCV, ISCHCV, ISSTPC, IRCFLU,
     &         CDTVAR, NITMAX, EPSILO)
C
C     Options numériques globales
        RELAXP = -999.D0
        EXTRAP = 0.D0
        CALL CSNUM2 (IVISSE, RELAXP, IPUCOU, EXTRAP, IMRGRA)
C       ===========
        IPHAS = 1
        EXTRAG(IPR(IPHAS)) = EXTRAP
        IF (IDTVAR.GE.0) RELAXV(IPR(IPHAS)) = RELAXP
C
C     Gravite, prop. phys
        CALL CSPHYS
C       ===========
     &             (NMODPP,
     &              IROVAR, IVIVAR,
     &              GX, GY, GZ,
     &              RO0, VISCL0, CP0, T0, P0)
C
C     Scamin, scamax
        CALL CSSCA2(ISCAVR, SCAMIN, SCAMAX)
C       ===========
C
C     Diffusivites
        CALL CSSCA3(ISCALT, ISCAVR, VISLS0, T0, P0)
C       ===========
C
C     Init turb (uref, almax) si necessaire (modele RANS)
        IPHAS = 1
        IF (ITYTUR(IPHAS).EQ.2 .OR. ITYTUR(IPHAS).EQ.3 .OR.
     &      ITYTUR(IPHAS).EQ.5 .OR. ITYTUR(IPHAS).EQ.6 ) THEN
          CALL CSTINI(UREF, ALMAX)
C         ===========
        ENDIF
C
        IAPPEL = 0
C
        CALL UIPROP
C       ===========
     &            (IROM, IVISCL, IVISCT, IVISLS, ICOUR, IFOUR,
     &             ISMAGO, IALE, ICP, ISCALT, ISCAVR,
     &             IPRTOT, IPPPRO, IPPROC, ICMOME,
     &             IPPTX, IPPTY, IPPTZ, IPPDT,
     &             IVISMA, IAPPEL)
C
        CALL UIMOYT
C       ===========
     &            (NDGMOX, ISCA,
     &             IPPPRO, IPPROC,  ICMOME,
     &             NTDMOM, IMOOLD, IDFMOM)
C
      ENDIF
C
#endif
C
C
C   - Sous-programme utilisateur
C     ==========================
C
      CALL USIPSU(NMODPP , IVERIF)
C     ===========
C
C
C --- Varpos
C      4ieme passage
      CALL VARPOS(IIIIII)
C     ===========
C
C
C=======================================================================
C 5. INITIALISATION DE PARAMETRES UTILISATEUR (entree sorties)
C=======================================================================
C
C --- Entree-sorties
C
C
#if defined(_CS_HAVE_XML)
C
C   - Interface Code_Saturne
C     ======================
C
C
      IF(IIHMPR.EQ.1) THEN
C
        IF(NBMOMT.GT.0 .OR. IPUCOU.EQ.1) THEN
          IAPPEL = 1
C
          CALL UIPROP
C         ===========
     &            (IROM, IVISCL, IVISCT, IVISLS, ICOUR, IFOUR,
     &             ISMAGO, IALE, ICP, ISCALT, ISCAVR,
     &             IPRTOT, IPPPRO, IPPROC, ICMOME,
     &             IPPTX, IPPTY, IPPTZ, IPPDT,
     &             IVISMA, IAPPEL)
C
        ENDIF
C
        DO II = 1,NVPPMX
          CALL FCNMVA (NOMVAR(II), LEN(NOMVAR(II)), II)
C         ===========
        ENDDO
C
        CALL CSENSO
C       ===========
     &     ( NVPPMX, NCAPT,  NTHIST, NTLIST,
     &       ICHRVL, ICHRBO, ICHRSY, ICHRMD,
     &       FMTCHR, LEN(FMTCHR), OPTCHR, LEN(OPTCHR),
     &       NTCHR,  IECAUX,
     &       ICHRVR, ILISVR, IHISVR, ISCA, ISCAPP,
     &       IPPRTP, IPPPRO, IPPROC,
     &       XYZCAP )
C
        DO II = 1,NVPPMX
          CALL CFNMVA(NOMVAR(II), LEN(NOMVAR(II)), II)
C         ===========
        ENDDO
C
        CALL UIUSAR(ICOFTU)
C       ===========
C
      ENDIF
C
#endif
C
C   - Sous-programme utilisateur
C     ==========================
C
      CALL USIPES(NMODPP, IVERIF)
C     ===========
C
C----
C FORMATS
C----
C
C
      RETURN
      END
c@z
