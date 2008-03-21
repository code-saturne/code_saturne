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
                        SUBROUTINE LAGOPT
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C
C***********************************************************************
C  FONCTION  :
C  ---------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC   1) Initialisation par defaut du parametrage du module
CFONC      lagrangien
CFONC
CFONC   2) Lecture du parametrage utilisateur
CFONC
CFONC   3) Verifications du parametrage utilisateur et
CFONC      controles de coherence
CFONC
CFONC   4) Initialisation des variables en COMMON et des pointeurs
CFONC      sur les tableaux lies aux particules, aux statistiques,
CFONC      aux conditions aux limites, aux variables parietales,
CFONC      aux donnees pour le couplage retour.
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
      INCLUDE "entsor.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagdim.h"
      INCLUDE "lagran.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "cpincl.h"
      INCLUDE "radiat.h"
C
C***********************************************************************
C
      INTEGER          II , IP , IRF , IPHAS , ICHA , I1 , I2 , I3, IOK
C
C***********************************************************************
C=======================================================================
C 0. INITIALISATION
C=======================================================================
C
C  Vérification que la norme IEEE 754 est respectée par l'architecture.
C  Dans le cas contraire la trajectographie des particules du module
C  lagrangien est susceptible de ne pas fonctionner correctement.
C
      CALL CSIEEE
C     ===========
C
C=======================================================================
C 1. INITIALISATIONS PAR DEFAUT DU MODULE LAGRANGIEN
C                        ^^^^^^
C    TOUTES CES INITIALISATIONS PRENNENT DES VALEURS NON ADMISSIBLES
C    A LA VERIFICATION POUR ETRE SUR QUE L'UTILISATEUR LES A MODIFIEES
C=======================================================================
C
C     ILPHAS : numero de la phase continue qui est la phase porteuse
C              des particules
C
      ILPHAS = 1
C
C
C     IILAGR = 0 : PAS DE CALCUL LAGRANGIEN
C            = 1 : DIPHASIQUE LAGRANGIEN SANS COUPLAGE RETOUR
C            = 2 : DIPHASIQUE LAGRANGIEN AVEC COUPLAGE RETOUR
C            = 3 : DIPHASIQUE LAGRANGIEN SUR CHAMP FIGE
C
      IILAGR = -1
C
C     ISUILA = 0 : PAS DE SUITE LAGRANGIENNE
C            = 1 : SUITE LAGRANGIENNE
C
      ISUILA = -1
C
C     Suite statistiques et TS couplage retour si =1
C
      ISUIST = -1
C
C     IPHYLA = 0 : PUREMENT DYNAMIQUE
C            = 1 : EQUATION SUR TP, DP et MP
C            = 2 : CHARBON
C
      IPHYLA = -1
C
C     EQUATION SUR LE DIAMETRE (O : NON , 1 : OUI)
C
      IDPVAR = -1
C
C     EQUATION SUR LA TEMPERATURE (O : NON , 1 : OUI)
C
      ITPVAR = -1
C
C     EQUATION SUR LA MASSE (O : NON , 1 : OUI)
C
      IMPVAR = -1
C
C     TEMPERATURE D'INITIALISATION DES PARTICULES DEJA PRESENTES
C
      TPART = -999.D0
C
C     CHALEUR MASSIQUE D'INITIALISATION DES PARTICULES DEJA PRESENTES
C
      CPPART = -999.D0
C
C     ENCRASSEMENT
C     IENCRA = 0 pas d'encrassement
C            = 1 encrassement
C     TPRENC : TEMPERATURE MINIMUM POUR l'ENCRASSEMENT
C     VISREF : VISCOSITE CRITIQUE
C
      IENCRA = -1
C
      DO ICHA = 1 , NCHARM2
        TPRENC(ICHA) = -999.D0
        VISREF(ICHA) = -999.D0
      ENDDO
C
C     NOMBRE DE PARTICULES MAXIMAL AUTORISE DANS LE DOMAINE
C
      NBPMAX = -1
C
C     NOMBRE DE VARIABLES SUPPLEMENTAIRE SUR LES PARTICULES
C
      NVLS = -1
C
C     CARACTERE STATIONNAIRE DE L'ECOULEMENT DE LA PHASE CONTINUE
C
      ISTTIO = -1
C
C     Nombre de passages absolus (i.e. suite comprise) avant de faire
C     une moyenne en temps (CALCUL STATIONNAIRE) des termes sources
C     de couplage retour (seulement si ISTTIO = 1)
C
      NSTITS = -1
C
C     COUPLAGE RETOUR SUR LA DYNAMIQUE (Vitesse + Turbulence)
C
      LTSDYN = -1
C
C     COUPLAGE UTILE SUR LA MASSE
C
      LTSMAS = -1
C
C     COUPLAGE UTILE SUR LA THERMIQUE OU LES VARIABLES CHARBON
C
      LTSTHE = -1
C
C     STATISTIQUES
C
C     Calcul statistiques si  = 1 ; si 0 pas de stat ; si -1 stop
C
      ISTALA = -1
C
C     Nombre de variables statistiques supplementaires
C
      NVLSTS = -1
C
C     Numero du pas de temps pour debut statistque
C
      IDSTNT = -1
C
C     Debut calcul stationnaire
C
      NSTIST = -1
C
C     Seuil en POIDS STAT de particules pour les stats
C
      SEUIL = -1.D0
C
C     Nom des variables
C
      DO II = 1,NVPLMX
C
        WRITE(NOMLAG(II),'(A6,I4.4)') 'MoyLag',II
C
        WRITE(NOMLAV(II),'(A6,I4.4)') 'VarLag',II
C
      ENDDO
C
C     Historique pour la moyenne et la variance
C
      DO II = 1,NVPLMX
        IHSLAG(II)  = 0
      ENDDO
C
C     INJECTION CONTINUE
C
      INJCON = -1
C
C     ROULETTE RUSSE
C
      IROULE = -1
C
C     ORDRE D'INTEGRATION
C
      NORDRE = -1
C
C     DISPERSION TURBULENTE
C
      IDISTU = -1
C
C     DIFFUSION TURBULENTE
C
      IDIFFL = -1
C
C     DISPERSION TURBULENTE
C
      MODCPL = -1
C
C     DISPERSION TURBULENTE
C
      IDIRLA = 0
C
C     CORRECTION DE PRESSION
C
      ILAPOI = -1
C
C     POSTPROCESSING MODE : TRAJECTOIRES
C
      IENSI1 = -1
C
C     POSTPROCESSING MODE : DEPLACEMENTS
C
      IENSI2 = -1
C
C     NOMBRE DE PARTICULES A VISUALISER MAXIMUM=NLISTE
C
      NBVIS  = -1
C
C     FREQUENCE D'AQUISITION DES DONNEES A VISUALISER
C
      NVISLA = -1
C
C     INITIALISATION  : PAR DEFAUT ON NE VISUALISE AUCUNE PARTICULE
C
      DO II = 1, NLISTE
        LISTE(II) = -1
        LIST0(II) = -1
        NPLIST(II) = 0
      ENDDO
C
C     POSTPROCESSING VARIABLE : VITESSE DU FLUIDE VU
C
      IVISV1  = -1
C
C     POSTPROCESSING VARIABLE : VITESSE DE LA PARTICULE
C
      IVISV2  = -1
C
C     POSTPROCESSING VARIABLE : TEMPS DE SEJOUR
C
      IVISTP  = -1
C
C     POSTPROCESSING VARIABLE : DIAMETRE
C
      IVISDM  = -1
C
C     POSTPROCESSING VARIABLE : MASSE
C
      IVISMP  = -1
C
C     POSTPROCESSING VARIABLE : TEMPERATURE
C
      IVISTE  = -1
C
C     POSTPROCESSING VARIABLE : TEMPERATURE CHARBON
C
      IVISHP  = -1
C
C     POSTPROCESSING VARIABLE : DIAMETRE DU COEUR RETRECISSANT
C
      IVISDK  = -1
C
C     POSTPROCESSING VARIABLE : MASSE CHARBON REACTIF
C
      IVISCH  = -1
C
C     POSTPROCESSING VARIABLE : MASSE COKE
C
      IVISCK  = -1
C
C     POSTPROCESSING MODE : INTERACTIONS PARTICULES/FRONTIERES
C
      IENSI3 = -1
C
C     POSTPROCESSING : STAT PARIETALES STATIONNAIRES
C
      NSTBOR = -1
C
C     INFORMATIONS A ENREGISTRER
C
      INBRBD = -1
      IFLMBD = -1
      IANGBD = -1
      IVITBD = -1
      IENCBD = -1
C
      NUSBOR = -1
C
C     Type de moyenne applicable pour affichage et post-processing
C
      DO II = 1,NUSBRD+10
        IMOYBR(II) = -1
      ENDDO
C
C=======================================================================
C 2.1 INITIALISATIONS UTILISATEUR DU MODULE LAGRANGIEN
C                     ^^^^^^^^^^^
C=======================================================================
C
      CALL USLAG1
C     ===========
C
      IF (IILAGR.EQ.0) RETURN
C
C=======================================================================
C 2.2 VERIFICATION DES INITIALISATIONS UTILISATEUR DU MODULE LAGRANGIEN
C=======================================================================
C
C on doit verifier toutes les options entrees par l'utilisateur
c  qui est inventif
C  et si une valeur ne correspond pas, on ne corrige pas, on s'arrete.
C
      IOK = 0
C
C     IILAGR
C
      IF (IILAGR.LT.0 .OR. IILAGR.GT.3) THEN
        WRITE(NFECRA,1010) IILAGR
        IOK = IOK + 1
      ENDIF
C
      IF (IILAGR.GT.0
     &    .AND. (NNOD.EQ.0 .OR. (LNDFAC.EQ.0 .AND. LNDFBR.EQ.0))) THEN
        WRITE(NFECRA,1011) NNOD, LNDFAC, LNDFBR
        IOK = IOK + 1
      ENDIF
C
C    CALCUL SUR CHAMP FIGE : SUITE OBLIGATOIRE
C    ATTENTION : LE CHAMP FIGE LAGRANGIEN N"EST PAS LE MEME QUE CELUI
C    DE L'EULERIEN. POUR LE LAGRANGIEN VITESSE PRESSION ET SCALAIRES
C    SONT CONSTANTS.
C
      IF (IILAGR.EQ.3 .AND. ISUITE.NE.1) THEN
        WRITE(NFECRA,1012) IILAGR, ISUITE
        IOK = IOK + 1
      ENDIF
C
      IF (IILAGR.EQ.3) ICCVFG = 1
C
      IF (IILAGR.NE.2 .AND. IPPMOD(ICPL3C).GE.1) THEN
        WRITE(NFECRA,1013) IILAGR, IPPMOD(ICPL3C)
        IOK = IOK + 1
      ENDIF
C
      IF (IILAGR.GT.0 .AND. IDTVAR.EQ.2) THEN
        WRITE(NFECRA,1014) IILAGR, IDTVAR
        IOK = IOK + 1
      ENDIF
C
C     ISUILA ISUIST
C
      IF (ISUILA.LT.0 .OR. ISUILA.GT.1) THEN
        WRITE(NFECRA,1020) ISUILA
        IOK = IOK + 1
      ENDIF
C
      IF (ISUILA.EQ.1 .AND. ISUITE.EQ.0) THEN
        WRITE(NFECRA,1021)
        IOK = IOK + 1
      ENDIF
C
      IF (ISUILA.EQ.1) THEN
        IF (ISUIST.LT.0 .OR. ISUIST.GT.1) THEN
          WRITE(NFECRA,1022) ISUIST
          IOK = IOK + 1
        ENDIF
      ELSE
        ISUIST = 0
      ENDIF
C
C     IPHYLA
C
      IF (IPHYLA.LT.0 .OR. IPHYLA.GT.2) THEN
        WRITE(NFECRA,1030) IPHYLA
        IOK = IOK + 1
      ENDIF
C
      IF (IOK.NE.0) CALL CSEXIT (1)
C                   ===========
C
C     IDPVAR ITPVAR IMPVAR
C
      IF (IPHYLA.EQ.1) THEN
        IF (IDPVAR.LT.0 .OR. IDPVAR.GT.1) THEN
          WRITE(NFECRA,1031) IDPVAR
          IOK = IOK + 1
        ENDIF
        IF (ITPVAR.LT.0 .OR. ITPVAR.GT.1) THEN
          WRITE(NFECRA,1032) ITPVAR
          IOK = IOK + 1
        ENDIF
        IF (IMPVAR.LT.0 .OR. IMPVAR.GT.1) THEN
          WRITE(NFECRA,1033) IMPVAR
          IOK = IOK + 1
        ENDIF
        II = ILPHAS
        IF (ITPVAR.EQ.1 .AND. ISCALT(II).EQ.-1) THEN
          WRITE(NFECRA,1034) ITPVAR, II, II, ISCALT(II)
          IOK = IOK + 1
        ENDIF
      ELSE
        ITPVAR = 0
        IMPVAR = 0
        IDPVAR = 0
      ENDIF
C
      IF (ISUILA.EQ.1 .AND. IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) THEN
        IF (CPPART.LT.0) THEN
          WRITE(NFECRA,1036) CPPART
          IOK = IOK + 1
        ENDIF
        IF (TPART.LT.TKELVN) THEN
          WRITE(NFECRA,1037) TKELVN, TPART
          IOK = IOK + 1
        ENDIF
      ENDIF
C
      IF (IOK.NE.0) CALL CSEXIT (1)
C                   ===========
C
C   APPEL A VARPOS SI RAYONNEMENT
C
      IF (IIRAYO.EQ.1 .AND. IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) THEN
        CALL VARPOS(IOK)
C       ===========
      ENDIF
C
C     IENCRA TPRENC VISREF
C
      IF (IPHYLA.EQ.2) THEN
        IF (IENCRA.LT.0 .OR. IENCRA.GT.1) THEN
          WRITE(NFECRA,1040) IENCRA
          IOK = IOK + 1
        ENDIF
C
        DO ICHA = 1 , NCHARB
          IF (IENCRA.EQ.1 .AND. VISREF(ICHA).LT.0) THEN
            WRITE(NFECRA,1041) IENCRA, VISREF(ICHA), ICHA
            IOK = IOK + 1
          ENDIF
          IF (IENCRA.EQ.1 .AND. TPRENC(ICHA).LT.TKELVN) THEN
            WRITE(NFECRA,1042) IENCRA, TKELVN, TPRENC(ICHA), ICHA
            IOK = IOK + 1
          ENDIF
        ENDDO
C
      ELSE
        IENCRA = 0
      ENDIF
C
      IF (IPHYLA.NE.2 .AND. IPPMOD(ICPL3C).GE.0) THEN
        WRITE(NFECRA,1043) IPHYLA, IPPMOD(ICPL3C)
        IOK = IOK + 1
      ENDIF
C
      IF ( IPHYLA.EQ.2 .AND. (IPPMOD(ICPL3C).LT.0 .AND.
     &                        IPPMOD(ICP3PL).LT.0) ) THEN
        WRITE(NFECRA,1044) IPHYLA, IPPMOD(ICPL3C), IPPMOD(ICP3PL)
        IOK = IOK + 1
      ENDIF
C
      IF (IOK.NE.0) CALL CSEXIT (1)
C                   ===========
C
C     NBPMAX
C
      IF (NBPMAX.LT.0) THEN
        WRITE(NFECRA,1050) NBPMAX
        IOK = IOK + 1
      ENDIF
C
C     NVLS
C
      IF (NVLS.LT.0 .OR. NVLS.GT.NUSVAR) THEN
        WRITE(NFECRA,1060) NUSVAR, NVLS
        IOK = IOK + 1
      ENDIF
C
      IF (IOK.NE.0) CALL CSEXIT (1)
C                   ===========
C
C     ISTTIO NSTITS LTSDYN LTSMAS LTSTHE
C
C     Si champs figes alors forcement en stationnaire
      IF (IILAGR.EQ.3) ISTTIO = 1
C
      IF (ISTTIO.LT.0 .OR. ISTTIO.GT.1) THEN
        WRITE(NFECRA,1061) ISTTIO
        IOK = IOK + 1
      ENDIF
C
      IF (IILAGR.EQ.2) THEN
        IF (ISTTIO.EQ.1 .AND. NSTITS.LT.1) THEN
          WRITE(NFECRA,1062) NSTITS
          IOK = IOK + 1
        ENDIF
        IF (LTSDYN.LT.0 .OR. LTSDYN.GT.1) THEN
          WRITE(NFECRA,1063) LTSDYN
          IOK = IOK + 1
        ENDIF
        IF (IPHYLA.EQ.1 .AND. (IMPVAR.EQ.1 .OR. IDPVAR.EQ.1)) THEN
          IF (LTSMAS.LT.0 .OR. LTSMAS.GT.1) THEN
            WRITE(NFECRA,1064) LTSMAS
            IOK = IOK + 1
          ENDIF
        ELSE
          LTSMAS = 0
        ENDIF
        IF ((IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) .OR. IPHYLA.EQ.2) THEN
          IF (LTSTHE.LT.0 .OR. LTSTHE.GT.1) THEN
            WRITE(NFECRA,1065) LTSTHE
            IOK = IOK + 1
          ENDIF
        ELSE
          LTSTHE = 0
        ENDIF
        IF (LTSDYN.EQ.1 .AND. ICCVFG.EQ.1) THEN
          WRITE(NFECRA,1066) LTSDYN, ICCVFG
          IOK = IOK + 1
        ENDIF
        IF (LTSDYN.NE.1 .AND. LTSTHE.NE.1 .AND. LTSMAS.NE.1) THEN
          WRITE(NFECRA,1067) IILAGR, LTSDYN, LTSTHE, LTSMAS
          IOK = IOK + 1
        ENDIF
      ELSE
        LTSDYN = 0
        LTSMAS = 0
        LTSTHE = 0
      ENDIF
C
      IF (IOK.NE.0) CALL CSEXIT (1)
C                   ===========
C
C     ISTALA SEUIL IDSTNT NSTIST NVLSTS
C
      IF (ISTALA.LT.0 .OR. ISTALA.GT.1) THEN
        WRITE(NFECRA,1070) ISTALA
        IOK = IOK + 1
      ENDIF
C
      IF (ISTALA.EQ.1) THEN
        IF (SEUIL.LT.0.D0) THEN
          WRITE(NFECRA,1071) SEUIL
          IOK = IOK + 1
        ENDIF
        IF (IDSTNT.LT.1) THEN
          WRITE(NFECRA,1072) IDSTNT
          IOK = IOK + 1
        ENDIF
        IF (ISTTIO.EQ.1) THEN
          IF (NSTIST.LT.IDSTNT) THEN
            WRITE(NFECRA,1073) IDSTNT, NSTIST
            IOK = IOK + 1
          ENDIF
        ENDIF
        IF (NVLSTS.LT.0 .OR. NVLSTS.GT.NUSSTA) THEN
          WRITE(NFECRA,1074) NUSSTA, NVLSTS
          IOK = IOK + 1
        ENDIF
        IF ( NBCLST .LT. 0 .OR. NBCLST.GT. NCLSTM ) THEN
          WRITE(NFECRA,1075) NCLSTM,NBCLST
          IOK = IOK + 1
        ENDIF
C
      ENDIF
C
      IF (IOK.NE.0) CALL CSEXIT (1)
C                   ===========
C
      IF (ISTALA.EQ.0) THEN
        ISUIST = 0
        SEUIL  = 0.D0
        IDSTNT = 0
        NSTIST = 0
      ENDIF
C
C     INJCON
C
      IF (INJCON.LT.0 .OR. INJCON.GT.1) THEN
        WRITE(NFECRA,1080) INJCON
        IOK = IOK + 1
      ENDIF
C
C     IROULE
C
      IF (IROULE.LT.0 .OR. IROULE.GT.2) THEN
        WRITE(NFECRA,1090) IROULE
        IOK = IOK + 1
      ENDIF
C
C     NORDRE
C
      IF (NORDRE.NE.1 .AND. NORDRE.NE.2) THEN
        WRITE(NFECRA,2000) NORDRE
        IOK = IOK + 1
      ENDIF
C
C     IDISTU
C
      IF (IDISTU.LT.0 .OR. IDISTU.GT.1) THEN
        WRITE(NFECRA,2010) IDISTU
        IOK = IOK + 1
      ENDIF
C
      II = ILPHAS
      IF (IDISTU.EQ.1 .AND. ITYTUR(II).NE.2 .AND. ITYTUR(II).NE.3
     &     .AND.  ITURB(II).NE.50 .AND. ITURB(II).NE.60 ) THEN
        WRITE(NFECRA,2011) IILAGR, IDISTU, II, ITURB(II)
        IOK = IOK + 1
      ELSE IF (IDISTU.EQ.0 .AND. ITURB(II).NE.0 .AND.
     &         ITYTUR(II).NE.2 .AND. ITYTUR(II).NE.3
     &     .AND.  ITURB(II).NE.50 .AND. ITURB(II).NE.60) THEN
        WRITE(NFECRA,2012) IILAGR, IDISTU, II, ITURB(II)
        IOK = IOK + 1
      ENDIF
C
C     IDISTU
C
      IF (IDIFFL.LT.0 .OR. IDIFFL.GT.1) THEN
        WRITE(NFECRA,2013) IDIFFL
        IOK = IOK + 1
      ENDIF
C
C     MODCPL IDIRLA
C
      IF (MODCPL.LT.0) THEN
        WRITE(NFECRA,2014) MODCPL
        IOK = IOK + 1
      ENDIF
      IF (MODCPL.GT.0) THEN
        IF (MODCPL.LT.IDSTNT) THEN
          WRITE(NFECRA,2015) MODCPL, IDSTNT
          IOK = IOK + 1
        ENDIF
        IF (ISTALA.EQ.0) THEN
          WRITE(NFECRA,2018) MODCPL, ISTALA
          IOK = IOK + 1
        ENDIF
        IF (IDIRLA.NE.1 .AND. IDIRLA.NE.2 .AND. IDIRLA.NE.3) THEN
          WRITE(NFECRA,2016) IDIRLA
          IOK = IOK + 1
        ENDIF
      ENDIF
C
C     ILAPOI
C
      IF (ILAPOI.LT.0 .OR. ILAPOI.GT.1) THEN
        WRITE(NFECRA,2017) ILAPOI
        IOK = IOK + 1
      ENDIF
C
      IF (IOK.NE.0) CALL CSEXIT (1)
C                   ===========
C
C     IENSI1 IENSI2
C
      IF (IENSI1.LT.0 .OR. IENSI1.GT.1) THEN
        WRITE(NFECRA,2030) IENSI1
        IOK = IOK + 1
      ENDIF
      IF (IENSI2.LT.0 .OR. IENSI2.GT.1) THEN
        WRITE(NFECRA,2031) IENSI2
        IOK = IOK + 1
      ENDIF
C
C     NBVIS NVISLA
C
      IF (IENSI1.EQ.1 .OR. IENSI2.EQ.1) THEN
        IF (NBVIS.GT.NBPMAX .OR. NBVIS.GT.NLISTE .OR.  NBVIS.LT.0) THEN
          WRITE(NFECRA,2032) NBPMAX, NLISTE, NBVIS
          IOK = IOK + 1
        ENDIF
        IF (NVISLA.LE.0 ) THEN
          WRITE(NFECRA,2033) NVISLA
          IOK = IOK + 1
        ENDIF
C
C     IVISV1 IVISV2 IVISTP IVISDM IVISTE
C
        IF (IVISV1.LT.0 .OR. IVISV1.GT.1) THEN
          WRITE(NFECRA,2040) IVISV1
          IOK = IOK + 1
        ENDIF
        IF (IVISV2.LT.0 .OR. IVISV2.GT.1) THEN
          WRITE(NFECRA,2041) IVISV2
          IOK = IOK + 1
        ENDIF
        IF (IVISTP.LT.0 .OR. IVISTP.GT.1) THEN
          WRITE(NFECRA,2042) IVISTP
          IOK = IOK + 1
        ENDIF
        IF (IVISDM.LT.0 .OR. IVISDM.GT.1) THEN
          WRITE(NFECRA,2043) IVISDM
          IOK = IOK + 1
        ENDIF
        IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) THEN
          IF (IVISTE.LT.0 .OR. IVISTE.GT.1) THEN
            WRITE(NFECRA,2044) IVISTE
            IOK = IOK + 1
           ENDIF
        ELSE
          IVISTE = 0
        ENDIF
C
C       IVISHP IVISDK IVISCH IVISCK
C
        IF (IPHYLA.EQ.2) THEN
          IF (IVISHP.LT.0 .OR. IVISHP.GT.1) THEN
            WRITE(NFECRA,2045) IVISHP
            IOK = IOK + 1
          ENDIF
          IF (IVISDK.LT.0 .OR. IVISDK.GT.1) THEN
            WRITE(NFECRA,2046) IVISDK
            IOK = IOK + 1
          ENDIF
          IF (IVISCH.LT.0 .OR. IVISCH.GT.1) THEN
            WRITE(NFECRA,2047) IVISCH
            IOK = IOK + 1
          ENDIF
          IF (IVISCK.LT.0 .OR. IVISCK.GT.1) THEN
            WRITE(NFECRA,2048) IVISCK
            IOK = IOK + 1
          ENDIF
        ELSE
          IVISHP = 0
          IVISDK = 0
          IVISCH = 0
          IVISCK = 0
        ENDIF
C
      ENDIF
C
C     IENSI3 NSTBOR
C
      IF (IENSI3.LT.0 .OR. IENSI3.GT.1) THEN
        WRITE(NFECRA,2050) IENSI3
        IOK = IOK + 1
      ENDIF
      IF (IENSI3.EQ.1 .AND. ISTTIO.EQ.1) THEN
        IF (NSTBOR.LT.1) THEN
          WRITE(NFECRA,2057) NSTBOR
          IOK = IOK + 1
        ENDIF
      ELSE
        NSTBOR = 1
      ENDIF
      IF (IENSI3.EQ.1) THEN
        IF (SEUILF.LT.0.D0) THEN
          WRITE(NFECRA,2058) SEUILF
          IOK = IOK + 1
        ENDIF
      ELSE
        SEUILF = 0.D0
      ENDIF
C
C     INBRBD IFLMBD IANGBD IVITBD IENCBD NUSBOR
C
      IF (IENSI3.EQ.1) THEN
C
        IF (IPHYLA.EQ.2 .AND. IENCRA.EQ.1) THEN
          IF (IENCBD.LT.0 .OR. IENCBD.GT.1) THEN
            WRITE(NFECRA,2051)
            IOK = IOK + 1
          ENDIF
        ELSE
          IENCBD = 0
        ENDIF
C
        IF (INBRBD.LT.0 .OR. INBRBD.GT.1) THEN
          WRITE(NFECRA,2052) INBRBD
          IOK = IOK + 1
        ENDIF
        IF (IFLMBD.LT.0 .OR. IFLMBD.GT.1) THEN
          WRITE(NFECRA,2053) IFLMBD
          IOK = IOK + 1
        ENDIF
        IF (IANGBD.LT.0 .OR. IANGBD.GT.1) THEN
          WRITE(NFECRA,2054) IANGBD
          IOK = IOK + 1
        ENDIF
        IF (IVITBD.LT.0 .OR. IVITBD.GT.1) THEN
          WRITE(NFECRA,2055) IVITBD
          IOK = IOK + 1
        ENDIF
        IF (NUSBOR.LT.0 .OR. NUSBOR.GT.NUSBRD) THEN
          WRITE(NFECRA,2056) NUSBRD, NUSBOR
          IOK = IOK + 1
        ENDIF
C
      IF (IOK.NE.0) CALL CSEXIT (1)
C                   ===========
C
        IRF = 0
C
        IF (INBRBD.EQ.1) THEN
          IRF = IRF + 1
          IF (IMOYBR(IRF).EQ.2) THEN
            WRITE(NFECRA,2060) INBRBD, IMOYBR(IRF)
          ENDIF
        ENDIF
        IF (IFLMBD.EQ.1) THEN
          IRF = IRF + 1
          IF (IMOYBR(IRF).EQ.2 .AND. INBRBD.EQ.0) THEN
            IOK = IOK + 1
            WRITE(NFECRA,2061) NOMBRD(IRF), 'IFLMBD',  IFLMBD,
     &                         IMOYBR(IRF), INBRBD
          ENDIF
        ENDIF
        IF (IANGBD.EQ.1) THEN
          IRF = IRF + 1
          IF (IMOYBR(IRF).EQ.2 .AND. INBRBD.EQ.0) THEN
            IOK = IOK + 1
            WRITE(NFECRA,2061) NOMBRD(IRF), 'IANGBD',  IANGBD,
     &                         IMOYBR(IRF), INBRBD
          ENDIF
        ENDIF
        IF (IVITBD.EQ.1) THEN
          IRF = IRF + 1
          IF (IMOYBR(IRF).EQ.2 .AND. INBRBD.EQ.0) THEN
            IOK = IOK + 1
            WRITE(NFECRA,2061) NOMBRD(IRF), 'IVITBD',  IVITBD,
     &                         IMOYBR(IRF), INBRBD
          ENDIF
        ENDIF
        IF (IPHYLA.EQ.2 .AND. IENCRA.EQ.1 .AND. IENCBD.EQ.1) THEN
          IRF = IRF + 1
          IF (IMOYBR(IRF).EQ.2 .AND. INBRBD.EQ.0) THEN
            IOK = IOK + 1
            WRITE(NFECRA,2061) NOMBRD(IRF), 'IENCBD',  IENCBD,
     &                         IMOYBR(IRF), INBRBD
          ENDIF
        ENDIF
        IF (NUSBOR.GT.0) THEN
          DO II = 1,NUSBOR
            IRF = IRF + 1
            IF (IMOYBR(IRF).EQ.2 .AND. INBRBD.EQ.0) THEN
              IOK = IOK + 1
              WRITE(NFECRA,2062) NOMBRD(IRF), II, IMOYBR(IRF), INBRBD
            ENDIF
          ENDDO
        ENDIF
C
        DO II = 1, IRF
          IF (IMOYBR(II).NE.0 .AND.
     &        IMOYBR(II).NE.1 .AND. IMOYBR(II).NE.2) THEN
          IOK = IOK + 1
            WRITE(NFECRA,2063) IMOYBR(IRF), NOMBRD(IRF)
        ENDIF
        ENDDO
C
      ENDIF
C
      IF (IOK.NE.0) CALL CSEXIT (1)
C                   ===========
C
C=======================================================================
C 3. INITIALISATIONS DES VARIABLES EN COMMON
C
C                            ATTENTION :
C                            ^^^^^^^^^^^
C
C    CES INITIALISATIONS NE DOIVENT ETRE MODIFIEES PAR L'UTILISATEUR
C
C=======================================================================
C
C
C ------------------------------------------------------------------
C 3.1 GENERALITES (D'AUTRES INITIALISATIONS SONT FAITES DANS LAGLEC)
C ------------------------------------------------------------------
C
C     NOMBRE DE PASSAGES ABSOLUS DANS LE MODULE LAGRANGIEN
C
      IPLAS = 0
C
C     NOMBRE DE PASSAGES RELATIFS DANS LE MODULE LAGRANGIEN
C
      IPLAR = 0
C
C     PAS DE TEMPS LAGRANGIEN (LAGUNE) : Par defaut le pas de temps
C                                 de reference de la phase continue
      DTP = DTREF
C
C     TEMPS COURANT PHYSIQUE LAGRANGIEN
C
      TTCLAG = 0.D0
C
C     INDICATEUR D'ERREUR (LAGCEL)
C
      IERR = 0
C
C    NBPART/DNBPAR : NOMBRE DE PARTICULES PRESENTES DANS LE DOMAINE
C                      DE CALCUL A CHAQUE ITERATION
C
      NBPART = 0
      DNBPAR = 0.D0
C
C     NBPERR/DNBPER : NOMBRE DE PARTICULES ELIMINES EN ERREUR
C
      NBPERR = 0
      DNBPER = 0.D0
C
C     NBPERT : NOMBRE DE PARTICULES ELIMINEES EN ERREUR DANS
C              LE CALCUL DEPUIS LE DEBUT SUITE COMPRISE
C
      NBPERT = 0
C
C     NBPTOT : NOMBRE DE PARTICULES DU CALCUL (SUITES COMPRISES)
C
      NBPTOT = 0
C
C     NBPOUT/DNBPOU : Contient les particules sorties de facon normale,
C                       plus les particules sorties en erreur de reperage.
C
      NBPOUT = 0
      DNBPOU = 0.D0
C
C     NDEPOT : Nombre de particules deposees definitivement
C               dont on garde une trace en memoire pour le
C               post-processing en mode deplacement.
C
      NDEPOT = 0
C
C     NPCLON/DNPCLO : NOMBRE DE NOUVELLES PARTICULES PAR CLONNAGE
C     NPKILL/DNPCSU : NOMBRE DE PARTICULES VICTIMES DE LA ROULETTE RUSSE
C     NPCSUP/DNPKIL : NOMBRE DE PARTICULES QUI ON SUBI LE CLONNAGE
C
      NPCLON = 0
      DNPCLO = 0.D0
C
      NPCSUP = 0
      DNPCSU = 0.D0
C
      NPKILL = 0
      DNPKIL = 0.D0
C
C     NPENCR/DNPENC : nombre de grains de charbon "encrasses"
C
      NPENCR = 0
      DNPENC = 0.D0
C
C
C     CONDITIONS AUX LIMITES
C
      DO II = 1,NFLAGM
        ILFLAG(II) = 0
        IUSNCL(II) = 0
        IUSCLB(II) = 0
        IUSMOY(II) = 0
        DEBLAG(II) = 0
      ENDDO
      DO I1 = 1,NDLAIM
        DO I2 = 1,NFLAGM
          DO I3 = 1,NCLAGM
            IUSLAG(I3,I2,I1) = 0
          ENDDO
        ENDDO
      ENDDO
      DO I1 = 1,NDLAGM
        DO I2 = 1,NFLAGM
          DO I3 = 1,NCLAGM
            RUSLAG(I3,I2,I1) = 0
          ENDDO
        ENDDO
      ENDDO
C
C
C     STATISTIQUES VOLUMIQUES
C
C     Nombre de pas de temps DEPUIS LE DEBUT DU CALCUL STATIONNAIRES
C     des stats
C
      NPST =  0
C
C     Nombre de pas de temps total des stats depuis le debut
C     du calcul, partie instationnaire comprise
C
      NPSTT = 0
C
C     Temps physique des stats
C
      TSTAT = 0.D0
C
C     STATISTIQUES AUX FRONTIERES
C
C     Nombre de pas de temps DEPUIS LE DEBUT DU CALCUL STATIONNAIRES
C     des stats aux frontieres
C
      NPSTF =  0
C
C     Nombre de pas de temps total des stats aux frontieres
C     depuis le debut du calcul, partie instationnaire comprise
C
      NPSTFT = 0
C
C     Temps physique des stats aux frontieres
C
      TSTATP = 0.D0
C
C     COUPLAGE RETOUR
C
C     Nombre de pas de temps DEPUIS LE DEBUT DU CALCUL STATIONNAIRES
C     des termes sources pour le couplage retour
C
      NPTS =  0
C
C     Initialisation du sous-pas
C
      NOR = 0
C
C     NOMBRE D'ENREGISTREMENT POUR LE POST DEPLACEMENT (ENSWAF)
C
      ITLAG = 0
C
C     TEMPS PHYSIQUE LAGRANGIEN POUR LE POST (ENSWAF)
C
      DO II = 1,9999
        TIMLAG(II) = 0.D0
      ENDDO
C
C     FINALISATION DE LA LISTE DE PARTICULES A VISUALISER
C
      DO II = NBVIS+1, NLISTE
        LISTE (II) = -1
      ENDDO
C
C    les trous, les repetitions dans le tableau LISTE seront
C    suprimes, et les numeros seront ranges par ordre croissant.
C
      CALL LAGTRI
C     ===========
C
C ------------------------------------------------
C 3.2 DIMENSIONS DES TABLEAUX LIEES AUX PARTICULES
C ------------------------------------------------
C
C-->  NOMBRE MINIMAL DE VARIABLES LIEES AUX PARTICULES
C
C     NVP   : Variables sur les particules avec equation (ETTP et ETTPA)
C     NVEP  : Variables d'etat (reels)   sur les particules (TEPA)
C     NIVEP : Variables d'etat (entiers) sur les particules (ITEPA)
C
      NVP   = 11
      NVEP  = 2
      NIVEP = 1
C
      IF (NBCLST.GT.0) THEN
C
C       --> Statistique par classe :
C          1 VARIABLE D'ETAT ENTIERE SUPPLEMENTAIRE : Numero de la
C          classe statique a laquelle appartient la particule
C
        NIVEP = NIVEP +1
C
      ENDIF
C
C
      IF (IPHYLA.EQ.1) THEN
C
C       --> EQUATION SUR LA TEMPERATURE :
C           3 VARIABLES SUPPLEMENTAIRES Tp, Tf, Cp dans ETTP et ETTPA
C
        IF (ITPVAR.EQ.1) NVP = NVP + 3
C
C       --> EQUATION SUR LA TEMPERATURE ET RAYONNEMENT : EMISSIVITE
C
        IF (ITPVAR.EQ.1 .AND. IIRAYO.EQ.1) NVEP = NVEP + 1
C
      ELSE IF (IPHYLA.EQ.2) THEN
C
C    --> CHARBON :
C
C     ETTP et ETTPA :
C     -------------
C       5 VARIABLES SUPPLEMENTAIRES Tp, Tf, Mch, Mck, Cp
        NVP = NVP + 5
C
C     TEPA :
C     ----
C       3 VARIABLES D'ETATS REELLES SUPPLEMENTAIRES : Dck, D0P, R0P
C
        NVEP = NVEP + 3
C
C     ITEPA :
C     -----
C       1 VARIABLE D'ETAT ENTIERE SUPPLEMENTAIRE : Numero du charbon
        NIVEP = NIVEP + 1
C
      ENDIF
C
C-->  VARIABLES UTILISATEURS SUPPLEMENTAIRES : NVLS
C
      IF (NVLS.GT.0) NVP = NVP + NVLS
C
C-->  NVP1 represente le nombre de variables sur les particules en
C     enlevant position, vitesse particule et vitesse fluides (TSVAR)
C
      NVP1 = NVP - 9
C
C -------------------------------------------------------------------
C 3.3 DEFINITION DES POINTEURS SUR LES VARIABLES LIEES AUX PARTICULES
C -------------------------------------------------------------------
C
C   3.3.1 TABLEAU ETTP
C   ~~~~~~~~~~~~~~~~~~
C
C   Attention il faut que
C   JMP, JDP, JXP, JYP, JZP, JUP, JVP, JWP, JUF, JVF, JWF
C   soient les derniers pointeurs pour les tableaux ETTP et ETTPA
C   a cause du remplissage et de la dimension du tableau
C   TSVAR(NBPMAX,NVP1) et de son mode de lecture/ecriture
C   (cf. LAGITG).
C
C    JXP,JYP,JZP  : COORDONNES DE LA POSITION DE LA PARTICULE NPT
C    JUP,JVP,JWP  : COMPOSANTES DE LA VITESSE ABSOLUE
C    JUF,JVF,JWF  : COMPOSANTES DE LA VITESSE DU FLUIDE VU
C
C    JMP,JDP      : MASSE, DIAMETRE
C    JTP,JTF,JCP  : TEMPERATURE PARTICULE ET FLUIDE ET CHALEUR SPECIFIQUE
C    JVLS(NUSVAR) : VARIABLE SUPPLEMENTAIRES
C
C   Charbon
C   -------
C    JHP  : Temperature en degres Celsius grain de charbon
C    JMCH : MASSE DE CHARBON REACTIF
C    JMCK : MASSE DE COKE
C
      JTP  = 0
      JTF  = 0
      JCP  = 0
      JHP  = 0
      JMCH = 0
      JMCK = 0
      DO II = 1,NUSVAR
        JVLS(II) = 0
      ENDDO
C
      IRF  = 0
C
      IF (IPHYLA.EQ.1) THEN
C
        IF (ITPVAR.EQ.1) THEN
          JTP = IRF + 1
          JTF = JTP + 1
          JCP = JTF + 1
          IRF = JCP
        ENDIF
C
      ELSE IF (IPHYLA.EQ.2) THEN
C
        JHP  = IRF  + 1
        JTF  = JHP  + 1
        JMCH = JTF  + 1
        JMCK = JMCH + 1
        JCP  = JMCK + 1
        IRF  = JCP
C
      ENDIF
C
      IF (NVLS.GT.0) THEN
        DO II = 1,NVLS
          IRF = IRF + 1
          JVLS(II) = IRF
        ENDDO
        IRF = JVLS(NVLS)
      ENDIF
C
      JMP = IRF + 1
      JDP = JMP + 1
      JXP = JDP + 1
      JYP = JXP + 1
      JZP = JYP + 1
      JUP = JZP + 1
      JVP = JUP + 1
      JWP = JVP + 1
      JUF = JWP + 1
      JVF = JUF + 1
      JWF = JVF + 1
      IRF = JWF
C
      IF (IRF.GT.NVP) THEN
        WRITE(NFECRA,3004) IRF, NVP
        CALL CSEXIT(1)
      ENDIF
C
C   3.3.2 TABLEAU TEPA
C   ~~~~~~~~~~~~~~~~~~
C
C     JRTSP       : TEMPS DE SEJOUR DES PARTICULES
C     JRPOI       : POIDS DES PARTICULES
C     JREPS       : EMISSIVITE DES PARTICULES
C
C   Charbon
C   -------
C     JRDCK       : DIAMETRE DU COEUR RETRECISSANT
C     JRD0P       : DIAMETRE INITIAL DES PARTICULES
C     JRR0P       : MASSE VOLUMIQUE INITIALE DES PARTICULES
C
      JREPS = 0
      JRDCK = 0
      JRR0P = 0
      JRR0P = 0
C
      JRTSP = 1
      JRPOI = 2
      IRF   = JRPOI
C
      IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1 .AND. IIRAYO.EQ.1) THEN
        JREPS = IRF   + 1
        IRF   = JREPS
      ENDIF
C
      IF (IPHYLA.EQ.2) THEN
        JRDCK = IRF   + 1
        JRD0P = JRDCK + 1
        JRR0P = JRD0P + 1
        IRF   = JRD0P
      ENDIF
C
      IF (IRF.GT.NVEP) THEN
        WRITE(NFECRA,3005) IRF, NVEP
        CALL CSEXIT(1)
      ENDIF
C
C   3.3.3 TABLEAU ITEPA
C   ~~~~~~~~~~~~~~~~~~~
C
C     JISOR       : MAILLE D'ARRIVEE
C
C   Statistique par classe
C   ----------------------
C
C     JCLST       : classe (statique) a laquelle la particule appartient
C
C   Charbon
C   -------
C     JINCH       : NUMERO DU CHARBON DE LA PARTICULE
C
      JINCH = 0
C
      JISOR = 1
      IRF   = JISOR
C
      IF (NBCLST .GT. 0) THEN
        JCLST = IRF + 1
        IRF = JCLST
      ENDIF
C
      IF (IPHYLA.EQ.2) THEN
        JINCH = IRF + 1
        IRF   = JINCH
      ENDIF
C
      IF (IRF.GT.NIVEP) THEN
        WRITE(NFECRA,3006) IRF, NIVEP
        CALL CSEXIT(1)
      ENDIF
C
C --------------------------------------------------------------
C 3.4 DEFINITION DES POINTEURS LIES A L'INJECTION DES PARTICULES
C --------------------------------------------------------------
C
C   3.4.1 TABLEAU RUSLAG (DONNEES D'ENTREE)
C   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C     IUNO  : Norme de la vitesse
C     IUPT  : U par classe et zones
C     IVPT  : V par classe et zones
C     IWPT  : W par classe et zones
C     IDEBT : Debit
C     IPOIT : Poids de la particule
C     IDPT  : Diametre
C     IVDPT : Ecart-type du diametre
C     ITPT  : Temperature
C     ICPT  : Cp
C     IEPSI : Emissivite
C     IROPT : Masse volumique
C     IHPT  : Temperature
C     IMCHT : Masse de charbon reactif
C     IMCKT : Masse de coke
C     IDCKT : Diametre du coeur retrecissant
C
      IUNO   = 1
      IUPT   = IUNO  + 1
      IVPT   = IUPT  + 1
      IWPT   = IVPT  + 1
      ITPT   = IWPT  + 1
      IDPT   = ITPT  + 1
      IVDPT  = IDPT  + 1
      IROPT  = IVDPT + 1
      ICPT   = IROPT + 1
      IEPSI  = ICPT  + 1
      IPOIT  = IEPSI + 1
      IDEBT  = IPOIT + 1
      IRF    = IDEBT
C
C Specifique Charbon
C
      IHPT   = IRF   + 1
      IMCHT  = IHPT  + 1
      IMCKT  = IMCHT + 1
      IRF    = IMCKT
C
      IF (IRF.GT.NDLAGM) THEN
        WRITE(NFECRA,3001) IRF, NDLAGM
        CALL CSEXIT(1)
      ENDIF
C
C   3.4.2 TABLEAU IUSLAG (DONNEES D'ENTREE)
C   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C     IJNBP  : nbr de part par classe et zones frontieres
C     IJFRE  : frequence d'injection
C     ICLST : numero de groupe auquel appartient la particule
C             (uniquement si on souhaite des statistiques par groupe)
C     IJUVW  : type de condition vitesse
C     IJPRTP : profil de temperature pour les particules
C     IJPRPD : profil de poids statistiques
C     IJPRDP : profil de diametre pour les particules
C     INUCHL : numero du charbon de la particule (si IPHYLA=2)
C
      IJNBP  = 1
      IJFRE  = IJNBP  + 1
      ICLST  = IJFRE  + 1
      IJUVW  = ICLST  + 1
      IJPRTP = IJUVW  + 1
      IJPRDP = IJPRTP + 1
      IJPRPD = IJPRDP + 1
      IRF    = IJPRPD
C
C Specifique Charbon
C
      INUCHL = IRF + 1
      IRF    = INUCHL
C
      IF (IRF.GT.NDLAIM) THEN
        WRITE(NFECRA,3002) IRF, NDLAIM
        CALL CSEXIT(1)
      ENDIF
C
C -----------------------------------------------
C 3.5 DIEMSIONS DU TABLEAU LIEES AUX STATISTIQUES
C -----------------------------------------------
C
C     NVLSTA   : NOMBRE DE VARIABLES LIEES AUX STATISTIQUES
C                DIMENSION DU TABLEAU STATIS
C     NVLSTS   : VARIABLES LIEES AUX STATISTIQUES SUPPLEMENTAIRES
C
      IF (ISTALA.EQ.0) THEN
C
        NVLSTA = 0
C
      ELSE IF(ISTALA.EQ.1) THEN
C
        NVLSTA = 6 + NVLSTS
C
        IF (IPHYLA .EQ. 1) THEN
C
          IF (IDPVAR.EQ.1) NVLSTA = NVLSTA + 1
C
          IF (IMPVAR.EQ.1) NVLSTA = NVLSTA + 1
C
          IF (ITPVAR.EQ.1) NVLSTA = NVLSTA + 1
C
        ELSE IF (IPHYLA .EQ. 2) THEN
C
          NVLSTA = NVLSTA + 4
C
          IF (IENCRA.EQ.1) THEN
            NVLSTA = NVLSTA + 1
          ENDIF
C
        ENDIF
C
      ENDIF
C
      IF(NVLSTA.GT.NVPLMX) THEN
        WRITE(NFECRA,3003) NVLSTA, NVPLMX
        CALL CSEXIT(1)
      ENDIF
C
C --------------------------------------------------
C 3.6 DEFINITION DES POINTEURS LIES AUX STATISTIQUES
C --------------------------------------------------
C
C     ILVX,ILVY,ILVZ    : Vitesse
C     ILFV              : Concentrations volumiques
C     ILPD              : Somme des poids statistiques
C
C     ILTP              : Temperature
C     ILDP              : Diametre
C     ILMP              : Masse
C
C     ILHP              : Temperature
C     ILMCH             : Masse de charbon reactif
C     ILMCK             : Masse de coke
C     ILMDK             : Diametre du coeur retrecissant
C
C     ILVU(NUSSTA)      : Statistiques utilisateur
C
      IP  = 0
      ILVX  = IP    + 1
      ILVY  = ILVX  + 1
      ILVZ  = ILVY  + 1
      ILFV  = ILVZ  + 1
      ILTS  = ILFV  + 1
      IRF   = ILTS
C
      ILTP  = 0
      ILDP  = 0
      ILMP  = 0
      ILHP  = 0
      ILMCH = 0
      ILMCK = 0
      ILDCK = 0
C
      DO II = 1,NUSSTA
        ILVU(II) = 0
      ENDDO
C
      IF (IPHYLA.EQ.1) THEN
C
        IF (ITPVAR.EQ.1) THEN
          ILTP = IRF + 1
          IRF  = ILTP
        ENDIF
C
        IF (IDPVAR.EQ.1) THEN
          ILDP = IRF + 1
          IRF  = ILDP
        ENDIF
C
        IF (IMPVAR.EQ.1) THEN
          ILMP = IRF + 1
          IRF  = ILMP
        ENDIF
C
      ELSE IF (IPHYLA.EQ.2) THEN
C
        ILHP  = IRF   + 1
        ILMCH = ILHP  + 1
        ILMCK = ILMCH + 1
        ILDCK = ILMCK + 1
        IRF   = ILDCK
C
      ENDIF
C
      IF (NVLSTS.GT.0) THEN
        DO II = 1,NVLSTS
          ILVU(II) = IRF + II
        ENDDO
        IRF = IRF + NVLSTS
      ENDIF
C
      ILPD  = IRF  + 1
C
C -----------------------------------------------------------------
C 3.7 DEFINITION DES POINTEURS LIES AUX STATISTIQUES AUX FRONTIERES
C -----------------------------------------------------------------
C
C     INBRBD : NOMBRE D'INTERACTIONS PARTICULES/FRONTIERES
C     IFLMBD : FLUX DE MASSE PARTICULAIRE
C     IANGBD : ANGLE VITESSE
C     IVITBD : VITESSE DE LA PARTICULE
C     IENCBD : MASSE DE GRAINS DE CHARBON ENCRASSES
C     NUSBOR : INFORMATIONS UTILISATEUR SUPPLEMENTAIRES
C     NVISBR : NOMBRE TOTAL D'INTERACTIONS A ENREGISTRER
C
      IF (IENSI3.EQ.1) THEN
C
        IRF = 0
C
        IF (INBRBD.EQ.1) THEN
          IRF = IRF + 1
          INBR = IRF
        ENDIF
C
        IF (IFLMBD.EQ.1) THEN
          IRF = IRF + 1
          IFLM = IRF
        ENDIF
C
        IF (IANGBD.EQ.1) THEN
          IRF = IRF + 1
          IANG = IRF
        ENDIF
C
        IF (IVITBD.EQ.1) THEN
          IRF = IRF + 1
          IVIT = IRF
        ENDIF
C
        IF (IPHYLA.EQ.2 .AND. IENCRA.EQ.1 .AND. IENCBD.EQ.1) THEN
          IRF = IRF + 1
          IENC = IRF
        ENDIF
C
        IF (NUSBOR.GT.0) THEN
          DO II = 1,NUSBOR
            IRF = IRF + 1
            IUSB(II) = IRF
          ENDDO
        ENDIF
C
        NVISBR = IRF
C
      ELSE
C
        NVISBR = 0
C
      ENDIF
C
C ---------------------------------------------------------------
C 3.8 DEFINITION DES POINTEURS LIES AUX TERMES SOURCES LAGRANGIEN
C     POUR COUPLAGE RETOUR
C ---------------------------------------------------------------
C
C     Couplage-retour uniquement vers la phase continue numero ILPHAS
C
      IPHAS  = ILPHAS
C
C     Nombre de termes sources de couplage-retour
C
      NTERSL = 0
C
      IRF = 0
      ITSVX  = 0
      ITSVY  = 0
      ITSVZ  = 0
      ITSLI  = 0
      ITSKE  = 0
      ITSR11 = 0
      ITSR12 = 0
      ITSR13 = 0
      ITSR22 = 0
      ITSR23 = 0
      ITSR33 = 0
      ITSMAS = 0
      ITSTE  = 0
      ITSTI  = 0
      DO ICHA = 1,NCHARM2
        ITSMV1(ICHA) = 0
        ITSMV2(ICHA) = 0
      ENDDO
      ITSCO  = 0
      ITSFP4 = 0
C
C Dynamique : Vitesse + Turbulence
C
      IF (LTSDYN.EQ.1) THEN
C
        NTERSL = NTERSL + 4
C
        ITSVX  = IRF   + 1
        ITSVY  = ITSVX + 1
        ITSVZ  = ITSVY + 1
        ITSLI  = ITSVZ + 1
        IRF    = ITSLI
C
        IF (ITYTUR(IPHAS).EQ.2 .OR. ITURB(IPHAS).EQ.50
     &       .OR. ITURB(IPHAS).EQ.60) THEN
C K-eps, v2f et k-omega
          NTERSL = NTERSL + 1
C
          ITSKE  = IRF    + 1
          IRF    = ITSKE
C
        ELSE IF (ITYTUR(IPHAS).EQ.3) THEN
C RIJ
          NTERSL = NTERSL + 6
C
          ITSR11 = IRF    + 1
          ITSR12 = ITSR11 + 1
          ITSR13 = ITSR12 + 1
          ITSR22 = ITSR13 + 1
          ITSR23 = ITSR22 + 1
          ITSR33 = ITSR23 + 1
          IRF    = ITSR33
        ELSE
          WRITE(NFECRA,3010) IILAGR, LTSDYN, IPHAS, ITURB(IPHAS)
          CALL CSEXIT (1)
C         ===========
        ENDIF
C
      ENDIF
C
C Masse
C
      IF (LTSMAS.EQ.1) THEN
C
        NTERSL = NTERSL + 1
C
        ITSMAS = IRF    + 1
        IRF    = ITSMAS
C
      ENDIF
C
C Thermique
C
      IF (LTSTHE.EQ.1) THEN
C
        IF (IPHYLA.EQ.1) THEN
C
C     Temperature
C
          IF (ITPVAR.EQ.1) THEN
C
            NTERSL = NTERSL + 2
C
            ITSTE  = IRF    + 1
            ITSTI  = ITSTE  + 1
            IRF    = ITSTI
C
          ENDIF
C
C     Charbon
C
        ELSE IF (IPHYLA.EQ.2) THEN
C
          NTERSL = NTERSL + 4 + 2*NCHARB
C
          ITSTE  = IRF    + 1
          ITSTI  = ITSTE  + 1
C
          DO ICHA = 1,NCHARB
            ITSMV1(ICHA) = ITSTI + ICHA
          ENDDO
C
          DO ICHA = 1,NCHARB
            ITSMV2(ICHA) = ITSMV1(NCHARB) + ICHA
          ENDDO

          ITSCO  = ITSMV2(NCHARB) + 1
          ITSFP4 = ITSCO          + 1
          IRF    = ITSFP4
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C
C--------
C FORMATS
C--------
C
 1010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE MODULE LAGRANGIEN A UNE VALEUR     ',/,
     &'@       NON PERMISE (LAGOPT).                                ',/,
     &'@                                                            ',/,
     &'@    IILAGR DEVRAIT ETRE UN ENTIER EGAL A 0, 1, 2 OU 3.      ',/,
     &'@       IL VAUT ICI IILAGR = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LES CONNECTIVITES FACES -> SOMMETS OU LES COORDONNEES   ',/,
     &'@    DES NOEUDS SOMMETS NE SONT PAS FOURNIS     (LAGOPT).    ',/,
     &'@                                                            ',/,
     &'@  Le module lagrangien est active mais les noeuds sommets   ',/,
     &'@   ou les connectivites faces -> sommets n''ont pas ete     ',/,
     &'@   fournis au noyau.                                        ',/,
     &'@                                                            ',/,
     &'@    NNDOD  = ', I10                                          ,/,
     &'@    LNDFAC = ', I10                                          ,/,
     &'@    LNDFBR = ', I10                                          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier si le fichier geomet contient bien les           ',/,
     &'@    connectivites faces -> sommets si IFOENV = 0.           ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1012 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE MODULE LAGRANGIEN IILAGR ET        ',/,
     &'@      L''INDICATEUR DE SUITE ONT DES VALEURS INCOMPATIBLES  ',/,
     &'@      (LAGOPT).                                             ',/,
     &'@                                                            ',/,
     &'@       IILAGR = ', I10                                       ,/,
     &'@       ISUITE = ', I10                                       ,/,
     &'@                                                            ',/,
     &'@  Le module lagrangien est active en mode champs figes,     ',/,
     &'@   alors que le calcul de la phase continue n''est pas      ',/,
     &'@   une suite.                                               ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR dans la subroutine USLAG1 et ',/,
     &'@  verifier la valeur de ISUITE dans la subroutine USINI1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1013 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LA PHYSIQUE PARTICULIERE COMBUTION CHARBON PULVERISE    ',/,
     &'@      COUPLE AU TRANSPORT LAGRANGIEN DES PARTICULES         ',/,
     &'@      DE CHARBON EST ACTIVEE (USPPMO), ALORS QUE LE COUPLAGE',/,
     &'@      RETOUR DE LA PHASE DISPERSEE SUR LE PHASE CONTINUE    ',/,
     &'@      N''EST PAS ENCLENCHE (LAGOPT).                        ',/,
     &'@                                                            ',/,
     &'@       IILAGR = ', I10                                       ,/,
     &'@       IPPMOD(ICPL3C) = ', I10                               ,/,
     &'@                                                            ',/,
     &'@  Le module lagrangien doit etre active en mode couplage    ',/,
     &'@   retour pour etre couple avec la combustion d''une        ',/,
     &'@   flamme de charbon pulverise en phase continue.           ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR dans la subroutine USLAG1 et ',/,
     &'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
 1014 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE MODULE LAGRANGIEN IILAGR ET        ',/,
     &'@      LE CHOIX DU TYPE DE PAS DE TEMPS IDTVAR               ',/,
     &'@      ONT DES VALEURS INCOMPATIBLES (LAGOPT).               ',/,
     &'@                                                            ',/,
     &'@       IILAGR = ', I10                                       ,/,
     &'@       IDTVAR = ', I10                                       ,/,
     &'@                                                            ',/,
     &'@  Le module lagrangien ne peut pas etre active avec un pas  ',/,
     &'@   de temps variable en temps et en espace. Seuls les pas   ',/,
     &'@   de temps uniforme et constant, et variable en temps et   ',/,
     &'@   uniforme en espace sont possibles.                       ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR dans la subroutine USLAG1 et ',/,
     &'@  verifier la valeur de IDTVAR dans la subroutine USINI1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1020 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE SUITE DU MODULE LAGRANGIEN A UNE       ',/,
     &'@       VALEUR NON PERMISE (LAGOPT).                         ',/,
     &'@                                                            ',/,
     &'@    ISUILA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI ISUILA = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1021 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ALERTE A L''EXECUTION DU MODULE LAGRANGIEN  ',/,
     &'@    *********                                               ',/,
     &'@                                                  (LAGOPT). ',/,
     &'@                                                            ',/,
     &'@  Le module lagrangien est active en suite de calcul,       ',/,
     &'@   alors que le calcul de la phase continue n''est pas      ',/,
     &'@   une suite.                                               ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de ISUILA dans la subroutine USLAG1.   ',/,
     &'@  Verifier la valeur de ISUITE dans la subroutine USINI1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1022 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE SUITE DE CALCUL SUR LES STATISTIQUES   ',/,
     &'@       VOLUMIQUE ET AUX FRONTIERES, AINSI QUE SUR LES       ',/,
     &'@       TERMES SOURCES DE COUPLAGES RETOUR                   ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    ISUIST DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI ISUIST = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de ISUIST dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DES MODELES PHYSIQUES LIES AUX PARTICULES ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IPHYLA DEVRAIT ETRE UN ENTIER EGAL A 0 1 OU 2           ',/,
     &'@       IL VAUT ICI IPHYLA = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1031 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''EQUATION DU DIAMETRE DES           ',/,
     &'@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).        ',/,
     &'@                                                            ',/,
     &'@     IDPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1            ',/,
     &'@       IL VAUT ICI IDPVAR = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IDPVAR dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1032 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''EQUATION DE LA TEMPERATURE DES     ',/,
     &'@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).        ',/,
     &'@                                                            ',/,
     &'@     ITPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1            ',/,
     &'@       IL VAUT ICI ITPVAR = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de ITPVAR dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1033 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''EQUATION DE LA MASSE DES           ',/,
     &'@       PARTICULES A UNE VALEUR NON PERMISE (LAGOPT).        ',/,
     &'@                                                            ',/,
     &'@     IMPVAR DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1            ',/,
     &'@       IL VAUT ICI IMPVAR = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IMPVAR dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1034 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''EQUATION DE LA TEMPERATURE DES     ',/,
     &'@       PARTICULES EST ACTIVE (ITPVAR = ',I10,'),'            ,/,
     &'@       ALORS QU''AUCUN SCALAIRE THERMIQUE N''EST DISPONIBLE ',/,
     &'@       POUR LA PHASE CONTINUE IPHAS = ',I10,' (LAGOPT).     ',/,
     &'@                                                            ',/,
     &'@     ISCALT DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL 1      ',/,
     &'@       IL VAUT ICI ISCALT(',I10   ,') = ',I10                ,/,
     &'@                                                            ',/,
     &'@  La valeur de ISCALT est renseignee automatiquement        ',/,
     &'@    si une physique particuliere est activee dans USPPMO.   ',/,
     &'@                                                            ',/,
     &'@  L''utilisateur doit renseigner ISCALT dans USINI1 si      ',/,
     &'@    aucune physique particuliere n''est activee dans USPPMO.',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de ITPVAR dans la subroutine USLAG1,   ',/,
     &'@  verifier la valeur de ISCALT dans la subroutine USINI1 et ',/,
     &'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1036 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LA CHALEUR MASSIQUE D''INITIALISATION DES PARTICULES    ',/,
     &'@       DEJA PRESENTE DANS LE DOMAINE DE CALCUL              ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@     CPPART DEVRAIT ETRE UN REEL STRICTEMENT POSITIF        ',/,
     &'@       IL VAUT ICI CPPART = ', E14.5                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de CPPART dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1037 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LA TEMPERATURE D''INITIALISATION DES PARTICULES         ',/,
     &'@       DEJA PRESENTE DANS LE DOMAINE DE CALCUL              ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@     TPART DEVRAIT ETRE UN REEL SUPERIEUR A ',E14.5          ,/,
     &'@       (EN DEGRES CELSIUS)                                  ',/,
     &'@       IL VAUT ICI TPART = ', E14.5                          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de TPART dans la subroutine USLAG1.    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''ENCRASSEMENT DES PARTICULES        ',/,
     &'@       DE CHARBON A UNE VALEUR NON PERMISE (LAGOPT).        ',/,
     &'@                                                            ',/,
     &'@     IENCRA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1            ',/,
     &'@       IL VAUT ICI IENCRA = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IENCRA dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1041 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''ENCRASSEMENT DES PARTICULES        ',/,
     &'@       DE CHARBON EST ACTIVE (IENCRA = ',I10,')'             ,/,
     &'@       AVEC UNE VALEUR DE VISCOSITE CRITIQUE                ',/,
     &'@       NON PERMISE (LAGOPT).                                ',/,
     &'@                                                            ',/,
     &'@     VISREF DEVRAIT ETRE UN REEL STRICTEMENT POSITIF (Pa.s) ',/,
     &'@       IL VAUT ICI VISREF = ', E14.5                         ,/,
     &'@       POUR LE CHARBON :' , I10                              ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de VISREF dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1042 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''ENCRASSEMENT DES PARTICULES        ',/,
     &'@       DE CHARBON EST ACTIVE (IENCRA = ',I10,')'             ,/,
     &'@       AVEC UNE VALEUR DE TEMPERATURE SEUIL                 ',/,
     &'@       NON PERMISE (LAGOPT).                                ',/,
     &'@                                                            ',/,
     &'@     TPRENC DEVRAIT ETRE UN REEL SUPERIEUR A ',E14.5         ,/,
     &'@       (EN DEGRES CELSIUS)                                  ',/,
     &'@       IL VAUT ICI TPRENC = ', E14.5                         ,/,
     &'@       POUR LE CHARBON :' , I10                              ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de TPRENC dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1043 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LA PHYSIQUE PARTICULIERE COMBUTION CHARBON PULVERISE    ',/,
     &'@      COUPLE AU TRANSPORT LAGRANGIEN DES PARTICULES         ',/,
     &'@      DE CHARBON EST ACTIVEE (USPPMO), ALORS QUE L''OPTION  ',/,
     &'@      TRANSPORT DE PARTICULE DE CHARBON                     ',/,
     &'@      N''EST PAS ENCLENCHEE (LAGOPT).                       ',/,
     &'@                                                            ',/,
     &'@       IPHYLA = ', I10                                       ,/,
     &'@       IPPMOD(ICPL3C) = ', I10                               ,/,
     &'@                                                            ',/,
     &'@  Le module lagrangien doit etre active en mode transport   ',/,
     &'@   de particules de charbon pour etre couple avec la        ',/,
     &'@   combustion d''une flamme de charbon pulverise en phase   ',/,
     &'@   continue.                                                ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1 et ',/,
     &'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1044 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE TRANSPORT LAGRANGIEN DE PARTICULES DE CHARBON        ',/,
     &'@      EST ACTIVE (LAGOPT), ALORS QU''AUCUNE PHYSIQUE        ',/,
     &'@      PARTICULIERE SUR LA COMBUSTION DU CHABON PULVERISE    ',/,
     &'@      N''EST PAS ENCLENCHE (USPPMO).                        ',/,
     &'@                                                            ',/,
     &'@       IPHYLA = ', I10                                       ,/,
     &'@       IPPMOD(ICPL3C) = ', I10                               ,/,
     &'@       IPPMOD(ICP3PL) = ', I10                               ,/,
     &'@                                                            ',/,
     &'@  Le transport lagrangien de particule de charbon doit      ',/,
     &'@   etre couple avec la combustion d''une flamme de charbon  ',/,
     &'@   pulverise en phase continue.                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1 et ',/,
     &'@  verifier la valeur de IPPMOD dans la subroutine USPPMO.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1050 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE NOMBRE DE PARTIUCULES A TRAITER    ',/,
     &'@      PAR LE MODULE LAGRANGIEN A UNE VALEUR                 ',/,
     &'@       NON PERMISE (LAGOPT).                                ',/,
     &'@                                                            ',/,
     &'@    NBPMAX DEVRAIT ETRE UN ENTIER STRICTEMENT POSITIF       ',/,
     &'@       IL VAUT ICI NBPMAX = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NBPMAX dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1060 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE NOMBRE DE VARIABLES                ',/,
     &'@       SUPPLEMENTAIRES LIEES AUX PARTICULES                 ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@     NVLS DEVRAIT ETRE UN ENTIER ENTRE 0 ET ',I10            ,/,
     &'@       IL VAUT ICI NVLS = ', I10                             ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NVLS dans la subroutine USLAG1.     ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1061 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE CARACTERE STATIONNAIRE DE          ',/,
     &'@       L''ECOULEMENT DE LA PHASE CONTINUE                   ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    ISTTIO DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI ISTTIO = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de ISTTIO dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1062 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE DECLENCHEMENT DU CALCUL            ',/,
     &'@       STATIONNAIRE DES STATISTIQUES POUR UN COUPLAGE RETOUR',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    NSTITS DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1     ',/,
     &'@       IL VAUT ICI NSTITS = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NSTITS dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1063 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE COUPLAGE RETOUR SUR LA DYNAMIQUE   ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    LTSDYN DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI LTSDYN = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de LTSDYN dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1064 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE COUPLAGE RETOUR SUR LA MASSE       ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    LTSMAS DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI LTSMAS = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de LTSMAS dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1065 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE COUPLAGE RETOUR SUR LA THERMIQUE   ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    LTSTHE DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI LTSTHE = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de LTSTHE dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1066 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE COUPLAGE RETOUR SUR LA DYNAMIQUE   ',/,
     &'@       EST ACTIVE (LTSDYN = ',I10,') (LAGOPT),              ',/,
     &'@       ALORS QUE LA PHASE PORTEUSE EST CALCULEE AVEC        ',/,
     &'@       L''OPTION CHAMP FIGE  (ICCVFG = ',I10,') (USINI1).   ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de LTSDYN dans la subroutine USLAG1 et ',/,
     &'@  verifier la valeur de ICCVFG dans la subroutine USINI1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1067 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********   (LAGOPT).                                   ',/,
     &'@                                                            ',/,
     &'@    L''INDICATEUR SUR LE COUPLAGE RETOUR EST ACTIVE         ',/,
     &'@        IILAGR = ',I10                                       ,/,
     &'@      ALORS QU''AUCUN COUPLAGE RETOUR N''EST ENCLENCHE      ',/,
     &'@        DYNAMIQUE : LTSDYN = ',I10                           ,/,
     &'@        THERMIQUE : LTSTHE = ',I10                           ,/,
     &'@        MASSIQUE  : LTSMAS = ',I10                           ,/,
     &'@                                                            ',/,
     &'@    LES COUPLAGES RETOUR SUR LA THERMIQUE ET SUR LA MASSE   ',/,
     &'@      NECESSITENT L''ACTIVATION D''UNE PHYSIQUE ADEQUATE    ',/,
     &'@      ASSOCIEE AUX PARTICULES.                              ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR dans la subroutine USLAG1.   ',/,
     &'@  Verifier la valeur de IPHYLA dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1070 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@   L''INDICATEUR SUR LE LANCEMENT DU CALCUL DES STATISTIQUES',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    ISTALA DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI ISTALA = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de ISTALA dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1071 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE SEUIL POUR LE CALCUL DES STATISTIQUES  ',/,
     &'@       VOLUMIQUES A UNE VALEUR NON PERMISE (LAGOPT).        ',/,
     &'@                                                            ',/,
     &'@    SEUIL DEVRAIT ETRE UN REEL SUPERIEUR OU EGAL A 0.D0     ',/,
     &'@       IL VAUT ICI SEUIL = ', E14.5                          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de SEUIL dans la subroutine USLAG1.    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1072 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE SEUIL POUR LE CALCUL DES STATISTIQUES  ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IDSTNT DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1     ',/,
     &'@       IL VAUT ICI IDSTNT = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IDSTNT dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1073 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE CALCUL STATIONNAIRE DES STATISTIQUES   ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    NSTIST DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL         ',/,
     &'@       A IDSTNT  = ', I10                                    ,/,
     &'@       IL VAUT ICI NSTIST = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NSTIST dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1074 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE NOMBRE DE VARIABLES                ',/,
     &'@       SUPPLEMENTAIRES LIEES AUX STATISTIQUES               ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@     NVLSTS DEVRAIT ETRE UN ENTIER ENTRE 0 ET ',I10          ,/,
     &'@       IL VAUT ICI NVLSTS = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NVLSTS dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
1075  FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE NOMBRE DE CLASSE DE STATISTIQUE    ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@     NBCLST DEVRAIT ETRE UN ENTIER ENTRE 1 ET ',I10          ,/,
     &'@       IL VAUT ICI NBCLST = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NBCLST dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1080 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''UTILISATION DE L''INJECTION        ',/,
     &'@       CONTINUE A UNE VALEUR NON PERMISE (LAGOPT).          ',/,
     &'@                                                            ',/,
     &'@    INJCON DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI INJCON = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de INJCON dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 1090 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''UTILISATION DE LA METHODE          ',/,
     &'@       DE CLONAGE/FUSION DES PARTICULES                     ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IROULE DEVRAIT ETRE UN ENTIER COMPRIS ENTRE 0 ET 2      ',/,
     &'@       IL VAUT ICI IROULE = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IROULE dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''ORDRE D''INTEGRATION               ',/,
     &'@       DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES          ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    NORDRE DEVRAIT ETRE UN ENTIER EGAL A 1 OU 2             ',/,
     &'@       IL VAUT ICI NORDRE = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NORDRE dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LA PRISE EN COMPTE DE LA DISPERSION   ',/,
     &'@       TURBULENTE A UNE VALEUR NON PERMISE (LAGOPT).        ',/,
     &'@                                                            ',/,
     &'@    IDISTU DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IDISTU = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IDISTU dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2011 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,
     &'@    DE TURBULENCE SELECTIONNE (LAGOPT).                     ',/,
     &'@                                                            ',/,
     &'@   Le module lagrangien a ete active avec IILAGR = ',I10     ,/,
     &'@     et la dispersion turbulente est prise en compte        ',/,
     &'@                                     avec IDISTU = ',I10     ,/,
     &'@   Le modele de turbulence active pour la phase ',I6         ,/,
     &'@     correspond a ITURB  = ',I10                             ,/,
     &'@   Or, les seuls traitements de la turbulence compatibles   ',/,
     &'@     avec le module Lagrangien et la dispersion turbulente  ',/,
     &'@     sont k-epsilon et Rij-epsilon, v2f et k-omega          ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,
     &'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,
     &'@  USINI1.                                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2012 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,
     &'@    DE TURBULENCE SELECTIONNE (LAGOPT).                     ',/,
     &'@                                                            ',/,
     &'@   Le module lagrangien a ete active avec IILAGR = ',I10     ,/,
     &'@     et la dispersion turbulente n''est pas prise en compte ',/,
     &'@                                     avec IDISTU = ',I10     ,/,
     &'@   Le modele de turbulence active pour la phase ',I6         ,/,
     &'@     correspond a ITURB  = ',I10                             ,/,
     &'@   Or, les seuls traitements de la turbulence compatibles   ',/,
     &'@     avec le module lagrangien sont : calcul laminaire,     ',/,
     &'@     k-epsilon, Rij-epsilon, v2f et k-omega.                ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,
     &'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,
     &'@  USINI1.                                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2013 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LA PRISE EN COMPTE DE LA DIFFUSION    ',/,
     &'@       TURBULENTE A UNE VALEUR NON PERMISE (LAGOPT).        ',/,
     &'@                                                            ',/,
     &'@    IDIFFL DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IDIFFL = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IDIFFL dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2014 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE CHOIX DU MODELE DE DISPERSION      ',/,
     &'@       TURBULENTE A UNE VALEUR NON PERMISE (LAGOPT).        ',/,
     &'@                                                            ',/,
     &'@    MODCPL DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 0     ',/,
     &'@       IL VAUT ICI MODCPL = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de MODCPL dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2015 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE CHOIX DU MODELE DE DISPERSION      ',/,
     &'@       TURBULENTE EST INCOMPATIBLE AVEC CELUI DU CALCUL     ',/,
     &'@       DES STATISTIQUES (LAGOPT).                           ',/,
     &'@                                                            ',/,
     &'@    LE MODELE COMPLET DE DISPERSION TURBULENTE EST ACTIVE   ',/,
     &'@      (MODCPL = ',I10,')'                                    ,/,
     &'@      AVANT LE DEBUT DU CALCUL DES STATISTIQUES             ',/,
     &'@      (IDSTNT = ',I10,')'                                    ,/,
     &'@                                                            ',/,
     &'@  Il est necessaire d''avoir calcule des statistiques       ',/,
     &'@  pour declencher le modele de dispersion turbulent complet.',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de MODCPL dans la subroutine USLAG1.   ',/,
     &'@  Verifier la valeur de IDSTNT dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2016 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE CHOIX DE LA DIRECTION DU MODELE COMPLET              ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IDIRLA DEVRAIT ETRE UN ENTIER EGAL A 1, 2 OU 3          ',/,
     &'@       (LA VALEUR 1 POUR UN ECOULEMENT SELON L''AXE X,      ',/,
     &'@        LA VALEUR 2 POUR UN ECOULEMENT SELON L''AXE Y,      ',/,
     &'@        LA VALEUR 3 POUR UN ECOULEMENT SELON L''AXE Z)      ',/,
     &'@       IL VAUT ICI IDIRLA = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IDIRLA dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2017 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE RESOLUTION DE L''EQUATION DE POISSON   ',/,
     &'@      POUR LES VITESSE MOYENNES ET DE CORRECTION DES        ',/,
     &'@      VITESSES INSTANTANNEES DES PARTICULES                 ',/,
     &'@      A UNE VALEUR NON PERMISE (LAGOPT).                    ',/,
     &'@                                                            ',/,
     &'@    ILAPOI DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI ILAPOI = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de ILAPOI dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2018 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR LE CHOIX DU MODELE DE DISPERSION      ',/,
     &'@      TURBULENTE A UNE VALEUR NON COMPATIBLE AVEC           ',/,
     &'@       CELUI DE L''ENCLENCHEMENT DES STATISTIQUES (LAGOPT). ',/,
     &'@                                                            ',/,
     &'@    LE MODELE COMPLET DE DISPERSION TURBULENTE EST ACTIVE   ',/,
     &'@      (MODCPL = ',I10,')'                                    ,/,
     &'@      AVANT L''ENCLENCHEMENT DU CALCUL DES STATISTIQUES     ',/,
     &'@      (ISTALA = ',I10,')'                                    ,/,
     &'@                                                            ',/,
     &'@  Il est necessaire d''activer le calcul des statistiques   ',/,
     &'@  pour declencher le modele de dispersion turbulent complet.',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de MODCPL dans la subroutine USLAG1.   ',/,
     &'@  Verifier la valeur de ISTALA dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2030 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING EN MODE TRAJECTOIRES   ',/,
     &'@      A UNE VALEUR NON PERMISE (LAGOPT).                    ',/,
     &'@                                                            ',/,
     &'@    IENSI1 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IENSI1 = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IENSI1 dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2031 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING EN MODE DEPLACEMENTS   ',/,
     &'@      PARTICULAIRES A UNE VALEUR NON PERMISE (LAGOPT).      ',/,
     &'@                                                            ',/,
     &'@    IENSI2 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IENSI2 = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IENSI2 dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2032 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE NOMBRE DE PARTICULE A VISUALISER EN POST-PROCESSING  ',/,
     &'@      A UNE VALEUR NON PERMISE (LAGOPT).                    ',/,
     &'@                                                            ',/,
     &'@    NBVIS DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 0      ',/,
     &'@       INFERIEUR AU NOMBRE MAX DE PARTICULES NBPMAX = ',I10  ,/,
     &'@       INFERIEUR AU PARAMETRE NLISTE = ',I10                 ,/,
     &'@                                                            ',/,
     &'@       IL VAUT ICI NBVIS = ', I10                            ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NBVIS dans la subroutine USLAG1.    ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2033 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE LA FREQUENCE D''ACQUISITION DES DONNES ',/,
     &'@       POUR LES SORTIES DE POST-PROCESSING EN               ',/,
     &'@       MODE TRAJECTOIRES OU MODE DEPLACEMENTS               ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    NVISLA DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1     ',/,
     &'@       IL VAUT ICI NVISLA = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NVISLA dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2040 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,
     &'@       "VITESSE DU FLUIDE VU"                               ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IVISV1 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVISV1 = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVISV1 dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2041 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,
     &'@       "VITESSE DES PARTICULES"                             ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IVISV2 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVISV2 = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVISV2 dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2042 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,
     &'@       "TEMPS DE SEJOUR DES PARTICULES"                     ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IVISTP DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVISTP = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVISTP dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2043 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,
     &'@       "DIAMETRE DES PARTICULES"                            ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IVISDM DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVISDM = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVISDM dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2044 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,
     &'@       "TEMPERATURE DES PARTICULES"                         ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IVISTE DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVISTE = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVISTE dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2045 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,
     &'@       "TEMPERATURE OU ENTHALPIE DES PARTICULES DE CHARBON" ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IVISHP DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVISHP = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVISHP dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2046 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,
     &'@       "DIAMETRE DU COEUR RETRECISSANT DES PARTICULES       ',/,
     &'@        DE CHARBON" A UNE VALEUR NON PERMISE (LAGOPT).      ',/,
     &'@                                                            ',/,
     &'@    IVISDK DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVISDK = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVISDK dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2047 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,
     &'@       "MASSE CHARBON REACTIF DES PARTICULES DE CHARBON"    ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IVISCH DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVISCH = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVISCH dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2048 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE POST-PROCESSING SUR LA VARIABLE        ',/,
     &'@       "MASSE DE COKE DES PARTICULES DE CHARBON"            ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IVISCK DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVISCK = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVISCK dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2050 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR SUR L''ACTIVATION DES STATISTIQUES        ',/,
     &'@       SUR LES INTERACTIONS PARTICULES/FRONTIERES           ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IENSI3 DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IENSI3 = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IENSI3 dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2051 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,
     &'@       "MASSE DE GRAINS DE CHARBON ENCRASSES"               ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IENCBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IENCBD = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IENCBD dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2052 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,
     &'@       "NOMBRE D''INTERACTIONS PARTICULES/FRONTIERES"       ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    INBRBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI INBRBD = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de INBRBD dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2053 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,
     &'@       "FLUX DE MASSE PARTICULAIRE LIE AUX INTERACTIONS"    ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IFLMBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IFLMBD = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IFLMBD dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2054 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,
     &'@       "ANGLE ENTRE LA VITESSE DE LA PARTICULE ET LE PLAN   ',/,
     &'@        DE LA FACE FRONTIERE"                               ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IANGBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IANGBD = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IANGBD dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2055 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE CALCUL DE LA STATISTIQUE AUX FRONTIERES',/,
     &'@       "VITESSE DE LA PARTICULE AU MOMENT DE L''INTERACTION"',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    IVITBD DEVRAIT ETRE UN ENTIER EGAL A 0 OU 1             ',/,
     &'@       IL VAUT ICI IVITBD = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IVITBD dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2056 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE NOMBRE D''INFORMATIONS SUPPLEMENTAIRES POUR LE       ',/,
     &'@       CALCUL DES STATISTIQUES AUX FRONTIERES               ',/,
     &'@       A UNE VALEUR NON PERMISE (LAGOPT).                   ',/,
     &'@                                                            ',/,
     &'@    NUSBOR DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 0     ',/,
     &'@       ET INFERIEUR OU EGAL A NUSBRD = ',I10                 ,/,
     &'@       IL VAUT ICI NUSBOR = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NUSBOR dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2057 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE CALCUL STATIONNAIRE DES STATISTIQUES   ',/,
     &'@       AUX FRONTIERES A UNE VALEUR NON PERMISE (LAGOPT).    ',/,
     &'@                                                            ',/,
     &'@    NSTBOR DEVRAIT ETRE UN ENTIER SUPERIEUR OU EGAL A 1     ',/,
     &'@       IL VAUT ICI NSTBOR = ', I10                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de NSTBOR dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2058 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE SEUIL POUR LE CALCUL DES STATISTIQUES  ',/,
     &'@       AUX FRONTIERES A UNE VALEUR NON PERMISE (LAGOPT).    ',/,
     &'@                                                            ',/,
     &'@    SEUILF DEVRAIT ETRE UN REEL SUPERIEUR OU EGAL A 0.D0    ',/,
     &'@       IL VAUT ICI SEUILF = ', E14.5                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de SEUILF dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2060 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@    LA STATISTIQUE AUX FRONTIERES                           ',/,
     &'@     "NOMBRE D''INTERACTIONS PARTICULES/FRONTIERES"         ',/,
     &'@     EST ACTIVE AVEC APPLICATION DE LA MOYENNE              ',/,
     &'@     PARTICULAIRE (LAGOPT).                                 ',/,
     &'@                                                            ',/,
     &'@    LES INDICATEURS DE CALCUL DES STATISTIQUES VALENT :     ',/,
     &'@       INBRBD       = ',I10                                  ,/,
     &'@       IMOYBR(INBR) = ',I10                                  ,/,
     &'@                                                            ',/,
     &'@  Le calcul continue mais risque de donner des affichages   ',/,
     &'@    et des sorties graphiques incoherentes. Une suite       ',/,
     &'@    de calcul est possible sans pertes de donnees.          ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IMOYBR(INBR) dans la  USLAG1.       ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2061 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@    LA STATISTIQUE AUX FRONTIERES :                         ',/,
     &'@ ',A60                                                       ,/,
     &'@     EST ACTIVE AVEC APPLICATION DE LA MOYENNE              ',/,
     &'@     PARTICULAIRE, ALORS QUE LE COMPTE DU NOMBRE            ',/,
     &'@     D''INTERACTIONS PARTICULES/FRONTIERES N''EST PAS       ',/,
     &'@     ENCLENCHE (LAGOPT).                                    ',/,
     &'@                                                            ',/,
     &'@    LES INDICATEURS DE CALCUL DES STATISTIQUES VALENT :     ',/,
     &'@       ',A12,     ' = ',I10                                  ,/,
     &'@             IMOYBR = ',I10                                  ,/,
     &'@                                                            ',/,
     &'@    INBRBD DEVRAIT ETRE A 1, IL VAUT ICI INBRBD = ',I10      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de INBRBD dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2062 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,
     &'@    *********                                               ',/,
     &'@    LA STATISTIQUE AUX FRONTIERES UTILISATEUR :             ',/,
     &'@ ',A60                                                       ,/,
     &'@     EST ACTIVE AVEC APPLICATION DE LA MOYENNE              ',/,
     &'@     PARTICULAIRE, ALORS QUE LE COMPTE DU NOMBRE            ',/,
     &'@     D''INTERACTIONS PARTICULES/FRONTIERES N''EST PAS       ',/,
     &'@     ENCLENCHE (LAGOPT).                                    ',/,
     &'@                                                            ',/,
     &'@    LES INDICATEURS DE CALCUL DES STATISTIQUES VALENT :     ',/,
     &'@       STAT UTILISATEUR NUMERO ',I10                         ,/,
     &'@       IMOYBR                = ',I10                         ,/,
     &'@                                                            ',/,
     &'@    INBRBD DEVRAIT ETRE A 1, IL VAUT ICI INBRBD = ',I10      ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de INBRBD dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 2063 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    L''INDICATEUR DE MOYENNE POUR LE CALCUL DES STATISTIQUES',/,
     &'@       AUX FRONTIERES A UNE VALEUR NON PERMISE (LAGOPT).    ',/,
     &'@                                                            ',/,
     &'@    IMOYBR DEVRAIT ETRE UN ENTIER EGAL A 0, 1 OU 2          ',/,
     &'@       IL VAUT ICI IMOYBR = ', I10                           ,/,
     &'@       POUR LA STATISTIQUE AUX FRONTERES :                  ',/,
     &'@ ',A60                                                       ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IMOYBR dans la subroutine USLAG1.   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3001 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE NOMBRE DE DONNEES REELLES SUR L''INJECTION           ',/,
     &'@     DES PARTICULES DANS LE DOMAINE DE CALCUL EST           ',/,
     &'@     SUPERIEUR AU NOMBRE MAXIMAL PREVU PAR DEFAUT.          ',/,
     &'@    (LAGOPT).                                               ',/,
     &'@                                                            ',/,
     &'@    LE NOMBRE DE DONNEES REELLES DEMANDEES VAUT : ',I10      ,/,
     &'@                                                            ',/,
     &'@    LE NOMBRE DE DONNEES MAXIMAL PREVU PAR DEFAUT VAUT      ',/,
     &'@      NDLAGM = ',I10                                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3002 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE NOMBRE DE DONNEES ENTIERES SUR L''INJECTION          ',/,
     &'@     DES PARTICULES DANS LE DOMAINE DE CALCUL EST           ',/,
     &'@     SUPERIEUR AU NOMBRE MAXIMAL PREVU PAR DEFAUT.          ',/,
     &'@    (LAGOPT).                                               ',/,
     &'@                                                            ',/,
     &'@    LE NOMBRE DE DONNEES ENTIERES DEMANDEES VAUT : ',I10     ,/,
     &'@                                                            ',/,
     &'@    LE NOMBRE DE DONNEES MAXIMAL PREVU PAR DEFAUT VAUT      ',/,
     &'@      NDLAIM = ',I10                                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3003 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                       (LAGOPT)',/,
     &'@    LE NOMBRE DE STATISTIQUES VOLUMIQUES LAGRANGIENNES      ',/,
     &'@      NVLSTA = ', I10                                        ,/,
     &'@    EST SUPERIEUR AU NOMBRE MAXIMAL PREVU                   ',/,
     &'@      NVPLMX = ', I10                                        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3004 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                       (LAGOPT)',/,
     &'@    LE NOMBRE DE POINTEURS DANS LES TABLEAUX DE             ',/,
     &'@      VARIABLES PARTICULAIRES ETTP ET ETTPA                 ',/,
     &'@      IRF = ', I10                                           ,/,
     &'@    EST SUPERIEUR AU NOMBRE MAXIMAL CALCULE                 ',/,
     &'@      NVP = ', I10                                           ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3005 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                       (LAGOPT)',/,
     &'@    LE NOMBRE DE POINTEURS DANS LE TABLEAUX DE              ',/,
     &'@      PARAMETRES PARTICULAIRES TEPA                         ',/,
     &'@      IRF = ', I10                                           ,/,
     &'@    EST SUPERIEUR AU NOMBRE MAXIMAL CALCULE                 ',/,
     &'@      NVEP = ', I10                                          ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3006 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                       (LAGOPT)',/,
     &'@    LE NOMBRE DE POINTEURS DANS LE TABLEAUX DE              ',/,
     &'@      PARAMETRES PARTICULAIRES ITEPA                        ',/,
     &'@      IRF = ', I10                                           ,/,
     &'@    EST SUPERIEUR AU NOMBRE MAXIMAL CALCULE                 ',/,
     &'@      NIVEP = ', I10                                         ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
 3010 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,
     &'@    *********                                               ',/,
     &'@    LE MODULE LAGRANGIEN EST INCOMPATIBLE AVEC LE MODELE    ',/,
     &'@    DE TURBULENCE SELECTIONNE (LAGOPT).                     ',/,
     &'@                                                            ',/,
     &'@   Le module lagrangien a ete active avec IILAGR = ',I10     ,/,
     &'@     et le couplage inverse sur la dynamique est pris en    ',/,
     &'@                              compte avec LTSDYN = ',I10     ,/,
     &'@   Le modele de turbulence active pour la phase ',I6         ,/,
     &'@     correspond a ITURB  = ',I10                             ,/,
     &'@   Or, les seuls traitements de la turbulence compatibles   ',/,
     &'@     avec le module Lagrangien et le couplage inverse sur   ',/,
     &'@     la dynamique sont k-epsilon, Rij-epsilon, v2f          ',/,
     &'@     et k-omega                                             ',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@  Verifier la valeur de IILAGR et IDISTU dans la subroutine ',/,
     &'@  USLAG1 et verifier la valeur de ITURB  dans la subroutine ',/,
     &'@  USINI1.                                                   ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
      RETURN
C
      END
c@z
