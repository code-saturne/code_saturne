c@a
                        SUBROUTINE USLAG1
C                       *****************
C     -------------------------------------------------------------
C     -------------------------------------------------------------
C
C***********************************************************************
c@versb
CVERS    CODE_SATURNE  VERSION 1.1            COPYRIGHT EDF 1999, 2003 *
c@verse
C***********************************************************************
C  FONCTION  : 
C  ---------
c@foncb
CFONC
CFONC    SOUS-PROGRAMME DU MODULE LAGRANGIEN :
CFONC    -------------------------------------
CFONC
CFONC        SOUS-PROGRAMME UTILISATEUR (INTERVENTION OBLIGATOIRE)
CFONC
CFONC        ROUTINE UTILISATEUR D'INITIALISATION DE CERTAINS PARAMETRES
CFONC        DU MODULE LAGRANGIEN. ILS CONCERNENT LES MODELES PHYSIQUES,
CFONC        NUMERIQUES, ET LES OPTIONS DE POSTPROCESSING.
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
      INCLUDE "entsor.h"
      INCLUDE "lagdim.h"
      INCLUDE "lagpar.h"
      INCLUDE "lagran.h"
C
C***********************************************************************
C VARIABLES LOCALES
C
      INTEGER          II , IPV
      DOUBLE PRECISION SIO2 , AL2O3 , FE2O3 , CAO
C
C***********************************************************************
C
C
C=======================================================================
C 1. ENCLENCHEMENT ET TYPE D'UTILISATION MODULE LAGRANGIEN
C=======================================================================
C
C     IILAGR = 0 : PAS DE CALCUL LAGRANGIEN (PAR DEFAUT)
C
C            = 1 : DIPHASIQUE LAGRANGIEN SANS COUPLAGE RETOUR
C
C            = 2 : DIPHASIQUE LAGRANGIEN AVEC COUPLAGE RETOUR
C                  ATTENTION, LE COUPLAGE RETOUR N'EST PRIS EN
C                  COMPTE QUE POUR LA PHASE CONTINUE NUMERO 1
C                  (on peut choisir de coupler la dynamique,
C                   la thermique et la masse independamment)
C
C            = 3 : DIPHASIQUE LAGRANGIEN SUR CHAMPS FIGES
C                  (cette option necessite une suite de calcul de la
C                   phase continue ISUITE = 1 ; les champs figes sont
C                   la totalite des champs euleriens, cette option
C                   n'est donc pas redondante avec celle du mot-cle
C                   ICCVFG qui permet un calul sur champs de
C                   vitesse et turbulents figes ; lorsque IILAGR = 3,
C                   on impose automatiquement ICCVFG = 1)
C
      IILAGR = 2
C
C=======================================================================
C 2. SUITES LAGRANGIENNES
C=======================================================================
C
C     ISUILA = 0 : PAS DE SUITE LAGRANGIENNE (PAR DEFAUT)
C            = 1 : SUITE LAGRANGIENNE
C                  (cette option necessite une suite de calcul de la
C                   phase continue ISUITE = 1)
C
      ISUILA = 0
C
C     SUITE DE CALCUL SUR DES STATISTIQUES VOLUMIQUES ET AUX FRONTIERES,
C     AINSI QUE LES TERMES SOURCES DE COUPLAGE RETOUR
C     (UTILE SI ISUILA = 1)
C     (DEFAUT NON : 0 ; OUI : 1)
C   
      IF (ISUILA.EQ.1) ISUIST = 0
C
C=======================================================================
C 3. MODELES PHYSIQUES LIES AUX PARTICULES
C=======================================================================
C
C     IPHYLA = 0 : PUREMENT DYNAMIQUE (PAR DEFAUT)
C            = 1 : EQUATIONS SUR Temperature (en degres Celsius), 
C                  Diametre et Masse
C            = 2 : CHARBON (Les particules  sont des grains de
C                  charbon pulverise. Cette option n'est disponible que 
C                  si la phase porteuse represente une flamme de charbon
C                  pulverise.)
C
      IPHYLA = 0
C
C
C
C=======================================================================
C 4. NOMBRE DE PARTICULES DU DOMAINE
C=======================================================================
C
C     NOMBRE DE PARTICULES MAXIMAL AUTORISE DANS LE DOMAINE
C     (PAR DEFAUT : NBPMAX = 1000)
C     * Attention, la mémoire est réservée en conséquence, et le nombre
C       max de particules traitees au cours d'une iteration Lagrangienne
C       est limite a NBPMAX.
C     * L'ordre de grandeur de NBPMAX doit etre evalue au mieux lors
C       d'injection a frequence non nulle en tenant compte :
C       - de la quantite de particules injectees (nombre par classe ou
C         en fonction du debit massique)
C       - de nombre de particules clonees ou detruites par la
C         technique de reduction de variance (IROULE = 1)
C
      NBPMAX = 30000
C
C=======================================================================
C 5. OPTIONS SUR LE TRAITEMENT DE LA PHASE DISPERSEE
C=======================================================================
C
C
C 5.1 VARIABLES SUPPLEMENTAIRES LIEES AUX PARTICULES
C --------------------------------------------------
C
C     * ces variables supplementaires sont stokees dans les tableaux
C       ETTP et ETTPA,
C     * on renseigne ici le nombre NVLS de variables supplementaires,
C     * la limite superieure de ce nombre est NUSVAR = 10
C       (fixee dans lagpar.h),
C     * on accede aux variables supplementaires dans ETTP et ETTPA
C       en utilisant le pointeur JVLS de la maniere suivante :
C       
C               etape courante -> ETTP(NBPT,JVLS(NVUS))
C               etape precedente -> ETTPA(NBPT,JVLS(NVUS))
C
C           NBPT est le numero de la particule traitee
C             (entier compris entre 1 et NBPART),
C           NVUS est le numero de la variable supplementaire
C             (entier compris entre 1 et NVLS),
C     * l'integration des Equations Differentielles Stochastiques
C       associees a ces nouvelles variables necessite une intervention
C       dans le sous-programme utilisateur USLAED.
C
C
      NVLS = 0
C
C
C 5.2 CARACTERE STATIONNAIRE DE L'ECOULEMENT DE LA PHASE CONTINUE
C ---------------------------------------------------------------
C
C     * si calcul stationnaire   : ISTTIO = 1
C     * si calcul instationnaire : ISTTIO = 0
C     * Lorsque les champs de la phase porteuse sont figes (IILAGR = 3)
C       l'indicateur est force en stationnaire (ISTTIO = 1), sinon
C       l'utilisateur doit le renseigner en prenant en compte la nature
C       de l'ecoulement de la phase continue.
C
C     REMARQUE : si ISTTIO = 0, les moyennes statistiques calculees
C                sont REMISE A ZERO a chaque iteration Lagrangienne.
C
      IF (IILAGR.NE.3) ISTTIO = 1
C
C
C 5.3 COUPLAGE RETOUR : INFLUENCE DE LA PHASE DISPERSEE SUR LA PHASE
C     CONTINUE (UTILE UNIQUEMENT SI IILAGR = 2)
C     ATTENTION : LE COUPLAGE RETOUR N'EST PRIS EN COMPTE QUE POUR
C                 LA PHASE CONTINUE NUMERO 1
C ------------------------------------------------------------------
C
      IF (IILAGR.EQ.2) THEN
C
C     * Nombre d'iterations Lagrangiennes absolues (i.e. suites comprises)
C       a partir duquel une moyenne en temps des termes sources de
C       couplage retour est calculee.
C     * Utile si CALCUL STATIONNAIRE, i.e. si ISTTIO = 1.
C     * Si le nombre d'iterations Lagrangiennes absolues est strictement
C       inferieur a NSTITS, les termes sources transmis sont
C       instationnaires (i.e. ils sont REMIS A ZERO a chaque iteration
C       Lagrangienne).
C     * La valeur minimale admissible pour NSTITS est 1.
C
        NSTITS = 1500 

C     COUPLAGE RETOUR SUR LA DYNAMIQUE (Vitesse + Turbulence)
C     (DEFAUT NON : 0 ; OUI : 1)
C     (UTILE UNIQUEMENT SI ICCVFG = 0)
C
        LTSDYN = 1
C
C     COUPLAGE RETOUR SUR LA MASSE (SI IPHYLA = 1 ET IMPVAR = 1)
C     (DEFAUT NON : 0 ; OUI : 1)
C
        IF(IPHYLA.EQ.1 .AND. IMPVAR.EQ.1) LTSMAS = 0
C
C     COUPLAGE RETOUR SUR LA THERMIQUE (SI IPHYLA = 1 ET ITPVAR = 1) OU 
C     LES VARIABLES CHARBON (SI IPHYLA = 2)
C     (DEFAUT NON : 0 ; OUI : 1)
C
        IF((IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) .OR. IPHYLA.EQ.2) LTSTHE = 0

      ENDIF

C
C 5.4 CALCUL DES STATISTIQUES VOLUMIQUES
C --------------------------------------
C
C   5.4.1 PARAMETRES GENERAUX
C   ~~~~~~~~~~~~~~~~~~~~~~~~~
C
C     CALCUL DES STATISTIQUES VOLUMIQUES
C     (DEFAUT NON : 0 ; OUI : 1)
C
      ISTALA = 1
C
C
      IF (ISTALA.EQ.1) THEN
C
C     SEUIL POUR LA PRISE EN COMPTE DES STATISTIQUES VOLUMIQUES
C     * La valeur de SEUIL est un poids statistique.
C     * Chaque cellule du maillage contient une certaine
C       quantite de particules en terme de poids statistique
C       (somme des poids statistiques de toutes les particules
C       contenues dans la cellule) ;
C       SEUIL est la valeur minimale a partir de laquelle la
C       contribution en poids statistique d'une cellule n'est
C       plus prise en compte dans le modele complet de dispersion
C       turbulente, dans la resolution de l'equation de poisson
C       de correction des vitesses moyennes, et dans les sorties
C       listing et post-processing
C     (DEFAUT : SEUIL = 0.D0)
C
        SEUIL = 0.D0
C
C     CALCUL DES STATISTIQUES VOLUMIQUES A PARTIR DE L'ITERATION
C     LAGRANGIENNE ABSOLUE
C     * IDSTNT est un nombre d'iterations Lagrangiennes absolues
C       (i.e. suites comprises).
C     (DEFAUT : IDSTNT = 1)
C
        IDSTNT = 1
C
C     CALCUL STATIONNAIRE DES STATISTIQUES VOLUMIQUES A PARTIR DE
C     L'ITERATION LAGRANGIENNE ABSOLUE NSTIST
C     * NSTIST est un nombre d'iterations Lagrangiennes absolues
C       (i.e. suites comprises) a partir du quel les statistiques sont
C       moyennees en temps.
C     * Utile si CALCUL STATIONNAIRE, i.e. si ISTTIO = 1.
C     * Si le nombre d'iterations Lagrangiennes absolues est strictement

C       inferieur a NSTIST, les statistiques sont instationnaires
C       (i.e. elles sont REMISES A ZERO a chaque iteration
C       Lagrangienne).
C     * La valeur minimale admissible pour NSTIST est IDSTNT.
C     (DEFAUT : NSTIST = IDSTNT)
C
        NSTIST = 2000 
C
C   5.4.2 NOMBRE DE VARIABLES STATISTIQUES VOLUMIQUES,
C         NOM DES VARIABLES POUR AFFICHAGE
C   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C    * Attention ce nom est utlise pour retrouver les informations
C      dans le fichier suite de calcul, si un nom est modifie
C      entre deux calculs, le resultat concerne est perdu.
C
C    * A priori l'utilisateur n'intervient que dans la rubrique 4)
C
C   1) Par defaut, statistiques TOUJOURS calculees :
C         Vitesse moyenne (3) 
C         Fluctuation de vitesse (3)
C         Taux de presence (i.e. concentration volumique) (1)
C         Somme des poids statistiques (1)
C
        NOMLAG(1)  = 'VitPartX'
        NOMLAG(2)  = 'VitPartY'
        NOMLAG(3)  = 'VitPartZ'
        NOMLAG(4)  = 'FVitParX'
        NOMLAG(5)  = 'FVitParY'
        NOMLAG(6)  = 'FVitParZ'
        NOMLAG(7)  = 'TauxVol'
        NOMLAG(8)  = 'SomPoids'
        IPV = 8

C
C   4) VARIABLES STATISTIQUES VOLUMIQUES SUPPLEMENTAIRES : 
C       * Si l'utilisateur souhaite des calculs de statistiques
C         autres que ceux fournis en standard, il doit
C         1) donner leur nombre NVLSTS,  2) renseigner leur nom,
C         et 3) intervenir dans les sous-programmes utilisateur 
C         USLAST.F et USLAEN.F pour programmer ses nouvelles 
C         statistiques (voir les exemples).
C       * Nombre de statistiques supplementaires maximal : 20
C         (sinon changer le parametre NUSSTA dans l'include lagpar.h)
C
        NVLSTS = 0
C
        IF (NVLSTS.GT.0) THEN
          DO II = 1,NVLSTS
            IPV = IPV + 1
            WRITE(NOMLAG(IPV),'(A6,I4.4)') 'VarLag',II
          ENDDO    
        ENDIF 
C
      ENDIF 
C
C
C=======================================================================
C 6. OPTIONS SUR L'ENTREE DES PARTICULES
C=======================================================================
C 
C     INJECTION EN CONTINUE DES PARTICLES PENDANT LE PAS DE TEMPS
C     (ET PAS UNIQUEMENT AU DEBUT DU PAS DE TEMPS, CETTE OPTION
C      EVITE LES PAQUETS DE PARTICULES EN FACE LES ZONES D'INJECTION)
C     (DEFAUT NON : 0 ; OUI : 1)
C
      INJCON = 0
C
C
C=======================================================================
C 7. TECHNIQUE DE REDUCTION DE VARIANCE : CLONAGE/FUSION DES PARTICULES
C=======================================================================
C
C     UTILISATION DE LA ROULETTE RUSSE (DEFAUT NON : 0 ; OUI : 1)
C 
      IROULE = 0
C
C
C=======================================================================
C 8. OPTIONS SUR LE TRAITEMENT NUMERIQUE DE LA PHASE DISPERSEE
C=======================================================================
C
C     ORDRE D'INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
C     PAR DEFAUT 2 (VALEURS ADMISSIBLES : 1 OU 2)
C
C       NORDRE = 1 : INTEGRATION DES EDS PAR UN SCHEMA D'ORDRE 1
C       NORDRE = 2 : INTEGRATION DES EDS PAR UN SCHEMA D'ORDRE 2
C 
      NORDRE = 2
C
C
C     RESOLUTION DE L'EQUATION DE POISSON POUR LES VITESSE MOYENNES
C     DES PARTICULES ET CORRECTION DES VITESSES INSTANTANNEES
C     DES PARTICULES
C      = 0 : pas de correction des vitesses (VALEUR PAR DEFAUT)
C      = 1 : correction des vitesses instantannes
C
C     ATTENTION : OPTION STRICTEMENT DEVELOPPEUR, LAISSEZ LA VALEUR PAR
C     =========   DEFAUT !
C
      ILAPOI = 0
C
C
C=======================================================================
C 9. OPTIONS SUR LE TRAITEMENT DE LA DISPERSION TURBULENTE
C=======================================================================
C
C     ATTENTION : DANS CETTE VERSION, LA DISPERSION TURBULENTE NE
C     ^^^^^^^^^^  FONCTIONNE QUE SI LA PHASE CONTINUE EST CALCULEE
C                 AVEC UN MODELE k-eps OU Rij-eps
C
C
C-->  PRISE EN COMPTE DE LA DISPERSION TURBULENTE
C     (DEFAUT OUI : 1 ; NON : 0)
C
      IDISTU = 1
C
C
C-->  DISPERSION TURBULENTE IMPOSEE A CELLE DU FLUIDE
C
C     SI ACTIF, ALORS LA DISPERSION TURBULENTE DES PARTICULES EST CELLE
C     DES PARTICULES FLUIDES. ON SE PLACE DONC DANS UN CAS DE DIFFUSION
C     TURBULENTE (ON SUPRIME LES EFFETS DE CROISEMENT DE
C     TRAJECTOIRES). SI LES PARTICULES SIMULEES ONT LA MASSE
C     VOLUMIQUE DU FLUIDE ALORS, ON SIMULE LE DEPLACEMENT DE
C     PARTICULES FLUIDES LAGRANGIENNES.
C     (DEFAUT NON : 0 ; OUI : 1)
C
      IDIFFL = 0
C
C     MODCPL : 
C          = 0 pour le modele incomplet (VALEUR PAR DEFAUT)
C          > 0 pour le modele complet, est egal au numero de l'iteration
C          Lagrangienne absolue (i.e. suites comprises) a partir de
C          laquelle mise le modele complet est active
C          MODCPL ne doit pas etre inferieur a IDSTNT
C
      MODCPL = 1
C
C     IDIRLA  =1 ou 2 ou 3 : 1ere, 2eme ou 3eme direction
C       du modele complet. Correspond à la direction principale
C       de l'écoulement. Permet de calculer un echelle de temps
C       Lagrangienne non isotrope.
C       (DEFAUT IDIRLA = 1)
C
      IF (MODCPL.GT.0) IDIRLA = 3
C
C
C=======================================================================
C 10. POST-PROCESSING
C=======================================================================
C
C 10.1 POST-PROCESSING DES TRAJECTOIRES ET DEPLACEMENTS PARTICULAIRES
C -------------------------------------------------------------------
C
C     ATTENTION : COUTEUX EN TEMPS DE CALCUL
C     ^^^^^^^^^^  
C
C   10.1.1 PARAMETRES GENERAUX
C   ~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C     MODE TRAJECTOIRES
C     (DEFAUT NON : 0 ; OUI : 1)
C
      IENSI1 = 0
C
C     MODE DEPLACEMENTS PARTICULAIRES
C     (DEFAUT NON : 0 ; OUI : 1)
C
      IENSI2 = 0
C
C
C     NOMBRE DE PARTICULES A VISUALISER MAXIMUM=NLISTE
C     (PAR DEFAUT NBVIS = NLISTE)
C     ATTENTION : NBVIS NE DOIT NI ETRE SUPERIEUR A NBPMAX
C                 NI A NLISTE (parametre fixe a 500 dans lagpar.h)
C
      NBVIS = NLISTE
C
C     PERIODE   D'AQUISITION DES DONNEES A VISUALISER
C     (PAR DEFAUT NVISLA = 1)
C
      NVISLA = 1
C
C     LISTE contient les numeros des particules que l'on souhaite
C     visualiser
C     (PAR DEFAUT LISTE(...) = -1, AUCUNE PARTICULE A VISUALISER)
C     ATTENTION : si le numero est negatif, il n'y a pas de visualisation.
C
C   > EXEMPLE 1 :
C     je veux voir les NBVIS 1eres
C
      DO II = 1, NBVIS
        LISTE(II) = II
      ENDDO
C
C   > REMARQUE : les trous, les repetitions dans le tableau LISTE seront
C                suprimes, et les numeros seront ranges par ordre
C                croissant.
C
C   10.1.2 VARIABLES A VISUALISER SUR LES TRAJECTOIRES OU LES PARTICULES
C   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C     VARIABLE : VITESSE DU FLUIDE VU
C     (DEFAUT NON : 0 ; OUI : 1)
C
      IVISV1  = 0
C
C     VARIABLE : VITESSE DE LA PARTICULE
C     (DEFAUT NON : 0 ; OUI : 1)
C
      IVISV2  = 0
C
C     VARIABLE : TEMPS DE SEJOUR
C     (DEFAUT NON : 0 ; OUI : 1)
C
      IVISTP  = 0
C
C     VARIABLE : DIAMETRE
C     (DEFAUT NON : 0 ; OUI : 1)
C
      IVISDM  = 0
C
C     VARIABLE : TEMPERATURE
C     (DEFAUT NON : 0 ; OUI : 1)
      IF (IPHYLA.EQ.1 .AND. ITPVAR.EQ.1) IVISTE  = 0
C
C     VARIABLE : MASSE
C     (DEFAUT NON : 0 ; OUI : 1)
C
      IVISMP  = 0
C
C
C 10.2 STATISTIQUES AUX FRONTIERES : VISUALISATION DES
C      INTERACTIONS PARTICULES/FRONTIERES
C ------------------------------------------------
C
C   10.2.1 PARAMETRES GENERAUX
C   ~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C     MODE INTERACTION PARTICULES/FRONTIERES
C     (DEFAUT NON : 0 ; OUI : 1)
C
      IENSI3 = 0
C
C     CALCUL STATIONNAIRE DES STATISTIQUES AUX FRONTIERES A PARTIR DE
C     L'ITERATION LAGRANGIENNE ABSOLUE NSTBOR
C     * NSTBOR est un nombre d'iterations Lagrangiennes absolues
C       (i.e. suites comprises) a partir duquel les statistiques sont
C       moyennees (en temps ou par nombre d'interactions).
C     * Utile si CALCUL STATIONNAIRE, i.e. si ISTTIO = 1.
C     * Si le nombre d'iterations Lagrangiennes absolues est strictement
C       inferieur a NSTBOR, les statistiques sont instationnaires
C       (i.e. elles sont REMISES A ZERO a chaque iteration
C       Lagrangienne).
C     * La valeur minimale admissible pour NSTBOR est 1.
C     (DEFAUT : NSTBOR = 1)
C
      NSTBOR = 1
C
C     SEUILF POUR LA PRISE EN COMPTE DES STATISTIQUES AUX FRONTIERES
C     * La valeur de SEUILF est un poids statistique.
C     * Chaque face du maillage a subi un certain nombre d'interactions
C       avec des particules en terme de poids statistique
C       (somme des poids statistiques de toutes les particules
C       qui ont interagi avec la face de bord) ;
C       SEUILF est la valeur minimale a partir de laquelle la
C       contribution en terme statistique d'une face n'est
C       plus prise en compte dans les sorties listing
C       et post-processing.
C     (DEFAUT : SEUILF = 0.D0)
C
      SEUILF = 0.D0
C
C
C   10.2.2 INFORMATIONS A ENREGISTRER
C   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C    * Un certain nombre d'informations succeptibles d'interesser
C      l'utilisateur sont d'ores et deja ecrites dans le sous-programme
C      USLABO. Pour les activer il faut mettre le mot-cle correspondant
C      a 1 ci-dessous.
C    * La selection des types d'interaction (IREBOL, IDEPO1,... voir
C      le sous-programme USLAG2) qu'enclenchent l'enregistrement
C      des informations est faite dans USLABO. La selection par defaut
C      doit etre validee ou modifiee par l'utilisateur.
C    * Par defaut, dans un meme enregistrement, on place les informations
C      demandees pour toutes les interactions particules/frontieres
C      selectionnees. Modifier ce comportement necessite une 
C      intervention dans le sous-programme utilisateur USLABO.
C    * La statistique aux frontieres :
C      "NOMBRE D'INTERACTIONS PARTICULES/FRONTIERES" doit etre
C      selectionnee obligatoirement pour utiliser la moyenne
C      particulaire IMOYBR(...) = 2.
C
C
C     NOMBRE D'INTERACTIONS PARTICULES/FRONTIERES
C     (DEFAUT NON : 0 ; OUI : 1)
      INBRBD = 1
C
C     FLUX DE MASSE PARTICULAIRE LIE AUX INTERACTIONS
C     PARTICULES/FRONTIERES
C     (DEFAUT NON : 0 ; OUI : 1)
      IF (INBRBD.EQ.1) IFLMBD = 0
C
C     ANGLE ENTRE LA VITESSE DE LA PARTICULE
C     ET LE PLAN DE LA FACE FRONTIERE
C     (DEFAUT NON : 0 ; OUI : 1)
      IANGBD = 0
C
C     NORME DE LA VITESSE DE LA PARTICULE AU MOMENT 
C     DE SON INTERACTION AVEC LA FACE FRONTIERE
C     (DEFAUT NON : 0 ; OUI : 1)
      IVITBD = 0
C
C     MASSE DE GRAINS DE CHARBON ENCRASSES
C     (DEFAUT NON : 0 ; OUI : 1)
       IF (IPHYLA.EQ.2 .AND. IENCRA.EQ.1) IENCBD = 0
C
C     INFORMATIONS UTILISATEUR SUPPLEMENTAIRES A ENREGISTREES
C     (PAR EXEMPLE TAUX D'EROSION, TEMPERATURE...)
C     * ces enregistrements supplementaires sont stokees dans
C       le tableau PARBOR,
C     * on renseigne ici le nombre NUSBOR d'enregistrements
C       supplementaires,
C     * la limite supperieure de ce nombre est NUSBRD = 10
C       (fixee dans lagpar.h),
C     * voir un exemple de code d'un enregistrement supplementaire
C       dans le sous-programme USLABO.
C
      NUSBOR = 0
C    
C
C   10.2.3 NOM DES ENREGISTREMENTS POUR AFFICHAGE,
C          MOYENNE EN TEMPS OU MOYENNE PARTICULAIRE
C          DES STATISTIQUES AUX FRONTIERES
C   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C    * A priori l'utilisateur n'intervient que dans les infos
C      utilisateur supplementaires a enregistrer : il doit
C      donner le nom de l'enregistrement ainsi que le type
C      de moyenne que l'on souhaite lui appliquer pour les affichages
C      listing et en fin de calcul pour le post-processing.
C
C    * ATTENTION ce nom est utilise pour retrouver les informations
C      dans le fichier suite de calcul, si un nom est modifie
C      entre deux calculs, le resultat concerne n'est pas relu 
C      (et est donc perdu).
C
C    * La moyenne appliquee est donnee par l'intermediaire du
C      tableau IMOYBR :
C       - si IMOYBR(IUSB(II)) = 0 -> pas de moyenne appliquee
C       - si IMOYBR(IUSB(II)) = 1 -> on applique une moyenne temporelle,
C               c'est a dire que la statistique est divisee par le
C               dernier pas de temps dans le cas d'un enregistrement
C               instationnaire, ou stationnaire mais avec un nombre
C               d'iterations inferieur a NSTBOR ; ou que la statistique
C               est divisee par le temps d'enregistrement
C               stationnaire dans le cas stationnaire.
C       - si IMOYBR(IUSB(II)) = 2 -> on applique une moyenne
C               particulaire, c'est a dire qu'on divise la statistique
C               par le nombre d'interactions particules/frontieres
C               enregistrees (en terme de poids statistique)
C               dans PARBOR(NFABOR,INBR) (cf. USLABO).
C               Pour utiliser cette moyenne il faut que INBRBD=1.
C    * Les sauvegardes dans le fichier suite de calcul sont faites
C      sans application de cette moyenne.
C    * La moyenne est appliquee si le nombre d'interactions (en poids
C      statistique) de la face de bord consideree est strictement 
C      superieur a SEUILF, sinon la moyenne est mise a zero.
C
      IPV = 0
C
      IF (IENSI3.EQ.1) THEN
C
        IF (INBRBD.EQ.1) THEN
          IPV = IPV + 1
          INBR = IPV
          NOMBRD(INBR) = 'nombreImpact'
          IMOYBR(INBR) = 0
        ENDIF
C
        IF (IFLMBD.EQ.1) THEN
          IPV = IPV + 1
          IFLM = IPV
          NOMBRD(IFLM) = 'fluxDeMasse'
          IMOYBR(IFLM) = 2
        ENDIF
C
        IF (IANGBD.EQ.1) THEN
          IPV = IPV + 1
          IANG = IPV
          NOMBRD(IANG) = 'angleImpact'
          IMOYBR(IANG) = 1
        ENDIF
C
        IF (IVITBD.EQ.1) THEN
          IPV = IPV + 1
          IVIT = IPV
          NOMBRD(IVIT) = 'normeVitImpact'
          IMOYBR(IVIT) = 1
        ENDIF
C
        IF (IENCBD.EQ.1) THEN
          IPV = IPV + 1
          IENC = IPV
          NOMBRD(IENC) = 'masseEncras'
          IMOYBR(IENC) = 0
        ENDIF
C
        IF (NUSBOR.GT.0) THEN
          DO II = 1,NUSBOR
            IPV = IPV + 1
            IUSB(II) = IPV
            WRITE(NOMBRD(IUSB(II)),'(A8,I4.4)') 'enrSupp',II
            IMOYBR(IUSB(II)) = 0
          ENDDO
        ENDIF
C
      ENDIF
C
C
C=======================================================================
C
      RETURN
C
      END
c@z
