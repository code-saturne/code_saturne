!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

                  subroutine uslag1
!================



!===============================================================================
!  FONCTION  :
!  ---------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!       SOUS-PROGRAMME UTILISATEUR (INTERVENTION OBLIGATOIRE)

!       ROUTINE UTILISATEUR D'INITIALISATION DE CERTAINS PARAMETRES
!       DU MODULE LAGRANGIEN. ILS CONCERNENT LES MODELES PHYSIQUES,
!       NUMERIQUES, ET LES OPTIONS DE POSTPROCESSING.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "entsor.h"
include "lagdim.h"
include "lagpar.h"
include "lagran.h"

!===============================================================================

! VARIABLES LOCALES

integer          ii , ipv , icha
double precision sio2 , al2o3 , fe2o3 , cao

!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

iilagr = 0

if(1.eq.1) then
  return
endif

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================
! 1. ENCLENCHEMENT ET TYPE D'UTILISATION MODULE LAGRANGIEN
!===============================================================================

!     IILAGR = 0 : PAS DE CALCUL LAGRANGIEN (PAR DEFAUT)

!            = 1 : DIPHASIQUE LAGRANGIEN SANS COUPLAGE RETOUR

!            = 2 : DIPHASIQUE LAGRANGIEN AVEC COUPLAGE RETOUR
!                  ATTENTION, LE COUPLAGE RETOUR N'EST PRIS EN
!                  COMPTE QUE POUR LA PHASE CONTINUE NUMERO 1
!                  (on peut choisir de coupler la dynamique,
!                   la thermique et la masse independamment)

!            = 3 : DIPHASIQUE LAGRANGIEN SUR CHAMPS FIGES
!                  (cette option necessite une suite de calcul de la
!                   phase continue ISUITE = 1 ; les champs figes sont
!                   la totalite des champs euleriens, cette option
!                   n'est donc pas redondante avec celle du mot-cle
!                   ICCVFG qui permet un calul sur champs de
!                   vitesse et turbulents figes ; lorsque IILAGR = 3,
!                   on impose automatiquement ICCVFG = 1)

iilagr = 1

!===============================================================================
! 2. SUITES LAGRANGIENNES
!===============================================================================

!     ISUILA = 0 : PAS DE SUITE LAGRANGIENNE (PAR DEFAUT)
!            = 1 : SUITE LAGRANGIENNE
!                  (cette option necessite une suite de calcul de la
!                   phase continue ISUITE = 1)

isuila = 0

!     SUITE DE CALCUL SUR DES STATISTIQUES VOLUMIQUES ET AUX FRONTIERES,
!     AINSI QUE LES TERMES SOURCES DE COUPLAGE RETOUR
!     (UTILE SI ISUILA = 1)
!     (DEFAUT NON : 0 ; OUI : 1)

if (isuila.eq.1) isuist = 0

!===============================================================================
! 3. MODELES PHYSIQUES LIES AUX PARTICULES
!===============================================================================

!     IPHYLA = 0 : PUREMENT DYNAMIQUE (PAR DEFAUT)
!            = 1 : EQUATIONS SUR Temperature (en degres Celsius),
!                  Diametre et Masse
!            = 2 : CHARBON (Les particules  sont des grains de
!                  charbon pulverise. Cette option n'est disponible que
!                  si la phase porteuse represente une flamme de charbon
!                  pulverise.)

iphyla = 0


! 3.1 OPTION POUR EQUATION SUR TP, DP et MP (UNIQUEMENT SI IPHYLA = 1)


if (iphyla.eq.1) then

!      EQUATION SUR LE DIAMETRE
!      (DEFAUT NON : 0 ; OUI : 1)

  idpvar = 0

!      EQUATION SUR LA TEMPERATURE
!      (DEFAUT NON : 0 ; OUI : 1)
!      Cette option impose l'existance d'une variable thermique
!      sur la phase continue (physique paticuliere ou non).
!      La variable resolue est la temperature en degres Celsius.

  itpvar = 0

!      EQUATION SUR LA MASSE
!      (DEFAUT NON : 0 ; OUI : 1)

  impvar = 0

endif

!     * Dans le cas ou une equation sur la temperature des particules
!       est enclenchee entre deux suites de calcul (IPHYLA= 1 et ITPVAR= 1)
!       il faut fournir une temperature (en degres Celsius) et une
!       chaleur massique (J/kg/K) d'initialisation des particules
!       deja presentes dans le domaine de calcul.

if (isuila.eq.1 .and. iphyla.eq.1 .and. itpvar.eq.1) then
  tpart = 700.d0
  cppart = 5200.d0
endif

! 3.2 OPTIONS POUR l'ENCRASSEMENT (CHARBON UNIQUEMENT i.e. IPHYLA = 2)


!     RAPPORTS EDF/R&D DE REFERENCE : HI-81/00/030/A
!                                     HI-81/01/033/A

!     D'une maniere generale, un calcul d'encrassement est effectue
!     en deux etapes : 1) on realise un calcul "classique" de chaudiere
!     a charbon pulverise par une approche homogene ; 2) on fait un
!     calcul Lagrangien de grains de charbon sur les champs
!     figes du calcul precedent. L'option d'encrassement est donc
!     utilise en "post-processing" sur un calcul charbon pulverise
!     homogene.


!     On utilise la probabilite qu'une particule a la temperature Tp
!     colle sur sa surface d'impact. Cette probabilite est le rapport
!     d'une viscosite critique sur la viscosite des cendres.

!              VISREF
!     P(Tp) = --------   pour VISCEN >= VISREF
!              VISCEN

!           = 1 sinon

!     Pour evaluer la viscosite des cendres VISCEN, on utilise
!     l'expression de J.D. Watt et T.Fereday (J.Inst.Fuel-Vol42-p99)

!                         ENC1 * 1.0D+7
!     Log  (10*VISCEN) = --------------- + ENC2
!        10                            2
!                        (Tp(°C) - 150)

!     La viscosite critique VISREF est donnée par la litterature et
!     peut varier entre 8 Pa.set 1.D7 Pa.s           !
!     En general on prend 1.0D+4 Pa.s...



if (iphyla.eq.2) then

!       IENCRA = 0 pas d'encrassement (PAR DEFAUT)
!              = 1 encrassement

!       * Il faut preciser les frontieres du domaine sur lesquelles
!         les grains de charbon peuvent s'encrasser (voir USLAG2).
!       * Le traitement de l'encrassement se deroule dans USLABO.
!       * Pour visualiser la cartographie de la masse de particules
!         encrassees, il faut que IENSI3 = 1 et que IENCBD = 1
!         (voir la rubrique 10.2 du present sous-programme).

  iencra = 0

!     DEFINITION DES CRITERES D'ENCRASSEMENT POUR LES DIFFERENTS CHARBONS
!     ATTENTION A BIEN LES DEFINIR POUR CHACUN DES NCHARB CHARBONS INJECTES
!                                       ======

!     Premier charbon ICHA = 1

  icha = 1

!       TPRENC : TEMPERATURE SEUIL EN DESSOUS DE LAQUELLE LES GRAINS
!                DE CHARBON NE S'ENCRASSENT PAS (en degres Celsius)

  tprenc(icha) = 600.d0

!       VISREF : VISCOSITE CRITIQUE (en Pa.s)

  visref(icha) = 10000.d0

!    > Exemple de composition de charbon en matieres minerales :
!       (avec SiO2 + Al2O3 + Fe2O3 + CaO + MgO = 100% en masse)

  sio2   =  36.0d0
  al2o3  =  20.8d0
  fe2o3  =   4.9d0
  cao    =  13.3d0

!       ENC1 et ENC2 : COEFFICIENTS DE L'EXPRESSION DE Watt et Fereday

  enc1(icha) = 0.00835d0 * sio2 + 0.00601d0 * al2o3 - 0.109d0

  enc2(icha) = 0.0415d0 * sio2  + 0.0192d0 * al2o3                &
       + 0.0276d0 * fe2o3 + 0.016 * cao - 3.92d0

endif

!===============================================================================
! 4. NOMBRE DE PARTICULES DU DOMAINE
!===============================================================================

!     NOMBRE DE PARTICULES MAXIMAL AUTORISE DANS LE DOMAINE
!     (PAR DEFAUT : NBPMAX = 1000)
!     * Attention, la mémoire est réservée en conséquence, et le nombre
!       max de particules traitees au cours d'une iteration Lagrangienne
!       est limite a NBPMAX.
!     * L'ordre de grandeur de NBPMAX doit etre evalue au mieux lors
!       d'injection a frequence non nulle en tenant compte :
!       - de la quantite de particules injectees (nombre par classe ou
!         en fonction du debit massique)
!       - de nombre de particules clonees ou detruites par la
!         technique de reduction de variance (IROULE = 1 ou IROULE =2)

nbpmax = 1000

!===============================================================================
! 5. OPTIONS SUR LE TRAITEMENT DE LA PHASE DISPERSEE
!===============================================================================


! 5.1 VARIABLES SUPPLEMENTAIRES LIEES AUX PARTICULES
! --------------------------------------------------

!     * ces variables supplementaires sont stokees dans les tableaux
!       ETTP et ETTPA,
!     * on renseigne ici le nombre NVLS de variables supplementaires,
!     * la limite superieure de ce nombre est NUSVAR = 10
!       (fixee dans lagpar.h),
!     * on accede aux variables supplementaires dans ETTP et ETTPA
!       en utilisant le pointeur JVLS de la maniere suivante :

!               etape courante -> ETTP(NBPT,JVLS(NVUS))
!               etape precedente -> ETTPA(NBPT,JVLS(NVUS))

!           NBPT est le numero de la particule traitee
!             (entier compris entre 1 et NBPART),
!           NVUS est le numero de la variable supplementaire
!             (entier compris entre 1 et NVLS),
!     * l'integration des Equations Differentielles Stochastiques
!       associees a ces nouvelles variables necessite une intervention
!       dans le sous-programme utilisateur USLAED.


nvls = 0


! 5.2 CARACTERE STATIONNAIRE DE L'ECOULEMENT DE LA PHASE CONTINUE


!     * si calcul stationnaire   : ISTTIO = 1
!     * si calcul instationnaire : ISTTIO = 0
!     * Lorsque les champs de la phase porteuse sont figes (IILAGR = 3)
!       l'indicateur est force en stationnaire (ISTTIO = 1), sinon
!       l'utilisateur doit le renseigner en prenant en compte la nature
!       de l'ecoulement de la phase continue.

!     REMARQUE : si ISTTIO = 0, les moyennes statistiques calculees
!                sont REMISE A ZERO a chaque iteration Lagrangienne.

if (iilagr.ne.3) isttio = 0


! 5.3 COUPLAGE RETOUR : INFLUENCE DE LA PHASE DISPERSEE SUR LA PHASE
!     CONTINUE (UTILE UNIQUEMENT SI IILAGR = 2)
!     ATTENTION : LE COUPLAGE RETOUR N'EST PRIS EN COMPTE QUE POUR
!                 LA PHASE CONTINUE NUMERO 1


if (iilagr.eq.2) then

!     * Nombre d'iterations Lagrangiennes absolues (i.e. suites comprises)
!       a partir duquel une moyenne en temps des termes sources de
!       couplage retour est calculee.
!     * Utile si CALCUL STATIONNAIRE, i.e. si ISTTIO = 1.
!     * Si le nombre d'iterations Lagrangiennes absolues est strictement
!       inferieur a NSTITS, les termes sources transmis sont
!       instationnaires (i.e. ils sont REMIS A ZERO a chaque iteration
!       Lagrangienne).
!     * La valeur minimale admissible pour NSTITS est 1.

  nstits = 1

!     COUPLAGE RETOUR SUR LA DYNAMIQUE (Vitesse + Turbulence)
!     (DEFAUT NON : 0 ; OUI : 1)
!     (UTILE UNIQUEMENT SI ICCVFG = 0)

  ltsdyn = 0

!     COUPLAGE RETOUR SUR LA MASSE (SI IPHYLA = 1 ET IMPVAR = 1)
!     (DEFAUT NON : 0 ; OUI : 1)

  if(iphyla.eq.1 .and. (impvar.eq.1 .or. idpvar.eq.1)) ltsmas = 0

!     COUPLAGE RETOUR SUR LA THERMIQUE (SI IPHYLA = 1 ET ITPVAR = 1) OU
!     LES VARIABLES CHARBON (SI IPHYLA = 2)
!     (DEFAUT NON : 0 ; OUI : 1)

  if((iphyla.eq.1 .and. itpvar.eq.1) .or. iphyla.eq.2) ltsthe = 0

endif


! 5.4 CALCUL DES STATISTIQUES VOLUMIQUES
! --------------------------------------

!   5.4.1 PARAMETRES GENERAUX
!   ~~~~~~~~~~~~~~~~~~~~~~~~~

!     CALCUL DES STATISTIQUES VOLUMIQUES
!     (DEFAUT NON : 0 ; OUI : 1)

istala = 0


if (istala.eq.1) then

!     SEUIL POUR LA PRISE EN COMPTE DES STATISTIQUES VOLUMIQUES
!     * La valeur de SEUIL est un poids statistique.
!     * Chaque cellule du maillage contient une certaine
!       quantite de particules en terme de poids statistique
!       (somme des poids statistiques de toutes les particules
!       contenues dans la cellule) ;
!       SEUIL est la valeur minimale a partir de laquelle la
!       contribution en poids statistique d'une cellule n'est
!       plus prise en compte dans le modele complet de dispersion
!       turbulente, dans la resolution de l'equation de poisson
!       de correction des vitesses moyennes, et dans les sorties
!       listing et post-processing
!     (DEFAUT : SEUIL = 0.D0)

  seuil = 0.d0

!     CALCUL DES STATISTIQUES VOLUMIQUES A PARTIR DE L'ITERATION
!     LAGRANGIENNE ABSOLUE
!     * IDSTNT est un nombre d'iterations Lagrangiennes absolues
!       (i.e. suites comprises).
!     (DEFAUT : IDSTNT = 1)

  idstnt = 1

!     CALCUL STATIONNAIRE DES STATISTIQUES VOLUMIQUES A PARTIR DE
!     L'ITERATION LAGRANGIENNE ABSOLUE NSTIST
!     * NSTIST est un nombre d'iterations Lagrangiennes absolues
!       (i.e. suites comprises) a partir du quel les statistiques sont
!       moyennees en temps.
!     * Utile si CALCUL STATIONNAIRE, i.e. si ISTTIO = 1.
!     * Si le nombre d'iterations Lagrangiennes absolues est strictement
!       inferieur a NSTIST, les statistiques sont instationnaires
!       (i.e. elles sont REMISES A ZERO a chaque iteration
!       Lagrangienne).
!     (DEFAUT : NSTIST = IDSTNT)

  nstist = idstnt


!   5.4.2 NOMBRE DE VARIABLES STATISTIQUES VOLUMIQUES,
!         NOM DES VARIABLES POUR AFFICHAGE
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!    NOMLAG : nom de la moyenne
!    NOMLAV : nom de la variance
!    IHSLAG : sortie des enregistrements sur les capteurs definis
!             dans USINI1

!    * ATTENTION : respecter l'ordre d'apparition des statistiques.

!    * Attention ces noms sont utlises pour retrouver les informations
!      dans le fichier suite de calcul ; donc si un nom est modifie
!      entre deux calculs, la statistique concernee est perdue.

!    * A priori, l'utilisateur n'intervient que dans les rubriques 4) et 6)

!   1) Par defaut, statistiques TOUJOURS calculees :
!         Moyenne et variance des composantes de la vitesse
!         Moyenne et variance du taux de presence
!            (i    .e. concentration volumique)
!         Moyenne et variance du temps de séjour

  ipv  =  1
  NOMLAG(IPV)  = 'MoVitPtX'
  NOMLAV(IPV)  = 'VaVitPtX'
  ihslag(ipv)  = 2

  ipv  = ipv  + 1
  NOMLAG(IPV)  = 'MoVitPtY'
  NOMLAV(IPV)  = 'VaVitPtY'
  ihslag(ipv)  = 2

  ipv  = ipv  + 1
  NOMLAG(IPV)  = 'MoVitPtZ'
  NOMLAV(IPV)  = 'VaVitPtZ'
  ihslag(ipv)  = 2

  ipv  = ipv  + 1
  NOMLAG(IPV)  = 'MoTauVol'
  NOMLAV(IPV)  = 'VaTauVol'
  ihslag(ipv)  = 2

  ipv  = ipv  + 1
  NOMLAG(IPV)  = 'MoTpsSej'
  NOMLAV(IPV)  = 'VaTpsSej'
  ihslag(ipv)  = 2

!   2) Physiques pariculieres (IPHYLA = 1) suivant les options choisies:
!         Moyenne et variance de la temperature
!         Moyenne et variance du diametre
!         Moyenne et variance de la masse

  if (iphyla.eq.1) then

    if (itpvar.eq.1) then
      ipv  = ipv  + 1
      NOMLAG(IPV) = 'MoTempPt'
      NOMLAV(IPV) = 'VaTempPt'
      ihslag(ipv)  = 2
    endif
    if (idpvar.eq.1) then
      ipv  = ipv  + 1
      NOMLAG(IPV) = 'MoDiamPt'
      NOMLAV(IPV) = 'VaDiamPt'
      ihslag(ipv)  = 2
    endif
    if (impvar.eq.1) then
      ipv  = ipv  + 1
      NOMLAG(IPV) = 'MoMassPt'
      NOMLAV(IPV) = 'VaMassPt'
      ihslag(ipv)  = 2
    endif

  else if (iphyla.eq.2) then

!   3) Charbon pulverise (IPHYLA = 2) :
!         Moyenne et variance de la temperature
!         Moyenne et variance de la masse de charbon reactif
!         Moyenne et variance de la masse de coke
!         Moyenne et variance du diametre du coeur retrecissant

    ipv  = ipv  + 1
    NOMLAG(IPV) = 'MoTempPt'
    NOMLAV(IPV) = 'VaTempPt'
    ihslag(ipv)  = 2

    ipv  = ipv  + 1
    NOMLAG(IPV) = 'MoMchPt'
    NOMLAV(IPV) = 'VaMchPt'
    ihslag(ipv)  = 2

    ipv  = ipv  + 1
    NOMLAG(IPV) = 'MoMcokPt'
    NOMLAV(IPV) = 'VaMcokPt'
    ihslag(ipv)  = 2

    ipv  = ipv  + 1
    NOMLAG(IPV) = 'MoDckPt'
    NOMLAV(IPV) = 'VaDckPt'
    ihslag(ipv)  = 2

  endif

!   4) VARIABLES STATISTIQUES VOLUMIQUES SUPPLEMENTAIRES :
!       * Si l'utilisateur souhaite des calculs de statistiques
!         autres que ceux fournis en standard, il doit
!         1) donner leur nombre NVLSTS,  2) renseigner leur nom,
!         3) renseigner IHSLAG
!         et 4) intervenir dans les sous-programmes utilisateur
!         USLAST.F et USLAEN.F pour programmer ses nouvelles
!         statistiques (voir les exemples).
!       * Nombre maximal de statistiques supplementaires : 20
!         (sinon changer le parametre NUSSTA dans l'include lagpar.h)

  nvlsts = 0

  if (nvlsts.gt.0) then
    do ii = 1,nvlsts
      ilvu(ii) = ipv + ii
      WRITE(NOMLAG(ILVU(II)),'(A6,I4.4)') 'MoyLag',II
      WRITE(NOMLAV(ILVU(II)),'(A6,I4.4)') 'VarLag',II
      ihslag(ilvu(ii))  = 1
    enddo
    ipv = ipv + nvlsts
  endif

!   5) Par defaut, statistique TOUJOURS calculee :
!         Somme du poids statistiques associé aux particules
!            (i    .e. nombre de particules par cellules)

  ipv  = ipv  + 1
  NOMLAG(IPV)  = 'SomPoids'
  ihslag(ipv)  = 1

!   6) STATISTIQUES PAR GROUPE :
!       * Si l'utilisateur souhaite des calculs de statistiques
!         par groupe de particule (par defaut aucune statistique par groupe),
!         il doit :
!         1) donner NBCLST le nombre de groupe (limite a 100)
!         2) preciser dans USLAG2 le groupe auquel appartient chaque particule
!              a l'aide du tableau IUSLAG

!      * ATTENTION : NBCLST ne peut pas etre change lors d'une suite de
!        calcul (ISUILA=1) meme si le calcul des statistiques n'est pas
!        encore enclenche (ISTALA=0)

  nbclst = 0

endif


!===============================================================================
! 6. OPTIONS SUR L'ENTREE DES PARTICULES
!===============================================================================

!     INJECTION EN CONTINUE DES PARTICLES PENDANT LE PAS DE TEMPS
!     (ET PAS UNIQUEMENT AU DEBUT DU PAS DE TEMPS, CETTE OPTION
!      EVITE LES PAQUETS DE PARTICULES EN FACE LES ZONES D'INJECTION)
!     (DEFAUT NON : 0 ; OUI : 1)

injcon = 0


!===============================================================================
! 7. TECHNIQUE DE REDUCTION DE VARIANCE : CLONAGE/FUSION DES PARTICULES
!===============================================================================

!     UTILISATION DE LA ROULETTE RUSSE
!                    DEFAUT NON : 0
!                           OUI : 1 sans calcul de Y+
!                                 2 avec calcul de Y+
iroule = 0


!===============================================================================
! 8. OPTIONS SUR LE TRAITEMENT NUMERIQUE DE LA PHASE DISPERSEE
!===============================================================================

!     ORDRE D'INTEGRATION DES EQUATIONS DIFFERENTIELLES STOCHASTIQUES
!     PAR DEFAUT 2 (VALEURS ADMISSIBLES : 1 OU 2)

!       NORDRE = 1 : INTEGRATION DES EDS PAR UN SCHEMA D'ORDRE 1
!       NORDRE = 2 : INTEGRATION DES EDS PAR UN SCHEMA D'ORDRE 2

nordre = 2


!     RESOLUTION DE L'EQUATION DE POISSON POUR LES VITESSE MOYENNES
!     DES PARTICULES ET CORRECTION DES VITESSES INSTANTANNEES
!     DES PARTICULES
!      = 0 : pas de correction des vitesses (VALEUR PAR DEFAUT)
!      = 1 : correction des vitesses instantanees

!     ATTENTION : OPTION STRICTEMENT DEVELOPPEUR, LAISSEZ LA VALEUR PAR
!     =========   DEFAUT           !

ilapoi = 0


!===============================================================================
! 9. OPTIONS SUR LE TRAITEMENT DE LA DISPERSION TURBULENTE
!===============================================================================

!     ATTENTION : DANS CETTE VERSION, LA DISPERSION TURBULENTE NE
!     ^^^^^^^^^^  FONCTIONNE QUE SI LA PHASE CONTINUE EST CALCULEE
!                 AVEC UN MODELE k-eps OU Rij-eps


!-->  PRISE EN COMPTE DE LA DISPERSION TURBULENTE
!     (DEFAUT OUI : 1 ; NON : 0)

idistu = 1


!-->  DISPERSION TURBULENTE IMPOSEE A CELLE DU FLUIDE

!     SI ACTIF, ALORS LA DISPERSION TURBULENTE DES PARTICULES EST CELLE
!     DES PARTICULES FLUIDES. ON SE PLACE DONC DANS UN CAS DE DIFFUSION
!     TURBULENTE (ON SUPRIME LES EFFETS DE CROISEMENT DE
!     TRAJECTOIRES). SI LES PARTICULES SIMULEES ONT LA MASSE
!     VOLUMIQUE DU FLUIDE ALORS, ON SIMULE LE DEPLACEMENT DE
!     PARTICULES FLUIDES LAGRANGIENNES.
!     (DEFAUT NON : 0 ; OUI : 1)

idiffl = 0

!     MODCPL :
!          = 0 pour le modele incomplet (VALEUR PAR DEFAUT)
!          > 0 pour le modele complet, est egal au numero de l'iteration
!          Lagrangienne absolue (i.e. suites comprises) a partir de
!          laquelle mise le modele complet est active
!          MODCPL ne doit pas etre inferieur a IDSTNT

modcpl = 0

!     IDIRLA  =1 ou 2 ou 3 : 1ere, 2eme ou 3eme direction
!       du modele complet. Correspond à la direction principale
!       de l'écoulement. Permet de calculer un echelle de temps
!       Lagrangienne non isotrope.
!       (DEFAUT IDIRLA = 1)

if (modcpl.gt.0) idirla = 1


!===============================================================================
! 10. OPTIONS SUR LE TRAITEMENT DES FORCES PARTICULIERES
!===============================================================================

!--> ACTIVATION DES FORCES :

!      - de Van der Waals
!      - Electrostatiques

!     (DEFAUT NON : 0 ; OUI : 1)

!     ATTENTION : OPTION DEVELOPPEUR UNIQUEMENT
!     =========

ladlvo = 0

!-->  Constante pour les forces de Van der Waals

!    Constante d'Hamaker Particule/Eau/substrat :

cstham = 6.d-20

!-->  Constante pour les forces Electrostatiques

!    Constante de Faradet (C/mol)

cstfar = 9.648d4

!    Constante dielectrique du vide :

epsvid = 8.854d-12

!    Constante dielectrique de l'eau

epseau = 80.10d0

!    Potentiel solide 1 (Volt)

phi1 = 50.d-3

!    Potentiel solide 2 (Volt)

phi2 = -50.d-3

!    FORCE IONIQUE (mol/l)

fion = 1.d-2

!    DISTANCE DE COUPURE
!      dans la litterature elle est egale à : 1.58D-10
!                                        ou   1.65D-10

dcoup = 1.58d-10

!    DENSITE DE CHARGE

sigch = 0.d0

!   Distance minimum entre la particule et la paroi

dparmn = 1.d-10


!===============================================================================
! 11. ACTIVATION DU MOUVEMENT BROWNIEN
!===============================================================================


!--> ACTIVATION DU MOUVEMENT BROWNIEN :

!     (DEFAUT NON : 0 ; OUI : 1)

!     ATTENTION : OPTION DEVELOPPEUR UNIQUEMENT
!     =========

lamvbr = 0


!===============================================================================
! 12. POST-PROCESSING
!===============================================================================

! 12.1 POST-PROCESSING DES TRAJECTOIRES ET DEPLACEMENTS PARTICULAIRES


!     ATTENTION : COUTEUX EN TEMPS DE CALCUL
!     ^^^^^^^^^^

!   12.1.1 PARAMETRES GENERAUX
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~

!     MODE TRAJECTOIRES
!     (DEFAUT NON : 0 ; OUI : 1)

iensi1 = 0

!     MODE DEPLACEMENTS PARTICULAIRES
!     (DEFAUT NON : 0 ; OUI : 1)

iensi2 = 0


!     NOMBRE DE PARTICULES A VISUALISER MAXIMUM=NLISTE
!     (PAR DEFAUT NBVIS = NLISTE)
!     ATTENTION : NBVIS NE DOIT NI ETRE SUPERIEUR A NBPMAX
!                 NI A NLISTE (parametre fixe a 500 dans lagpar.h)

nbvis = nliste

!     PERIODE   D'AQUISITION DES DONNEES A VISUALISER
!     (PAR DEFAUT NVISLA = 1)

nvisla = 1

!     LISTE contient les numeros des particules que l'on souhaite
!     visualiser
!     (PAR DEFAUT LISTE(...) = -1, AUCUNE PARTICULE A VISUALISER)
!     ATTENTION : si le numero est negatif, il n'y a pas de visualisation.

!   > EXEMPLE 1 :
!     je veux voir les NBVIS 1eres

do ii = 1, nbvis
  liste(ii) = ii
enddo

!   > EXEMPLE 2 :
!     je veux voir la 3 5 67 23 1 76 35 36 ...etc..

!     LISTE(1) = 3
!     LISTE(2) = 5
!     LISTE(3) = 67
!     LISTE(4) = 23
!     LISTE(5) = 1
!     LISTE(6) = 76
!     LISTE(7) = 35
!     LISTE(8) = 36
!     ...etc...

!   > REMARQUE : les trous, les repetitions dans le tableau LISTE seront
!                suprimes, et les numeros seront ranges par ordre
!                croissant.

!   12.1.2 VARIABLES A VISUALISER SUR LES TRAJECTOIRES OU LES PARTICULES
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!     VARIABLE : VITESSE DU FLUIDE VU
!     (DEFAUT NON : 0 ; OUI : 1)

ivisv1  = 0

!     VARIABLE : VITESSE DE LA PARTICULE
!     (DEFAUT NON : 0 ; OUI : 1)

ivisv2  = 0

!     VARIABLE : TEMPS DE SEJOUR
!     (DEFAUT NON : 0 ; OUI : 1)

ivistp  = 0

!     VARIABLE : DIAMETRE
!     (DEFAUT NON : 0 ; OUI : 1)

ivisdm  = 0

!     VARIABLE : TEMPERATURE
!     (DEFAUT NON : 0 ; OUI : 1)
if (iphyla.eq.1 .and. itpvar.eq.1) iviste  = 0

!     VARIABLE : MASSE
!     (DEFAUT NON : 0 ; OUI : 1)

ivismp  = 0


if (iphyla.eq.2) then

!       VARIABLE CHARBON : TEMPERATURE EN DEGRES CELSIUS
!       (DEFAUT NON : 0 ; OUI : 1)

  ivishp = 0

!       VARIABLE CHARBON : DIAMETRE DU COEUR RETRECISSANT
!       (DEFAUT NON : 0 ; OUI : 1)

  ivisdk  = 0

!       VARIABLE CHARBON : MASSE DE CHARBON REACTIF
!       (DEFAUT NON : 0 ; OUI : 1)

  ivisch  = 0

!       VARIABLE CHARBON : MASSE DE COKE
!       (DEFAUT NON : 0 ; OUI : 1)

  ivisck  = 0

endif


! 12.2 STATISTIQUES AUX FRONTIERES : VISUALISATION DES
!      INTERACTIONS PARTICULES/FRONTIERES
! ------------------------------------------------

!   12.2.1 PARAMETRES GENERAUX
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~

!     MODE INTERACTION PARTICULES/FRONTIERES
!     (DEFAUT NON : 0 ; OUI : 1)

iensi3 = 0

!     CALCUL STATIONNAIRE DES STATISTIQUES AUX FRONTIERES A PARTIR DE
!     L'ITERATION LAGRANGIENNE ABSOLUE NSTBOR
!     * NSTBOR est un nombre d'iterations Lagrangiennes absolues
!       (i.e. suites comprises) a partir duquel les statistiques sont
!       moyennees (en temps ou par nombre d'interactions).
!     * Utile si CALCUL STATIONNAIRE, i.e. si ISTTIO = 1.
!     * Si le nombre d'iterations Lagrangiennes absolues est strictement
!       inferieur a NSTBOR, les statistiques sont instationnaires
!       (i.e. elles sont REMISES A ZERO a chaque iteration
!       Lagrangienne).
!     * La valeur minimale admissible pour NSTBOR est 1.
!     (DEFAUT : NSTBOR = 1)

nstbor = 1

!     SEUILF POUR LA PRISE EN COMPTE DES STATISTIQUES AUX FRONTIERES
!     * La valeur de SEUILF est un poids statistique.
!     * Chaque face du maillage a subi un certain nombre d'interactions
!       avec des particules en terme de poids statistique
!       (somme des poids statistiques de toutes les particules
!       qui ont interagi avec la face de bord) ;
!       SEUILF est la valeur minimale a partir de laquelle la
!       contribution en terme statistique d'une face n'est
!       plus prise en compte dans les sorties listing
!       et post-processing.
!     (DEFAUT : SEUILF = 0.D0)

seuilf = 0.d0


!   12.2.2 INFORMATIONS A ENREGISTRER
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!    * Un certain nombre d'informations succeptibles d'interesser
!      l'utilisateur sont d'ores et deja ecrites dans le sous-programme
!      USLABO. Pour les activer il faut mettre le mot-cle correspondant
!      a 1 ci-dessous.
!    * La selection des types d'interaction (IREBOL, IDEPO1,... voir
!      le sous-programme USLAG2) qu'enclenchent l'enregistrement
!      des informations est faite dans USLABO. La selection par defaut
!      doit etre validee ou modifiee par l'utilisateur.
!    * Par defaut, dans un meme enregistrement, on place les informations
!      demandees pour toutes les interactions particules/frontieres
!      selectionnees. Modifier ce comportement necessite une
!      intervention dans le sous-programme utilisateur USLABO.
!    * La statistique aux frontieres :
!      "NOMBRE D'INTERACTIONS PARTICULES/FRONTIERES" doit etre
!      selectionnee obligatoirement pour utiliser la moyenne
!      particulaire IMOYBR(...) = 2.


!     NOMBRE D'INTERACTIONS PARTICULES/FRONTIERES
!     (DEFAUT NON : 0 ; OUI : 1)
inbrbd = 1

!     FLUX DE MASSE PARTICULAIRE LIE AUX INTERACTIONS
!     PARTICULES/FRONTIERES
!     (DEFAUT NON : 0 ; OUI : 1)
iflmbd = 1

!     ANGLE ENTRE LA VITESSE DE LA PARTICULE
!     ET LE PLAN DE LA FACE FRONTIERE
!     (DEFAUT NON : 0 ; OUI : 1)
iangbd = 0

!     NORME DE LA VITESSE DE LA PARTICULE AU MOMENT
!     DE SON INTERACTION AVEC LA FACE FRONTIERE
!     (DEFAUT NON : 0 ; OUI : 1)
ivitbd = 0

!     MASSE DE GRAINS DE CHARBON ENCRASSES
!     (DEFAUT NON : 0 ; OUI : 1)
 if (iphyla.eq.2 .and. iencra.eq.1) iencbd = 0

!     INFORMATIONS UTILISATEUR SUPPLEMENTAIRES A ENREGISTREES
!     (PAR EXEMPLE TAUX D'EROSION, TEMPERATURE...)
!     * ces enregistrements supplementaires sont stokees dans
!       le tableau PARBOR,
!     * on renseigne ici le nombre NUSBOR d'enregistrements
!       supplementaires,
!     * la limite supperieure de ce nombre est NUSBRD = 10
!       (fixee dans lagpar.h),
!     * voir un exemple de code d'un enregistrement supplementaire
!       dans le sous-programme USLABO.

nusbor = 0


!   12.2.3 NOM DES ENREGISTREMENTS POUR AFFICHAGE,
!          MOYENNE EN TEMPS OU MOYENNE PARTICULAIRE
!          DES STATISTIQUES AUX FRONTIERES
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!    * A priori l'utilisateur n'intervient que dans les infos
!      utilisateur supplementaires a enregistrer : il doit
!      donner le nom de l'enregistrement ainsi que le type
!      de moyenne que l'on souhaite lui appliquer pour les affichages
!      listing et en fin de calcul pour le post-processing.

!    * ATTENTION ce nom est utilise pour retrouver les informations
!      dans le fichier suite de calcul, si un nom est modifie
!      entre deux calculs, le resultat concerne n'est pas relu
!      (et est donc perdu).

!    * La moyenne appliquee est donnee par l'intermediaire du
!      tableau IMOYBR :
!       - si IMOYBR(IUSB(II)) = 0 -> pas de moyenne appliquee
!       - si IMOYBR(IUSB(II)) = 1 -> on applique une moyenne temporelle,
!               c'est a dire que la statistique est divisee par le
!               dernier pas de temps dans le cas d'un enregistrement
!               instationnaire, ou stationnaire mais avec un nombre
!               d'iterations inferieur a NSTBOR ; ou que la statistique
!               est divisee par le temps d'enregistrement
!               stationnaire dans le cas stationnaire.
!       - si IMOYBR(IUSB(II)) = 2 -> on applique une moyenne
!               particulaire, c'est a dire qu'on divise la statistique
!               par le nombre d'interactions particules/frontieres
!               enregistrees (en terme de poids statistique)
!               dans PARBOR(NFABOR,INBR) (cf. USLABO).
!               Pour utiliser cette moyenne il faut que INBRBD=1.
!    * Les sauvegardes dans le fichier suite de calcul sont faites
!      sans application de cette moyenne.
!    * La moyenne est appliquee si le nombre d'interactions (en poids
!      statistique) de la face de bord consideree est strictement
!      superieur a SEUILF, sinon la moyenne est mise a zero.

ipv = 0

if (iensi3.eq.1) then

  if (inbrbd.eq.1) then
    ipv = ipv + 1
    inbr = ipv
    NOMBRD(INBR) = 'nombreImpact'
    imoybr(inbr) = 0
  endif

  if (iflmbd.eq.1) then
    ipv = ipv + 1
    iflm = ipv
    NOMBRD(IFLM) = 'fluxDeMasse'
    imoybr(iflm) = 1
  endif

  if (iangbd.eq.1) then
    ipv = ipv + 1
    iang = ipv
    NOMBRD(IANG) = 'angleImpact'
    imoybr(iang) = 1
  endif

  if (ivitbd.eq.1) then
    ipv = ipv + 1
    ivit = ipv
    NOMBRD(IVIT) = 'normeVitImpact'
    imoybr(ivit) = 1
  endif

  if (iencbd.eq.1) then
    ipv = ipv + 1
    ienc = ipv
    NOMBRD(IENC) = 'masseEncras'
    imoybr(ienc) = 0
  endif

  if (nusbor.gt.0) then
    do ii = 1,nusbor
      ipv = ipv + 1
      iusb(ii) = ipv
      WRITE(NOMBRD(IUSB(II)),'(A8,I4.4)') 'enrSupp',II
      imoybr(iusb(ii)) = 0
    enddo
  endif

endif

!===============================================================================
! 13. LISTING LAGRANGIEN
!===============================================================================
! Periode de sortie du listing lagrangien

ntlal = 1

!===============================================================================

return

end
