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
C                              lagran.h
C***********************************************************************
C
C=======================================================================
C
C     Include pour le module Lagrangien
C
C         Trois fichiers complementaires
C                            lagran.h qui porte les non dimensions
C                            lagdim.h qui porte les dimensions variables
C                            lagpar.h qui porte les parametres
C
C=======================================================================
C  1. Base
C
C     IILAGR = 0 : PAS DE CALCUL LAGRANGIEN
C            = 1 : DIPHASIQUE LAGRANGIEN SANS COUPLAGE RETOUR
C            = 2 : DIPHASIQUE LAGRANGIEN AVEC COUPLAGE RETOUR
C
C     ISUILA = 0 : PAS SUITE LAGRANGIEN
C            = 1 :     SUITE LAGRANGIEN
C
C     ISTTIO = 0 : calcul instationnaire pour le lagrangien
C            = 1 : calcul stationnaire   pour le lagrangien
C
C     ILPHAS     : Numero de la phase continue sur laquelle le
C                  module prend ses infos et vers laquelle il
C                  renvoie les termes sources de couplage retour.
C
      INTEGER           IILAGR , ISUILA , ISTTIO , ILPHAS
      COMMON / ILAG11 / IILAGR , ISUILA , ISTTIO , ILPHAS
C
C=======================================================================
C 2. Compteurs de particules (sans et avec poids statistique)
C
C     NBPART/DNBPAR : NOMBRE DE PARTICULES PRESENTES DANS LE DOMAINE
C                        A CHAQUE ITERATION
C
C     NBPNEW/DNBPNW : NOMBRE DE NOUVELLES PARTICULES ENTRANTES
C
C     NBPERR/DNBPER : NOMBRE DE PARTICULES ELIMINEES EN ERREUR
C
C     NBPERT : NOMBRE DE PARTICULES ELIMINEES EN ERREUR DANS
C                LE CALCUL DEPUIS LE DEBUT, SUITES COMPRISES
C
C     NBPTOT : NOMBRE DE PARTICULES TOTAL INJECTE DANS
C                LE CALCUL DEPUIS LE DEBUT SUITE COMPRISE
C
C     NBPOUT/DNBPOU : Contient les particules sorties de facon normal,
C                       plus les particules sorties en erreur de reperage.
C
C     NDEPOT : Nombre de particules deposees definitivement
C               dont on garde une trace en memoire pour le
C               post-processing en mode deplacement.
C
C     NPCLON/DNPCLO : NOMBRE DE NOUVELLES PARTICULES PAR CLONNAGE
C
C     NPKILL/DNPKIL : NOMBRE DE PARTICULES VICTIMES DE LA ROULETTE RUSSE
C
C     NPCSUP/DNPCSU : NOMBRE DE PARTICULES QUI ON SUBIT LE CLONNAGE
C
C
      INTEGER           NBPART , NBPNEW , NBPERR , NBPTOT , NBPOUT ,
     &                  NBPERT , NDEPOT
      COMMON / ILAG21 / NBPART , NBPNEW , NBPERR , NBPTOT , NBPOUT ,
     &                  NBPERT , NDEPOT
C
      DOUBLE PRECISION  DNBPAR , DNBPNW , DNBPER ,          DNBPOU
      COMMON / RLAG21 / DNBPAR , DNBPNW , DNBPER ,          DNBPOU
C
      INTEGER           NPCLON , NPKILL , NPCSUP
      COMMON / ILAG22 / NPCLON , NPKILL , NPCSUP
C
      DOUBLE PRECISION  DNPCLO , DNPKIL , DNPCSU
      COMMON / RLAG22 / DNPCLO , DNPKIL , DNPCSU
C
C=======================================================================
C 4. Physiques particulieres
C
C       SI IPHYLA = 1 ALORS
C
C          ITPVAR : EQUATION SUR LA TEMPERATURE
C          IDPVAR : EQUATION SUR LE DIAMETRE
C          IMPVAR : EQUATION SUR LA MASSE
C          ILUMN  : POINTEUR QUI PERMET DE REPERER L INTEGRALE DE LA
C                   LUMINANCE DANS LA TABLEAU PROPCE
C
C
      INTEGER           IPHYLA, ITPVAR, IDPVAR, IMPVAR, ILUMN
      COMMON / ILAG41 / IPHYLA, ITPVAR, IDPVAR, IMPVAR, ILUMN
C
C
C       SI SUITE ET ENCLENCHEMENT ITPVAR =1 EN COUR DE CALCUL
C
C          TPART  : Temperature d initialisation en degres Celsius
C          CPPART : Chaleur massique specifique (J/kg/K)
C
      DOUBLE PRECISION  TPART , CPPART
      COMMON / RLAG41 / TPART , CPPART
C
C
C=======================================================================
C 5. Pas de temps Lagrangien
C
C    IPLAS : NOMBRE DE PASSAGES ABSOLUS DANS LE MODULE LAGRANGIEN
C    IPLAR : NOMBRE DE PASSAGES RELATIFS DANS LE MODULE LAGRANGIEN
C
      INTEGER           IPLAS , IPLAR
      COMMON / ILAG51 / IPLAS , IPLAR
C
C    DTP :  duree d une iteration lagrangienne
C    TTCLAG : temps courant physique lagrangien
C
      DOUBLE PRECISION  DTP , TTCLAG
      COMMON / RLAG51 / DTP , TTCLAG
C
C=======================================================================
C 6. Indicateur d erreur
C
      INTEGER           IERR
      COMMON / ILAG61 / IERR
C
C=======================================================================
C 3. Pointeurs particules
C
C   Tableau ETTP
C   ^^^^^^^^^^^^
C
C    JXP,JYP,JZP  : COORDONNES DE LA POSITION DE LA PARTICULE
C    JUP,JVP,JWP  : COMPOSANTES DE LA VITESSE ABSOLUE
C    JUF,JVF,JWF  : COMPOSANTES DE LA VITESSE DU FLUIDE VU
C
C    JMP,JDP      : MASSE, DIAMETRE
C    JTP,JTF,JCP  : TEMPERATURE PARTICULE ET FLUIDE ET CHALEUR SPECIFIQUE
C    JVLS(NUSVAR) : VARIABLE SUPPLEMENTAIRES
C
C   Charbon
C   -------
C    JHP          : TEMPERATURE DES GRAINS DE CHARBON
C    JMCH         : MASSE DE CHARBON REACTIF
C    JMCK         : MASSE DE COKE
C
C
      INTEGER           JXP , JYP , JZP ,
     &                  JUP , JVP , JWP ,
     &                  JUF , JVF , JWF ,
     &                  JMP , JDP , JTP , JTF , JCP ,
     &                  JHP , JMCH, JMCK,
     &                  JVLS(NUSVAR)
C
      COMMON / ILAG31 / JXP , JYP , JZP ,
     &                  JUP , JVP , JWP ,
     &                  JUF , JVF , JWF ,
     &                  JMP , JDP , JTP , JTF , JCP ,
     &                  JHP , JMCH, JMCK,
     &                  JVLS
C
C   Tableau TEPA
C   ^^^^^^^^^^^^
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
      INTEGER           JRTSP, JRPOI, JREPS, JRD0P, JRR0P, JRDCK
      COMMON / ILAG32 / JRTSP, JRPOI, JREPS, JRD0P, JRR0P, JRDCK
C
C   Tableau ITEPA
C   ^^^^^^^^^^^^^
C
C     JISOR       : MAILLE D ARRIVEE
C
C   Statistique par classe
C   ----------------------
C
C     JCLST       : classe (statique ) a laquelle la particule appartient
C
C   Charbon
C   -------
C     JINCH       : NUMERO DU CHARBON DE LA PARTICULE
C
      INTEGER           JISOR, JINCH , JCLST
      COMMON / ILAG33 / JISOR, JINCH , JCLST
C
C    NVLS         : NOMBRE DE VARIABLES UTILISATEUR SUPPLEMENTAIRES
C                   (DEJA CONTENU DANS NVP et NVP1)
C
      INTEGER           NVLS
      COMMON / ILAG34 / NVLS
C
C=======================================================================
C 7. Conditions aux limites
C
C     TABLEAUX POUR LES CONDITIONS AUX LIMITES
C     ----------------------------------------
C
C     NFRLAG  : nbr de zones frontieres
C     INJCON  : INJECTION CONTINUE OU NON
C     ILFLAG  : liste des numeros des zones frontieres
C     IUSNCL  : nbr de classes par zones
C     IUSCLB  : conditions au bord pour les particules
C          = IENTRL
C          = ISORTL -> particule sortie du domaine par une face fluide
C          = IREBOL -> rebond elastique
C          = IDEPO1 -> deposition definitive (particule eliminee de la memoire)
C          = IDEPO2 -> deposition definitive (part. non eliminee de la memoire)
C          = IDEPO3 -> deposition temporaire (remise en suspension possible)
C          = IENCRL -> encrassement (Charbon uniquement IPHYLA = 2)
C          = JBORD1 -> interactions utilisateur
C          = JBORD2 -> interactions utilisateur
C          = JBORD3 -> interactions utilisateur
C          = JBORD4 -> interactions utilisateur
C          = JBORD5 -> interactions utilisateur
C     IUSMOY  : tableau si on fait une moyenne par zone sur la zone considere
C     IUSLAG  : tableau d info par classe et par frontieres
C     DEBLAG  : debit massique par zone
C
C
      INTEGER           NFRLAG, INJCON,
     &                  ILFLAG(NFLAGM),
     &                  IUSNCL(NFLAGM),
     &                  IUSCLB(NFLAGM),
     &                  IUSMOY(NFLAGM),
     &                  IUSLAG(NCLAGM, NFLAGM, NDLAIM)
C
      COMMON / ILAG71 / NFRLAG, INJCON,
     &                  ILFLAG, IUSNCL,
     &                  IUSCLB,
     &                  IUSMOY, IUSLAG
C
      DOUBLE PRECISION  DEBLAG(NFLAGM)
      COMMON / RLAG71 / DEBLAG
C
C
C
C     IJNBP  : nbr de part par classe et zones frontieres
C     IJFRE  : frequence d injection
C               (si < 0 : on ne rentre des particles qu a la 1ere iter)
C     IJUVW  : type de condition vitesse
C          = -1 vitesse fluide imposee
C          =  0 vitesse imposee selon la direction normale
C               a la face de bord et de norme IUNO
C          =  1 vitesse imposee : on donne IUPT IVPT IWPT
C          =  2 profil de vitesse donne par l'utilisateur
C     IJPRPD = 1 distribution uniforme
C            = 2 profil de taux de presence donne par l'utilisateur
C     IJPRTP = 1 profil plat de temperature donne par la valeur dans uslag2
C            = 2 profil de temperature donne par l'utilisateur
C     IJPRDP = 1 profil plat de diametre donne par la valeur dans uslag2
C            = 2 profil dediametre donne par l'utilisateur
C     INUCHL : numero du charbon de la particule (si IPHYLA=2)
C     ICLST  : numero du groupe de statistiques
C
C
      INTEGER           IJNBP, IJFRE, IJUVW, IJPRTP, IJPRDP, IJPRPD
      INTEGER           INUCHL, ICLST
      COMMON / ILAG72 / IJNBP, IJFRE, IJUVW, IJPRTP, IJPRDP, IJPRPD,
     &                  INUCHL, ICLST
C
C
C     RUSLAG  : tableau d info par classe et par frontieres
C
C
      DOUBLE PRECISION  RUSLAG(NCLAGM, NFLAGM, NDLAGM)
      COMMON / RLAG73 / RUSLAG
C
C
C     IUNO  : Norme de la vitesse
C     IUPT  : U par classe et zones
C     IVPT  : V par classe et zones
C     IWPT  : W par classe et zones
C     IDEBT : Debit
C     IPOIT : Poids de la particule
C     IDPT  : Diametre
C     IVDPT : Variance du diametre
C     ITPT  : Temperature
C     ICPT  : Cp
C     IEPSI : Emissivite des particules
C     IROPT : Masse volumique
C     IHPT  : Temperature
C     IMCHT : Masse de charbon reactif
C     IMCKT : Masse de coke
C     IDCKT : Diametre du coeur retrecissant
C
C
      INTEGER           IUNO, IUPT, IVPT, IWPT,
     &                  ITPT, IDPT, IVDPT, IROPT,
     &                  ICPT, IPOIT, IDEBT, IEPSI,
     &                  IHPT, IMCHT, IMCKT, IDCKT
C
      COMMON / ILAG74 / IUNO, IUPT, IVPT, IWPT,
     &                  ITPT, IDPT, IVDPT, IROPT,
     &                  ICPT, IPOIT, IDEBT, IEPSI,
     &                  IHPT, IMCHT, IMCKT, IDCKT
C
C=======================================================================
C 8. Statistiques
C
C     POINTEURS POUR LES STATISTIQUES
C     -------------------------------
C
C     ILVX,ILVY,ILVZ    : Vitesse
C     ILFV              : Concentration volumique
C     ILPD              : Somme des poids statistiques
C     ILTS              : Temps de sejour
C
C     ILTP              : Temperature
C     ILDP              : Diametre
C     ILMP              : Masse
C
C     ILHP              : Temperature
C     ILMCH             : Masse de charbon reactif
C     ILMCK             : Masse de coke
C     ILDCK             : Diametre du coeur retrecissant
C
C     ILVU(NUSSTA)      : Statistiques supplementaires utilisateur
C
      INTEGER           ILVX  , ILVY  , ILVZ  ,
     &                  ILPD  , ILFV  , ILTS  ,
     &                  ILTP  , ILDP  , ILMP  ,
     &                  ILHP  , ILMCH , ILMCK , ILDCK ,
     &                  ILVU(NUSSTA)
C
      COMMON / ILAGR7 / ILVX  , ILVY  , ILVZ  ,
     &                  ILPD  , ILFV  , ILTS  ,
     &                  ILTP  , ILDP  , ILMP  ,
     &                  ILHP  , ILMCH , ILMCK , ILDCK ,
     &                  ILVU
C
C
C     DONNEES POUR LES STATISTIQUES VOLUMIQUES
C     ----------------------------------------
C
C      ISTALA : Calcul statistiques       si  >= 1 sinon pas de stat
C      ISUIST : Suite calcul statistiques si  >= 1 sinon pas de stat
C      NVLSTS : NOMBRE DE VARIABLES STATISTIQUES SUPPLEMENTAIRES
C               UTILISATEUR (CONTENU DANS NVLSTA)
C      IDSTNT : Numero du pas de temps pour debut statistque
C      NSTIST : Debut calcul stationnaire
C      NPST   : Nombre de pas de temps pour le cumul des stats
C      NPSTT  : Nombre de pas de temps total des stats depuis le debut
C               du calcul, partie instationnaire comprise
C      TSTAT  : Temps physique des stats volumiques
C      SEUIL  : Seuil en POIDS STAT de particules pour les stats
C
C
      INTEGER           ISTALA , ISUIST , NVLSTS ,
     &                  IDSTNT , NSTIST ,
     &                  NPST   , NPSTT
C
      COMMON / ILASTA / ISTALA , ISUIST , NVLSTS ,
     &                  IDSTNT , NSTIST ,
     &                  NPST   , NPSTT
C
      DOUBLE PRECISION   TSTAT , SEUIL
      COMMON / RLASTA /  TSTAT , SEUIL
C
C
C     NOMS DES VARIABLES STATISTIQUES (MOYENNES ET VARIANCES)
C     -------------------------------------------------------
C     Taille limitee par le fait qu on utilise NOMBRD dans
C       l ecriture des fichiers suites (lagout)
C
      CHARACTER*50      NOMLAG(NVPLMX) , NOMLAV(NVPLMX)
      COMMON / ALASTS / NOMLAG         , NOMLAV
C
C     OPTION POUR LES HISTORIQUES SUR LES STATS
C     -----------------------------------------
C
      INTEGER           IHSLAG(NVPLMX)
      COMMON / ILOPHL / IHSLAG
C
C     STATISTIQUE PAR ZONE ET PAR CLASSE
C     ----------------------------------
C
      INTEGER           NBCLST
      COMMON / ILSTCL / NBCLST
C
C=======================================================================
C 9. Termes Sources
C
C     OPTION TERMES SOURCES
C     ---------------------
C       Dynamique
C       Masse
C       Thermique
C
      INTEGER          LTSDYN , LTSMAS , LTSTHE
      COMMON /ILOPTS / LTSDYN , LTSMAS , LTSTHE
C
C     POINTEURS POUR LES TERMES SOURCES
C     ---------------------------------
C
C    ITSVX,ITSVY,ITVZ    : Termes sources sur la vitesse
C    ITSLI               : Terme source implicite (vitesse+turbulence)
C    ITSKE               : Terme source sur la turbulence en k-eps
C    ITSR11,ITR12,ITSR13 : Termes sources sur la turbulence en Rij-Eps
C    ITSR22,ITR23,ITSR33
C    ITSTE, ITSTI        : Termes sources pour la thermique
C    ITSMAS              : Terme source pour la masse
C    ITSMV1              : Terme source sur F1 (MV legeres)
C    ITSMV2              : Terme source sur F2 (MV loudres)
C    ITSCO               : Terme source sur F3 (C sous forme de CO)
C    ITSFP4              : Variance du traceur relatif a l air
C
C
      INTEGER           ITSVX  , ITSVY  , ITSVZ  , ITSLI ,
     &                  ITSKE  ,
     &                  ITSR11 , ITSR12 , ITSR13 ,
     &                  ITSR22 , ITSR23 , ITSR33 ,
     &                  ITSTE  , ITSTI  ,
     &                  ITSMAS , ITSMV1(NCHARM2), ITSMV2(NCHARM2) ,
     &                  ITSCO  , ITSFP4
C
      COMMON / ILAG91 / ITSVX  , ITSVY , ITSVZ, ITSLI ,
     &                  ITSKE  ,
     &                  ITSR11 , ITSR12 , ITSR13 ,
     &                  ITSR22 , ITSR23 , ITSR33 ,
     &                  ITSTE  , ITSTI  ,
     &                  ITSMAS , ITSMV1 , ITSMV2 ,
     &                  ITSCO  , ITSFP4
C
C
C     DONNEES POUR LES TERMES SOURCES
C     -------------------------------
C
C     NSTITS : debut calcul terme source stationnaire
C     NPTS   : nombre de pas de temps pour le cumul des termes sources
C     NTXERR : nombre de cellules qui un taux vol > 0.8
C     VMAX   : taux volumique max atteint
C     TMAMAX : taux massique max atteint
C
      INTEGER           NSTITS , NPTS , NTXERR
      COMMON / ILAG92 / NSTITS , NPTS , NTXERR
C
      DOUBLE PRECISION  VMAX , TMAMAX
      COMMON / ILAG93 / VMAX , TMAMAX
C
C=======================================================================
C 10. Clonage/fusion des particules
C
C     INDICATEUR D ACTIVATION DE LA ROULETTE RUSSE
C
C
      INTEGER           IROULE
      COMMON / ILA101 / IROULE
C
C=======================================================================
C 11. Encrassement
C
C     DONNEES POUR L ENCRASSEMENT
C
      INTEGER           IENCRA , NPENCR
      COMMON / ILA111 / IENCRA , NPENCR
C
C
      DOUBLE PRECISION  ENC1(NCHARM2) , ENC2(NCHARM2) ,
     &                  TPRENC(NCHARM2) , VISREF(NCHARM2) , DNPENC
      COMMON / RLA112 / ENC1 , ENC2 , TPRENC , VISREF , DNPENC
C
C
C=======================================================================
C 12. Forces chimiques
C
C       1) FORCES DE VAN DER WAALS
C       2) FORCES ELECTROSTATIQUES
C
      INTEGER           LADLVO
      COMMON / ILADLV / LADLVO
C
C      CSTHAM : constante d'Hamaker
C      CSTFAR : constant de FARADET
C      EPSEAU : Constante dielectrique de l'eau
C      EPSEAU : Constante dielectrique du vide
C      PHI1   : potentiel solide 1
C      PHI1   : potentiel solide 2
C      FION   : force ionique
C      GAMASV : energie de surface
C      DPARMN : distance entre particule/paroi minimum
C
      DOUBLE PRECISION  CSTHAM , EPSEAU  , EPSVID , PHI1 , PHI2
      DOUBLE PRECISION  FION   , GAMASV  , DCOUP  , SIGCH
      DOUBLE PRECISION  CSTFAR , DPARMN
      COMMON / RLADLV / CSTHAM , EPSEAU  , EPSVID , PHI1 , PHI2 ,
     &                  FION   , GAMASV  , DCOUP  , SIGCH,
     &                  CSTFAR , DPARMN
C
C
C=======================================================================
C 13. Mouvement brownien
C
C     ACTIVATION DU MOUVEMENT BROWNIEN :
C
      INTEGER           LAMVBR
      COMMON / ILAMBR / LAMVBR
C
      DOUBLE PRECISION KBOLTZ
      PARAMETER          (KBOLTZ = 1.38D-23)
C
C
C=======================================================================
C 14. Schema en temps, dispersion turbulente et equation de poisson
C
C     NOR    : numero du sous-pas Lagrangien (1 ou 2)
C
C     NORDRE : ordre de la methode d integration (1 ou 2)
C
C     MODCPL : = 0 pour le modele incomplet
C              > 0 pour le modele complet, est egal au nombre de
C                 passages avant mise en route du modele complet
C
C     IDIRLA : = 1 ou 2 ou 3 direction du modele complet
C
C     IDISTU : = 0 pas de prise en compte de la dispersion turbulente (la
C                  vitesse instantanee est egale a la vitesse moyenne)
C              > 0 prise en compte de la dispersion turbulente (si k-eps
C                  ou Rij-eps)
C
C     IDIFFL : =1 la dispersion turbulente de la particule est celle de
C                 la particule fluide (=0 sinon)
C
C     ILAPOI : = 0 Pas de correction de pression
C              = 1 Correction de pression
C
C
      INTEGER           NOR , NORDRE , MODCPL , IDIRLA ,
     &                  IDISTU , IDIFFL , ILAPOI
C
      COMMON / ILA121 / NOR , NORDRE , MODCPL , IDIRLA ,
     &                  IDISTU , IDIFFL , ILAPOI
C
C
C=======================================================================
C 15. Traitement des statistiques interactions particules/frontieres
C
C     DONNEES POUR LES STATISTIQUES AUX FRONTIERES
C     --------------------------------------------
C
C      NUSBOR : NOMBRE DE VARIABLES A ENREGISTRER SUR LES FRONTIERES
C               SUPPLEMENTAIRES UTILISATEUR (CONTENU DANS NVISBR)
C      NSTBOR : debut calcul stationnaire
C      NPSTF  : nombre de pas de temps pour le cumul des stats
C      NPSTF  : nombre de pas de temps total des stats depuis le debut
C               du calcul, partie instationnaire comprise
C      TSTATP : Temps physique des stats aux frontieres stationnaires
C      SEUILF : Seuil en POIDS STAT de particules pour les stats
C      IMOYBR : Type de moyenne applicable pour affichage et
C               post-procesing
C
      INTEGER           NUSBOR , NSTBOR ,
     &                  NPSTF  , NPSTFT ,
     &                  INBRBD , IFLMBD , IANGBD , IVITBD , IENCBD ,
     &                  INBR   , IFLM   , IANG   , IVIT   , IENC   ,
     &                  IUSB(NUSBRD)    , IMOYBR(NUSBRD+10)
C
      COMMON / LAGBRD / NUSBOR , NSTBOR ,
     &                  NPSTF  , NPSTFT ,
     &                  INBRBD , IFLMBD , IANGBD , IVITBD , IENCBD ,
     &                  INBR   , IFLM   , IANG   , IVIT   , IENC   ,
     &                  IUSB   , IMOYBR
C
C
      DOUBLE PRECISION  TSTATP , SEUILF
      COMMON / RLABRD / TSTATP , SEUILF
C
C     NOMS DES VARIABLES STATISTIQUES
C     -------------------------------
C     Taille limitee par le fait qu on utilise NOMBRD dans
C       l ecriture des fichiers suites (lagout)
C
      CHARACTER*50      NOMBRD(NVPLMX)
      COMMON / ALABRD / NOMBRD
C
C IIFRLA Pointeur dans IA sur IFRLAG pour reperage des zones
C          frontieres associees aux faces de bord
C
      INTEGER           IIFRLA
      COMMON / IRLORD / IIFRLA
C
C=======================================================================
C 16. Visu
C
C... NBVIS  : nombre de particules a visualiser a l instant t
C    LISTE  : numero des particules a visualiser
C    LIST0  : sauvegarde de LISTE pour post-processing trajectoires
C    NPLIST : nombre d enregistrement par particule
C    NVISLA : periode d aquisition
C
      INTEGER           NBVIS, LISTE(NLISTE), LIST0(NLISTE),
     &                  NPLIST(NLISTE), NVISLA
      COMMON / IENLA1 / NBVIS, LISTE, LIST0, NPLIST,  NVISLA
C
C... Type de visualisation :
C    IENSI1 : trajectoires
C    IENSI2 : deplacements
C    IENSI3 : interaction particules/frontieres
C
C
      INTEGER           IENSI1 , IENSI2 , IENSI3
      COMMON / IENLA3 / IENSI1 , IENSI2 , IENSI3
C
C
C... Contenu des flichiers resultats
C
C
      INTEGER           IVISV1 , IVISV2 , IVISTP ,
     &                  IVISDM , IVISTE , IVISMP ,
     &                  IVISHP , IVISCH , IVISCK , IVISDK
      COMMON / IENLA5 / IVISV1 , IVISV2 , IVISTP ,
     &                  IVISDM , IVISTE , IVISMP ,
     &                  IVISHP , IVISCH , IVISCK , IVISDK
C
C
C... visualisation de type deplacement
C    ITLAG : nombre d enregistrement
C    TIMLAG : temps physiques lagrangien pour la visualisation
C
C
      INTEGER           ITLAG
      COMMON / IENLA6 / ITLAG
C
      DOUBLE PRECISION  TIMLAG(9999)
      COMMON / RENLA7 / TIMLAG
C
C=======================================================================
C
C FIN
c@z

