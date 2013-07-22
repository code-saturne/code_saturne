!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!> \file lagran.f90
!> Module for Lagrangian model.

module lagran

  !===========================================================================

  !         Trois modules complementaires
  !                            lagran qui porte les non dimensions
  !                            lagdim qui porte les dimensions variables
  !                            lagpar qui porte les parametres

  use lagpar

  !=============================================================================
  !  1. Base

  !     IILAGR = 0 : PAS DE CALCUL LAGRANGIEN
  !            = 1 : DIPHASIQUE LAGRANGIEN SANS COUPLAGE RETOUR
  !            = 2 : DIPHASIQUE LAGRANGIEN AVEC COUPLAGE RETOUR

  !     ISUILA = 0 : PAS SUITE LAGRANGIEN
  !            = 1 :     SUITE LAGRANGIEN

  !     ISTTIO = 0 : calcul instationnaire pour le lagrangien
  !            = 1 : calcul stationnaire   pour le lagrangien

  integer, save ::           iilagr , isuila , isttio

  !=============================================================================
  ! 2. Compteurs de particules (sans et avec poids statistique)

  !     NBPART/DNBPAR : NOMBRE DE PARTICULES PRESENTES DANS LE DOMAINE
  !                        A CHAQUE ITERATION

  !     NBPNEW/DNBPNW : NOMBRE DE NOUVELLES PARTICULES ENTRANTES

  !     NBPERR/DNBPER : NOMBRE DE PARTICULES ELIMINEES EN ERREUR

  !     NBPDEP/DNBDEP : NOMBRE DE PARTICULES DEPOSEES

  !     NBPRES/DNBRES : NOMBRE DE PARTICULES RESUSPENDUES

  !     NBPERT : NOMBRE DE PARTICULES ELIMINEES EN ERREUR DANS
  !                LE CALCUL DEPUIS LE DEBUT, SUITES COMPRISES

  !     NBPTOT : NOMBRE DE PARTICULES TOTAL INJECTE DANS
  !                LE CALCUL DEPUIS LE DEBUT SUITE COMPRISE

  !     NBPOUT/DNBPOU : Contient les particules sorties de facon normal,
  !                       plus les particules sorties en erreur de reperage.

  !     NDEPOT : Nombre de particules deposees definitivement
  !               dont on garde une trace en memoire pour le
  !               post-processing en mode deplacement.

  !     NPCLON/DNPCLO : NOMBRE DE NOUVELLES PARTICULES PAR CLONNAGE

  !     NPKILL/DNPKIL : NOMBRE DE PARTICULES VICTIMES DE LA ROULETTE RUSSE

  !     NPCSUP/DNPCSU : NOMBRE DE PARTICULES QUI ON SUBIT LE CLONNAGE


  integer, save ::           nbpart , nbpnew , nbperr , nbptot , nbpout ,   &
                             nbpert , ndepot , nbpdep , nbpres
  double precision, save ::  dnbpar , dnbpnw , dnbper , dnbpou, dnbdep,     &
                             dnbres

  integer, save ::           npclon , npkill , npcsup
  double precision, save ::  dnpclo , dnpkil , dnpcsu

  !=============================================================================
  ! 4. Physiques particulieres

  !       SI IPHYLA = 1 ALORS

  !          ITPVAR : EQUATION SUR LA TEMPERATURE
  !          IDPVAR : EQUATION SUR LE DIAMETRE
  !          IMPVAR : EQUATION SUR LA MASSE

  integer, save ::           iphyla, itpvar, idpvar, impvar

  !       SI SUITE ET ENCLENCHEMENT ITPVAR =1 EN COUR DE CALCUL

  !          TPART  : Temperature d initialisation en degres Celsius
  !          CPPART : Chaleur massique specifique (J/kg/K)

  double precision, save ::  tpart , cppart




  ! 4.1 Particules deposition submodel
  !==================================


  !     IDEPST = 0 : no deposition submodel activated
  !            = 1 : deposition submodel used

  integer, save ::     idepst

  !     IDLVO  = 0 : no DLVO conditions
  !            = 1 : DLVO conditions

  integer, save ::         idlvo

  !     NGEOL : geometric parameters stored

  integer ngeol
  parameter (ngeol = 13)

  ! Additional pointers in the ITEPA array
  ! ITEPA contains the particule state

  integer, save ::   jimark , jdiel  , jdfac , jdifel , jtraj , jptdet , jinjst
  integer, save ::   jryplu , jrinpf

  ! 4.2 Resuspension model
  !=======================

  !     IREENT = 0 : no resuspension model
  !            = 1 : resuspension model

  integer, save ::         ireent

  ! Additional pointers in the ITEPA and TEPA array
  ! ITEPA contains the particule state

  integer, save ::   jroll , jnbasg , jnbasp , jdepo
  integer, save ::   jfadh , jmfadh , jndisp

  ! Parameters of the particle resuspension model

   double precision, save :: espasg , denasp, modyeq , rayasp, rayasg

  !=============================================================================
  ! 5. Pas de temps Lagrangien

  !    IPLAS : NOMBRE DE PASSAGES ABSOLUS DANS LE MODULE LAGRANGIEN
  !    IPLAR : NOMBRE DE PASSAGES RELATIFS DANS LE MODULE LAGRANGIEN

  integer, save ::           iplas , iplar

  !    DTP :  duree d une iteration lagrangienne
  !    TTCLAG : temps courant physique lagrangien

  double precision, save ::  dtp , ttclag

  !=============================================================================
  ! 6. Indicateur d erreur

  integer, save :: ierr

  !=============================================================================
  ! 3. Pointeurs particules

  !   Tableau ETTP
  !   ^^^^^^^^^^^^

  !    JXP,JYP,JZP  : COORDONNES DE LA POSITION DE LA PARTICULE
  !    JUP,JVP,JWP  : COMPOSANTES DE LA VITESSE ABSOLUE
  !    JUF,JVF,JWF  : COMPOSANTES DE LA VITESSE DU FLUIDE VU

  !    JMP,JDP      : MASSE, DIAMETRE
  !    JTAUX        : POUR L'ORDRE 2, AUXILIAIRE DE CALCUL UTILE EN PARALLELE
  !    JTP,JTF,JCP  : TEMPERATURE PARTICULE ET FLUIDE ET CHALEUR SPECIFIQUE
  !    JVLS(NUSVAR) : VARIABLE SUPPLEMENTAIRES

  !   Charbon
  !   -------
  !    JHP          : TEMPERATURE DES GRAINS DE CHARBON
  !    JMWAT        : MASSE D EAU
  !    JMCH         : MASSE DE CHARBON REACTIF
  !    JMCK         : MASSE DE COKE

  integer, save ::           jxp , jyp , jzp ,                               &
                             jup , jvp , jwp ,                               &
                             juf , jvf , jwf ,                               &
                             jmp , jdp , jtp , jtf , jcp ,                   &
                             jhp , jmwat, jmch, jmck, jtaux,                 &
                             jvls(nusvar)

  !   Tableau TEPA
  !   ^^^^^^^^^^^^

  !     jrval       : random number associated with a particle
  !     jrtsp       : temps de sejour des particules
  !     jrpoi       : poids des particules
  !     jreps       : emissivite des particules

  !   Charbon
  !   -------
  !     jrdck       : diametre du coeur retrecissant
  !     jrd0p       : diametre initial des particules
  !     jrr0p       : masse volumique initiale des particules
  !     jrhock      : masse volumique du coke


  integer, save :: jrval, jrtsp, jrpoi, jreps, jrd0p, jrr0p, jrdck, jrhock

  !   Tableau ITEPA
  !   ^^^^^^^^^^^^^

  !     JISOR       : MAILLE D ARRIVEE

  !   Statistique par classe
  !   ----------------------

  !     JCLST       : classe (statique ) a laquelle la particule appartient

  !   Charbon
  !   -------
  !     JINCH       : NUMERO DU CHARBON DE LA PARTICULE

  integer, save ::           jisor, jinch , jclst

  !    NVLS         : NOMBRE DE VARIABLES UTILISATEUR SUPPLEMENTAIRES
  !                   (DEJA CONTENU DANS NVP et NVP1)

  integer, save ::           nvls

  !=============================================================================
  ! 7. Conditions aux limites

  !     TABLEAUX POUR LES CONDITIONS AUX LIMITES
  !     ----------------------------------------

  !     NFRLAG  : nbr de zones frontieres
  !     INJCON  : INJECTION CONTINUE OU NON
  !     ILFLAG  : liste des numeros des zones frontieres
  !     IUSNCL  : nbr de classes par zones
  !     IUSCLB  : conditions au bord pour les particules
  !          = IENTRL
  !          = ISORTL -> particule sortie du domaine par une face fluide
  !          = IREBOL -> rebond elastique
  !          = IDEPO1 -> deposition definitive (particule eliminee de la memoire)
  !          = IDEPO2 -> deposition definitive (part. non eliminee de la memoire)
  !          = IENCRL -> encrassement (Charbon uniquement IPHYLA = 2)
  !          = JBORD1 -> interactions utilisateur
  !          = JBORD2 -> interactions utilisateur
  !          = JBORD3 -> interactions utilisateur
  !          = JBORD4 -> interactions utilisateur
  !          = JBORD5 -> interactions utilisateur
  !     IUSMOY  : tableau si on fait une moyenne par zone sur la zone considere
  !     IUSLAG  : tableau d info par classe et par frontieres
  !     DEBLAG  : debit massique par zone

  integer, save ::           nfrlag, injcon,                                 &
                             ilflag(nflagm),                                 &
                             iusncl(nflagm),                                 &
                             iusclb(nflagm),                                 &
                             iusmoy(nflagm)

  integer, allocatable, dimension(:,:,:) :: iuslag

  double precision, save ::  deblag(nflagm)

  !     IJNBP  : nbr de part par classe et zones frontieres
  !     IJFRE  : frequence d injection
  !               (si < 0 : on ne rentre des particles qu a la 1ere iter)
  !     IJUVW  : type de condition vitesse
  !          = -1 vitesse fluide imposee
  !          =  0 vitesse imposee selon la direction normale
  !               a la face de bord et de norme IUNO
  !          =  1 vitesse imposee : on donne IUPT IVPT IWPT
  !          =  2 profil de vitesse donne par l'utilisateur
  !     IJPRPD = 1 distribution uniforme
  !            = 2 profil de taux de presence donne par l'utilisateur
  !     IJPRTP = 1 profil plat de temperature donne par la valeur dans uslag2
  !            = 2 profil de temperature donne par l'utilisateur
  !     IJPRDP = 1 profil plat de diametre donne par la valeur dans uslag2
  !            = 2 profil dediametre donne par l'utilisateur
  !     INUCHL : numero du charbon de la particule (si IPHYLA=2)
  !     ICLST  : numero du groupe de statistiques

  integer, save ::           ijnbp, ijfre, ijuvw, ijprtp, ijprdp, ijprpd
  integer, save ::           inuchl, iclst

  !     RUSLAG  : tableau d info par classe et par frontieres

  double precision, allocatable, dimension(:,:,:) :: ruslag

  !     IUNO  : Norme de la vitesse
  !     IUPT  : U par classe et zones
  !     IVPT  : V par classe et zones
  !     IWPT  : W par classe et zones
  !     IDEBT : Debit
  !     IPOIT : Poids de la particule
  !     IDPT  : Diametre
  !     IVDPT : Variance du diametre
  !     ITPT  : Temperature
  !     ICPT  : Cp
  !     IEPSI : Emissivite des particules
  !     IROPT : Masse volumique
  !     IHPT  : Temperature
  !     IMWAT : Masse d eau dans le charbon
  !     IMCHT : Masse de charbon reactif
  !     IMCKT : Masse de coke
  !     IDCKT : Diametre du coeur retrecissant

  integer, save ::           iuno, iupt, ivpt, iwpt,                         &
                             itpt, idpt, ivdpt, iropt,                       &
                             icpt, ipoit, idebt, iepsi,                      &
                             ihpt, imwat, imcht, imckt, idckt

  !=============================================================================
  ! 8. Statistiques

  !     POINTEURS POUR LES STATISTIQUES
  !     -------------------------------

  !     ILVX,ILVY,ILVZ    : Vitesse
  !     ILFV              : Concentration volumique
  !     ILPD              : Somme des poids statistiques
  !     ILTS              : Temps de sejour

  !     ILTP              : Temperature
  !     ILDP              : Diametre
  !     ILMP              : Masse

  !     ILHP              : Temperature
  !     ILMWAT            : Masse d eau
  !     ILMCH             : Masse de charbon reactif
  !     ILMCK             : Masse de coke
  !     ILDCK             : Diametre du coeur retrecissant

  !     ILVU(NUSSTA)      : Statistiques supplementaires utilisateur

  integer, save ::           ilvx  , ilvy  , ilvz  ,                         &
                             ilpd  , ilfv  , ilts  ,                         &
                             iltp  , ildp  , ilmp  ,                         &
                             ilhp  , ilmwat, ilmch ,                         &
                             ilmck , ildck , ilvu(nussta)

  integer, save ::           iactfv, iactvx, iactvy, iactvz, iactts

  !     DONNEES POUR LES STATISTIQUES VOLUMIQUES
  !     ----------------------------------------

  !      ISTALA : Calcul statistiques       si  >= 1 sinon pas de stat
  !      ISUIST : Suite calcul statistiques si  >= 1 sinon pas de stat
  !      NVLSTS : NOMBRE DE VARIABLES STATISTIQUES SUPPLEMENTAIRES
  !               UTILISATEUR (CONTENU DANS NVLSTA)
  !      IDSTNT : Numero du pas de temps pour debut statistque
  !      NSTIST : Debut calcul stationnaire
  !      NPST   : Nombre de pas de temps pour le cumul des stats
  !      NPSTT  : Nombre de pas de temps total des stats depuis le debut
  !               du calcul, partie instationnaire comprise
  !      TSTAT  : Temps physique des stats volumiques
  !      SEUIL  : Seuil en POIDS STAT de particules pour les stats

  integer, save ::           istala , isuist , nvlsts ,                      &
                             idstnt , nstist ,                               &
                             npst   , npstt

  double precision, save ::  tstat , seuil

  !     NOMS DES VARIABLES STATISTIQUES (MOYENNES ET VARIANCES)
  !     -------------------------------------------------------
  !     Taille limitee par le fait qu on utilise NOMBRD dans
  !       l ecriture des fichiers suites (lagout)

  character*32, save ::      nomlag(nvplmx) , nomlav(nvplmx)

  !     OPTION POUR LES HISTORIQUES SUR LES STATS
  !     -----------------------------------------

  integer, save ::           ihslag(nvplmx)

  !     STATISTIQUE PAR ZONE ET PAR CLASSE
  !     ----------------------------------

  integer, save ::           nbclst

  !===============================================================================
  ! 9. Termes Sources

  !     OPTION TERMES SOURCES
  !     ---------------------
  !       Dynamique
  !       Masse
  !       Thermique

  integer, save ::          ltsdyn , ltsmas , ltsthe

  !     POINTEURS POUR LES TERMES SOURCES
  !     ---------------------------------

  !    ITSVX,ITSVY,ITVZ    : Termes sources sur la vitesse
  !    ITSLI               : Terme source implicite (vitesse+turbulence)
  !    ITSKE               : Terme source sur la turbulence en k-eps
  !    ITSR11,ITR12,ITSR13 : Termes sources sur la turbulence en Rij-Eps
  !    ITSR22,ITR23,ITSR33
  !    ITSTE, ITSTI        : Termes sources pour la thermique
  !    ITSMAS              : Terme source pour la masse
  !    ITSMV1              : Terme source sur F1 (MV legeres)
  !    ITSMV2              : Terme source sur F2 (MV lourdes)
  !    ITSCO               : Terme source sur F3 (C sous forme de CO)
  !    ITSFP4              : Variance du traceur relatif a l air

  integer, save ::           itsvx  , itsvy  , itsvz  , itsli ,              &
                             itske  ,                                        &
                             itsr11 , itsr12 , itsr13 ,                      &
                             itsr22 , itsr23 , itsr33 ,                      &
                             itste  , itsti  ,                               &
                             itsmas , itsmv1(ncharm2), itsmv2(ncharm2) ,     &
                             itsco  , itsfp4

  !     DONNEES POUR LES TERMES SOURCES
  !     -------------------------------

  !     NSTITS : debut calcul terme source stationnaire
  !     NPTS   : nombre de pas de temps pour le cumul des termes sources
  !     NTXERR : nombre de cellules qui un taux vol > 0.8
  !     VMAX   : taux volumique max atteint
  !     TMAMAX : taux massique max atteint

  integer, save ::           nstits , npts , ntxerr

  double precision, save ::  vmax , tmamax

  !=============================================================================
  ! 10. Clonage/fusion des particules

  !     INDICATEUR D ACTIVATION DE LA ROULETTE RUSSE

  integer, save ::           iroule

  !=============================================================================
  ! 11. Encrassement

  !     DONNEES POUR L ENCRASSEMENT

  integer, save ::           iencra , npencr

  double precision, save ::  enc1(ncharm2) , enc2(ncharm2) ,                 &
                             tprenc(ncharm2) , visref(ncharm2) , dnpenc

  !=============================================================================
  ! 12. Physico-chemical (DLVO) parameters

  !      cstham : Hamaker constant for the particle/fluid/substrate system
  !      epseau : Dielectric constant of the fluid
  !      phi1   : Electrokinetic potential of the first solid
  !      phi2   : Electrokinetic potential of the second solid
  !      fion   : Ionic force

  double precision, save ::  cstham , epseau,  phi1 , phi2, fion


!    Faraday constant (C/mol)

  double precision cstfar
  parameter(cstfar = 9.648d4)

!    Vacuum permittivity (F/m):

  double precision epsvid
  parameter(epsvid = 8.854d-12)

!    Boltzmann constant (J/K):

  double precision kboltz
  parameter(kboltz = 1.38d-23)

!    Cut-off distance for adhesion forces (assumed to be the Born distance) (m)

  double precision dcutof
  parameter(dcutof = 1.65d-10)

!    Characteristic retardation wavelength (m) for Hamaker constant

  double precision lambwl
  parameter(lambwl = 1000.d-9)

  !=============================================================================
  ! 13. Mouvement brownien

  !     ACTIVATION DU MOUVEMENT BROWNIEN :

  integer, save :: lamvbr


  !=============================================================================
  ! 14. Schema en temps, dispersion turbulente et equation de poisson

  !     NOR    : numero du sous-pas Lagrangien (1 ou 2)

  !     NORDRE : ordre de la methode d integration (1 ou 2)

  !     MODCPL : = 0 pour le modele incomplet
  !              > 0 pour le modele complet, est egal au nombre de
  !                 passages avant mise en route du modele complet

  !     IDIRLA : = 1 ou 2 ou 3 direction du modele complet

  !     IDISTU : = 0 pas de prise en compte de la dispersion turbulente (la
  !                  vitesse instantanee est egale a la vitesse moyenne)
  !              > 0 prise en compte de la dispersion turbulente (si k-eps
  !                  ou Rij-eps)

  !     IDIFFL : =1 la dispersion turbulente de la particule est celle de
  !                 la particule fluide (=0 sinon)

  !     ILAPOI : = 0 Pas de correction de pression
  !              = 1 Correction de pression

  integer, save ::           nor , nordre , modcpl , idirla ,                &
                             idistu , idiffl , ilapoi

  !=============================================================================
  ! 15. Traitement des statistiques interactions particules/frontieres

  !     DONNEES POUR LES STATISTIQUES AUX FRONTIERES
  !     --------------------------------------------

  !      NUSBOR : NOMBRE DE VARIABLES A ENREGISTRER SUR LES FRONTIERES
  !               SUPPLEMENTAIRES UTILISATEUR (CONTENU DANS NVISBR)
  !      NSTBOR : debut calcul stationnaire
  !      NPSTF  : nombre de pas de temps pour le cumul des stats
  !      NPSTF  : nombre de pas de temps total des stats depuis le debut
  !               du calcul, partie instationnaire comprise
  !      TSTATP : Temps physique des stats aux frontieres stationnaires
  !      SEUILF : Seuil en POIDS STAT de particules pour les stats
  !      IMOYBR : Type de moyenne applicable pour affichage et
  !               post-procesing

  integer, save ::           nusbor   , nstbor   ,                         &
                             npstf    , npstft   ,                         &
                             inbrbd   , iflmbd   , iangbd   , ivitbd   ,   &
                             iencnbbd , iencmabd , iencdibd , iencckbd ,   &
                             inbr     , iflm     , iang     , ivit     ,   &
                             iencnb   , iencma   , iencdi   , iencck   ,   &
                             iusb(nusbrd)    , imoybr(nusbrd+10)

  double precision, save ::  tstatp , seuilf

  !     NOMS DES VARIABLES STATISTIQUES
  !     -------------------------------
  !     Taille limitee par le fait qu on utilise NOMBRD dans
  !       l ecriture des fichiers suites (lagout)

  character*50, save ::      nombrd(nvplmx)

  ! IIFRLA Pointeur dans IA sur IFRLAG pour reperage des zones
  !          frontieres associees aux faces de bord

  integer, save ::           iifrla

  !=============================================================================
  ! 16. Visu

  !... Type de visualisation :
  !    IENSI3 : interaction particules/frontieres

  integer, save ::           iensi3

  !... Contenu des fichiers resultats


  integer, save ::           ivisv1 , ivisv2 , ivistp ,                      &
                             ivisdm , iviste , ivismp ,                      &
                             iviswat, ivisch , ivisck ,                      &
                             ivisdk

  !=============================================================================

contains

  !=============================================================================

  ! Local function to initialize Lagrangian module parameters for
  ! a given zone ii and class jj

  subroutine init_zone_class_param(ii, jj)

    use cstnum
    implicit none

    ! Arguments

    integer :: ii, jj

    ! Local variables

    integer :: kk

    ! define defaults (impossible values the user should override)

    do kk = 1, ndlaim
      iuslag(ii, jj, kk) = 0
    enddo
    iuslag(ii,jj,ijuvw) = -2
    iuslag(ii,jj,ijprtp) = -2
    iuslag(ii,jj,ijprdp) = -2
    iuslag(ii,jj,ijprpd) = -2
    do kk = 1, ndlagm
      ruslag(ii, jj, kk) = 0.d0
    enddo
    ruslag(ii,jj,iuno)  = -grand
    ruslag(ii,jj,iupt)  = -grand
    ruslag(ii,jj,ivpt)  = -grand
    ruslag(ii,jj,iwpt)  = -grand
    ruslag(ii,jj,ipoit) = -grand
    ruslag(ii,jj,idpt)  = -grand
    ruslag(ii,jj,ivdpt) = -grand
    ruslag(ii,jj,iropt) = -grand
    if (iphyla.eq.1) then
      if (itpvar.eq.1) then
        ruslag(ii,jj,itpt)  = -grand
        ruslag(ii,jj,icpt)  = -grand
        ruslag(ii,jj,iepsi) = -grand
      endif
    else if ( iphyla .eq. 2 ) then
      ruslag(ii,jj,ihpt)   = -grand
      ruslag(ii,jj,imwat)  = -grand
      ruslag(ii,jj,imcht)  = -grand
      ruslag(ii,jj,imckt)  = -grand
      ruslag(ii,jj,icpt)   = -grand
    endif

  end subroutine init_zone_class_param

  !=============================================================================

  !> \brief Initialize Lagrangian module parameters for a given zone and class

  !> \param[in]  i_cz_params  integer parameters for this class and zone
  !> \param[in]  r_cz_params  real parameters for this class and zone

  subroutine lagr_init_zone_class_param(i_cz_params, r_cz_params)  &
    bind(C, name='cs_lagr_init_zone_class_param')

    use, intrinsic :: iso_c_binding
    use cstnum
    implicit none

    ! Arguments

    integer, dimension(ndlaim)          :: i_cz_params
    double precision, dimension(ndlagm) :: r_cz_params

    ! Local variables

    integer :: ii

    ! define defaults (impossible values the user should override)

    do ii = 1, ndlaim
      i_cz_params(ii) = 0
    enddo
    i_cz_params(ijuvw) = -2
    i_cz_params(ijprtp) = -2
    i_cz_params(ijprdp) = -2
    i_cz_params(ijprpd) = -2

    do ii = 1, ndlagm
      r_cz_params(ii) = 0.d0
    enddo
    if (iphyla.eq.1) then
      if (itpvar.eq.1) then
        r_cz_params(itpt)  = -grand
        r_cz_params(icpt)  = -grand
        r_cz_params(iepsi) = -grand
      endif
    else if (iphyla .eq. 2) then
      r_cz_params(ihpt)  = -grand
      r_cz_params(imwat) = -grand
      r_cz_params(imcht) = -grand
      r_cz_params(imckt) = -grand
      r_cz_params(icpt)  = -grand
    endif

  end subroutine lagr_init_zone_class_param

  !=============================================================================

  !> \brief Define Lagrangian module parameters for a given zone and class

  !> \param[in]     class_id     id of given particle class
  !> \param[in]     zone_id      id of given boundary zone
  !> \param[in]     i_cz_params  integer parameters for this class and zone
  !> \param[in]     r_cz_params  real parameters for this class and zone

  subroutine lagr_define_zone_class_param(class_id, zone_id,         &
                                          i_cz_params, r_cz_params)  &
    bind(C, name='cs_lagr_define_zone_class_param')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer, value                      :: class_id
    integer, value                      :: zone_id
    integer, dimension(ndlaim)          :: i_cz_params
    double precision, dimension(ndlagm) :: r_cz_params

    ! Local variables

    integer :: ncmax, nzmax, mcmxp, nzmxp, ii, jj, kk
    integer, allocatable, dimension(:,:,:) :: itmp
    integer, dimension(3) :: shpe
    double precision, allocatable, dimension(:,:,:) :: rtmp

    ! Allocate on first pass

    if (.not.allocated(iuslag) .or. .not.allocated(ruslag)) then
      allocate(iuslag(1,1,ndlaim))
      allocate(ruslag(1,1,ndlagm))
      call init_zone_class_param(1, 1)
    endif

    ! Reallocate arrays if required
    ! (use size margin to avoid reallocating too often, though this
    ! should only impact the first time step)

    shpe = shape(iuslag)
    ncmax = shpe(1)
    nzmax = shpe(2)

    if (class_id.gt.ncmax .or. zone_id.gt.nzmax) then

      mcmxp = ncmax
      nzmxp = nzmax

      ncmax = max(class_id, mcmxp+5)
      nzmax = max(zone_id, nzmxp+5)

      ! Save iuslag and ruslag arrays

      allocate(itmp(mcmxp,nzmxp,ndlaim))
      allocate(rtmp(mcmxp,nzmxp,ndlagm))

      do ii = 1, mcmxp
        do jj = 1, nzmxp
          do kk = 1, ndlaim
            itmp(ii,jj,kk) = iuslag(ii,jj,kk)
          enddo
          do kk = 1, ndlagm
            rtmp(ii,jj,kk) = ruslag(ii,jj,kk)
          enddo
        enddo
      enddo

      ! Reallocate iuslag and ruslag arrays

      deallocate(iuslag)
      deallocate(ruslag)
      allocate(iuslag(ncmax,nzmax,ndlaim))
      allocate(ruslag(ncmax,nzmax,ndlagm))

      ! Restore saved values, and initialize new entries

      do ii = 1, mcmxp
        do jj = 1, nzmxp
          do kk = 1, ndlaim
            iuslag(ii,jj,kk) = itmp(ii,jj,kk)
          enddo
          do kk = 1, ndlagm
            ruslag(ii,jj,kk) = rtmp(ii,jj,kk)
          enddo
        enddo
        do jj = nzmxp + 1, nzmax
          call init_zone_class_param(ii, jj)
        enddo
      enddo
      do ii = mcmxp + 1, ncmax
        do jj = 1, nzmxp
          call init_zone_class_param(ii, jj)
        enddo
      enddo

      deallocate(rtmp)
      deallocate(itmp)

    endif

    ! Now copy defined values

    do kk = 1, ndlaim
      iuslag(class_id, zone_id, kk) = i_cz_params(kk)
    enddo
    do kk = 1, ndlagm
      ruslag(class_id, zone_id, kk) = r_cz_params(kk)
    enddo

  end subroutine lagr_define_zone_class_param

  !=============================================================================

  !> \brief Return Lagrangian model status.

  !> \param[out]   model     0 without Lagrangian, 1 or 2 with Lagrangian
  !> \param[out]   restart   1 for Lagrangian restart, 0 otherwise
  !> \param[out]   frozen    1 for frozen Eulerian flow, 0 otherwise

  subroutine lagr_status(model, restart, frozen)  &
    bind(C, name='cs_lagr_status')

    use, intrinsic :: iso_c_binding
    implicit none

    ! Arguments

    integer(c_int), intent(out) :: model, restart, frozen

    model = iilagr
    restart = isuila
    frozen = isttio

    return

  end subroutine lagr_status

  !=============================================================================

end module lagran
