!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

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

!                              optcal.h
!===============================================================================

! DEFINITION DES EQUATIONS
!   ISTAT
!     = 1 PRISE EN COMPTE DU TERME INSTATIONNAIRE
!     = 0 PRISE EN COMPTE DU TERME INSTATIONNAIRE
!   ICONV
!     = 1 PRISE EN COMPTE DE LA CONVECTION
!     = 0 NON PRISE EN COMPTE DE LA CONVECTION
!   IDIFF
!     = 1 PRISE EN COMPTE DE LA DIFFUSION (MOLECULAIRE ET TURBULENTE)
!     = 0 NON PRISE EN COMPTE DE LA DIFFUSION (MOLECULAIRE ET TURBULENTE)
!   IDIFFT : SI IDIFF = 1
!     = 1 PRISE EN COMPTE DE LA DIFFUSION TURBULENTE
!     = 0 NON PRISE EN COMPTE DE LA DIFFUSION TURBULENTE

integer           istat (nvarmx),iconv (nvarmx),idiff (nvarmx),   &
                  idifft(nvarmx)
common / iequat / istat         ,iconv         ,idiff         ,   &
                  idifft


! PROPRIETES PHYSIQUES RHO ET VISCL CONSTANTES OU VARIABLES
!    =1 VARIABLE, =0 CONSTANT
!     SERT LORS DES LECTURES DE FICHIER SUITE POUR EVITER D'ECRASER
!     LA VALEUR FOURNIE PAR LA VALEUR DE L'ANCIEN CALCUL.
integer           irovar(nphsmx),ivivar(nphsmx)
common / iphvar / irovar        ,ivivar

! SCHEMA EN TEMPS

!  ISCHTP : INDICATEUR DE SCHEMA EN TEMPS
!     = 2 : ORDRE 2
!     = 1 : STANDARD
!  ISTMPF : INDICATEUR DE SCHEMA FLUX DE MASSE
!     = 2 THETA SCHEMA avec THETA > 0 (= 0.5 : ordre 2)
!     = 0 THETA SCHEMA avec THETA = 0 (explicite)
!     = 1 SCHEMA STANDARD V1.0
!  NTERUP : Nombre d'iteration sur navier-stokes pour couplage vitesse/
!           pression
!  ISNO2T : INDICATEUR D'EXTRAPOLATION DE TERMES SOURCES NAVIER STOKES
!           POUR LE SCHEMA EN TEMPS
!  ISTO2T : INDICATEUR D'EXTRAPOLATION DE TERMES SOURCES DES GRANDEURS
!           TURBULENTES POUR LE SCHEMA EN TEMPS
!  ISSO2T : INDICATEUR D'EXTRAPOLATION DE TERMES SOURCES DES SCALAIRES
!           POUR LE THETA SCHEMA EN TEMPS
!  IROEXT : INDICATEUR D'EXTRAPOLATION DE LA MASSE VOLUMIQUE
!           POUR LE SCHEMA EN TEMPS
!  IVIEXT : INDICATEUR D'EXTRAPOLATION DE LA VISCOSITE TOTALE
!           POUR LE SCHEMA EN TEMPS
!  IVSEXT : INDICATEUR D'EXTRAPOLATION DE LA DIFFUSIVITE SCALAIRE

!  INITVI : =1 SI VISCOSITE TOTALE RELUE DANS UN SUITE

!  INITRO : =1 SI MASSE VOLUMIQUE RELUE DANS UN SUITE

!  ICPEXT : INDICATEUR D'EXTRAPOLATION DE LA MASSE VOLUMIQUE
!           POUR LE SCHEMA EN TEMPS

!  INITCP : =1 SI  CHALEUR SPECIFIQUE RELUE DANS UN SUITE
!  INITVS : =1 SI  DIFFUSIVITE SCALAIRE RELUE DANS UN SUITE

!  THETAV : PONDERATION ENTRE LES PAS DE TEMPS N ET N+1 POUR LES
!           VARIABLE PRINCIPALES
!     = 1 : SCHEMA EULER IMPLICITE
!     =1/2: SCHEMA CENTRE EN TEMPS

!  THETSN : SCHEMA EN TEMPS POUR LES TERMES SOURCES DE NAVIER STOKES
!     = 0 : VISCOSITE SECONDAIRE EXPLICITE
!     =1/2: VISCOSITE SECONDAIRE EXTRAPOLEE EN N+1/2
!     = 1 : VISCOSITE SECONDAIRE EXTRAPOLEE EN N+1
!  THETST : SCHEMA EN TEMPS POUR LES TERMES SOURCES DES GRANDEURS TURBULENTES
!     = 0 : VISCOSITE SECONDAIRE EXPLICITE
!     =1/2: VISCOSITE SECONDAIRE EXTRAPOLEE EN N+1/2
!     = 1 : VISCOSITE SECONDAIRE EXTRAPOLEE EN N+1
!  THETSS : SCHEMA EN TEMPS POUR LES TERMES SOURCES DES SCALAIRES
!     = 0 : VISCOSITE SECONDAIRE EXPLICITE
!     =1/2: VISCOSITE SECONDAIRE EXTRAPOLEE EN N+1/2
!     = 1 : VISCOSITE SECONDAIRE EXTRAPOLEE EN N+1
!  THETFL : SCHEMA EN TEMPS POUR LE FLUX DE MASSE
!     = 0 : FLUX DE MASSE EXPLICITE
!     =1/2: FLUX DE MASSE EXTRAPOLE EN N+1/2
!     = 1 : FLUX DE MASSE EXTRAPOLE EN N+1
!  THETVI : SCHEMA EN TEMPS POUR LA VISCOSITE TOTALE
!     = 0 : VISCOSITE TOTALE EXPLICITE
!     =1/2: VISCOSITE TOTALE EXTRAPOLEE EN N+1/2
!     = 1 : VISCOSITE TOTALE EXTRAPOLEE EN N+1
!  THETRO : SCHEMA EN TEMPS POUR LA MASSE VOLUMIQUE
!     = 0 : MASSE VOLUMIQUE TOTALE EXPLICITE
!     =1/2: MASSE VOLUMIQUE TOTALE EXTRAPOLEE EN N+1/2
!     = 1 : MASSE VOLUMIQUE EXTRAPOLEE EN N+1
!  THETCP : SCHEMA EN TEMPS POUR LA MASSE VOLUMIQUE
!     = 0 : CHALEUR SPECIFIQUE TOTALE EXPLICITE
!     =1/2: CHALEUR SPECIFIQUE TOTALE EXTRAPOLEE EN N+1/2
!     = 1 : CHALEUR SPECIFIQUE EXTRAPOLEE EN N+1
!  EPSUP  : TESTS DE CONVERGENCE DU SYSTEME VITESSE/PRESSION QUAND CE
!           DERNIER EST RESOLU PAR SOUS-ITERATIONS (POINT FIXE)
!  XNRMU  : NORME DE U(k+1) - U(k)
!  XNRMU0 : NORME DE U(0)

integer           nterup,                                         &
                  ischtp(nphsmx), istmpf(nphsmx),                 &
                  isno2t(nphsmx), isto2t(nphsmx), isso2t(nscamx), &
                  iroext(nphsmx),                                 &
                  iviext(nphsmx), icpext(nphsmx), ivsext(nscamx), &
                  initro(nphsmx), initvi(nphsmx),                 &
                  initcp(nphsmx), initvs(nscamx)
common / ievtmp / nterup,                                         &
                  ischtp        , istmpf        ,                 &
                  isno2t        , isto2t        , isso2t        , &
                  iroext        ,                                 &
                  iviext        , icpext        , ivsext        , &
                  initro        , initvi        ,                 &
                  initcp        , initvs
double precision  thetav(nvarmx), thetsn(nphsmx), thetst(nphsmx), &
                  thetss(nscamx),                                 &
                  thetfl(nphsmx), thetro(nphsmx), thetvi(nphsmx), &
                  thetcp(nphsmx), thetvs(nscamx), epsup (nphsmx), &
                  xnrmu0(nphsmx), xnrmu (nphsmx)
common / revtmp / thetav        , thetsn        , thetst        , &
                  thetss        ,                                 &
                  thetfl        , thetro        , thetvi        , &
                  thetcp        , thetvs        , epsup         , &
                  xnrmu0        , xnrmu

! SCHEMA CONVECTIF

!  BLENCV : 100*(1-BLENCV) EST LE POURCENTAGE D'UPWIND
!     = 1 : PAS D'UPWIND EN DEHORS DU TEST DE PENTE
!     = 0 : UPWIND
!  ISCHCV : SCHEMA CONVECTIF CENTRE OU SECOND ORDER
!     = 1 : CENTRE
!     = 0 : SECOND ORDER
!  ISSTPC : INDICATEUR SANS OU AVEC TEST DE PENTE
!     = 1 : SANS TEST DE PENTE
!     = 0 : AVEC TEST DE PENTE

integer           ischcv(nvarmx), isstpc(nvarmx)
common / icnvsc / ischcv        , isstpc
double precision                  blencv(nvarmx)
common / rcnvsc /                 blencv


! RECONSTRUCTION DES GRADIENTS ET DES SECONDS MEMBRES
!   IMRGRA : METHODE DE RECONTRUCTION DES GRADIENTS
!     = 0  : RECONTRUCTION 97
!     = 1  : MOINDRES CARRES 99
!     = 2  : MOINDRES CARRES SUPPORT ETENDU COMPLET
!     = 3  : MOINDRES CARRES AVEC SELECTION DU SUPPORT ETENDU
!     = 4  : RECONSTRUCTION 97 AVEC INITIALISATION MOINDRES CARRES
!   ANOMAX : ANGLE DE NON ORTHOGONALITE DES FACES EN RADIAN AU DELA DUQUEL
!            ON RETIENT DANS LE SUPPORT ETENDU DES CELLULES VOISINES
!            DE LA FACE LES CELLULES DONT UN NOEUD EST SUR LA FACE
!   NSWRGR : NOMBRE DE SWEEPS DE RECONSTRUCTION DES GRADIENTS 97
!   NSWRSM : NOMBRE DE SWEEPS DE RECONSTRUCTION DES SECONDS MEMBRES
!   EPSRGR : PRECISION POUR LA   RECONSTRUCTION DES GRADIENTS 97
!   EPSRSM : PRECISION POUR LA   RECONSTRUCTION DES SECONDS MEMBRES
!   IMLIGR : LIMITATION DES GRADIENTS
!     < 0  : PAS DE LIMITATION DES GRADIENTS
!     = 0  : PREMIER ORDRE
!     = 1  : SECOND ORDRE
!   CLIMGR : FACTEUR DE LIMITATION (>=1, =1 : FORTE LIMITATION)
!   IRCFLU : RECONSTRUCTION DES FLUX AUX FACES
!     = 0  : NON
!     = 1  : OUI
!   EXTRAG : EXTRAPOLATION DES GRADIENTS AU BORD (0 <= EXTRAG <= 1)
!     = 0  : NON
!     = 1  : OUI

integer           imrgra, nswrgr(nvarmx), nswrsm(nvarmx),         &
                  imligr(nvarmx)        , ircflu(nvarmx)
common / irecgr / imrgra, nswrgr        , nswrsm        ,         &
                  imligr                , ircflu

double precision  anomax ,                                        &
                  epsrgr(nvarmx), epsrsm(nvarmx),                 &
                  climgr(nvarmx), extrag(nvarmx)
common / rrecgr / anomax ,                                        &
                  epsrgr        , epsrsm        ,                 &
                  climgr        , extrag


! SOLVEURS ITERATIFS
!   NITMAX : NOMBRE D'ITERATIONS MAX
!   EPSILO : PRECISION RELATIVE CHERCHEE
!   IRESOL
!     =-1 : CALCULE AUTOMATIQUEMENT (0 SI ICONV=0, 1 SINON)
!     = 0 : GRADIENT CONJUGUE
!     = 1 : JACOBI
!     = 2 : BI-CGSTAB
!    et ON AJOUTE IPOL*1000 OU IPOL EST LE DEGRE DU POLYNOME DE
!       PRECONDITIONNEMENT DE NEUMANN
!     En pratique, il semble que ce preconditonnement ne soit pas efficace
!        on gagne 10% CPU sur un cas, on perd 3% sur un autre avec IPOL=1
!        on perd avec IPOL=2
!        Ces valeurs ont ete obtenues sur de petits cas.
!   IDIRCL : DECALAGE DE LA DIAGONALE DE LA MATRICE S'IL N'Y A PAS DE DIRICHLET
!     = 0 : NON
!     = 1 : OUI
!     Le code calcule automatiquement pour chaque variable NDIRCL, nombre de
!        CL de Dirichlet, et en deduit s'il doit decaler ou pas la diagonale

integer           nitmax(nvarmx),iresol(nvarmx),idircl(nvarmx),   &
                  ndircl(nvarmx)
common / inivcv / nitmax        ,iresol        ,idircl        ,   &
                  ndircl

double precision                 epsilo(nvarmx)
common / rnivcv /                epsilo


! MULTIGRILLE
!   IMGR
!     = 0 PAS DE MULTIGRILLE
!     = 1        MULTIGRILLE ALGEBRIQUE
!   NCYMAX : NOMBRE MAX DE CYCLES
!   NITMGF : NOMBRE D'ITER SUR MAILLAGE FIN
!   RLXP1  :

integer           imgr(nvarmx), ncymax(nvarmx), nitmgf(nvarmx)
common / imultg / imgr        , ncymax        , nitmgf

double precision  rlxp1
common / rmultg / rlxp1


! GESTION DU CALCUL
!   ISUITE : SUITE DE CALCUL
!     = 0 POUR SFS
!     = 1 POUR SUITE DE CALCUL
!   ISCOLD : correspondance nouveaux-anciens scalaires
!   IECAUX : ecriture du suite auxiliaire
!   ILEAUX : lecture  du suite auxiliaire
!   ISUIT1 : suite du module thermique 1d en paroi
!   ISUICT : suite du module aerorefrigerant
!   ISUIVO : suite de la methode des vortex

integer           isuite , ileaux, iecaux, iscold(nscamx),        &
                  isuit1 , isuict, isuivo
common / istart / isuite , ileaux, iecaux, iscold, isuit1, isuict, isuivo


! GESTION DES PAS DE TEMPS
!   NTPABS : PAS DE TEMPS PRECEDENT ABSOLU
!   NTCABS : PAS DE TEMPS COURANT   ABSOLU
!   NTMABS : PAS DE TEMPS MAX       ABSOLU
!   TTPABS :        TEMPS PRECEDENT ABSOLU
!   TTCABS :        TEMPS COURANT   ABSOLU
!   TTMABS :        TEMPS MAX       ABSOLU
!   INPDT0 : INDICATEUR "ZERO PAS DE TEMPS"

!   NTMABS = numero absolu du dernier pas de temps desire
!            Si on a deja fait 10 pas de temps
!              et qu'on veut en faire 10 autres,
!              il faut affecter 10 + 10 = 20 a NTMABS
!   NTPABS = numero relu dans le fichier suite
!   NTCABS = incremente au debut du pas de temps
!              et donc initialise a NTPABS
!   INPDT0 = 1 pour ne faire aucun pas de temps (0 sinon)
!              Pour les calculs non suite :
!                on saute uniquement les resolutions (navier-stokes,
!                  turbulence, scalaires...)
!              Pour les calculs suite :
!                on saute les resolutions (navier-stokes,
!                  turbulence, scalaires...) et le calcul des proprietes
!                  physiques, les conditions aux limites (les grandeurs
!                  sont lues dans le fichier suite)

integer           ntpabs, ntcabs, ntmabs, inpdt0
common / itemps / ntpabs, ntcabs, ntmabs, inpdt0

double precision  ttpabs, ttcabs
common / rtemps / ttpabs, ttcabs


! OPTION PAS DE TEMPS
!   IDTVAR : PAS DE TEMPS VARIABLE
!     = 0 : PAS DE TEMPS CONSTANT
!     = 1 : PAS DE TEMPS UNIFORME EN ESPACE ET VARIABLE EN TEMPS
!     = 2 : PAS DE TEMPS VARIABLE EN ESPACE ET VARIABLE EN TEMPS
!   IPTLRO : LIMITATION DU PAS DE TEMPS LIEE AUX EFFETS DE DENSITE
!     = 0 : NON
!     = 1 : OUI
!   COUMAX : NOMBRE DE COURANT         MAXIMUM        (IDTVAR NON NUL)
!   FOUMAX : NOMBRE DE         FOURIER MAXIMUM        (IDTVAR NON NUL)
!   VARRDT : VARIATION RELATIVE PERMISE DE DT         (IDTVAR NON NUL)
!   DTMIN, DTMAX : VALEUR LIMITE MIN ET MAX DE DT     (IDTVAR NON NUL)
!       PRENDRE POUR DTMAX = MAX (Ld/Ud, SQRT(Lt/(gDelta rho/rho)), ...)
!   CDTVAR : COEF MULTIPLICATIF POUR LE PAS DE TEMPS DE CHAQUE VARIABLE
!         POUR U,V,W,P IL EST INUTILISE
!         POUR K,E    ON PREND LA MEME VALEUR : CELLE DE K
!         POUR RIJ, E ON PREND LA MEME VALEUR : CELLE DE R11
!   RELAXV : RELAXATION DES VARIABLES (1 PAS DE RELAX)
!   RELXST : COEFFICIENT DE RELAXATION DE BASE STATIONNAIRE

integer           idtvar,iptlro
common / iptvar / idtvar,iptlro

double precision  dtref,coumax,foumax,                            &
                  dtmin,dtmax ,varrdt,cdtvar(nvarmx),             &
                  relaxv(nvarmx), relxst
common / rptvar / dtref,coumax,foumax,                            &
                  dtmin,dtmax ,varrdt,cdtvar,relaxv,relxst


! TURBULENCE
!  ITURB
!    = 0  PAS DE TURBULENCE
!    = 10 LONGUEUR DE MELANGE
!    = 20, 21 K-EPSILON
!         * 20 MODELE STANDARD
!         * 21 MODELE A PRODUCTION LINEAIRE
!    = 30, 31 RIJ-EPSILON
!         * 30 MODELE STANDARD (LRR)
!         * 31 MODELE SSG
!    = 40, 41, 42 LES
!         * 40 MODELE DE SMAGORINSKY CONSTANT
!         * 41 MODELE DE SMAGORINSKY DYNAMIQUE "CLASSIQUE"
!         * 42 MODELE DE SMAGORINSKY DYNAMIQUE DE "PIOMELLI ET LIU"
!    = 50 v2f phi-model
!    = 60 K-OMEGA SST
!  ITYTUR
!    = INT(ITURB/10) POUR DISTINGUER RAPIDEMENT LES CLASSES DE MODELES
!  IDEUCH
!    = 0 UNE ECHELLE       (DEUX ECHELLES = FAUX)
!    = 1 DEUX ECHELLES     (DEUX ECHELLES = VRAI)
!    = 2 DEUX ECHELLES LIMITATION DE YPLUS A YPLULI (SCALABLE WALL FUNCTION)
!  ILOGPO
!    = 0 UNE ECHELLE  AVEC LOI EN PUISSANCE
!    = 1 UNE ECHELLES AVEC LOI LOG
!  ICLKEP
!    = 0 CLIPPING EN VALEUR ABSOLUE DE K ET EPSILON
!    = 1 CLIPPING COUPLE K-EPSILON BASE SUR DES RELATIONS PHYSIQUES
!  IGRHOK
!    = 1     PRISE EN COMPTE DE 2/3 RHO GRAD K DANS NAVIER STOKES
!    = 0 NON PRISE EN COMPTE DE 2/3 RHO GRAD K DANS NAVIER STOKES
!  IGRAKE
!    = 1 GRAVITE DANS K-EPSILON
!    = 0 SINON
!  IGRARI
!    = 1 GRAVITE DANS RIJ-EPSILON
!    = 0 SINON
!  ISCALT NUMERO DU SCALAIRE QUI TIENT LIEU DE TEMPERATURE
!    DONC VARIABLE ISCA(ISCALT)
!  IKECOU
!    = 1 K-EPSILON COUPLE EN INCREMENTS
!    = 0 SINON
!  IRIJNU
!         = 1 VISCOSITE DANS LA MATRICE EN INCREMENTS DE VITESSE (RIJ)
!         = 0 SINON
!  IRIJRB
!         = 1 TRAITEMENT PRECIS DE RIJ AU BORD, VOIR CONDLI      (RIJ)
!         = 0 SINON
!  IDIFRE
!         = 1 TRAITEMENT COMPLET DE LA DIAGONALE DU TENSEUR DE
!             DIFFUSION DE RIJ ET EPSILON (RIJ)
!         = 0 TRAITEMENT SIMPLIFIE
!  ICLSYR
!         = 1 IMPLICITATION PARTIELLE DE RIJ DANS LES CL DE SYMETRIE
!         = 0 PAS D'IMPLICITATION
!  ICLPTR
!         = 1 IMPLICITATION PARTIELLE DE RIJ ET EPSILON DANS LES CL
!             DE PAROI TURBULENTE
!         = 0 PAS D'IMPLICITATION
!  IDRIES : AMORTISSEMENT DE TYPE VAN DRIEST A LA PAROI
!         = 0 SANS AMORTISSEMENT
!         = 1 AVEC AMORTISSEMENT
!  IVRTEX : UTILISATION DE LA METHODE DES VORTEX
!         = 0 SANS METHODE DES VORTEX
!         = 1 AVEC METHODE DES VORTEX

integer           iturb(nphsmx) , itytur(nphsmx),                 &
                  ideuch(nphsmx), ilogpo(nphsmx), iclkep(nphsmx), &
                  igrhok(nphsmx), igrake(nphsmx),                 &
                  iscalt(nphsmx), ikecou(nphsmx),                 &
                  irijnu(nphsmx), irijrb(nphsmx), irijec(nphsmx), &
                  igrari(nphsmx), idifre(nphsmx), iclsyr(nphsmx), &
                  iclptr(nphsmx), idries(nphsmx), ivrtex
common / iturbu / iturb         , itytur        ,                 &
                  ideuch        , ilogpo        , iclkep        , &
                  igrhok        , igrake        ,                 &
                  iscalt        , ikecou        ,                 &
                  irijnu        , irijrb        , irijec        , &
                  igrari        , idifre        , iclsyr        , &
                  iclptr        , idries        , ivrtex



!   IVISSE PRISE EN COMPTE DE -2/3 GRAD(MU DIV(U)) + DIV(MU (GRAD_T(U)))

integer           ivisse(nphsmx)
common / ivisc2 / ivisse

! STOKES
!   IREVMC
!     = 2 POUR RECONSTRUCTION DES VITESSES DE TYPE RT0
!     = 1 POUR RECONSTRUCTION DES VITESSES AVEC GRADIENT DE L'INCREMENT
!           DE PRESSION PAR MOINDRES CARRES
!     = 0 SINON
!   IPRCO
!     = 0 POUR CALCUL SANS PRESSION CONTINUITE
!     = 1 POUR CALCUL AVEC PRESSION CONTINUITE
!   ARAK PROPORTION D'ARAKAWA (1 POUR ARAKAWA COMPLET)
!   RELAXV RELAXATION DES VARIABLES (1 PAS DE RELAX)
!   RNORMP NORMALISATION POUR LA CONVERGENCE DE RESOLP

integer           irevmc(nphsmx), iprco , irnpnw
common / istoke / irevmc        , iprco , irnpnw

double precision  rnormp(nphsmx), arak(nphsmx)
common / rstoke / rnormp        , arak


!   IPUCOU ALGORITHME COUPLAGE INSTATIONNAIRE VITESSE/PRESSION

integer           ipucou
common / coupup / ipucou


!   ICCVFG CALCUL A CHAMP DE VITESSE FIGE

integer           iccvfg
common / icfige / iccvfg

! CALCUL DE LA VISCOSITE

integer           imvisf
common / rvscfa / imvisf

!  TYPE DES CONDITIONS LIMITES ET INDEX MIN ET MAX
!                  DES SOUS LISTES DEFACES DE BORD

integer           idebty(ntypmx,nphsmx), ifinty(ntypmx,nphsmx)
common / itycli / idebty        , ifinty


!  ITRBRB = 1 TRAITEMENT PRECIS DE LA TEMPERATURE AU BORD, VOIR CONDLI
!             (UTILISE POUR COUPLAGE SYRTHES)
!         = 0 SINON
!  ICPSYR = 1 SI SCALAIRE COUPLE A SYRTHES
!    DONC POUR LE MOMENT VAUT 1 POUR ISCALT UNIQUEMENT

integer           itbrrb, icpsyr(nscamx)
common / couplb / itbrrb, icpsyr

!   PRISE EN COMPTE DE l'EQUILIBRE ENTRE LE GRADIENT DE PRESSION
!        ET LES TERMES SOURCES DE GRAVITE ET DE PERTE DE CHARGE

!     IPHYDR = 0 ALGORITHME SANS PRISE EN COMPTE DE L'EQUILIBRE
!            = 1 ALGORITHME AVEC PRISE EN COMPTE DE L'EQUILIBRE
!     ICALHY = 0 PAS DE CALCUL DE LA PRESSION HYDROSTATIQUE POUR LES
!                DIRICHLETS DE PRESSION EN SORTIE
!            = 1        CALCUL DE LA PRESSION HYDROSTATIQUE POUR LES
!                DIRICHLETS DE PRESSION EN SORTIE

integer           iphydr, icalhy
common / iprehy / iphydr, icalhy


!   CALCUL DES ESTIMATEURS

integer           iescal(nestmx,nphsmx)
common / icaest / iescal


!   CALCUL DES MOYENNES TEMPORELLES (CALCUL DES MOMENTS)

!  NBMOMT : NOMBRE DE MOYENNES DEMANDEES
!  NBDTCM : NOMBRE DE TABLEAUX NCEL POUR LE TEMPS CUMULE
!  NTDMOM : NUMERO DU PAS DE TEMPS INITIAL POUR LE CALCUL DU MOMENT
!  IMOOLD : NUMERO DE L'ANCIEN MOMENT CORRESPONDANT EN CAS DE SUITE
!  ICMOME : POINTEUR POUR LES MOMENTS (donne un numero de propriete)
!           s'utilise ainsi PROPCE(IEL,IPPROC(ICMOME(IMOM)))
!  IDTMOM : NUMERO DU TEMPS CUMULE ASSOCIE AUX MOMENTS
!           ce numero va de 1 a n pour les temps cumules non uniformes
!                     et de -1 a -p pour les temps cumules uniformes
!           s'utilise ainsi :
!              si IDTMOM(IMOM) > 0 PROPCE(IEL,IPROPC(ICDTMO(IDTMOM(IMOM))))
!              si IDTMOM(IMOM) < 0 DTCMOM(-IDTMOM(IMOM))
!  IDFMOM : NUMERO DES VARIABLES COMPOSANT LE MOMENT IDFMOM(JJ,IMOM)
!  IDGMOM : DEGRE DU MOMENT
!  ICDTMO : NUMERO DE PROPRIETE DU TEMPS CUMULE (voir IDTMOM)
!  IPPMOM : REPERE POUR LE POST SI ON DOIT DIVISER LA VARIABLE
!           PAR UN TEMPS CUMULE (voir memtri et useevo)
!  DTCMOM : VALEUR DU PAS DE TEMPS CUMULE QUAND IL EST UNIFORME (voir IDTMOM).

integer           nbmomt, nbdtcm,                                 &
                  ntdmom(nbmomx), imoold(nbmomx),                 &
                  icmome(nbmomx), idtmom(nbmomx),                 &
                  idfmom(ndgmox,nbmomx),          idgmom(nbmomx), &
                  icdtmo(nbmomx), ippmom(nvppmx)
common / imomen / nbmomt, nbdtcm,                                 &
                  ntdmom        , imoold        ,                 &
                  icmome        , idtmom        ,                 &
                  idfmom                        , idgmom        , &
                  icdtmo        , ippmom
double precision  dtcmom(nbmomx)
common / rmomen / dtcmom


!   INDICATEUR PERTES DE CHARGE GLOBAL (IE SOMME SUR LES PROCESSEURS
!       DE NCEPDC)

integer           ncpdct(nphsmx)
common / icpdct / ncpdct

!   INDICATEUR MODULE THERMIQUE 1D GLOBAL (IE SOMME SUR LES PROCESSEURS
!       DE NFPT1D)

integer           nfpt1t
common / ict1dt / nfpt1t

!   INDICATEUR TERMES SOURCES DE MASSE GLOBAL (IE SOMME SUR LES PROCESSEURS
!       DE NCETSM)

integer           nctsmt(nphsmx)
common / ictsmt / nctsmt

!   INDICATEUR DE PASSAGE DANS L'INITIALISATION DES
!                         VARIABLES PAR L'UTILISATEUR
!          IUSINI = 1 PASSAGE DANS USINIV OU PPINIV
!                   0 PAS DE PASSAGE (NI IUSINI NI PPINIV)
!          IUSCFP = 1 PASSAGE DANS USCFPV
!                   0 PAS DE PASSAGE

integer           iusini, iuscfp
common / iusspg / iusini, iuscfp

! PARAMETRES NUMERIQUES POUR LE CALCUL DE LA DISTANCE A LA PAROI

! INEEDY : = 1 DISTANCE A LA PAROI EST NECESSAIRE POUR LE CALCUL
!          = 0 DISTANCE A LA PAROI N'EST PAS NECESSAIRE
! IMAJDY : = 1 DISTANCE A LA PAROI A ETE MISE A JOUR
!          = 0 DISTANCE A LA PAROI N'A PAS ETE MISE A JOUR
! ICDPAR : = 1 CALCUL STANDARD (ET RELECTURE EN SUITE DE CALCUL)
!          = 2 CALCUL ANCIEN   (ET RELECTURE EN SUITE DE CALCUL)
!          =-1 FORCER LE RECALCUL EN SUITE (PAR CALCUL STANDARD)
!          =-2 FORCER LE RECALCUL EN SUITE (PAR CALCUL ANCIEN)
! NITMAY : NOMBRE MAX D'ITERATIONS POUR LES RESOLUTIONS ITERATIVES
! NSWRSY : NOMBRE DE SWEEP POUR RECONSTRUCTION DES S.M.
! NSWRGY : NOMBRE DE SWEEP POUR RECONSTRUCTION DES GRADIENTS
! IMLIGY : METHODE DE LIMITATION DU GRADIENT
! IRCFLY : INDICATEUR POUR RECONSTRUCTION DES FLUX
! ISCHCY : INDICATEUR DU SCHEMA EN ESPACE
! ISSTPY : INDICATEUR POUR TEST DE PENTE
! IMGRPY : MULTIGRILLE
! IWARNY : NIVEAU D'IMPRESSION
! NTCMXY : NOMBRE MAX D'ITERATION POUR LA CONVECTION DE Y

integer           ineedy       , imajdy        , icdpar      ,    &
                  nitmay       , nswrsy        , nswrgy      ,    &
                  imligy       , ircfly        , ischcy      ,    &
                  isstpy       , imgrpy        , iwarny      ,    &
                  ntcmxy
common / idpopt / ineedy       , imajdy        , icdpar      ,    &
                  nitmay       , nswrsy        , nswrgy      ,    &
                  imligy       , ircfly        , ischcy      ,    &
                  isstpy       , imgrpy        , iwarny      ,    &
                  ntcmxy

! BLENCY : 1 - PROPORTION D'UPWIND
! EPSILY : PRECISION POUR RESOLUTION ITERATIVE
! EPSRSY : PRECISION POUR LA RECONSTRUCTION DU SECOND MEMBRE
! EPSRGY : PRECISION POUR LA RECONSTRUCTION DES GRADIENTS
! CLIMGY : COEF GRADIENT*DISTANCE/ECART
! EXTRAY : COEF D'EXTRAPOLATION DES GRADIENTS
! COUMXY : VALEUR MAX   DU COURANT POUR EQUATION CONVECTION
! EPSCVY : PRECISION POUR CONVERGENCE EQUATION CONVECTION STATIONNAIRE
! YPLMXY : VALEUR MAX   DE YPLUS AU DESSUS DE LAQUELLE L'AMORTISSEMENT DE
!          VAN DRIEST EST SANS EFFET ET DONC POUR LAQUELLE UN CALCUL DE
!          YPLUS MOINS PRECIS EST SUFFISANT

double precision  blency      , epsily         , epsrsy      ,    &
                  epsrgy      , climgy         , extray      ,    &
                  coumxy      , epscvy         , yplmxy
common / rdpopt / blency      , epsily         , epsrsy      ,    &
                  epsrgy      , climgy         , extray      ,    &
                  coumxy      , epscvy         , yplmxy


! PARAMETRES NUMERIQUES POUR LE CALCUL DES EFFORTS AUX BORDS

! INEEDF : = 1 ON CALCULE LES EFFORTS AUX PAROIS
!          = 0 ON NE CALCULE PAS LES EFFORTS AUX PAROIS
integer ineedf
common / iforbr / ineedf

! FIN

