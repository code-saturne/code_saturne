!-------------------------------------------------------------------------------

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

subroutine iniusi(iverif)
!================

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE APPELANT LES ROUTINES UTILISATEUR POUR L'ENTREE DES
!   PARAMETRES DE CALCUL : ON PASSE ICI POUR TOUT CALCUL

! CETTE ROUTINE PERMET DE CACHER A L'UTILISATEUR LES APPELS
!   A VARPOS ET AU LECTEUR XML DE L'IHM

! LE DECOUPAGE DE L'ANCIEN USINI1 PERMET EGALEMENT DE MIEUX
!   CONTROLER LES ZONES OU SONT INITIALISES LES VARIABLES (PAR
!   LE BIAIS DE PARAMETRES PASSES EN ARGUMENT)

! ON INITIALISE EGALEMENT ICI NPHAS ET ISCAPH QUI VALENT 1 EN
!   PRATIQUE DANS TOUS LES CALCULS (C'EST UN PARAMETRE UTILISATEUR,
!   MAIS SA VALEUR EST QUELQUE PEU CONTRAINTE ENCORE...)


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
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
include "cstnum.h"
include "dimens.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "vector.h"
include "albase.h"
include "parall.h"
include "period.h"
include "ihmpre.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"
include "ppcpfu.h"
include "radiat.h"

!===============================================================================

! Arguments

integer          iverif

! VARIABLES LOCALES

integer          ii, iphas , iscal , nmodpp
integer          nphmax, nscmax, nesmax, nphusi, nscusi
integer          ieepre, ieeder, ieecor, ieetot, iihmpu
integer          ialgce, imgrpr
integer          iappel
double precision relaxp, extrap

!===============================================================================

!===============================================================================
! 0. INITIALISATION DE L'INFORMATION "FICHIER XML (IHM) REQUIS & EXISTE"
!===============================================================================


!   - Interface Code_Saturne
!     ======================

!     Avec Xml, on regarde si le fichier a ete ouvert (requis et existe,
!       selon les tests realises dans cs_main)
!     IIHMPR a ete initialise a 0 juste avant (INIINI)

call csihmp(iihmpr)
!==========

if(iihmpr.eq.1) then

  call uiinit
  !==========

endif

!===============================================================================
! 1. INITIALISATION DE PARAMETRES DEPENDANT DU NOMBRE DE PHASES
!===============================================================================

! --- Nombre de phases

!     Egal a 1 en version 1.2
!     =======================

!     Le nombre de phases maximal est donne par NPHSMX dans paramx.h
!     On teste la valeur donnee par l'utilisateur avant de completer
!       les tableaux dimensionnes a NPHSMX (ex. ITURB(NPHSMX)).

nphas = 1


! --- Varpos
!     Verification du nombre de phases
!      1er passage
call varpos(nmodpp)
!==========

! --- Parametres dependant du nombre de phases

!     Turbulence
!     Chaleur massique variable ou non

!   - Interface Code_Saturne
!     ======================

if(iihmpr.eq.1) then

  call csturb(iturb, ideuch, igrake, igrari, xlomlg)
  !==========

  call cscpva(icp)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

nphmax = nphsmx
nphusi = nphas
iihmpu = iihmpr
call usipph                                                       &
!==========
 (nphmax , nphusi , iihmpu , nfecra , iturb  , icp , iverif)

!===============================================================================
! 2. INITIALISATION DE PARAMETRES DEPENDANT DU NOMBRE DE SCALAIRES
!===============================================================================

! --- Nombre de scalaires utilisateurs


!   - Interface Code_Saturne
!     ======================

if(iihmpr.eq.1) then

  call csnsca(nscaus)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

iihmpu = iihmpr
call usinsc(iihmpu , nfecra , nscaus , iverif)
!==========


! --- Dans le cas de physiques particulieres definies par des modules
!        specifiques du code tels que charbon, combustion, electrique
!        le sous-programme USPPMO doit etre complete imperativement
!        par l'utilisateur

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  call uippmo                                                     &
  !==========
 ( ippmod, icod3p, icodeq, icoebu, icobml,                        &
   icolwc, icp3pl, icpl3c, icfuel,                                &
   ieljou, ielarc, ielion, icompf, iatmos,                        &
   iaeros, indjon, ieqco2)

endif

!   - Sous-programme utilisateur
!     ==========================

call usppmo
!==========

! --- Varpos
!     Verification et construction de ISCAPP
!      2ieme passage
call varpos(nmodpp)
!==========


! --- Activation du module transferts radiatifs

!     Il est necessaire de connaitre l'activation du module transferts
!     radiatifs tres tot de maniere a pouvoir reserver less variables
!     necessaires dans certaines physiques particuliere

!   - Interface Code_Saturne
!     ======================

if(iihmpr.eq.1) then

  call uiray1(iirayo, isuird, ndirec, nfreqr, idiver, iimpar, iimlum)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

call usray1
!==========

! --- Parametres dependant du nombre de scalaires utilisateurs

!     Moyenne du carre des fluctuations d'un scalaire UTILISATEUR
!     Diffusivite variable ou non


!   - Interface Code_Saturne
!     ======================

if(iihmpr.eq.1) then

  call  csisca(iscavr)
!       ============

  call  csivis(iscavr, ivisls, iscalt, iscsth)
!       ============

endif

!   - Sous-programme utilisateur
!     ==========================

nscmax = nscamx
nscusi = nscaus
iihmpu = iihmpr
call usipsc                                                       &
!==========
 ( nscmax , nscusi , iihmpu , nfecra , iscavr , ivisls , iverif )


! --- Parametres dependant du nombre de scalaires : exception
!     IPHSCA indique le numero de la phase porteuse pour chaque
!                                                    scalaire UTILISATEUR

!       Dans le cas d'une seule phase, IPHSCA(ISCAL) = 1 : ne rien changer.
!       ==================================================================

!       Sinon noter bien que :

!       La phase porteuse des scalaires UTILISATEUR ISCAL qui
!         representent la moyenne du carre des fluctuations d'un
!         scalaire utilisateur K sera la meme que celle de ce scalaire
!         utilisateur K.
!         Donc, pour de tels scalaires ISCAL (reperes par ISCAVR(ISCAL)>0),
!                          on ne doit pas renseigner IPHSCA(ISCAL) ici.
!         C'est l'objet du test inclus dans l'exemple ci-dessous.

!       Pour les scalaires non utilisateur relatifs a des physiques
!         particulieres, (charbon, combustion, electrique : voir usppmo)
!         implicitement definis selon le modele,
!         les informations sont donnees automatiquement par ailleurs :
!                                         on ne modifie pas IPHSCA ici.

do iscal = 1, nscaus
  if(iscavr(iscal).le.0) then
    iphsca(iscal) = 1
  endif
enddo


!===============================================================================
! 3. INITIALISATION DE PARAMETRES "GLOBAUX"
!===============================================================================


! --- Parametres globaux

!     Pas de temps
!     Couplage vitesse/pression
!     Prise en compte de la pression hydrostatique
!     Estimateurs (pas encore dans l'IHM)


!   - Interface Code_Saturne
!     ======================

if(iihmpr.eq.1) then

  call csidtv(idtvar)
  !==========

  call csiphy(iphydr)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

nphmax = nphsmx
nesmax = nestmx
ieepre = iespre
ieeder = iesder
ieecor = iescor
ieetot = iestot
nphusi = nphas
iihmpu = iihmpr
!     IALGCE permet de remplir la variable cs_glob_maillage_grd_cdg_cel dans
!       cs_maillage_grd.c, a travers la routine ALGCEN.
!     cs_glob_maillage_grd_cdg_cel est initialise a 0 dans cs_maillage_grd.c,
!       et on ne change sa valeur ici que si elle a vraiment ete touchee par
!       l'utilisateur (pour garder l'initialisation en un seul endroit).
!     Le blindage en erreur est dans cs_maillage_grd.c (erreur si IALGCE>1,
!       cs_glob_maillage_grd_cdg_cel inchange si IALGCE<0)
ialgce = -999

call usipgl                                                       &
!==========
 ( nphmax , nesmax ,                                              &
   ieepre , ieeder , ieecor , ieetot ,                            &
   nphusi , iihmpu , nfecra ,                                     &
   idtvar , ipucou , iphydr , ialgce , iescal , iverif )

if (ialgce.ne.-999) call algcen(ialgce)

! --- Parametres de la methode ALE

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  call uialin (iale, nalinf, nalimx, epalim, iortvm)
  !==========
endif

!   - Sous-programme utilisateur
!     ==========================

call usalin
!==========

! --- Varpos
!     Positionnement de pointeurs
!     Verifications
!     Determination de IPR, IU ... ISCA, NVAR
!     Determination de IPP...

!      3ieme passage
call varpos(nmodpp)
!==========


!===============================================================================
! 4. INITIALISATION DE PARAMETRES UTILISATEUR SUPPLEMENTAIRES
!===============================================================================

! --- Format des fichiers aval (entsor.h)
! --- Options du calcul (optcal.h)
! --- Constantes physiques (cstphy.h)


!   - Interface Code_Saturne
!     ======================

if(iihmpr.eq.1) then

  call csvnum                                                     &
  !==========
            (nvar,                                                &
             iu, iv, iw, ipr,                                     &
             iturb, ik, iep,                                      &
             ir11, ir22, ir33,                                    &
             ir12, ir13, ir23,                                    &
             iomg, iphi, ifb,                                     &
             iale, iuma, ivma, iwma,                              &
             isca, iscapp)

!     Suite de calcul, relecture fichier auxiliaire, champ de vitesse figé
  call csisui(isuite, ileaux, iccvfg)
  !==========

!     Pas de temps (seulement NTMABS, DTREF, INPDT0)
  call cstime                                                     &
  !==========
             (inpdt0, iptlro, ntmabs, idtvar, dtref, dtmin,       &
              dtmax, coumax, foumax, varrdt, relxst)

!     Temperature ou enthalpie (hors physiques particulieres)
  if(nmodpp.eq.0) then
    call cssca1(iscalt, iscsth)
    !==========

  endif

!      Options numériques locales

  call uinum1                                                     &
  !==========
        (isca, iscapp, blencv, ischcv, isstpc, ircflu,            &
         cdtvar, nitmax, epsilo)

!     Options numériques globales
  relaxp = -999.d0
  extrap = 0.d0
  imgrpr = 0
  call csnum2 (ivisse, relaxp, ipucou, extrap, imrgra, imgrpr)
  !==========
  iphas = 1
  extrag(ipr(iphas)) = extrap
  if (idtvar.ge.0) relaxv(ipr(iphas)) = relaxp
  imgr(ipr(iphas)) = imgrpr

!     Gravite, prop. phys
  call csphys                                                     &
  !==========
             (nmodpp,                                             &
              irovar, ivivar, icorio,                             &
              gx, gy, gz, omegax, omegay, omegaz ,                &
              ro0, viscl0, cp0, t0, p0)

!     Scamin, scamax
  call cssca2(iscavr, scamin, scamax)
  !==========

!     Diffusivites
  call cssca3(iscalt, iscavr, visls0, t0, p0)
  !==========

!     Init turb (uref, almax) si necessaire (modele RANS)
  iphas = 1
  if (itytur(iphas).eq.2 .or. itytur(iphas).eq.3 .or.             &
      itytur(iphas).eq.5 .or. itytur(iphas).eq.6 ) then
    call cstini(uref, almax)
    !==========
  endif

  iappel = 0

  call uiprop                                                     &
  !==========
            (irom, iviscl, ivisct, ivisls, icour, ifour,          &
             ismago, iale, icp, iscalt, iscavr,                   &
             iprtot, ipppro, ipproc, icmome,                      &
             ipptx, ippty, ipptz, ippdt,                          &
             ivisma, iappel)

  call uimoyt (ndgmox, ntdmom, imoold, idfmom)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

call usipsu(nmodpp , iverif)
!==========


! --- Varpos
!      4ieme passage
call varpos(nmodpp)
!==========


!===============================================================================
! 5. INITIALISATION DE PARAMETRES UTILISATEUR (entree sorties)
!===============================================================================

! --- Entree-sorties


!   - Interface Code_Saturne
!     ======================


if(iihmpr.eq.1) then

  if(nbmomt.gt.0 .or. ipucou.eq.1) then
    iappel = 1

    call uiprop                                                   &
    !==========
            (irom, iviscl, ivisct, ivisls, icour, ifour,          &
             ismago, iale, icp, iscalt, iscavr,                   &
             iprtot, ipppro, ipproc, icmome,                      &
             ipptx, ippty, ipptz, ippdt,                          &
             ivisma, iappel)

  endif

  do ii = 1,nvppmx
    call fcnmva (nomvar(ii), len(nomvar(ii)), ii)
    !==========
  enddo

  call csenso                                                     &
  !==========
     ( nvppmx, ncapt,  nthist, ntlist,                            &
       ichrvl, ichrbo, ichrsy, ichrmd,                            &
       fmtchr, len(fmtchr), optchr, len(optchr),                  &
       ntchr,  iecaux,                                            &
       ipstdv, ipstyp, ipstcl, ipstft, ipstfo,                    &
       ichrvr, ilisvr, ihisvr, isca, iscapp,                      &
       ipprtp, xyzcap )

  do ii = 1,nvppmx
    call cfnmva(nomvar(ii), len(nomvar(ii)), ii)
    !==========
  enddo

  call nvamem
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

call usipes(nmodpp, iverif)
!==========

!----
! FORMATS
!----


return
end
