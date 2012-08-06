!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

subroutine iniusi
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


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use albase
use mltgrd
use parall
use period
use ihmpre
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use radiat
use cs_coal_incl

!===============================================================================

implicit none

! Arguments

! Local variables

integer          ii, iscal , nmodpp
integer          nscmax, nesmax, nscusi
integer          ieepre, ieeder, ieecor, ieetot, iihmpu
integer          ialgce, icwfps
integer          iappel, ioptit, ioplsq
double precision relaxp, extrap, cwfthr

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
! 1. INITIALISATION DE PARAMETRES POUR LA PHASE CONTINUE
!===============================================================================

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

iihmpu = iihmpr
call usipph(iihmpu , nfecra , iturb  , irccor , icp)
!==========

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
call usinsc(iihmpu , nfecra , nscaus)
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
   icolwc, iccoal, icpl3c, icfuel,                                &
   ieljou, ielarc, ielion, icompf, iatmos,                        &
   iaeros, ieos  , ieqco2)

endif

!   - Sous-programme utilisateur
!     ==========================

call usppmo
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

! --- Varpos
!     Verification et construction de ISCAPP
!      1ier passage
call varpos(nmodpp)
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
 ( nscmax , nscusi , iihmpu , nfecra , iscavr , ivisls )

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

if (iihmpr.eq.1) then

  call csidtv(idtvar)
  !==========

  call csiphy(iphydr)
  !==========

  ! Mesh related options

  call uicwf
  !=========

endif

!   - Sous-programme utilisateur
!     ==========================

nesmax = nestmx
ieepre = iespre
ieeder = iesder
ieecor = iescor
ieetot = iestot
iihmpu = iihmpr
!     IALGCE permet de remplir la variable cs_glob_maillage_grd_cdg_cel dans
!       cs_maillage_grd.c, a travers la routine ALGCEN.
!     cs_glob_maillage_grd_cdg_cel est initialise a 0 dans cs_maillage_grd.c,
!       et on ne change sa valeur ici que si elle a vraiment ete touchee par
!       l'utilisateur (pour garder l'initialisation en un seul endroit).
!     Le blindage en erreur est dans cs_maillage_grd.c (erreur si IALGCE>1,
!       cs_glob_maillage_grd_cdg_cel inchange si IALGCE<0)
ialgce = -999
icwfps = 0     ! Set to 1 to postprocess cutting of warped faces
cwfthr = -1.d0 ! Threshold (in degrees) to triangulate warped faces if positive

call usipgl                                                       &
!==========
 ( nesmax ,                                                       &
   ieepre , ieeder , ieecor , ieetot ,                            &
   iihmpu , nfecra ,                                              &
   idtvar , ipucou , idilat , iphydr , ialgce , iescal ,          &
   icwfps,  cwfthr )

if (ialgce.ne.-999) call algcen(ialgce)
if (cwfthr.ge.0.d0) call setcwf(icwfps, cwfthr)

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

!      2ieme passage
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
             iomg, iphi, ifb, ial,                                &
             inusa,                                               &
             iale, iuma, ivma, iwma,                              &
             isca, iscapp)

!     Suite de calcul, relecture fichier auxiliaire, champ de vitesse figé

  call csisui(ntsuit, ileaux, iccvfg)
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
         cdtvar, nitmax, epsilo, iresol, imgr, nswrsm)

!     Options numériques globales
  relaxp = -999.d0
  extrap = 0.d0
  call csnum2 (ivisse, relaxp, ipucou, extrap, imrgra, nterup)
  !==========
  extrag(ipr) = extrap
  if (idtvar.ge.0) relaxv(ipr) = relaxp

!     Gravite, prop. phys
  call csphys                                                     &
  !==========
             (nmodpp,                                             &
              irovar, ivivar, icorio,                             &
              gx, gy, gz, omegax, omegay, omegaz ,                &
              ro0, viscl0, viscv0, cp0, t0, p0, xmasmr)

!     Scamin, scamax
  call cssca2(iscavr, scamin, scamax)
  !==========

  ! Diffusivites
  call cssca3(iscalt, iscsth, iscavr, visls0, t0, p0)
  !==========

!     Init turb (uref, almax) si necessaire (modele RANS)
  if (itytur.eq.2 .or. itytur.eq.3 .or.             &
      itytur.eq.5 .or. itytur.eq.6 .or.             &
      itytur.eq.7) then
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
             ivisma, idtvar, ipucou, iappel)

  call uimoyt (ndgmox, ntdmom, imoold, idfmom)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

call usipsu(nmodpp)
!==========

call clmopt(mltmmn, mltmgl, mltmmr, mltmst, mlttyp)
!==========

call indsui(isuite)
!==========

! Choose if the 3x3 dimensionless matrix cocg is computed for the iterative
! algorithm and the Least square method for ivelco = 1.
if (ivelco.eq.1) then
  if (imrgra.eq.0) then
    ioptit = 1
    ioplsq = 0
  elseif (imrgra.eq.1.or.imrgra.eq.2.or.imrgra.eq.3) then
    ioptit = 0
    ioplsq = 1
  elseif (imrgra.eq.4) then
    ioptit = 1
    ioplsq = 1
  endif
else
  ioptit = 0
  ioplsq = 0
endif
call comcoc(ioptit, ioplsq)

! --- Varpos
!      3ieme passage
call varpos(nmodpp)
!==========


!===============================================================================
! 5. INITIALISATION DE PARAMETRES UTILISATEUR (entree sorties)
!===============================================================================

! --- Entree-sorties


!   - Interface Code_Saturne
!     ======================


if(iihmpr.eq.1) then

    iappel = 1

    call uiprop                                                   &
    !==========
            (irom, iviscl, ivisct, ivisls, icour, ifour,          &
             ismago, iale, icp, iscalt, iscavr,                   &
             iprtot, ipppro, ipproc, icmome,                      &
             ipptx, ippty, ipptz, ippdt,                          &
             ivisma, idtvar, ipucou, iappel)

  do ii = 1,nvppmx
    call fcnmva (nomvar(ii), len(nomvar(ii)), ii)
    !==========
  enddo

  call csenso                                                     &
  !==========
     ( nvppmx, ncapt,  nthist, frhist, ntlist, iecaux,            &
       ipstdv, ipstyp, ipstcl, ipstft, ipstfo,                    &
       ichrvr, ilisvr, ihisvr, tplfmt, isca, iscapp,              &
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

call usipes(nmodpp)
!==========

!----
! FORMATS
!----


return
end subroutine
