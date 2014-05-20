!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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
use cs_c_bindings
use field

!===============================================================================

implicit none

! Arguments

! Local variables

integer          n_fields, nmodpp
integer          nscmax, nesmax, nscusi
integer          ieepre, ieeder, ieecor, ieetot, iihmpu
integer          ialgce
integer          iappel
double precision relaxp, extrap

!===============================================================================

! Check for restart and read matching time steps

call parameters_read_restart_info

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

if (iihmpr.eq.1) then

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

if (iihmpr.eq.1) then

  call csther(itherm, itpscl)
  !==========

  call csturb(iturb, iturt, ideuch, igrake, igrari, xlomlg)
  !==========

  call cscpva(icp)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

iihmpu = iihmpr
call usipph(iihmpu , nfecra , iturb , irccor , idirsm , itherm, icp, icavit)

! --- ALE parameters

!   - Code_Saturne GUI
!     ================

if (iihmpr.eq.1) then
  call uialin (iale, nalinf, nalimx, epalim, iortvm)
endif

!   - User sub-routines
!     =================

call usalin


!===============================================================================
! 2. INITIALISATION DE PARAMETRES DEPENDANT DU NOMBRE DE SCALAIRES
!===============================================================================

! --- Nombre de scalaires utilisateurs

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

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

  call cfnmtd(ficfpp, len(ficfpp))
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

! Initialize specific physics modules not available at the moment

ippmod(icobml) = -1  ! premix model of Bray - Moss - Libby
ippmod(icodeq) = -1  ! diffusion flame with fast equilibrium chemistry
ippmod(ielion) = -1  ! ionic mobility

! User initialization

iihmpu = iihmpr
call usppmo(iihmpu)
!==========

! --- Activation du module transferts radiatifs

!     Il est necessaire de connaitre l'activation du module transferts
!     radiatifs tres tot de maniere a pouvoir reserver less variables
!     necessaires dans certaines physiques particuliere

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  call uiray1(iirayo, isuird, i_quadrature, ndirec, nfreqr, idiver, iimpar, iimlum)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

call usray1
!==========

!     Define fields for variables, check and build iscapp

call fldvar(nmodpp)
!==========

if (ippmod(icompf).ge.0) then

!     For compressible model, call to uscfx1 to get ieos.
!     With GUI, ieos has been read below in the call to uippmo.
  call uscfx1
  !==========

endif

! --- Parametres dependant du nombre de scalaires utilisateurs

!     Moyenne du carre des fluctuations d'un scalaire UTILISATEUR
!     Diffusivite variable ou non


!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  call csisca(iscavr, itherm)
  !==========

  call csivis(iscavr, ivisls, iscalt, itherm, itempk)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

nscmax = nscamx
nscusi = nscaus
iihmpu = iihmpr
call usipsc(nscmax , nscusi , iihmpu , nfecra , iscavr , ivisls)
!==========

if (ippmod(icompf).ge.0) then
!     For compressible model, call to uscfx1 to get ivisls(itempk) et iviscv.
!     With GUI, iviscv has been read below in the first call to varpos (csvvva)
!     and ivisl(itempk) below in the call to csivis.

  call uscfx1
  !==========
!     Dynamic viscosity of reference of the scalar total energy (ienerg).
  if(ivisls(itempk).gt.0 .or. icv.gt.0) then
    ivisls(ienerg) = 1
  else
    ivisls(ienerg) = 0
  endif
endif
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

call usipgl                                                       &
!==========
 ( nesmax ,                                                       &
   ieepre , ieeder , ieecor , ieetot ,                            &
   iihmpu , nfecra ,                                              &
   idtvar , ipucou , idilat , iphydr , ialgce , iescal )

! If time step is local or variable, pass information to C layer, as it
! may ne needed for some field (or moment) definitions.
if (idtvar.ne.0) then
  call time_step_define_variable(1)
endif
if (idtvar.eq.2.or.idtvar.eq.-1) then
  call time_step_define_local(1)
endif

if (ialgce.ne.-999) call algcen(ialgce)

! Define main properties (pointers, checks, ipp)

call fldprp
!==========

!===============================================================================
! 4. INITIALISATION DE PARAMETRES UTILISATEUR SUPPLEMENTAIRES
!===============================================================================

! --- Format des fichiers aval (entsor.h)
! --- Options du calcul (optcal.h)
! --- Constantes physiques (cstphy.h)


!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

!     Suite de calcul, relecture fichier auxiliaire, champ de vitesse figé

  call csisui(ntsuit, ileaux, iccvfg)
  !==========

!     Pas de temps (seulement NTMABS, DTREF, INPDT0)
  call cstime                                                     &
  !==========
             (inpdt0, iptlro, ntmabs, idtvar, dtref, dtmin,       &
              dtmax, coumax, foumax, varrdt, relxst)

!      Options numériques locales

  call uinum1                                                     &
  !==========
        (blencv, ischcv, isstpc, ircflu,                          &
         cdtvar, nitmax, epsilo, iresol, imgr, nswrsm)

!     Options numériques globales
  relaxp = -999.d0
  extrap = 0.d0
  call csnum2 (ivisse, relaxp, ipucou, extrap, imrgra, nterup)
  !==========
  extrag(ipr) = extrap
  if (idtvar.ge.0) relaxv(ipr) = relaxp

!     Gravite, prop. phys
  call csphys                                                         &
  !==========
             (nmodpp,                                                 &
              irovar, ivivar, icorio,                                 &
              gx, gy, gz, omegax, omegay, omegaz ,                    &
              ro0, viscl0, viscv0, visls0, cp0, t0,                   &
              p0, xmasmr, itempk, itherm, itpscl)

!     Scamin, scamax, turbulent flux model
  call cssca2(iscavr, itytur, iturt)
  !==========

  ! Diffusivites
  call cssca3(itherm, iscalt, iscavr, visls0, t0, p0, cp0)
  !==========

!     Init turb (uref, almax) si necessaire (modele RANS)
  if (itytur.eq.2 .or. itytur.eq.3 .or.             &
      itytur.eq.5 .or. itytur.eq.6 .or.             &
      itytur.eq.7) then
    call cstini(uref, almax)
    !==========
  endif

  iappel = 0

  call uiprop (ivisls, ismago, iale, icp, iscavr)
  !==========

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


if (ippmod(icompf).ge.0) then
!      For compressible model, call to uscfx2 to get visls0(itempk), viscv0,
!      xmasmr and ivivar
!      With GUI, visls0(itempk), viscv0, xmasmr and ivivar have been read
!      below in the call to csphys.
  call uscfx2
  !==========
endif

! Choose which 3x3 cocg matrixes are computed for gradient algorithms.
call comcoc(imrgra)

! --- Varpos
!      1er passage
call varpos
!==========

!===============================================================================
! 5. INITIALISATION DE PARAMETRES UTILISATEUR (entree sorties)
!===============================================================================

! --- Entree-sorties


!   - Interface Code_Saturne
!     ======================

call field_get_n_fields(n_fields)

if (iihmpr.eq.1) then

  iappel = 1

  call uiprop (ivisls, ismago, iale, icp, iscavr)
  !==========

  call csenso                                                   &
  !==========
     ( nvppmx, ncapt,  nthist, frhist, ntlist, iecaux,          &
       ipstdv, ihisvr, tplfmt, isca, iscapp,                    &
       ipprtp, xyzcap )

endif

!----
! Formats
!----


return
end subroutine
