!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

integer          nmodpp, ifcvsl
integer          nscmax, nscusi
integer          iihmpu
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

  call csturb(iturb, iwallf, igrake, igrari, xlomlg)
  !==========

  call cscpva(icp)
  !==========

endif

! ALE parameters
!---------------

! GUI

if (iihmpr.eq.1) then
  call uialin (iale, nalinf, nalimx, epalim, iortvm)
endif

! User sub-routines
! =================

iihmpu = iihmpr
call usipph(iihmpu, iturb, itherm, iale, icavit)

! Other model parameters, including user-defined scalars
!-------------------------------------------------------

! GUI

if (iihmpr.eq.1) then
  call cs_gui_user_variables
endif

! User sub-routines

call cs_user_model

!===============================================================================
! 2. Initialize parameters for specific physics
!===============================================================================

! GUI

if (iihmpr.eq.1) then

  call uippmo                                                     &
  !==========
 ( ippmod, icod3p, icodeq, icoebu, icobml,                        &
   icolwc, iccoal, icpl3c, icfuel,                                &
   ieljou, ielarc, ielion, icompf, iatmos,                        &
   iaeros, ieos  , ieqco2, idarcy)

  call cfnmtd(ficfpp, len(ficfpp))
  !==========

endif

! --- Activation du module transferts radiatifs

!     Il est necessaire de connaitre l'activation du module transferts
!     radiatifs tres tot de maniere a pouvoir reserver les variables
!     necessaires dans certaines physiques particuliere

!   - Interface Code_Saturne
!     ======================

if (iihmpr.eq.1) then

  call uiray1(iirayo, isuird, i_quadrature, ndirec, nfreqr, idiver, iimpar, iimlum)

endif

! User subroutine

! Initialize specific physics modules not available at the moment

ippmod(icobml) = -1  ! premix model of Bray - Moss - Libby
ippmod(icodeq) = -1  ! diffusion flame with fast equilibrium chemistry
ippmod(ielion) = -1  ! ionic mobility

! User initialization

iihmpu = iihmpr
call usppmo(iihmpu)
!==========

! Define fields for variables, check and build iscapp

call fldvar(nmodpp)
!==========

if (iihmpr.eq.1) then
  call csivis
  !==========
endif

nscmax = nscamx
nscusi = nscaus
iihmpu = iihmpr

if (ippmod(icompf).ge.0) then

  ! For the compressible model, call to uscfx1 to get ieos as well as
  ! scalar_diffusivity_id for itempk and iviscv.
  ! With the GUI, ieos has been read above in the call to uippmo, and
  ! iviscv has been read below in the call to fldvar (csvvva).

  call uscfx1
  !==========
  ! Dynamic viscosity of reference of the scalar total energy (ienerg).
  call field_get_key_int(ivarfl(isca(itempk)), kivisl, ifcvsl)
  if (ifcvsl.ge.0 .or. icv.gt.0) then
    call field_set_key_int(ivarfl(isca(ienerg)), kivisl, 0)
  else
    call field_set_key_int(ivarfl(isca(ienerg)), kivisl, -1)
  endif
endif


! ---> Physique particuliere : darcy

if (ippmod(idarcy).ge.0) then
  call daini1
  !==========
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
         cdtvar, epsilo, nswrsm)

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
              gx, gy, gz,                                             &
              ro0, viscl0, viscv0, visls0, cp0, t0, p0,               &
              xmasmr, itempk)

!     Scamin, scamax, turbulent flux model
  call cssca2(itytur, iturt)
  !==========

  ! Diffusivites
  call cssca3(visls0, t0, p0, cp0)
  !==========

!     Init turb (uref, almax) si necessaire (modele RANS)
  if (itytur.eq.2 .or. itytur.eq.3 .or.             &
      itytur.eq.5 .or. itytur.eq.6 .or.             &
      itytur.eq.7) then
    call cstini(uref, almax)
    !==========
  endif

  call uiipsu(iporos)
  !==========

endif

!   - Sous-programme utilisateur
!     ==========================

call usipsu(nmodpp)
!==========

! If time step is local or variable, pass information to C layer, as it
! may be needed for some field (or moment) definitions.
if (idtvar.ne.0) then
  call time_step_define_variable(1)
endif
if (idtvar.eq.2.or.idtvar.eq.-1) then
  call time_step_define_local(1)
endif

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

! Choose the porous model
call compor(iporos)

! --- Varpos
call varpos
!==========

!----
! Formats
!----

return
end subroutine
