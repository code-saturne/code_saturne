!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
!> \file atmcls.f90
!> \brief Compute u*, q0, e0, (momentum, sensible heat and latent heat fluxes)
!>   for a non neutral atmospheric surface layer using the explicit formula
!>   developed for the ECMWF by Louis (1982)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
! mode               name       role                                           !
!______________________________________________________________________________!
!> \param[in]     ifac      treated boundary face
!> \param[in]     iel       boundary cell
!> \param[in]     utau      tangentiel mean
!> \param[in]     yplus     adim distance to he boundary faces
!> \param[out]    uet       friction velocity
!> \param[out]    gredu     reduced gravity for non horizontal wall
!> \param[out]    cfnnu     non neutral correction coefficients for profiles of wind
!> \param[out]    cfnns     non neutral correction coefficients for profiles of scalar
!> \param[out]    cfnnk     non neutral correction coefficients for profiles of k
!> \param[out]    cfnne     non neutral correction coefficients for profiles of eps
!> \param[in]     icodcl        code for boundary conditions at boundary faces
!>                              (nfabor,nvar)
!>-                           = 1   -> dirichlet
!>-                           = 3   -> densite de flux
!>-                           = 4   -> glissemt et u.n=0 (vitesse)
!>-                           = 5   -> frottemt et u.n=0 (vitesse)
!>-                           = 6   -> rugosite et u.n=0 (vitesse)
!>-                           = 9   -> entree/sortie libre (vitesse
!>                                      entrante eventuelle     bloquee
!> \param[in]     rcodcl         valeur des conditions aux limites
!>                                    (nfabor,nvar) aux faces de bord
!>-                           rcodcl(1) = valeur du dirichlet
!>-                           rcodcl(2) = valeur du coef. d'echange
!>                              ext. (infinie si pas d'echange)
!>-                           rcodcl(3) = valeur de la densite de
!>           flux (negatif si gain) w/m2 ou hauteur de rugosite (m) si icodcl=6
!>--                   pour les vitesses      (vistl+visct)*gradu
!>--                   pour la pression       dt*gradp
!>--                   pour les scalaires     cp*(viscls+visct/turb_schmidt)*gradt
! ______________________________________________________________________________!

subroutine atmcls &
 ( ifac   , iel    ,                                              &
   utau   , yplus  ,                                              &
   uet    ,                                                       &
   gredu  ,                                                       &
   cfnnu ,  cfnns  , cfnnk  , cfnne  ,                            &
   icodcl ,                                                       &
   rcodcl )

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use dimens, only: nvar
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use atincl

!===============================================================================

implicit none

! Arguments

integer          ifac   , iel

integer          icodcl(nfabor,nvar)

double precision utau, yplus, uet
double precision gredu
double precision cfnnu, cfnns, cfnnk,cfnne

double precision rcodcl(nfabor,nvar,3)

! Local variables

double precision rib, lmo, q0, e0
double precision tpot1,tpot2,tpotv1,tpotv2
double precision rscp1,rscp2
double precision actu,actt,b,c,d
double precision fm,fh,fmden1,fmden2,fhden
double precision rugd,rugt,distbf
double precision, dimension(:), pointer :: cvar_totwt, cvar_t

!===============================================================================

!===============================================================================
! 1. Initialisations
!===============================================================================

b = 5.d0
c = 5.d0
d = 5.d0

rugd = rcodcl(ifac,iu,3)
distbf = yplus*rugd
rugt = rcodcl(ifac,iv,3)
! 1/U+
actu = xkappa/log((distbf+rugd)/rugd)
actt = xkappa/log((distbf+rugt)/rugt)

if (ippmod(iatmos).eq.2) then
  call field_get_val_s(ivarfl(isca(itotwt)), cvar_totwt)
endif
call field_get_val_s(ivarfl(isca(iscalt)), cvar_t)

! prise en compte de l'humidite dans le rapport r/cp
if (ippmod(iatmos).eq.2) then
  rscp1 = (rair/cp0)*(1.d0 + (rvsra-cpvcpa)*rcodcl(ifac, isca(itotwt),1))
  rscp2 = (rair/cp0)*(1.d0 + (rvsra-cpvcpa)*cvar_totwt(iel))
else
  rscp1 = rair/cp0
  rscp2 = rair/cp0
endif

tpot1 = rcodcl(ifac,isca(iscalt),1)
tpot2 = cvar_t(iel)

!     .........    ............................................
!     3.2 - compute virtual potential temperature at two levels
!     .........    ............................................

if (ippmod(iatmos).eq.2) then
  tpotv1 = tpot1*(1.d0 + (rvsra - 1.d0)* rcodcl(ifac, isca(itotwt),1))
  tpotv2 = tpot2*(1.d0 + (rvsra - 1.d0)* cvar_totwt(iel))
else
  tpotv1 = tpot1
  tpotv2 = tpot2
endif

! Patch for the initial time step when thermal field is not initalized
if (ntcabs.eq.1) tpotv2 = tpotv1

!     .........    .....................................
!     3.3 - compute layer average Richardson number
!     .........    .....................................

! NB: rib = 0 if thermal flux conditions are imposed and tpot1 not defined
if (abs(utau).le.epzero.or.icodcl(ifac,isca(iscalt)).eq.3) then
 rib = 0.d0
else
 rib = 2.d0*gredu*distbf*(tpotv2 - tpotv1)/(tpotv1 + tpotv2)/utau/utau
endif

!     .........    ..................................................
!     3.4 - compute correction factors based on ECMWF parametrisation
!         Louis (1982)
!     ...............................................................

if (rib.ge.epzero) then
  ! Stable case
  fm = 1.d0/(1.d0 + 2.d0*b*rib/sqrt(1.d0 + d*rib))
  fh = 1.d0/(1.d0 + 3.d0*b*rib*sqrt(1.d0 + d*rib))

else
  ! Unstable case
  fmden1 = (distbf + rugt)*abs(rib)/rugt
  fmden2 = 1.d0 + 3.d0*b*c*actu*actt*sqrt(fmden1)
  fm = 1.d0 - 2.d0*b*rib/fmden2
  fhden = 3.d0*b*c*actu*actt*sqrt((distbf + rugt)/rugt)
  fh = 1.d0 - (3.d0*b*rib)/(1.d0 + fhden*sqrt(abs(rib)))
endif

if (fm.le.epzero) fm = epzero
if (abs(fh).le.epzero) fh = epzero

cfnnu = 1.d0/sqrt(fm)
cfnns = fh/sqrt(fm)
if ((1.d0-rib).gt.epzero)then
  cfnnk = sqrt(1.d0 - rib)  ! +correction with turbulent Prandtl
  cfnne = (1.d0-rib)/sqrt(fm)
else
  cfnnk = 1.d0
  cfnne = 1.d0
endif

!     ------------------------------------
!     4 - compute friction velocity  uet
!     ------------------------------------

uet = actu * utau * sqrt(fm)

!     -----------------------------------------
!     5 - compute surface sensible heat flux q0
!     -----------------------------------------

q0 = (tpot1-tpot2) * uet * actt * fh / sqrt(fm)

!     -----------------------------------------------
!     6 - compute Monin-Obukhov length (Garratt p 38)
!     -----------------------------------------------

if (abs(gredu*q0).le.epzero) then
  lmo = -99999.d0
else
  lmo = -uet**3*(t0 + tkelvi)/(xkappa*abs(gredu)*q0)
endif

!     ---------------------------------------
!     7 - compute surface latent heat flux e0
!     ---------------------------------------
!
!  if ( ippmod(iatmos).eq.2 ) then
!    e0 = (qvs(ifac)-qv(iel)) * uet * actt * fh / sqrt(fm)
!  endif

!----
! End
!----

return
end subroutine
