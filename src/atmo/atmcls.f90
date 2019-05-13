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
!> \brief Compute friction velocity u* and surface sensible heat flux q0
!> for a non neutral atmospheric surface layer using the explicit formula
!> developed for the ECMWF by Louis (1982)
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
! mode               name       role                                           !
!______________________________________________________________________________!
!> \param[in]     ifac          treated boundary face
!> \param[in]     utau          tangential mean
!> \param[in]     yplus         adim distance to he boundary faces
!> \param[out]    uet           friction velocity
!> \param[out]    gredu         reduced gravity for non horizontal wall
!> \param[out]    cfnnu         non neutral correction coefficients for profiles of wind
!> \param[out]    cfnns         non neutral correction coefficients for profiles of scalar
!> \param[out]    cfnnk         non neutral correction coefficients for profiles of k
!> \param[out]    cfnne         non neutral correction coefficients for profiles of eps
!> \param[in]     temp          potential temperature in boundary cell
!> \param[in]     totwt         total water content in boundary cell
!> \param[in]     liqwt         liquid water content in boundary cell
!> \param[in]     icodcl        face boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 11 Boundary value related to the next cell
!>                                 value by an affine function
!>                               - 13 Dirichlet for the advection operator and
!>                                 Neumann for the diffusion operator
!> \param[in]     rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
! ______________________________________________________________________________!

subroutine atmcls &
 ( ifac   ,                                                       &
   utau   , yplus  ,                                              &
   uet    ,                                                       &
   gredu  ,                                                       &
   cfnnu ,  cfnns  , cfnnk  , cfnne  ,                            &
   temp  ,  totwt  , liqwt  ,                                     &
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

integer          ifac

integer          icodcl(nfabor,nvar)

double precision utau, yplus, uet
double precision gredu
double precision cfnnu, cfnns, cfnnk,cfnne
double precision temp, totwt, liqwt
double precision rcodcl(nfabor,nvar,3)

! Local variables

double precision rib, q0
double precision tpot1, tpot2, tpotv1, tpotv2
double precision rscp1, rscp2
double precision actu, actt, b, c, d
double precision fm, fh, fmden1, fmden2, fhden
double precision rugd, rugt, distbf

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

! Take into account humidity in ratio r/cp
if (ippmod(iatmos).eq.2) then
  rscp1 = (rair/cp0)*(1.d0 + (rvsra-cpvcpa) * rcodcl(ifac, isca(itotwt),1))
  ! Bouzerau PhD
  rscp2 = (rair/cp0)*(1.d0 + (rvsra-cpvcpa) * (totwt-liqwt))
else
  rscp1 = rair/cp0
  rscp2 = rair/cp0
endif

tpot1 = rcodcl(ifac,isca(iscalt),1)
tpot2 = temp

! Compute virtual potential temperature at two levels

if (ippmod(iatmos).eq.2) then
  tpotv1 = tpot1*(1.d0 + (rvsra - 1.d0) * rcodcl(ifac, isca(itotwt),1))
  ! Bouzerau PhD
  tpotv2 = tpot2*(1.d0 + (rvsra - 1.d0) * (totwt-liqwt))
else
  tpotv1 = tpot1
  tpotv2 = tpot2
endif

! Patch for the initial time step when thermal field is not initalized
if (ntcabs.eq.1) tpotv2 = tpotv1

! Compute layer average Richardson number

! NB: rib = 0 if thermal flux conditions are imposed and tpot1 not defined
if (abs(utau).le.epzero.or.icodcl(ifac,isca(iscalt)).eq.3) then
 rib = 0.d0
else
 rib = 2.d0*gredu*distbf*(tpotv2 - tpotv1)/(tpotv1 + tpotv2)/utau/utau
endif

! Compute correction factors based on ECMWF parametrisation
! Louis (1982)

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

! Compute friction velocity  uet
uet = actu * utau * sqrt(fm)

! Compute surface sensible heat flux q0
q0 = (tpot1-tpot2) * uet * actt * fh / sqrt(fm)

!----
! End
!----

return
end subroutine
