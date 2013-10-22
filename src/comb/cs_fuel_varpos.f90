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

subroutine cs_fuel_varpos
!========================
!===============================================================================
!  FONCTION  :
!  ---------

!       INIT DES POSITIONS DES VARIABLES TRANSPORTEES POUR
!                COMBUSTION FUEL

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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use cs_fuel_incl
use ppincl
use ppcpfu

!===============================================================================

implicit none

integer        is, isc, icla

!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================

! ---> Variables propres a la phase dispersee

do icla = 1, nclafu
  is          = 1 + icla
  ing(icla)   = iscapp(is)
  is          = 1 + 1*nclafu+icla
  iyfol(icla) = iscapp(is)
  is          = 1 + 2*nclafu+icla
  ih2(icla)   = iscapp(is)
enddo

! ---> Variables propres a la phase gaz

is       = is + 1
ifvap    = iscapp(is)
! Oxydant 2
if ( noxyd .ge. 2 ) then
  is = is+1
  if4m  = iscapp(is)
endif
! Oxydant 3
if ( noxyd .ge. 3 ) then
  is = is+1
  if5m  = iscapp(is)
endif
! combustion heterogene
is       = is + 1
if7m   = iscapp(is)
! Variance
is       = is + 1
ifvp2m  = iscapp(is)
! Transport du CO2
if ( ieqco2 .ge. 1 ) then
  is = is+1
  iyco2 = iscapp(is)
endif
!  Transport du NOx : HCN, NOx et Hox
if ( ieqnox .eq. 1 ) then
  is = is+1
  iyhcn = iscapp(is)
  is = is+1
  iyno  = iscapp(is)
  is = is+1
  ihox = iscapp(is)
endif

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!    - PROPRES AUX SCALAIRES   : IVISLS, ISCAVR
!      Rq : pas de variance associee a un scalaire dans notre cas
!    - PROPRES A LA SUSPENSION : ICP
!===============================================================================

do isc = 1, nscapp

  if ( iscavr(iscapp(isc)) .le. 0 ) then

! ---- Viscosite dynamique de reference relative au scalaire
!      ISCAPP(ISC)
    ivisls(iscapp(isc)) = 0

  endif

enddo

! ---- Bien que l on soit en enthalpie on conserve un CP constant

icp    = 0

!----
! End
!----

return
end subroutine
