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

subroutine cs_coal_varpos
!========================

!===============================================================================
!  FONCTION  :
!  ---------
!       INIT DES POSITIONS DES VARIABLES TRANSPORTEES POUR
!                COMBUSTION CHARBON PULVERISE
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
use ppincl
use ppcpfu
use ihmpre
use cs_coal_incl

!===============================================================================

implicit none

integer        icla,  is, icha, isc , is1

!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================

! ---> Variables propres a la suspension gaz - particules
is = 1
ihm   = iscapp(is)

! ---> Variables propres a la phase dispersee

do icla = 1, nclacp
  is = 1+icla
  inp(icla) = iscapp(is)
  is = 1+1*nclacp+icla
  ixch(icla)= iscapp(is)
  is = 1+2*nclacp+icla
  ixck(icla) = iscapp(is)
  if ( ippmod(iccoal) .eq. 1 ) then
    is = 1+3*nclacp+icla
    ixwt(icla) = iscapp(is)
    is = 1+4*nclacp+icla
    ih2(icla) = iscapp(is)
  else
    is = 1+3*nclacp+icla
    ih2(icla) = iscapp(is)
  endif
enddo
! ---> Variables propres a la phase continue
!
! Matiere volatiles legeres (F8) et lourdes (F9)
is1 = is
do icha = 1, ncharb
  is          = is1+icha
  if1m(icha)  = iscapp(is)
  is          = is1+ncharb+icha
  if2m(icha)  = iscapp(is)
enddo
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
! Humidite
if (ippmod(iccoal) .eq. 1 ) then
  is = is+1
  if6m  = iscapp(is)
endif
! Produits de la combustion du coke par O2
is = is+1
if7m  = iscapp(is)
! Produits de la combustion du coke par CO2
if ( ihtco2 .eq. 1 ) then
  is = is+1
  if8m  = iscapp(is)
endif
! Produits de la combustion du coke par H2O
if ( ihth2o .eq. 1 ) then
  is = is+1
  if9m  = iscapp(is)
endif
! variance
is = is+1
ifvp2m  = iscapp(is)
! Transport du CO ou du CO2
if ( ieqco2 .ge. 1 ) then
  is    = is+1
  iyco2 = iscapp(is)
endif
! Transport du NOx : HCN, NOx et Tair
if ( ieqnox .eq. 1 ) then
  is     = is+1
  iyhcn  = iscapp(is)
  is     = is+1
  iyno   = iscapp(is)
  is     = is+1
  ihox   = iscapp(is)
endif
!   - Interface Code_Saturne
!     ======================
!     Construction de l'indirection entre la numerotation du noyau et XML

if (iihmpr.eq.1) then
!
  call uicpsc (ncharb, nclacp, noxyd, ippmod,               &
               iccoal, ieqnox, ieqco2, ihtco2,              &
               ihth2o, ihm, inp, ixch, ixck, ixwt, ih2,     &
               if1m, if2m, if4m, if5m, if6m,                &
               if7m, if8m, ifvp2m, iyco2, if9m,     &
               iyhcn, iyno, ihox )
endif

!===============================================================================
! 2. PROPRIETES PHYSIQUES
!    A RENSEIGNER OBLIGATOIREMENT (sinon pb dans varpos)
!    - PROPRES AUX SCALAIRES   : IVISLS, ISCAVR
!      Rq : pas de variance associee a un scalaire dans notre cas
!    - PROPRES A LA SUSPENSION : ICP
!===============================================================================

do isc = 1, nscapp

  if ( iscavr(iscapp(isc)).le.0 ) then

! ---- Viscosite dynamique de reference relative au scalaire
!      ISCAPP(ISC)
    ivisls(iscapp(isc)) = 0

  endif

enddo

! ---- Bien que l'on soit en enthalpie on conserve un CP constant

icp    = 0

return
end subroutine
