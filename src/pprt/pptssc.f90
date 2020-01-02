!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!===============================================================================
! Function:
! ---------

!> \file pptssc.f90
!>
!> \brief This subroutine defines the source terms for scalars which are part of
!> specific physics models. Source terms are defined over one time step.
!>
!> Warning: source terms are treated differently from the way they are in ustssc.
!> rovsdt*d(var) = smbrs is solved. rovsdt and smbrs already hold possible user
!> source terms values and thus have to be incremented (and not overwritten).
!>
!> For stability reasons, only positive terms are added to rovsdt, while there
!> are no such constrains on values to be added to smbrs.
!>
!> In the case of a source term of the form cexp + cimp*var, the source term
!> should be implemented as follows:
!> \f[
!>   smbrs  = smbrs  + cexp + cimp*var
!> \f]
!> \f[
!>   rovsdt = rovsdt + max(-cimp,0)
!> \f]
!>
!> rovsdt and smbrs are provided here respectively in kg/s and in kg/s*[scalar].
!> Examples:
!>   velocity \f$ kg m/s^2 \f$
!>   temperature \f$ kg K/s \f$
!>   enthalpy \f$ J/s \f$
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iscal         scalar number
!> \param[in,out] smbrs         explicit source term part
!> \param[in,out] rovsdt        implicite source term part
!> \param[in]     tslagr        coupling term for the Lagrangian module
!______________________________________________________________________________!


subroutine pptssc &
 ( iscal  ,                                                         &
   smbrs  , rovsdt , tslagr )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe, only: itypfb
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use lagran
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          iscal

double precision smbrs(ncelet), rovsdt(ncelet)
double precision tslagr(ncelet,*)

! Local variables

!===============================================================================

! Soot model

if (isoot.eq.1) then
  call sootsc(iscal, smbrs, rovsdt)
endif

! ---> Flamme de premelange : Modele EBU

if (ippmod(icoebu).ge.0) then
  call ebutss(iscal, smbrs, rovsdt)
endif

! ---> Flamme de premelange : Modele BML

!      IF ( IPPMOD(ICOBML).GE.0 )
!           CALL BMLTSS

! ---> Flamme de premelange : Modele LWC

if (ippmod(icolwc).ge.0) then
  call lwctss(iscal, smbrs, rovsdt)
endif

! ---> Flamme charbon pulverise

if ( ippmod(iccoal).ge.0 ) then
  call cs_coal_scast(iscal, smbrs, rovsdt)
endif

! ---> Flamme charbon pulverise couplee Transport Lagrangien
!      des particules de charbon

if ( ippmod(icpl3c).ge.0 .and. iilagr.eq.2 ) then
  call cpltss(iscal, itypfb, smbrs, rovsdt, tslagr)
endif

! ---> Flamme fuel

if ( ippmod(icfuel).ge.0 ) then
  call cs_fuel_scast(iscal, smbrs, rovsdt)
endif

! ---> Versions electriques :
!             Effet Joule
!             Arc Electrique
!             Conduction ionique

if (ippmod(ieljou).ge.1 .or.                                      &
    ippmod(ielarc).ge.1       ) then
  call eltssc(iscal, smbrs)
endif

! ---> Version atmospherique :

if (ippmod(iatmos).ge.0) then
  call attssc(iscal,smbrs)
endif


! ---> Version aerorefrigerant :

if (ippmod(iaeros).ge.0) then
  call cs_ctwr_source_term(ivarfl(isca(iscal)), p0, molmass_rat, smbrs, rovsdt)
endif

!----
! End
!----

return

end subroutine
