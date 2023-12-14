!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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
!> Warning: source terms are treated differently from the way they are
!>          in cs_user_source_terms.
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
!______________________________________________________________________________!

subroutine pptssc     &
 ( iscal  ,           &
   smbrs  , rovsdt )  &
  bind(C, name='cs_physical_model_source_terms')

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
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

integer(c_int), value :: iscal

real(c_double) :: smbrs(ncelet), rovsdt(ncelet)

! Local variables

integer           ivar, f_id

!===============================================================================

! Modele de la flamme de diffusion: steady laminar flamelet

if (ippmod(islfm).ge.0) then
  call cs_steady_laminar_flamelet_source_terms(iscal, smbrs, rovsdt)
endif

! Soot model

if (isoot.ge.1) then
  ivar = isca(iscal)
  if (ivar.eq.isca(ifsm).or.ivar.eq.isca(inpm)) then
    ! Scalar f_id
    f_id = ivarfl(ivar)
    call cs_soot_production(f_id, smbrs, rovsdt)
  endif
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

if (ippmod(iccoal).ge.0) then
  call cs_coal_scast(iscal, smbrs, rovsdt)
endif

! ---> Atmospheric version

if (ippmod(iatmos).ge.0) then
  call attssc(iscal,smbrs)
endif

!----
! End
!----

return

end subroutine
