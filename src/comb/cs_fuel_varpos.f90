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
use field

!===============================================================================

implicit none

integer        isc, icla

integer          keyccl, kscmin, kscmax, f_id
character*80     f_label, f_name

!===============================================================================
! 0. Definitions for fields
!===============================================================================

! Key id of the fuel scalar class
call field_get_key_id("scalar_class", keyccl)

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

!===============================================================================
! 1. Definition of fields
!===============================================================================

! Thermal model

itherm = 2
call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
iscalt = ihm

! Set min and max clipping
f_id = ivarfl(isca(iscalt))
call field_set_key_double(f_id, kscmin, -grand)
call field_set_key_double(f_id, kscmax, grand)

! ---> Variables propres a la phase dispersee

do icla = 1, nclafu

  write(f_name,'(a8,i2.2)') 'nd_fuel_', icla
  write(f_label,'(a6,i2.2)') 'NG_FOL', icla
  call add_model_scalar_field(f_name, f_label, ing(icla))
  f_id = ivarfl(isca(inp(icla)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, rinfin)

enddo

do icla = 1, nclafu

  write(f_name,'(a10,i2.2)') 'yfol_fuel_', icla
  write(f_label,'(a8,i2.2)') 'YFOL_FOL', icla
  call add_model_scalar_field(f_name, f_label, iyfol(icla))
  f_id = ivarfl(isca(iyfol(icla)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! Set min and max clipping
  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 4.d-1)

enddo

do icla = 1, nclafu

  write(f_name,'(a9,i2.2)') 'hlf_fuel_', icla
  write(f_label,'(a7,i2.2)') 'HLF_FOL', icla
  call add_model_scalar_field(f_name, f_label, ih2(icla))
  f_id = ivarfl(isca(ih2(icla)))

  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)
  ! min clipping is set later, as h02fol is known after reading fuel data
  ! call field_set_key_double(f_id, kscmin, h02fol)
  call field_set_key_double(f_id, kscmax, +grand)

enddo

! Variables specific to gas phase

call add_model_scalar_field('vap_fraction', 'Fr_VAP', ifvap)
f_id = ivarfl(isca(ifvap))

call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 1.d0)

! Oxydants

if (noxyd .ge. 2) then

  call add_model_scalar_field('oxyd2_fraction', 'Fr_OXYD2', if4m)
  f_id = ivarfl(isca(if4m))

  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

endif

if (noxyd .ge. 3) then

  call add_model_scalar_field('oxyd3_fraction', 'Fr_OXYD3', if5m)
  f_id = ivarfl(isca(if5m))

  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

endif

! Combustion heterogene

call add_model_scalar_field('het_fraction', 'Fr_HET', if7m)
f_id = ivarfl(isca(if7m))

call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 1.d0)

! Variance

call add_model_scalar_field('cb_variance', 'Var_CB', ifvp2m)
f_id = ivarfl(isca(ifvp2m))

call field_set_key_double(f_id, kscmin, 0.d0)
call field_set_key_double(f_id, kscmax, 0.25d0)

! CO2 transport
if (ieqco2 .ge. 1) then

  call add_model_scalar_field('co2_fraction', 'FR_CO2', iyco2)
  f_id = ivarfl(isca(iyco2))

  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

endif

!  NOx transport: HCN, NOx and Hox

if (ieqnox .eq. 1) then

  call add_model_scalar_field('hcn_fraction', 'FR_HCN', iyhcn)
  f_id = ivarfl(isca(iyhcn))

  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  call add_model_scalar_field('no_fraction', 'FR_NO', iyno)
  f_id = ivarfl(isca(iyno))

  call field_set_key_double(f_id, kscmin, 0.d0)
  call field_set_key_double(f_id, kscmax, 1.d0)

  call add_model_scalar_field('ox_enthalpy', 'Enth_Ox', ihox)
  f_id = ivarfl(isca(ihox))

  call field_set_key_double(f_id, kscmin, -grand)
  call field_set_key_double(f_id, kscmax, +grand)

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
