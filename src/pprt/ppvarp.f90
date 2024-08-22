!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine ppvarp

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES POSITIONS DES VARIABLES SELON
!   LE TYPE DE PHYSIQUE PARTICULIERE
! REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use atchem
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
use cs_coal_incl
use ppcpfu
use atincl
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer          f_id, scal_id
integer          kscmin, kscmax

double precision scmaxp, scminp

!===============================================================================
! Interfaces
!===============================================================================

interface

  subroutine cs_cf_add_variable_fields() &
    bind(C, name='cs_cf_add_variable_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_cf_add_variable_fields

  subroutine cs_elec_add_variable_fields()  &
    bind(C, name='cs_elec_add_variable_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_elec_add_variable_fields

  subroutine cs_field_pointer_map_gas_mix()  &
    bind(C, name='cs_field_pointer_map_gas_mix')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_field_pointer_map_gas_mix

  subroutine cs_ctwr_add_variable_fields()  &
    bind(C, name='cs_ctwr_add_variable_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_ctwr_add_variable_fields

  subroutine cs_rad_transfer_add_variable_fields()  &
    bind(C, name='cs_rad_transfer_add_variable_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_rad_transfer_add_variable_fields

  subroutine cs_atmo_add_variable_fields() &
    bind(C, name='cs_atmo_add_variable_fields')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_atmo_add_variable_fields

end interface

!===============================================================================

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)

! 1. Gas combustion
!------------------

if (     ippmod(icod3p).ge.0   &
    .or. ippmod(islfm ).ge.0   &
    .or. ippmod(icoebu).ge.0   &
    .or. ippmod(icolwc).ge.0) then
  call covarp
endif

! Number of Diracs for LWC model

if (     ippmod(icolwc).eq.0  &
    .or. ippmod(icolwc).eq.1) then
  ndirac = 2

else if (      ippmod(icolwc).eq.2  &
         .or.  ippmod(icolwc).eq.3) then
  ndirac = 3

else if (     ippmod(icolwc).eq.4  &
         .or. ippmod(icolwc).eq.5) then
  ndirac = 4
endif

! 2. Coal combustion
!-------------------

if (ippmod(iccoal).ge.0) then
  call cs_coal_varpos
endif

! 3. Compressible model
!----------------------

if (ippmod(icompf).ge.0) then
  call cs_cf_add_variable_fields

  ! Fortran field mappings
  call map_variable_field_try('total_energy', ienerg)
  call map_variable_field_try('temperature', itempk)
endif

! 4. Electric arcs model
!-----------------------

if (ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1) then
  call cs_elec_add_variable_fields
endif

! 6. Atmospheric model
!---------------------

if (ippmod(iatmos).ge.0) then
  call cs_atmo_add_variable_fields()

  call init_chemistry_pointers()
  ! Update scalar id in Fortran no need in c version
  if (ippmod(iatmos).eq.2) then
    call field_get_id('ym_water', f_id)
    call field_get_key_int_by_name(f_id, "scalar_id", iymw)

    call field_get_id('number_of_droplets', f_id)
    call field_get_key_int_by_name(f_id, "scalar_id", intdrp)
  end if
endif

! 7. Cooling towers model
!------------------------

if (ippmod(iaeros).ge.0) then
  call cs_ctwr_add_variable_fields
endif

! 8. Gas mixtures modelling
!--------------------------

if (ippmod(igmix).ge.0) then

  if (ippmod(icompf).lt.0) then
    itherm = 2
    call add_model_scalar_field('enthalpy', 'Enthalpy', ihm)
    iscalt = ihm
  else
    call field_set_key_int(ivarfl(isca(itempk)), kivisl, 0)
  endif

  call field_set_key_int(ivarfl(isca(iscalt)), kivisl, 0)

  ! Clipping of mass fractions
  scminp = 0.d0
  scmaxp = 1.d0

  if (ippmod(igmix).lt.5) then
    call add_model_scalar_field('y_o2', 'Y_O2', scal_id)
    f_id = ivarfl(isca(scal_id))
    call gas_mix_add_species(f_id)
    call field_set_key_int(f_id, kivisl, 0)
    call field_set_key_double(f_id, kscmin, scminp)
    call field_set_key_double(f_id, kscmax, scmaxp)

    call add_model_scalar_field('y_n2', 'Y_N2', scal_id)
    f_id = ivarfl(isca(scal_id))
    call gas_mix_add_species(f_id)
    call field_set_key_int(f_id, kivisl, 0)
    call field_set_key_double(f_id, kscmin, scminp)
    call field_set_key_double(f_id, kscmax, scmaxp)

    if (ippmod(igmix).eq.3) then
      call add_model_scalar_field('y_he', 'Y_He', scal_id)
      f_id = ivarfl(isca(scal_id))
      call gas_mix_add_species(f_id)
      call field_set_key_int(f_id, kivisl, 0)
      call field_set_key_double(f_id, kscmin, scminp)
      call field_set_key_double(f_id, kscmax, scmaxp)

    elseif (ippmod(igmix).eq.4) then
      call add_model_scalar_field('y_h2', 'Y_H2', scal_id)
      f_id = ivarfl(isca(scal_id))
      call gas_mix_add_species(f_id)
      call field_set_key_int(f_id, kivisl, 0)
      call field_set_key_double(f_id, kscmin, scminp)
      call field_set_key_double(f_id, kscmax, scmaxp)

    endif
  else ! ippmod(igmix).eq.5

    call add_model_scalar_field('y_n2', 'Y_N2', scal_id)
    f_id = ivarfl(isca(scal_id))
    call gas_mix_add_species(f_id)
    call field_set_key_int(f_id, kivisl, 0)
    call field_set_key_double(f_id, kscmin, scminp)
    call field_set_key_double(f_id, kscmax, scmaxp)

    call add_model_scalar_field('y_he', 'Y_He', scal_id)
    f_id = ivarfl(isca(scal_id))
    call gas_mix_add_species(f_id)
    call field_set_key_int(f_id, kivisl, 0)
    call field_set_key_double(f_id, kscmin, scminp)
    call field_set_key_double(f_id, kscmax, scmaxp)

  endif

  ! MAP to C API
  call cs_field_pointer_map_gas_mix
endif

! Radiative transfer
call cs_rad_transfer_add_variable_fields()

return
end subroutine
