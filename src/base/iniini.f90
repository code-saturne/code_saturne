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

!-------------------------------------------------------------------------------

!> \file iniini.f90
!> \brief Commons default initialization before handing over the user.
!>
!------------------------------------------------------------------------------

subroutine iniini () &
  bind(C, name='cs_f_iniini')

!===============================================================================
! Module files
!===============================================================================

use atincl
use paramx
use cstnum
use dimens
use numvar
use optcal
use cstphy
use entsor
use pointe
use parall
use period
use ppincl
use ppcpfu
use mesh
use field
use vof
use radiat
use ctincl
use cfpoin
use vof
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer        ii, iscal

procedure() :: ledevi, tstjpe

!===============================================================================

!===============================================================================
! 0. Global field keys
!===============================================================================

call field_get_key_id("label", keylbl)
call field_get_key_id('log', keylog)
call field_get_key_id('post_vis', keyvis)

call field_get_key_id("inner_mass_flux_id", kimasf)
call field_get_key_id("boundary_mass_flux_id", kbmasf)

call field_get_key_id("diffusivity_id", kivisl)
call field_get_key_id("diffusivity_ref", kvisl0)

call field_get_key_id("is_temperature", kscacp)

call field_get_key_id("gradient_weighting_id", kwgrec)

call field_get_key_id("source_term_prev_id", kstprv)
call field_get_key_id("source_term_id", kst)

call field_get_key_id("turbulent_schmidt", ksigmas)

icrom = -1
ibrom = -1

!===============================================================================
! Map Fortran pointers to C global data
!===============================================================================

call atmo_init
call time_step_init
call time_step_options_init
call thermal_model_init
call turb_model_init
call turb_rans_model_init
call turb_les_model_init
call turb_hybrid_model_init
call turb_model_constants_init
call wall_functions_init
call physical_constants_init
call fluid_properties_init
call space_disc_options_init
call time_scheme_options_init
call velocity_pressure_options_init
call restart_auxiliary_options_init
call turb_reference_values_init
call radiat_init
call ctwr_properties_init
call cf_model_init
call vof_model_init

!===============================================================================
! Get mesh metadata.
!===============================================================================

call ledevi(iperio)
call tstjpe(iperio)

!===============================================================================
! Position of variables in numvar.f90
!===============================================================================

! Initialize mappings of field ids

do ii = 1, nvarmx
  ivarfl(ii) = -1
enddo

! Scalar to variable mappings

do iscal = 1, nscamx
  isca  (iscal) = 0
  iscapp(iscal) = 0
enddo

return
end subroutine
