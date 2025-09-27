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

!-------------------------------------------------------------------------------

!> \file ppinii.f90
!> \brief Default initialization of specific modules
!> (only non-map fortran common variables of modules)
!>
!------------------------------------------------------------------------------

subroutine ppinii () &
  bind(C, name='cs_f_ppinii')

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use cstphy
use ppincl
use atincl
use field

!===============================================================================

implicit none

!===============================================================================

! Mappings to C

call pp_models_init

!===============================================================================
! 3. Global variables for atmospheric flows (module atincl.f90)
!===============================================================================

! ------------------------------------
! 1d radiative transfer model:
! ----------------------------

! iatra1 -->  flag for the use of the 1d atmo radiative model

iatra1 = 0
! key id for optimal interpolation

call field_get_key_id("opt_interp_id", kopint)

! --> Initialisation for the aerosol chemistry model:

! Default values (climatic ones) for radiative transfer and
! aerosols
aod_o3_tot=0.20d0
aod_h2o_tot=0.10d0
gaero_o3=0.66d0
gaero_h2o=0.64d0
piaero_o3=0.84d0
piaero_h2o=0.84d0
black_carbon_frac=0.d0
zaero = 6000d0

return
end subroutine ppinii
