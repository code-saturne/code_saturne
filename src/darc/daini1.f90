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

!===============================================================================
! Function:
! ---------

!> \file daini1.f90
!>
!> \brief Initialize global settings for darcy module.
!>

subroutine daini1
!================

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
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
use ppincl
use ihmpre
use darcy_module
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

type(var_cal_opt) :: vcopt
integer iscal, key_decay
double precision decay_rate

!=============================================================================

interface

  !---------------------------------------------------------------------------

  ! Interface to C function defining field keys for Ground water flow module.

  subroutine cs_gwf_parameters_define_field_keys()  &
    bind(C, name='cs_gwf_parameters_define_field_keys')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gwf_parameters_define_field_keys

  !---------------------------------------------------------------------------

end interface

!===============================================================================

!===============================================================================
! 1. Set options
!===============================================================================

call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

! Standard gradient calculation option
vcopt%iwgrec = 0

call field_set_key_struct_var_cal_opt(ivarfl(ipr), vcopt)

imrgra = 1
imvisf = 1

! Mass flux reconstruction option
irecmf = 1

! Reference density and pressure
ro0 = 1.d0
p0 = 0.d0

! Set permeability and dispersion to isotropic
darcy_anisotropic_permeability = 0
darcy_anisotropic_dispersion = 0

! Steady flow
darcy_unsteady = 0

! Convergence criteron of the Newton scheme over pressure
darcy_convergence_criterion = 0

! By default, gravity is not needed
darcy_gravity = 0

! Default gravity direction is z
darcy_gravity_x = 0.d0
darcy_gravity_y = 0.d0
darcy_gravity_z = 1.d0

! Unsaturated zone is taken into account by default
darcy_unsaturated = 1

! Definition of sorption parameters
call cs_gwf_parameters_define_field_keys

if (iihmpr.eq.1) then
  call uidai1(darcy_anisotropic_permeability,     &
              darcy_anisotropic_dispersion,       &
              darcy_unsteady,                     &
              darcy_gravity,                      &
              darcy_unsaturated)
endif

!===============================================================================
! 2. Define user settings
!===============================================================================

call user_darcy_ini1
!===================

return
end subroutine
