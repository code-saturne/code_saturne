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
!> \file atprop.f90
!> \brief Add if needed the variables fields for temperature and liquid water
!>
!> \brief Add if needed the variable field for temperature and liquid water \n
!>    Nota : ippmod(iatmos) = 1 --> Dry atmosphere, = 2 --> Humid atmosphere
subroutine atprop

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
use field
use atincl

!===============================================================================

implicit none

! Local variables

integer          f_id, idim1, ityloc, itycat

!===============================================================================

!===============================================================================
! 1. POSITIONNEMENT DES PROPRIETES
!    Atmospheric modules:  dry and humid atmosphere
!===============================================================================

! Temperature (ippmod(iatmos) = 1 or 2)

if (ippmod(iatmos).ge.1) then
  call add_property_field_1d('real_temperature', 'RealTemp', itempc)

  call add_boundary_property_field_owner('non_neutral_scalar_correction',  &
                                         'Non Neutral Scalar Correction', f_id)
  call field_set_key_int(f_id, keylog, 0)
endif

! Liquid water content (ippmod(iatmos) = 2)

if (ippmod(iatmos).eq.2) then
  call add_property_field_1d('liquid_water', 'LiqWater', iliqwt)

  ! sedimentation and deposition model enabled
  if (modsedi.eq.1.and.moddep.gt.0) then
    idim1  = 1
    itycat = FIELD_INTENSIVE + FIELD_PROPERTY
    ityloc = 3 ! boundary faces

    ! wall friction velocity if not already created
    call field_find_or_create('ustar', itycat, ityloc, idim1, f_id)

    ! boundary roughness
    call add_boundary_property_field_owner('boundary_roughness',  &
                                           'Boundary Roughness', f_id)

    ! boundary thermal roughness
    call add_boundary_property_field_owner('boundary_thermal_roughness',  &
                                           'Boundary Thermal Roughness', f_id)
  endif
endif

return
end subroutine atprop
