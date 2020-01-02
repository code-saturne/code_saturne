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
! Function :
! ----------
!> \file clpalp.f90
!> \brief Clipping of alpha in the framwork of the Rij-EBRSM model.
!>        Also called for alpha of scalars for EB-DFM.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           nom           role
!______________________________________________________________________________!
!> \param[in]     f_id          field id of alpha variable
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     alpha_min     minimum acceptable value for alpha
!______________________________________________________________________________!


subroutine clpalp &
 ( f_id, ncelet , ncel  , alpha_min)

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use entsor
use numvar
use cstnum
use parall
use cs_c_bindings
use field

!===============================================================================

implicit none

! Arguments

integer          f_id, ncelet, ncel
double precision alpha_min(ncelet)

! VARIABLES LOCALES

integer          iel
integer          iclpmn(1), iclpmx(1)
integer          kclipp, clip_a_id
double precision vmin(1), vmax(1), var

double precision, dimension(:), pointer :: cvar_al
double precision, dimension(:), pointer :: cpro_a_clipped

!===============================================================================

call field_get_val_s(f_id, cvar_al)

call field_get_key_id("clipping_id", kclipp)

! Postprocess clippings?
call field_get_key_int(f_id, kclipp, clip_a_id)
if (clip_a_id.ge.0) then
  call field_get_val_s(clip_a_id, cpro_a_clipped)
endif

!===============================================================================
!  ---> Stockage Min et Max pour log
!===============================================================================

vmin(1) =  grand
vmax(1) = -grand
do iel = 1, ncel
  var = cvar_al(iel)
  vmin(1) = min(vmin(1),var)
  vmax(1) = max(vmax(1),var)
enddo

do iel = 1, ncel
  if (clip_a_id.ge.0) &
    cpro_a_clipped(iel) = 0.d0
enddo

! ---> Clipping (modif pour eviter les valeurs exactement nulles)

iclpmn(1) = 0
iclpmx(1) = 0
do iel = 1, ncel
  if (cvar_al(iel).lt.alpha_min(iel)) then
    if (clip_a_id.ge.0) &
      cpro_a_clipped(iel) = alpha_min(iel) - cvar_al(iel)
    iclpmn(1) = iclpmn(1) + 1
    cvar_al(iel) = alpha_min(iel)
  elseif(cvar_al(iel).gt.1.d0) then
    if (clip_a_id.ge.0) &
      cpro_a_clipped(iel) = cvar_al(iel) - 1.d0
    iclpmx(1) = iclpmx(1) + 1
    cvar_al(iel) = 1.d0
  endif
enddo

call log_iteration_clipping_field(f_id, iclpmn(1), iclpmx(1), vmin, vmax,iclpmn(1), iclpmx(1))

return

end
