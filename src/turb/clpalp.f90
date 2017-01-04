!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
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
double precision vmin(1), vmax(1), var

double precision, dimension(:), pointer :: cvar_al

!===============================================================================

call field_get_val_s(f_id, cvar_al)

!===============================================================================
!  ---> Stockage Min et Max pour listing
!===============================================================================

vmin(1) =  grand
vmax(1) = -grand
do iel = 1, ncel
  var = cvar_al(iel)
  vmin(1) = min(vmin(1),var)
  vmax(1) = max(vmax(1),var)
enddo

! ---> Clipping (modif pour eviter les valeurs exactement nulles)

iclpmn(1) = 0
iclpmx(1) = 0
do iel = 1, ncel
  if (cvar_al(iel).lt.alpha_min(iel)) then
    iclpmn(1) = iclpmn(1) + 1
    cvar_al(iel) = alpha_min(iel)
  elseif(cvar_al(iel).gt.1.d0) then
    iclpmx(1) = iclpmx(1) + 1
    cvar_al(iel) = 1.d0
  endif
enddo

call log_iteration_clipping_field(f_id, iclpmn(1), iclpmx(1), vmin, vmax,iclpmn(1), iclpmx(1))

return

end
