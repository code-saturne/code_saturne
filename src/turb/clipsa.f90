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

subroutine clipsa &
!================

 ( ncel )

!===============================================================================
! Purpose:
! --------

! Clipping of nusa for the Spalart-Allmaras model

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncel             ! i  ! <-- ! number of cells                                !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstphy
use cstnum
use entsor
use field
use optcal
use parall
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          ncel

! Local variables

integer          iclpnu(1),iel, iclmx(1)
integer          kclipp, clip_nusa_id
double precision xnu, var, epz2
double precision vmin(1), vmax(1)

double precision, dimension(:), pointer :: cvar_nusa
double precision, dimension(:), pointer :: cpro_nusa_clipped

!===============================================================================

call field_get_val_s(ivarfl(inusa), cvar_nusa)

call field_get_key_id("clipping_id", kclipp)

! Postprocess clippings?
call field_get_key_int(ivarfl(inusa), kclipp, clip_nusa_id)
if (clip_nusa_id.ge.0) then
  call field_get_val_s(clip_nusa_id, cpro_nusa_clipped)
endif

! Une petite valeur pour eviter des valeurs exactement nulles.

epz2 = epzero**2
iclmx(1) = 0
!===============================================================================
! ---> Stockage Min et Max pour log
!===============================================================================

vmin(1) =  grand
vmax(1) = -grand
do iel = 1, ncel
  var = cvar_nusa(iel)
  vmin(1) = min(vmin(1),var)
  vmax(1) = max(vmax(1),var)
enddo

if (clip_nusa_id.ge.0) then
  do iel = 1, ncel
    cpro_nusa_clipped(iel) = 0.d0
  enddo
endif

!===============================================================================
! ---> Clipping "standard" NUSA>0
!===============================================================================

iclpnu(1) = 0
do iel = 1, ncel
  xnu = cvar_nusa(iel)
  if (xnu.lt.0.d0) then
    if (clip_nusa_id.ge.0) then
      cpro_nusa_clipped(iel) = - xnu
    endif
    iclpnu(1) = iclpnu(1) + 1
    cvar_nusa(iel) = 0.d0
  endif
enddo

call log_iteration_clipping_field(ivarfl(inusa), iclpnu(1), 0, vmin, vmax,iclpnu(1), iclmx(1))

return

end subroutine
