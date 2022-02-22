!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file tspdcv.f90
!>
!> \brief This subroutine computes the explicit contribution of
!> headlosses terms.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     vela          velocity at the previous time step
!> \param[in]     ckupdc        work array for the head loss
!> \param[out]    trav          right hand side
!_______________________________________________________________________________


subroutine tspdcv &
 ( ncepdp , icepdc ,                                              &
   vela   ,                                                       &
   ckupdc , trav   )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use mesh
use field
!===============================================================================

implicit none

! Arguments

integer          ncepdp
integer          icepdc(ncepdp)

double precision ckupdc(6,ncepdp)
double precision trav(3,ncelet)
double precision vela  (3  ,ncelet)

! Local variables

integer          iel   , ielpdc
integer          key_t_ext_id
integer          iroext
double precision romvom, vit1  , vit2  , vit3
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23
double precision, dimension(:), pointer :: crom, croma

!===============================================================================

! Time extrapolation?
call field_get_key_id("time_extrapolated", key_t_ext_id)

call field_get_val_s(icrom, crom)

call field_get_key_int(icrom, key_t_ext_id, iroext)
if (iroext.gt.0.and.isno2t.gt.0) then
  call field_get_val_prev_s(icrom, croma)
endif

!     La matrice est toujours "implicite"

do ielpdc = 1, ncepdp

  iel    = icepdc(ielpdc)
  romvom =-crom(iel)*cell_f_vol(iel)
  cpdc11 = ckupdc(1,ielpdc)
  cpdc22 = ckupdc(2,ielpdc)
  cpdc33 = ckupdc(3,ielpdc)
  cpdc12 = ckupdc(4,ielpdc)
  cpdc23 = ckupdc(5,ielpdc)
  cpdc13 = ckupdc(6,ielpdc)
  vit1   = vela(1,iel)
  vit2   = vela(2,iel)
  vit3   = vela(3,iel)

  trav(1,ielpdc) = romvom*(cpdc11*vit1 + cpdc12*vit2 + cpdc13*vit3)
  trav(2,ielpdc) = romvom*(cpdc12*vit1 + cpdc22*vit2 + cpdc23*vit3)
  trav(3,ielpdc) = romvom*(cpdc13*vit1 + cpdc23*vit2 + cpdc33*vit3)

enddo

return
end subroutine
