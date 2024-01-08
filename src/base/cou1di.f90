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

!> \brief Input data for 1D wall thermal coupling

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine cou1di()  &
 bind(C, name='cs_f_cou1di')

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use cstphy
use entsor
use pointe
use field
use radiat
use cs_c_bindings
use ppincl, only: icondb

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

! Arguments

! Local variables


integer          ii , ivar
integer          ifac
integer          icldef
logical          update_bnd_temp
double precision, dimension(:), pointer :: b_temp

integer, dimension(:), pointer :: ifpt1d
double precision, dimension(:), pointer :: tppt1d

integer, pointer, dimension(:,:) :: icodcl
double precision, pointer, dimension(:,:,:) :: rcodcl

!===============================================================================

! Get the 1D wall thermal module arrays
call cs_1d_wall_thermal_get_faces(ifpt1d)
call cs_1d_wall_thermal_get_temp(tppt1d)

! Update boundary temperature field for radiative transfer or wall condensation

update_bnd_temp = ((iirayo.ge.1).or.(icondb.ge.0)).and.(nfpt1d.gt.0)

if (update_bnd_temp) then

  call field_get_val_s(itempb, b_temp)

  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    if (itypfb(ifac).eq.iparoi.or.itypfb(ifac).eq.iparug) then
      if (itpscl.eq.2) then
        b_temp(ifac) = tppt1d(ii) - tkelvi
      else
        b_temp(ifac) = tppt1d(ii)
      endif
    endif
  enddo

endif

! Sans specification, une face couplee est une face de type paroi FIXME, ca pose
! probleme si on a une couleur pour plusieurs types de face, paroi + autre
! Ces conditions peuvent etre ecrasees sur la face est couplee au flux radiatif dans
! cs_user_radiative_transfer_bcs.

icldef = 5

ivar = isca(iscalt)

call field_build_bc_codes_all(icodcl, rcodcl) ! Get map

do ii = 1, nfpt1d

   ifac = ifpt1d(ii)

   if ((icodcl(ifac,ivar).ne.1 .and.                           &
        icodcl(ifac,ivar).ne.5 .and.                           &
        icodcl(ifac,ivar).ne.6).and.                           &
       (itypfb(ifac).eq.iparoi .or. itypfb(ifac).eq.iparug))   &
     icodcl(ifac,ivar) = icldef

   rcodcl(ifac,ivar,1) = tppt1d(ii)
   rcodcl(ifac,ivar,2) = rinfin
   rcodcl(ifac,ivar,3) = 0.d0

enddo

! Conversion eventuelle temperature -> enthalpie

if (itherm.eq.2) then

  do ii = 1, nfpt1d
    ifac = ifpt1d(ii)
    icodcl(ifac,ivar) = - icodcl(ifac,ivar)
  enddo

endif

end subroutine


