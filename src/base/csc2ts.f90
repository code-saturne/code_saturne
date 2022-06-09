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

!> \file csc2ts.f90
!> \brief Code-code coupling with source terms.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     ncecpl        number of coupling
!> \param[in]     lcecpl
!> \param[in]     f_id          field index
!> \param[in]     f_dim         field dimension
!> \param[in]     rvcpce        distant variable array
!> \param[in,out] crvexp        explicit source term
!______________________________________________________________________________

subroutine csc2ts &
 ( ncecpl,  lcecpl , f_id, f_dim , rvcpce , &
   crvexp )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use entsor
use optcal
use cstphy
use parall
use period
use cplsat
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          ncecpl
integer          lcecpl(ncecpl)
integer          f_id, f_dim

double precision crvexp(f_dim, ncelet)
double precision rvcpce(f_dim, ncecpl)

! Local variables

integer          isou
integer          ipt    , ielloc
double precision xdis   , xloc   , xtau   , rovtau
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: cvara_s
double precision, dimension(:,:), pointer :: cvara_v

!----------------------------------------------------------------------------------

call field_get_val_s(icrom, crom)

xtau = 100.d0*dtref

! For scalars
if (f_dim.eq.1) then
  call field_get_val_prev_s(f_id, cvara_s)

  do ipt = 1, ncecpl

    ielloc = lcecpl(ipt)

    rovtau = cell_f_vol(ielloc)*crom(ielloc)/xtau

    xdis = rvcpce(1, ipt)
    xloc = cvara_s(ielloc)
    crvexp(1,ielloc) = crvexp(1,ielloc) + rovtau*(xdis-xloc) ! TODO Should be implicit to remove constrainte on the time step...

  enddo

  ! For vectors and tensors
else
  call field_get_val_prev_v(f_id, cvara_v)

  do ipt = 1, ncecpl

    ielloc = lcecpl(ipt)

    rovtau = cell_f_vol(ielloc)*crom(ielloc)/xtau

    do isou = 1, f_dim
      xdis = rvcpce(isou,ipt)
      xloc = cvara_v(isou,ielloc)
      crvexp(isou,ielloc) = crvexp(isou,ielloc) + rovtau*(xdis-xloc) ! TODO Should be implicit to remove constrainte on the time step...
    enddo

  enddo
endif

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
