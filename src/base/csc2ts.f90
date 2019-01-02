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
!> \param[in]     vela          variable value at time step beginning
!> \param[out]    crvexp        explicit source term
!> \param[in]     rvcpce
!______________________________________________________________________________

subroutine csc2ts &
 ( ncecpl,  lcecpl ,                                              &
   vela   , crvexp , rvcpce )

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

double precision crvexp(3,ncelet)
double precision rvcpce(3,ncecpl)
double precision vela(3,ncelet)

! Local variables

integer          isou
integer          ipt    , ielloc
double precision xdis   , xloc   , xtau   , rovtau
double precision, dimension(:), pointer ::  crom
!----------------------------------------------------------------------------------


call field_get_val_s(icrom, crom)

xtau = 100.d0*dtref

do ipt = 1, ncecpl

  ielloc = lcecpl(ipt)

  rovtau = cell_f_vol(ielloc)*crom(ielloc)/xtau

  do isou = 1, 3
    xdis = rvcpce(isou,ipt)
    xloc = vela(isou,ielloc)
    crvexp(isou,ielloc) = crvexp(isou,ielloc) + rovtau*(xdis-xloc)
  enddo

enddo

!--------
! Formats
!--------

!----
! End
!----

return

end subroutine
