!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file cs_user_porosity.f90
!>
!> \brief This function computes the porosity (volume factor \f$ \epsilon \f$
!> when porosity module is activated (iporos = 1 in cs_user_parameters.f90).
!>
!> See \subpage cs_porosity for examples.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine usporo
!================

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstnum
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Local variables

!< [loc_var_dec]
integer          iel, ii, jj
double precision x, pormin, pormax, hc, ll, dhc

double precision, dimension(:), pointer :: cpro_porosi
!< [loc_var_dec]

!===============================================================================

!< [init]
! Retrieve porosity field
call field_get_val_s(ipori, cpro_porosi)
!< [init]

!< [example_1]
! Example: fixe a linear by part porosity profile

ll = 2.d0
dhc = 1.d0

do iel = 1, ncel
  x = xyzcen(1,iel)
  if (x.le.(ll/2.d0)) then
    hc = 1.d0 - 2.d0*dhc*x/ll
  else
    hc = 1.d0 - dhc + dhc*(2.d0*x-ll)/ll
  endif

  cpro_porosi(iel) = hc

  pormin = min(pormin,cpro_porosi(iel))
  pormax = max(pormax,cpro_porosi(iel))
enddo
!< [example_1]

!< [parallelism]
! Periodicity and parallelism treatment
if (irangp.ge.0) then
  call parmax (pormax)
  call parmin (pormin)
endif
!< [parallelism]

return
end subroutine usporo
