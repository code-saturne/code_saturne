!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file usporo.f90
!>
!> \brief This function computes the porosity (volume factor \f$ \epsilon \f$
!> when porosity module is activated (iporos = 1 in cs_user_parameters.f90).
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

integer          iel, ii, jj
double precision x, pormin, pormax, hc, ll, dhc

double precision, dimension(:), pointer :: porosi

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END
!===============================================================================

! Retrieve porosity field
call field_get_val_s(ipori, porosi)

! Example: fixe a linear by part porosity profile

do iel = 1, ncel
  x = xyzcen(1,iel)
  if (x.le.(ll/2.d0)) then
    hc = 1.d0 - 2.d0*dhc*x/ll
  else
    hc = 1.d0 - dhc + dhc*(2.d0*x-ll)/ll
  endif

  porosi(iel) = hc

  pormin = min(pormin,porosi(iel))
  pormax = max(pormax,porosi(iel))
enddo

! Periodicity and parallelism treatment
if (irangp.ge.0) then
  call parmax (pormax)
  call parmin (pormin)
endif

if (iperio.eq.1.or.irangp.ge.0) then
  call synsca(porosi)
endif

return
end subroutine usporo
